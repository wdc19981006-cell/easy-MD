#!/usr/bin/env python
"""Resumable easy-md protein-ligand docking pipeline."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import traceback
from pathlib import Path
from typing import Any

SCRIPTS = Path(__file__).resolve().parent
SKILL_ROOT = SCRIPTS.parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from analysis.plip_analysis import analyze  # noqa: E402
from common import (  # noqa: E402
    WORK_MARKER,
    atomic_write_json,
    configured_path,
    load_config,
    new_status,
    safe_slug,
    update_stage,
)
from docking.pose_rank import rank_poses  # noqa: E402
from ligand.ligand_search import download_ligand_sdf, search_ligand  # noqa: E402
from protein.normalize import normalize_protein  # noqa: E402
from protein.pdb_search import download_pdb, enrich_entry, search_structures  # noqa: E402
from pdb_selection import rank_pdb_entries  # noqa: E402
from visualization.plip_pymol_native import render_native  # noqa: E402
from visualization.pymol_generator import generate_pml, run_pymol  # noqa: E402
from visualization.visualization_validate import validate as validate_visualization  # noqa: E402


def _run(command: list[str], cwd: Path, log: Path | None = None) -> subprocess.CompletedProcess[str]:
    process = subprocess.run(command, cwd=str(cwd), capture_output=True, text=True)
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            f"COMMAND: {json.dumps(command)}\n\nSTDOUT\n{process.stdout}\n\nSTDERR\n{process.stderr}\n",
            encoding="utf-8",
        )
    if process.returncode != 0:
        raise RuntimeError(
            f"Command failed ({process.returncode}): {command[1] if len(command) > 1 else command[0]}\n"
            f"{process.stderr[-3000:]}\n{process.stdout[-3000:]}"
        )
    return process


def _copy(source: Path, target: Path) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def _layout(run_dir: Path) -> dict[str, Path]:
    return {
        "run": run_dir,
        "docking": run_dir / "docking",
        "structure": run_dir / "structure",
        "analysis": run_dir / "analysis",
        "visualization": run_dir / "visualization",
        "status": run_dir / "run_status.json",
        "metadata": run_dir / "metadata.json",
    }


def _create_run(output_root: Path, protein: str, ligand: str) -> dict[str, Path]:
    run_dir = output_root / f"{safe_slug(protein, 'protein')}_{safe_slug(ligand, 'ligand')}"
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty result directory: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)
    paths = _layout(run_dir)
    for key in ("docking", "structure", "analysis", "visualization"):
        paths[key].mkdir(parents=True, exist_ok=True)
    (run_dir / WORK_MARKER).write_text("easy-md managed work directory\n", encoding="utf-8")
    atomic_write_json(paths["status"], new_status())
    return paths


def _resolution(entry: dict[str, Any]) -> float | None:
    values = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    return float(values[0]) if values and values[0] is not None else None


def _verified_modified_site(pdb_path: Path, ligand_name: str, title: str) -> tuple[str, str] | None:
    """Return a curated covalent site marker only when the PDB itself verifies it."""
    ligand_key = ligand_name.casefold().replace("-", " ")
    title_key = title.casefold().replace("-", " ")
    expected = None
    if ligand_key in {"aspirin", "acetylsalicylic acid"} and "acetylated" in title_key:
        expected = "OAS"
    if not expected:
        return None
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if line.startswith("MODRES") and line[12:15].strip().upper() == expected:
            return expected, line[16:17].strip() or "A"
    return None


def prepare(args: argparse.Namespace) -> int:
    config_path = Path(args.config or SKILL_ROOT / "config.yaml").resolve()
    config = load_config(config_path)
    python = configured_path(config, "python") or sys.executable

    identity = normalize_protein(args.protein, args.species)
    if identity["ambiguous"]:
        print(json.dumps(identity, indent=2, ensure_ascii=False))
        print("Protein identity is ambiguous; choose a gene symbol and rerun.", file=sys.stderr)
        return 2
    if not identity["resolved"]:
        print("Protein identity could not be resolved.", file=sys.stderr)
        return 1
    selected_identity = identity["selected"]
    gene = selected_identity["symbol"]

    ligand_metadata: dict[str, Any]
    if args.smiles:
        ligand_metadata = {
            "input": args.ligand, "canonical_smiles": args.smiles,
            "isomeric_smiles": args.smiles, "source": "user SMILES",
        }
    else:
        ligand_metadata = search_ligand(args.ligand)

    paths = _create_run(Path(args.output_root).resolve(), gene, args.ligand)
    run_dir, status_path = paths["run"], paths["status"]
    update_stage(status_path, "acquisition", "RUNNING")
    try:
        _run([python, str(SCRIPTS / "preflight.py"), "--config", str(config_path), "--json"], run_dir)

        metadata_candidates: dict[str, dict[str, Any]] = {}
        if args.protein_pdb:
            pdb_id = (args.pdb_id or Path(args.protein_pdb).stem).upper()
            selected_entry = enrich_entry(pdb_id) if len(pdb_id) == 4 else {}
            original_pdb = run_dir / f"{pdb_id}.pdb"
            _copy(Path(args.protein_pdb).resolve(), original_pdb)
            ranked: list[dict[str, Any]] = []
            protein_source = "user/MCP supplied PDB"
        else:
            if args.pdb_id:
                pdb_id = args.pdb_id.upper()
                metadata_candidates = {pdb_id: enrich_entry(pdb_id)}
            else:
                metadata_candidates = search_structures(gene, int(config.get("acquisition", {}).get("pdb_candidates", 15)))
            usable = {key: value for key, value in metadata_candidates.items() if "error" not in value}
            exact_gene = {
                key: value for key, value in usable.items()
                if gene.upper() in value.get("easy_md", {}).get("gene_names", [])
                or gene.casefold() in json.dumps(
                    value.get("easy_md", {}).get("descriptions", []), ensure_ascii=False
                ).casefold()
            }
            if exact_gene:
                usable = exact_gene
            ligand_key = args.ligand.casefold().replace("-", " ")
            for entry in usable.values():
                quality = entry.setdefault("easy_md", {})
                title = str(entry.get("struct", {}).get("title", ""))
                normalized_title = title.casefold().replace("-", " ")
                if ligand_key and ligand_key in normalized_title:
                    quality["relevant_ligand_evidence"] = f"entry title: {title}"
                if quality.get("binding_site_defined") and any(
                    token in normalized_title for token in (" bound ", " complex ", " inhibitor")
                ):
                    quality["active_site_evidence"] = "resolved ligand-bound co-crystal structure"
            ranked = rank_pdb_entries(usable, config.get("pdb_scoring", {}))
            if not ranked:
                raise RuntimeError(f"No usable PDB structures found for {gene}")
            pdb_id = ranked[0]["pdb_id"]
            selected_entry = usable[pdb_id]
            if not selected_entry.get("easy_md", {}).get("is_human"):
                raise RuntimeError("No evidence that the recommended PDB is human; provide --pdb-id explicitly")
            original_pdb = run_dir / f"{pdb_id}.pdb"
            download_pdb(pdb_id, original_pdb)
            protein_source = "RCSB REST fallback"

        selected_title = str(selected_entry.get("struct", {}).get("title", ""))
        verified_site = _verified_modified_site(original_pdb, args.ligand, selected_title)

        if args.ligand_sdf:
            ligand_sdf = run_dir / "ligand_input.sdf"
            _copy(Path(args.ligand_sdf).resolve(), ligand_sdf)
            ligand_source = "user/MCP supplied SDF"
        elif args.smiles:
            ligand_sdf = None
            ligand_source = "user SMILES"
        else:
            ligand_sdf = run_dir / "ligand_input.sdf"
            download_ligand_sdf(ligand_metadata.get("cid") or args.ligand, ligand_sdf)
            ligand_source = "PubChem PUG REST fallback"

        update_stage(status_path, "acquisition", "SUCCESS", artifacts=[
            str(original_pdb), *( [str(ligand_sdf)] if ligand_sdf else [] )
        ])
    except Exception as exc:
        update_stage(status_path, "acquisition", "FAILED", error=str(exc))
        traceback.print_exc()
        return 1

    update_stage(status_path, "preparation", "RUNNING")
    try:
        protein_command = [
            python, str(SCRIPTS / "protein_prepare.py"), str(original_pdb),
            "--pdb-id", pdb_id, "--work-dir", str(run_dir), "--config", str(config_path),
        ]
        effective_chain = args.chain
        if verified_site:
            protein_command.extend(["--site-resname", verified_site[0]])
            effective_chain = effective_chain or verified_site[1]
        if effective_chain:
            protein_command.extend(["--chain", effective_chain])
        _run(protein_command, run_dir, run_dir / "protein_prepare.log")

        ligand_command = [
            python, str(SCRIPTS / "ligand_prepare.py"),
            "--name", args.ligand, "--work-dir", str(run_dir), "--config", str(config_path),
        ]
        if ligand_sdf:
            ligand_command.extend(["--sdf", str(ligand_sdf)])
        else:
            ligand_command.insert(2, args.smiles)
        _run(ligand_command, run_dir, run_dir / "ligand_prepare.log")

        box_command = [
            python, str(SCRIPTS / "box_prepare.py"), "--pdb", str(original_pdb),
            "--work-dir", str(run_dir), "--config", str(config_path), "--json",
        ]
        if args.center:
            box_command.extend(["--center", *[str(value) for value in args.center]])
            if args.size:
                box_command.extend(["--size", *[str(value) for value in args.size]])
        elif not (run_dir / "native_ligand_info.txt").is_file():
            raise RuntimeError(
                "No co-crystal ligand was found. Provide experimentally supported --center and --size; blind docking is not automatic."
            )
        box_result = _run(box_command, run_dir, run_dir / "box_prepare.log")
        box_json = json.loads(box_result.stdout[box_result.stdout.find("{"):])

        _copy(original_pdb, paths["structure"] / "original.pdb")
        _copy(run_dir / "protein_clean.pdb", paths["structure"] / "protein_clean.pdb")
        _copy(run_dir / "ligand.sdf", paths["structure"] / "ligand.sdf")
        _copy(run_dir / "receptor.pdbqt", paths["docking"] / "receptor.pdbqt")
        _copy(run_dir / "ligand.pdbqt", paths["docking"] / "ligand.pdbqt")
        _copy(run_dir / "config.txt", paths["docking"] / "config.txt")

        metadata = {
            "schema_version": 1,
            "skill": "easy-md",
            "protein": {
                "input": args.protein, "gene_symbol": gene,
                "official_name": selected_identity.get("name"), "pdb_id": pdb_id,
                "resolution": _resolution(selected_entry),
                "organism": selected_entry.get("easy_md", {}).get("organism"),
                "source": protein_source,
            },
            "ligand": {**ligand_metadata, "name": args.ligand, "source": ligand_source},
            "selection": {"ranking": ranked[:10], "candidate_count": len(metadata_candidates)},
            "site_evidence": {
                "type": "verified_modified_residue" if verified_site else "co_crystal_ligand",
                "marker": verified_site[0] if verified_site else box_json.get("source"),
                "chain": effective_chain,
            },
            "box": box_json,
            "docking": config.get("docking", {}),
            "redock": bool(args.redock),
            "scientific_warnings": ([
                "Aspirin covalently acetylates COX-2, but standard AutoDock Vina models only a noncovalent pose and score. "
                "The verified OAS 530 marker locates the pocket; this workflow does not model acyl transfer or a covalent bond."
            ] if verified_site and args.ligand.casefold() in {"aspirin", "acetylsalicylic acid"} else []),
            "config": str(config_path),
            "run_dir": str(run_dir),
        }
        atomic_write_json(paths["metadata"], metadata)
        update_stage(status_path, "preparation", "SUCCESS", artifacts=[
            str(paths["structure"] / "protein_clean.pdb"),
            str(paths["docking"] / "receptor.pdbqt"),
            str(paths["docking"] / "ligand.pdbqt"),
            str(paths["docking"] / "config.txt"),
        ])
        update_stage(status_path, "docking", "AWAITING_CONFIRMATION")
    except Exception as exc:
        update_stage(status_path, "preparation", "FAILED", error=str(exc))
        traceback.print_exc()
        return 1

    print("DOCKING CONFIRMATION")
    print(f"Protein: {gene} ({metadata['protein']['official_name']}) | Organism: {metadata['protein']['organism']}")
    print(f"PDB: {pdb_id} | Method: {ranked[0]['info']['method'] if ranked else 'user supplied'} | Resolution: {metadata['protein']['resolution']} | Chain: {effective_chain or 'auto'}")
    if ranked:
        print(f"Selection score: {ranked[0]['score']} | Evidence: {'; '.join(ranked[0]['positives'])}")
        if ranked[0]["negatives"]:
            print(f"Selection cautions: {'; '.join(ranked[0]['negatives'])}")
    print(f"Ligand: {args.ligand} | CID: {ligand_metadata.get('cid')} | Formula: {ligand_metadata.get('formula')} | SMILES: {ligand_metadata.get('isomeric_smiles') or ligand_metadata.get('canonical_smiles')}")
    print(f"Binding site: {box_json['source']} | Evidence: {metadata['site_evidence']['type']}")
    print(f"Box: center={box_json['center']} | size={box_json['size']}")
    print(f"Vina: rigid receptor / flexible ligand | exhaustiveness={metadata['docking'].get('exhaustiveness')} | num_modes={metadata['docking'].get('num_modes')} | energy_range={metadata['docking'].get('energy_range')} | runs={metadata['docking'].get('num_runs')} | seeds={metadata['docking'].get('seeds')} | redock={metadata['redock']}")
    print(f"Expected runtime: {metadata['docking'].get('expected_runtime', 'CPU and box dependent')}")
    for warning in metadata.get("scientific_warnings", []):
        print(f"SCIENTIFIC WARNING: {warning}")
    print(f"Output: {run_dir}")
    print(f"To continue after explicit approval: pipeline.py execute --run-dir \"{run_dir}\" --confirm")
    return 0


def execute(args: argparse.Namespace) -> int:
    run_dir = Path(args.run_dir).resolve()
    paths = _layout(run_dir)
    if not args.confirm:
        print("Refusing to run Vina without --confirm after user approval.", file=sys.stderr)
        return 2
    metadata = json.loads(paths["metadata"].read_text(encoding="utf-8"))
    config_path = Path(metadata["config"])
    config = load_config(config_path)
    python = configured_path(config, "python") or sys.executable
    seeds = [int(value) for value in config.get("docking", {}).get("seeds", [42, 2024, 8888])]
    obabel = configured_path(config, "obabel") or "obabel"

    update_stage(paths["status"], "docking", "RUNNING")
    try:
        _run([
            python, str(SCRIPTS / "docking_run.py"), "--config", "config.txt",
            "--work-dir", str(run_dir), "--skill-config", str(config_path), "--json",
        ], run_dir, paths["docking"] / "docking.log")
        ranking = rank_poses(run_dir, seeds, obabel)
        _copy(run_dir / "best_pose.pdbqt", paths["docking"] / "best_pose.pdbqt")
        _copy(run_dir / "scores.csv", paths["docking"] / "scores.csv")
        _copy(run_dir / "best_pose.pdb", paths["structure"] / "best_pose.pdb")
        runs_dir = paths["docking"] / "runs"
        for run in range(1, len(seeds) + 1):
            shutil.copytree(run_dir / f"run{run}", runs_dir / f"run{run}", dirs_exist_ok=True)
        if metadata.get("redock"):
            _run([
                python, str(SCRIPTS / "redocking.py"), "--pdb", str(paths["structure"] / "original.pdb"),
                "--work-dir", str(run_dir), "--config", str(config_path),
            ], run_dir, paths["docking"] / "redocking.log")
        metadata["pose_ranking"] = ranking
        atomic_write_json(paths["metadata"], metadata)
        update_stage(paths["status"], "docking", "SUCCESS", artifacts=[
            str(paths["docking"] / "best_pose.pdbqt"),
            str(paths["docking"] / "docking.log"), str(paths["docking"] / "scores.csv"),
        ])
    except Exception as exc:
        update_stage(paths["status"], "docking", "FAILED", error=str(exc))
        traceback.print_exc()
        return 1

    update_stage(paths["status"], "interaction", "RUNNING")
    try:
        best = metadata["pose_ranking"]["best"]
        analysis_result = analyze(
            paths["structure"] / "protein_clean.pdb",
            paths["docking"] / "best_pose.pdbqt",
            paths["structure"] / "ligand.sdf",
            metadata["protein"]["gene_symbol"], metadata["ligand"]["name"],
            float(best["affinity"]), paths["structure"], paths["analysis"], obabel,
        )
        two_d = subprocess.run([
            python, str(SCRIPTS / "interaction_2d.py"),
            "--json", analysis_result["summary_json"],
            "--sdf", str(paths["structure"] / "ligand.sdf"),
            "--complex", analysis_result["complex"],
            "--output-dir", str(paths["analysis"]),
            "--protein-name", metadata["protein"]["gene_symbol"],
            "--ligand-name", metadata["ligand"]["name"],
            "--ligand-resname", analysis_result["ligand_resname"],
            "--best-affinity", str(best["affinity"]),
        ], capture_output=True, text=True, cwd=str(run_dir))
        (paths["analysis"] / "interaction_2d.log").write_text(
            two_d.stdout + "\n" + two_d.stderr, encoding="utf-8"
        )
        analysis_result["two_d_status"] = "SUCCESS" if two_d.returncode == 0 else "FAILED"
        metadata["interaction"] = analysis_result
        atomic_write_json(paths["metadata"], metadata)
        update_stage(paths["status"], "interaction", "SUCCESS", artifacts=[
            analysis_result["master"], analysis_result["summary_txt"], analysis_result["xml"],
        ])
    except Exception as exc:
        update_stage(paths["status"], "interaction", "FAILED", error=str(exc))
        traceback.print_exc()
        return 1

    update_stage(paths["status"], "visualization", "RUNNING")
    visual_metadata: dict[str, Any] = {}
    try:
        master_path = Path(analysis_result["master"])
        pse = paths["visualization"] / "interaction_visualization.pse"
        pml = paths["visualization"] / "interaction_visualization.pml"
        generate_pml(master_path, Path(analysis_result["complex"]), pml, pse)
        pymol_exe = configured_path(config, "pymol")
        pymol_python = configured_path(config, "pymol_python")
        if not pymol_python and pymol_exe:
            pymol_python = str(Path(pymol_exe).with_name("python.exe"))
        if not pymol_exe or not pymol_python:
            raise RuntimeError("PyMOL executable and Python runtime are required")
        import plip
        plip_site = Path(plip.__file__).resolve().parent.parent
        native = render_native(
            Path(analysis_result["visualizer_data"]), master_path, pse,
            pymol_python, plip_site,
            int(config.get("interaction", {}).get("pymol_timeout", 120)),
        )
        source = "plip_native" if native.get("success") else "custom_pml_fallback"
        if not native.get("success"):
            fallback = run_pymol(pymol_exe, pml)
            if fallback.returncode != 0:
                raise RuntimeError(f"Native and fallback PyMOL failed: {native.get('stderr')} | {fallback.stderr}")
        validation = validate_visualization(
            master_path, pse, pml, pymol_python,
            int(config.get("validation", {}).get("min_pse_bytes", 1024)),
        )
        if not validation["success"] and source == "plip_native":
            fallback = run_pymol(pymol_exe, pml)
            source = "custom_pml_fallback_after_native_validation"
            validation = validate_visualization(master_path, pse, pml, pymol_python)
        visual_metadata = {
            "source": source, "native_attempt": native, "validation": validation,
            "interaction_count": analysis_result["interaction_count"],
        }
        atomic_write_json(paths["visualization"] / "visualization_metadata.json", visual_metadata)
        if not validation["success"]:
            raise RuntimeError(f"Visualization validation failed: {validation}")
        update_stage(paths["status"], "visualization", "SUCCESS", artifacts=[str(pse), str(pml)])
    except Exception as exc:
        visual_metadata.setdefault("error", str(exc))
        atomic_write_json(paths["visualization"] / "visualization_metadata.json", visual_metadata)
        update_stage(paths["status"], "visualization", "FAILED", error=str(exc))
        traceback.print_exc()
        return 1

    update_stage(paths["status"], "packaging", "RUNNING")
    report = validate_run(run_dir, pymol_python)
    metadata["validation"] = report
    atomic_write_json(paths["metadata"], metadata)
    if report["success"]:
        update_stage(paths["status"], "packaging", "SUCCESS", artifacts=report["required_files"])
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return 0
    update_stage(paths["status"], "packaging", "FAILED", error=json.dumps(report["missing"]))
    return 1


def validate_run(run_dir: Path, pymol_python: str | None = None) -> dict[str, Any]:
    paths = _layout(run_dir)
    required = [
        paths["docking"] / "receptor.pdbqt", paths["docking"] / "ligand.pdbqt",
        paths["docking"] / "best_pose.pdbqt", paths["docking"] / "docking.log",
        paths["analysis"] / "interaction_summary.txt", paths["analysis"] / "interaction_master.json",
        paths["visualization"] / "interaction_visualization.pse",
        paths["visualization"] / "interaction_visualization.pml",
    ]
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    return {"success": not missing, "missing": missing, "required_files": [str(path) for path in required]}


def status(args: argparse.Namespace) -> int:
    path = Path(args.run_dir).resolve() / "run_status.json"
    print(path.read_text(encoding="utf-8"))
    return 0


def validate_command(args: argparse.Namespace) -> int:
    result = validate_run(Path(args.run_dir).resolve())
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["success"] else 1


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    root.add_argument("--config")
    commands = root.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    prep.add_argument("--protein", required=True)
    prep.add_argument("--ligand", required=True)
    prep.add_argument("--output-root", required=True)
    prep.add_argument("--species", default="human")
    prep.add_argument("--pdb-id")
    prep.add_argument("--protein-pdb")
    prep.add_argument("--ligand-sdf")
    prep.add_argument("--smiles")
    prep.add_argument("--chain")
    prep.add_argument("--center", nargs=3, type=float)
    prep.add_argument("--size", nargs=3, type=float)
    prep.add_argument("--redock", action="store_true")
    prep.set_defaults(func=prepare)
    run = commands.add_parser("execute")
    run.add_argument("--run-dir", required=True)
    run.add_argument("--confirm", action="store_true")
    run.set_defaults(func=execute)
    show = commands.add_parser("status")
    show.add_argument("--run-dir", required=True)
    show.set_defaults(func=status)
    check = commands.add_parser("validate")
    check.add_argument("--run-dir", required=True)
    check.set_defaults(func=validate_command)
    return root


def main() -> int:
    args = parser().parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
