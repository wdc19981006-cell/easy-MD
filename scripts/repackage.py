#!/usr/bin/env python
"""Repackage completed easy-md runs without rerunning Vina or PLIP."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
import tempfile

from common import atomic_write_json, configured_path, load_config, update_stage
from pipeline import SKILL_ROOT, _create_run, _layout, render_interactions, validate_run, write_vina_input_paths


def _copy(source: Path, target: Path) -> None:
    if not source.is_file() or not source.stat().st_size:
        raise ValueError(f"Missing or empty required result: {source}")
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def _copy_optional(source: Path, target: Path) -> None:
    if source.is_file():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


def _coordinates(path: Path, pdbqt: bool) -> Counter:
    coords: Counter = Counter()
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line.split()[-1] if pdbqt else line[76:78].strip()
        if element in {"H", "HD", "HS"}:
            continue
        coords[tuple(float(line[i:i + 8]) for i in (30, 38, 46))] += 1
    return coords


def _relocate(value, temporary: Path, destination: Path):
    if isinstance(value, dict):
        return {key: _relocate(item, temporary, destination) for key, item in value.items()}
    if isinstance(value, list):
        return [_relocate(item, temporary, destination) for item in value]
    if isinstance(value, str) and (value == str(temporary) or value.startswith(str(temporary) + "\\")):
        return str(destination / Path(value).relative_to(temporary))
    if isinstance(value, str) and value.startswith(temporary.as_posix() + "/"):
        return str(destination / Path(value).relative_to(temporary))
    return value


def repackage(source: Path, output_root: Path, config_path: Path) -> Path:
    source = source.resolve()
    old = _layout(source)
    metadata = json.loads(old["metadata"].read_text(encoding="utf-8"))
    status = json.loads(old["status"].read_text(encoding="utf-8"))
    if any(status["stages"][stage]["status"] != "SUCCESS" for stage in ("docking", "interaction")):
        raise ValueError(f"Docking and PLIP must already be complete: {source}")
    final = output_root.resolve() / source.name
    if final.exists() or final.is_relative_to(source):
        raise FileExistsError(f"Refusing to overwrite or write inside the original run: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    config = load_config(config_path)

    with tempfile.TemporaryDirectory(prefix=".repackage-", dir=final.parent) as folder:
        # _create_run uses an internal safe name; the original result name is restored on delivery.
        paths = _create_run(Path(folder), "result", "work")
        work = paths["run"]
        for name in ("protein_clean.pdb", "ligand.sdf", "best_pose.pdb", "complex_for_visualization.pdb"):
            _copy(old["structure"] / name, paths["structure"] / name)
        for name in ("original.pdb", "original_chain_I.pdb", "ligand_input.sdf", "best_pose.sdf"):
            _copy_optional(old["structure"] / name, paths["structure"] / name)
        for name in ("receptor.pdbqt", "ligand.pdbqt", "best_pose.pdbqt"):
            candidate = old["structure"] / name
            _copy(candidate if candidate.is_file() else old["docking"] / name, paths["structure"] / name)
        best = _coordinates(paths["structure"] / "best_pose.pdbqt", pdbqt=True)
        if not best or best != _coordinates(paths["structure"] / "best_pose.pdb", pdbqt=False):
            raise ValueError(f"The analysis pose coordinates differ from the Vina pose: {source}")
        write_vina_input_paths(old["docking"] / "config.txt", paths["docking"] / "config.txt", "../structure/")
        for name in ("docking.log", "scores.csv"):
            _copy(old["docking"] / name, paths["docking"] / name)
        for name in ("docking_stats.json", "pose_ranking.json", "native_ligand_info.txt",
                     "protein_prepare.log", "ligand_prepare.log", "box_prepare.log", "redocking.log"):
            candidate = old["docking"] / name
            _copy_optional(candidate if candidate.is_file() else source / name, paths["docking"] / name)
        for name in ("redocking_info.json", "redocking_result.txt", "docked_native.pdbqt"):
            legacy = source / ("redock" if name == "docked_native.pdbqt" else "") / name
            candidate = old["docking"] / "redocking" / name
            _copy_optional(candidate if candidate.is_file() else legacy, paths["docking"] / "redocking" / name)
        docking = metadata["docking"]
        count = int(docking.get("num_runs", len(docking["seeds"])))
        for index in range(1, count + 1):
            folder_name = f"run{index}"
            for name in ("docked.pdbqt", "log.txt"):
                path = Path("runs") / folder_name / name
                candidate = old["docking"] / path
                duplicate = source / folder_name / name
                if duplicate.is_file() and duplicate.read_bytes() != candidate.read_bytes():
                    raise ValueError(f"Root and docking Vina outputs disagree: {candidate}")
                _copy(candidate, paths["docking"] / path)
        for name in ("plip_report.xml", "interaction_master.json", "interaction_summary.txt", "interaction_summary.json"):
            _copy(old["analysis"] / name, paths["analysis"] / name)
        # VisualizerData is a pickle. Only repackage trusted local easy-md runs.
        _copy_optional(old["analysis"] / "visualizer_data.pkl", paths["analysis"] / "visualizer_data.pkl")
        master = paths["analysis"] / "interaction_master.json"
        interaction = {
            **metadata["interaction"],
            "master": str(master),
            "summary_txt": str(paths["analysis"] / "interaction_summary.txt"),
            "summary_json": str(paths["analysis"] / "interaction_summary.json"),
            "xml": str(paths["analysis"] / "plip_report.xml"),
            "visualizer_data": str(paths["analysis"] / "visualizer_data.pkl"),
            "complex": str(paths["structure"] / "complex_for_visualization.pdb"),
            "ligand_sdf": str(paths["structure"] / "ligand.sdf"),
            "interaction_count": json.loads(master.read_text(encoding="utf-8"))["interaction_count"],
        }
        interaction.pop("two_d_status", None)
        metadata.update(run_dir=str(work), config=str(config_path), interaction=interaction)
        atomic_write_json(paths["metadata"], metadata)
        atomic_write_json(paths["status"], status)
        outputs = render_interactions(paths, interaction, config)
        report = validate_run(work)
        if not report["success"]:
            raise RuntimeError(f"Missing packaged files: {report['missing']}")
        stage_paths = {
            "acquisition": [p for p in paths["structure"].iterdir()
                            if p.name.startswith("original") or p.suffix == ".sdf" and p.name.startswith("ligand")],
            "preparation": [paths["structure"] / name for name in
                            ("protein_clean.pdb", "ligand.sdf", "receptor.pdbqt", "ligand.pdbqt")],
            "docking": [paths["docking"] / "scores.csv", paths["structure"] / "best_pose.pdbqt"],
            "interaction": [interaction[key] for key in ("master", "summary_txt", "xml")],
            "visualization": outputs,
            "packaging": report["required_files"],
        }
        for stage, files in stage_paths.items():
            update_stage(paths["status"], stage, "SUCCESS",
                         artifacts=[str(final / Path(path).relative_to(work)) for path in files])
        metadata["validation"] = report
        atomic_write_json(paths["metadata"], _relocate(metadata, work, final))
        visual = paths["visualization"] / "visualization_metadata.json"
        atomic_write_json(visual, _relocate(json.loads(visual.read_text(encoding="utf-8")), work, final))
        if final.exists():
            raise FileExistsError(f"Refusing to overwrite {final}")
        work.rename(final)
    return final


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path, help="Trusted local result directories")
    parser.add_argument("--config", type=Path, default=SKILL_ROOT / "config.yaml")
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    config_path = args.config.resolve()
    configured = configured_path(load_config(config_path), "output_root")
    if args.output_root is None and not configured:
        parser.error("Set paths.output_root or pass --output-root")
    root = args.output_root or Path(configured)
    for run in args.runs:
        print(repackage(run, root, config_path), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
