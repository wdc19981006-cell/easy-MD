#!/usr/bin/env python
"""Run PLIP once, emit XML, master JSON, reports, and VisualizerData."""

from __future__ import annotations

import argparse
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from analysis.interaction_parser import parse_interactions  # noqa: E402
from analysis.interaction_report import write_reports  # noqa: E402

FORBIDDEN = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "HOH", "WAT", "DOD", "UNK",
}


def ligand_resname(name: str) -> str:
    clean = re.sub(r"[^A-Za-z]", "", name).upper()
    for start in range(max(0, len(clean) - 2)):
        value = clean[start:start + 3]
        if len(value) == 3 and value not in FORBIDDEN:
            return value
    return (clean[:2] + "1") if len(clean) >= 2 else "LIG"


def convert_pose(source: Path, target: Path, obabel: str) -> None:
    result = subprocess.run(
        [obabel, str(source), "-O", str(target), "-f", "1", "-l", "1"],
        capture_output=True, text=True,
    )
    if result.returncode != 0 or not target.is_file():
        raise RuntimeError(f"OpenBabel pose conversion failed: {result.stderr}")


def build_complex(
    protein: Path,
    ligand: Path,
    output: Path,
    resname: str,
    chain: str = "L",
    residue_number: int = 1,
) -> None:
    protein_lines: list[str] = []
    max_serial = 0
    for line in protein.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            protein_lines.append(line)
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except ValueError:
                pass
        elif line.startswith("TER"):
            protein_lines.append(line)
    ligand_lines: list[str] = []
    serial = max_serial
    for line in ligand.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        serial += 1
        padded = line.ljust(80)
        atom_name = padded[12:16]
        altloc = padded[16:17]
        suffix = padded[26:]
        rebuilt = (
            f"HETATM{serial:5d} {atom_name}{altloc}{resname:>3s} {chain:1s}"
            f"{residue_number:4d}{suffix}"
        )
        ligand_lines.append(rebuilt[:80])
    if not ligand_lines:
        raise ValueError("Converted docking pose contains no ligand atoms")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "\n".join(protein_lines + ["TER"] + ligand_lines + ["TER", "END"]) + "\n",
        encoding="utf-8",
    )


def select_site(mol: Any, resname: str, chain: str, residue_number: int) -> tuple[Any, Any]:
    candidates: list[tuple[Any, Any]] = []
    for key, interaction_set in mol.interaction_sets.items():
        ligand = interaction_set.ligand
        if (
            str(ligand.hetid).upper() == resname.upper()
            and str(ligand.chain).strip() == chain
            and int(ligand.position) == int(residue_number)
        ):
            candidates.append((key, interaction_set))
    if len(candidates) != 1:
        found = [str(key) for key in mol.interaction_sets]
        raise ValueError(
            f"Expected exactly one PLIP site for {resname}:{chain}:{residue_number}; "
            f"found {len(candidates)} (all sites: {found})"
        )
    return candidates[0]


def analyze(
    protein: Path,
    best_pose_pdbqt: Path,
    ligand_sdf: Path,
    protein_name: str,
    ligand_name: str,
    affinity: float | None,
    structure_dir: Path,
    analysis_dir: Path,
    obabel: str,
) -> dict[str, Any]:
    from plip.basic.remote import VisualizerData
    from plip.exchange.report import StructureReport
    from plip.structure.preparation import PDBComplex

    structure_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)
    pose_pdb = structure_dir / "best_pose.pdb"
    convert_pose(best_pose_pdbqt, pose_pdb, obabel)
    resname = ligand_resname(ligand_name)
    complex_pdb = structure_dir / "complex_for_visualization.pdb"
    build_complex(protein, pose_pdb, complex_pdb, resname)

    xml_path = analysis_dir / "plip_report.xml"
    with tempfile.TemporaryDirectory(prefix=".plip-", dir=analysis_dir.resolve()) as scratch:
        mol = PDBComplex()
        mol.output_path = scratch + os.sep
        mol.load_pdb(str(complex_pdb.resolve()))
        mol.analyze()
        site, interaction_set = select_site(mol, resname, "L", 1)
        report = StructureReport(mol, outputprefix="plip_report")
        report.write_xml()
        shutil.copyfile(Path(scratch) / "plip_report.xml", xml_path)
        master = parse_interactions(interaction_set, complex_pdb)
        visualizer_data = VisualizerData(mol, site)
    master.update({
        "protein": protein_name,
        "ligand": ligand_name,
        "ligand_site": {"resname": resname, "chain": "L", "residue_number": 1},
    })
    master_path = analysis_dir / "interaction_master.json"
    master_path.write_text(
        json.dumps(master, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    summary_txt = analysis_dir / "interaction_summary.txt"
    summary_json = analysis_dir / "interaction_summary.json"
    write_reports(master, protein_name, ligand_name, affinity, summary_txt, summary_json)

    # PLIP derives pdbid from the full local path for file-based complexes.
    # PyMOL sanitizes that name while PyMOLVisualizer keeps the unsanitized
    # value, so later official methods cannot find the protein object.  Store a
    # short, PyMOL-safe identifier before serialization.
    visualizer_data.pdbid = re.sub(r"[^A-Za-z0-9_]", "_", protein_name)[:40] or "easy_md"
    visualizer_path = analysis_dir / "visualizer_data.pkl"
    with visualizer_path.open("wb") as handle:
        pickle.dump(visualizer_data, handle, protocol=4)

    return {
        "complex": str(complex_pdb),
        "master": str(master_path),
        "summary_txt": str(summary_txt),
        "summary_json": str(summary_json),
        "xml": str(xml_path),
        "visualizer_data": str(visualizer_path),
        "interaction_count": master["interaction_count"],
        "ligand_resname": resname,
        "ligand_sdf": str(ligand_sdf),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--best-pose-pdbqt", required=True)
    parser.add_argument("--ligand-sdf", required=True)
    parser.add_argument("--protein-name", required=True)
    parser.add_argument("--ligand-name", required=True)
    parser.add_argument("--affinity", type=float)
    parser.add_argument("--structure-dir", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--obabel", required=True)
    args = parser.parse_args()
    result = analyze(
        Path(args.protein).resolve(), Path(args.best_pose_pdbqt).resolve(),
        Path(args.ligand_sdf).resolve(), args.protein_name, args.ligand_name,
        args.affinity, Path(args.structure_dir).resolve(),
        Path(args.analysis_dir).resolve(), args.obabel,
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
