#!/usr/bin/env python
"""Generate a strict PyMOL fallback from PLIP coordinates only."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

COLORS = {
    "hydrogen_bond": "blue",
    "hydrophobic": "gray60",
    "salt_bridge": "yellow",
    "pi_stacking": "green",
    "pi_cation": "green",
    "halogen_bond": "cyan",
    "water_bridge": "marine",
    "metal_coordination": "hotpink",
}


def _endpoint(row: dict[str, Any], side: str) -> list[float] | None:
    direct = row.get(f"{side}_coordinate")
    if direct:
        return direct
    for key in ("ring_center", "charge_center"):
        value = row.get(key, {}).get(side)
        if value:
            return value
    return None


def generate_pml(
    master_path: Path,
    complex_pdb: Path,
    output_pml: Path,
    output_pse: Path,
) -> dict[str, Any]:
    master = json.loads(master_path.read_text(encoding="utf-8"))
    records = master.get("interactions", [])
    ligand = master.get("ligand_site", {"resname": "LIG", "chain": "L", "residue_number": 1})
    lines = [
        "reinitialize",
        f"load \"{complex_pdb.as_posix()}\", Complex",
        "select Protein, Complex and polymer.protein",
        f"select Ligand, Complex and resn {ligand['resname']} and chain {ligand['chain']} and resi {ligand['residue_number']}",
        "hide everything, all",
        "show cartoon, Protein",
        "color gray80, Protein",
        "show sticks, Ligand",
        "color cyan, Ligand",
    ]
    residue_selectors: list[str] = []
    interaction_objects: list[str] = []
    marker_names: list[str] = []
    rendered = 0
    for row in records:
        chain = row.get("protein_chain", "")
        resi = row.get("protein_residue_number", "")
        if resi:
            residue_selectors.append(f"(chain {chain} and resi {resi})" if chain else f"resi {resi}")
        protein = _endpoint(row, "protein")
        ligand_coord = _endpoint(row, "ligand")
        if not protein or not ligand_coord:
            continue
        iid = row["interaction_id"]
        pobj, lobj = f"_easy_md_p_{iid}", f"_easy_md_l_{iid}"
        marker = f"_easy_md_interaction_{iid}"
        distance = f"Interaction_{iid}"
        midpoint = [(protein[i] + ligand_coord[i]) / 2.0 for i in range(3)]
        lines.extend([
            f"pseudoatom {pobj}, pos=[{protein[0]},{protein[1]},{protein[2]}]",
            f"pseudoatom {lobj}, pos=[{ligand_coord[0]},{ligand_coord[1]},{ligand_coord[2]}]",
            f"distance {distance}, {pobj}, {lobj}",
            f"set dash_color, {COLORS.get(row['interaction_type'], 'white')}, {distance}",
            f"pseudoatom {marker}, pos=[{midpoint[0]},{midpoint[1]},{midpoint[2]}]",
            f"hide everything, {pobj} or {lobj} or {marker}",
        ])
        if row["interaction_type"] != "hydrogen_bond":
            lines.append(f"hide labels, {distance}")
        interaction_objects.append(distance)
        marker_names.append(marker)
        rendered += 1

    if residue_selectors:
        lines.extend([
            f"select BindingSite, Protein and ({' or '.join(sorted(set(residue_selectors)))})",
            "show sticks, BindingSite",
            "color yellow, BindingSite",
        ])
    else:
        lines.append("select BindingSite, none")
    if interaction_objects:
        lines.append(f"group Interactions, {' '.join(interaction_objects)}")
    if marker_names:
        lines.extend([
            f"group _easy_md_validation, {' '.join(marker_names)}",
            "disable _easy_md_validation",
        ])
    lines.extend([
        "group Structures, Protein Ligand BindingSite",
        "bg_color white",
        "set dash_gap, 0.35",
        "set dash_length, 0.25",
        "set dash_radius, 0.08",
        "orient Ligand",
        "zoom Ligand or BindingSite, 5",
        "python",
        f"cmd.save({output_pse.as_posix()!r})",
        "python end",
        "quit",
    ])
    output_pml.parent.mkdir(parents=True, exist_ok=True)
    output_pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "source": "custom_pml_from_plip_master",
        "interaction_count": len(records),
        "rendered_interaction_count": rendered,
        "pml": str(output_pml),
        "pse": str(output_pse),
    }


def run_pymol(pymol: str, pml: Path, timeout: int = 120) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [pymol, "-cq", str(pml)], capture_output=True, text=True, timeout=timeout,
        cwd=str(pml.parent),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master", required=True)
    parser.add_argument("--complex", required=True)
    parser.add_argument("--pml", required=True)
    parser.add_argument("--pse", required=True)
    parser.add_argument("--pymol")
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    metadata = generate_pml(
        Path(args.master).resolve(), Path(args.complex).resolve(),
        Path(args.pml).resolve(), Path(args.pse).resolve(),
    )
    if args.run:
        if not args.pymol:
            parser.error("--run requires --pymol")
        result = run_pymol(args.pymol, Path(args.pml).resolve())
        metadata.update({"returncode": result.returncode, "stdout": result.stdout[-2000:], "stderr": result.stderr[-2000:]})
    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    return 0 if metadata["rendered_interaction_count"] == metadata["interaction_count"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
