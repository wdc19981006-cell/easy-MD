#!/usr/bin/env python
"""Generate human-readable and compatibility JSON reports from PLIP master data."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


DISPLAY = {
    "hydrogen_bond": "Hydrogen Bonds",
    "hydrophobic": "Hydrophobic Contacts",
    "salt_bridge": "Salt Bridges",
    "pi_stacking": "Pi Stacking",
    "pi_cation": "Pi-Cation",
    "halogen_bond": "Halogen Bonds",
    "water_bridge": "Water Bridges",
    "metal_coordination": "Metal Coordination",
}


def write_reports(
    master: dict[str, Any],
    protein: str,
    ligand: str,
    affinity: float | None,
    text_path: str | Path,
    json_path: str | Path,
) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in master.get("interactions", []):
        grouped[row["interaction_type"]].append(row)
    lines = [
        "=" * 68,
        "Protein-Ligand Interaction Analysis",
        "Source: PLIP (sole interaction classifier)",
        "=" * 68,
        f"Protein: {protein}",
        f"Ligand: {ligand}",
        f"Docking affinity: {affinity if affinity is not None else 'N/A'} kcal/mol",
        "",
    ]
    for kind in DISPLAY:
        rows = grouped.get(kind, [])
        lines.append(f"{DISPLAY[kind]} ({len(rows)})")
        lines.append("-" * 68)
        if not rows:
            lines.append("None detected")
        for row in rows:
            distance = "N/A" if row.get("distance") is None else f"{row['distance']:.2f} A"
            angle = "" if row.get("angle") is None else f", angle {row['angle']:.1f} deg"
            lines.append(f"{row['interaction_id']}: {row['residue']} - {distance}{angle}")
        lines.append("")
    lines.extend([
        f"Total PLIP interactions: {master.get('interaction_count', 0)}",
        "PyMOL is used only to render these PLIP-defined interactions.",
        "=" * 68,
    ])
    Path(text_path).write_text("\n".join(lines) + "\n", encoding="utf-8")

    compatibility = {
        "protein": protein,
        "ligand": ligand,
        "best_affinity_kcal_mol": affinity,
        "source": "PLIP",
        "interaction_count": master.get("interaction_count", 0),
        "counts": master.get("counts", {}),
        "interactions": {
            "hydrogen_bonds": {"count": len(grouped["hydrogen_bond"]), "details": grouped["hydrogen_bond"]},
            "hydrophobic_contacts": {"count": len(grouped["hydrophobic"]), "details": grouped["hydrophobic"]},
            "salt_bridges": {"count": len(grouped["salt_bridge"]), "details": grouped["salt_bridge"]},
            "pi_stacking": {"count": len(grouped["pi_stacking"]), "details": grouped["pi_stacking"]},
            "pi_cation": {"count": len(grouped["pi_cation"]), "details": grouped["pi_cation"]},
            "halogen_bonds": {"count": len(grouped["halogen_bond"]), "details": grouped["halogen_bond"]},
            "water_bridges": {"count": len(grouped["water_bridge"]), "details": grouped["water_bridge"]},
            "metal_coordination": {"count": len(grouped["metal_coordination"]), "details": grouped["metal_coordination"]},
        },
    }
    Path(json_path).write_text(
        json.dumps(compatibility, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return compatibility


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master", required=True)
    parser.add_argument("--protein", required=True)
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--affinity", type=float)
    parser.add_argument("--text", required=True)
    parser.add_argument("--json", required=True)
    args = parser.parse_args()
    master = json.loads(Path(args.master).read_text(encoding="utf-8"))
    write_reports(master, args.protein, args.ligand, args.affinity, args.text, args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
