#!/usr/bin/env python
"""Validate the saved PyMOL objects and their actual PLIP distance geometry."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

RUNNER = r'''
import json
import math
import sys
from pathlib import Path
from pymol import cmd
scene_dir, master_path, pse_path, complex_path = sys.argv[1:]
sys.path.insert(0, scene_dir)
from scene import (DASH_SETTINGS, SCENE_OBJECTS, contact_residues,
                   interaction_segments, match_segments, measurement_segments, scene_atom_indices)
master = json.loads(Path(master_path).read_text(encoding="utf-8"))
cmd.load(pse_path)
objects = cmd.get_names("objects")
expected = interaction_segments(master)
segments = measurement_segments(cmd, "Interactions") if "Interactions" in objects else []
settings = {key: cmd.get_setting_float(key, "Interactions") for key in DASH_SETTINGS} if segments else {
    key: cmd.get_setting_float(key) for key in DASH_SETTINGS
}
checks = {
    "four_objects_only": set(objects) == set(SCENE_OBJECTS) and not cmd.get_names("selections"),
    "interactions_match_plip_geometry": match_segments(expected, segments),
    "dash_settings_match": all(math.isclose(settings[key], value, abs_tol=1e-6)
                               for key, value in DASH_SETTINGS.items()),
}
residues = []
if checks["four_objects_only"]:
    residues = sorted({(atom.chain, atom.resi, atom.resn)
                       for atom in cmd.get_model("Contact_Residues").atom})
    checks["structure_objects_valid"] = (
        cmd.get_type("Protein") == cmd.get_type("Ligand") == "object:molecule"
        and cmd.count_atoms("Protein") > 0 and cmd.count_atoms("Ligand") > 0
        and all(any(chain == c and resi == r and (not resn or resn == n)
                    for c, r, n in residues) for chain, resi, resn in contact_residues(master))
    )
    if complex_path:
        cmd.load(complex_path, "_reference")
        ligand, contacts = scene_atom_indices(cmd, master, "_reference")
        checks["structure_objects_valid"] &= (
            cmd.count_atoms("Ligand") == len(ligand)
            and cmd.count_atoms("Contact_Residues") == len(contacts)
            and cmd.count_atoms("Protein") + len(ligand) == cmd.count_atoms("_reference")
        )
else:
    checks["structure_objects_valid"] = False
print(json.dumps({"objects": objects, "segment_count": len(segments), "checks": checks,
                  "dash_settings": settings, "contact_residues": residues}))
'''


def validate(
    master_path: Path,
    pse_path: Path,
    pml_path: Path,
    pymol_python: str,
    min_pse_bytes: int = 1024,
    timeout: int = 120,
    *,
    complex_pdb: Path | None = None,
) -> dict[str, object]:
    master = json.loads(master_path.read_text(encoding="utf-8"))
    expected = int(master.get("interaction_count", 0))
    details: dict[str, object] = {}
    error = None
    try:
        process = subprocess.run(
            [pymol_python, "-B", "-c", RUNNER, str(Path(__file__).resolve().parent),
             str(master_path.resolve()), str(pse_path.resolve()),
             str(complex_pdb.resolve()) if complex_pdb else ""], capture_output=True, text=True,
            timeout=timeout, cwd=str(pse_path.parent),
        )
        if process.returncode != 0:
            error = process.stderr[-2000:] or f"PyMOL validation exited {process.returncode}"
        else:
            details = json.loads(process.stdout.splitlines()[-1])
    except Exception as exc:
        error = str(exc)
    pse_size = pse_path.stat().st_size if pse_path.is_file() else 0
    checks = {
        "pml_exists_nonempty": pml_path.is_file() and pml_path.stat().st_size > 0,
        "pse_size_valid": pse_size >= min_pse_bytes,
        "master_count_valid": expected == len(master["interactions"]),
        "session_readable": bool(details) and error is None,
        **details.get("checks", {}),
    }
    return {
        "success": all(checks.values()),
        "expected_interaction_count": expected,
        "rendered_segment_count": details.get("segment_count"),
        "dash_settings": details.get("dash_settings", {}),
        "contact_residues": details.get("contact_residues", []),
        "pse_bytes": pse_size,
        "checks": checks,
        "pymol_objects": details.get("objects", []),
        "error": error,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master", required=True)
    parser.add_argument("--pse", required=True)
    parser.add_argument("--pml", required=True)
    parser.add_argument("--pymol-python", required=True)
    parser.add_argument("--min-pse-bytes", type=int, default=1024)
    parser.add_argument("--metadata")
    args = parser.parse_args()
    result = validate(
        Path(args.master).resolve(), Path(args.pse).resolve(),
        Path(args.pml).resolve(), args.pymol_python, args.min_pse_bytes,
    )
    if args.metadata:
        Path(args.metadata).write_text(
            json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
