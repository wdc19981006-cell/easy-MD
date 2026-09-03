#!/usr/bin/env python
"""Validate PyMOL artifacts and one-to-one PLIP interaction markers."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

RUNNER = r'''
import json
from pathlib import Path
import pymol
pymol.finish_launching(["pymol", "-cq"])
from pymol import cmd
cmd.load(r"{pse}")
objects = cmd.get_names("objects")
markers = [name for name in objects if name.startswith("_easy_md_interaction_")]
Path(r"{result}").write_text(json.dumps({{"objects": objects, "marker_count": len(markers)}}, indent=2), encoding="utf-8")
cmd.quit()
'''


def validate(
    master_path: Path,
    pse_path: Path,
    pml_path: Path,
    pymol_python: str,
    min_pse_bytes: int = 1024,
    timeout: int = 120,
) -> dict[str, object]:
    master = json.loads(master_path.read_text(encoding="utf-8"))
    expected = int(master.get("interaction_count", 0))
    result_file = pse_path.parent / "_visualization_validation.json"
    runner = pse_path.parent / "_visualization_validation.py"
    runner.write_text(
        RUNNER.format(pse=str(pse_path), result=str(result_file)), encoding="utf-8"
    )
    details: dict[str, object] = {}
    error = None
    try:
        process = subprocess.run(
            [pymol_python, str(runner)], capture_output=True, text=True,
            timeout=timeout, cwd=str(pse_path.parent),
        )
        if result_file.exists():
            details = json.loads(result_file.read_text(encoding="utf-8"))
        if process.returncode != 0:
            error = process.stderr[-2000:] or f"PyMOL validation exited {process.returncode}"
    except Exception as exc:
        error = str(exc)
    finally:
        runner.unlink(missing_ok=True)
        result_file.unlink(missing_ok=True)
    pse_size = pse_path.stat().st_size if pse_path.is_file() else 0
    checks = {
        "pml_exists_nonempty": pml_path.is_file() and pml_path.stat().st_size > 0,
        "pse_size_valid": pse_size >= min_pse_bytes,
        "interaction_marker_count_matches": details.get("marker_count") == expected,
    }
    return {
        "success": all(checks.values()),
        "expected_interaction_count": expected,
        "rendered_marker_count": details.get("marker_count"),
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
