#!/usr/bin/env python
"""Generate a strict PyMOL fallback from PLIP coordinates only."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path
from typing import Any

if __package__:
    from .scene import interaction_segments
else:
    from scene import interaction_segments


def generate_pml(
    master_path: Path,
    complex_pdb: Path,
    output_pml: Path,
    output_pse: Path,
) -> dict[str, Any]:
    master = json.loads(master_path.read_text(encoding="utf-8"))
    interaction_segments(master)
    output_pml.parent.mkdir(parents=True, exist_ok=True)
    base = output_pml.parent.resolve()
    relative = lambda path: Path(os.path.relpath(path.resolve(), base)).as_posix()
    lines = [
        "# Run from this script's directory. Open the PSE for the saved session.",
        "python",
        "import json",
        "from pathlib import Path",
        "from pymol import cmd",
        Path(__file__).with_name("scene.py").read_text(encoding="utf-8"),
        f"master = json.loads(Path({relative(master_path)!r}).read_text(encoding='utf-8'))",
        f"create_scene(cmd, master, {relative(complex_pdb)!r}, interaction_segments(master))",
        f"cmd.save({relative(output_pse)!r})",
        "python end",
    ]
    output_pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "source": "custom_pml_from_plip_master",
        "interaction_count": len(master["interactions"]),
        "pml": str(output_pml),
        "pse": str(output_pse),
    }


def run_pymol(pymol_python: str, pml: Path, timeout: int = 120) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [pymol_python, "-B", "-c", "from pymol import cmd; import sys; cmd.do('@' + sys.argv[1], echo=0)",
         pml.name], capture_output=True, text=True, timeout=timeout,
        cwd=str(pml.parent),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master", required=True)
    parser.add_argument("--complex", required=True)
    parser.add_argument("--pml", required=True)
    parser.add_argument("--pse", required=True)
    parser.add_argument("--pymol-python")
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    metadata = generate_pml(
        Path(args.master).resolve(), Path(args.complex).resolve(),
        Path(args.pml).resolve(), Path(args.pse).resolve(),
    )
    if args.run:
        if not args.pymol_python:
            parser.error("--run requires --pymol-python")
        result = run_pymol(args.pymol_python, Path(args.pml).resolve())
        metadata.update({"returncode": result.returncode, "stdout": result.stdout[-2000:], "stderr": result.stderr[-2000:]})
    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    return int(args.run and (result.returncode != 0 or not Path(args.pse).is_file()))


if __name__ == "__main__":
    raise SystemExit(main())
