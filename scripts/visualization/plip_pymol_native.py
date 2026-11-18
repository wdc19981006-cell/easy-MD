#!/usr/bin/env python
"""Render PLIP VisualizerData with PLIP's official PyMOLVisualizer."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


RUNNER = r'''
import json
import pickle
import sys
from pathlib import Path

import numpy
import pymol

pymol.finish_launching(["pymol", "-cq"])
sys.path.append(r"{plip_site}")

from pymol import cmd
sys.path.insert(0, {scene_dir!r})
from scene import create_scene, interaction_segments, match_segments, measurement_segments
from plip.basic import config
from plip.visualization.pymol import PyMOLVisualizer

with open(r"{visualizer_data}", "rb") as handle:
    data = pickle.load(handle)
master = json.loads(Path(r"{master}").read_text(encoding="utf-8"))

config.MODEL = 1
vis = PyMOLVisualizer(data)
vis.set_initial_representations()
cmd.read_pdbstr(data.source_pdb_file_content, f"{{data.pdbid}}_orig")
cmd.frame(1)
current = cmd.get_object_list(selection="(all)")[0]
cmd.set_name(current, data.pdbid)
cmd.hide("everything", "all")
cmd.select(vis.ligname, "resn %s and chain %s and resi %s*" % (data.hetid, data.chain, data.position))
cmd.show("sticks", vis.ligname)
cmd.color("myblue")
cmd.color("myorange", vis.ligname)
cmd.util.cnc("all")
cmd.deselect()

vis.make_initial_selections()
vis.show_hydrophobic()
vis.show_hbonds()
vis.show_halogen()
vis.show_stacking()
vis.show_cationpi()
vis.show_sbridges()
vis.show_wbridges()
vis.show_metal()
segments = [segment for name in cmd.get_names("objects")
            for segment in measurement_segments(cmd, name)]
if not match_segments(interaction_segments(master), segments):
    raise ValueError("Official PLIP PyMOL distances do not match the master geometry")
create_scene(cmd, master, {complex_pdb!r}, segments)
cmd.save(r"{pse}")
Path(r"{result}").write_text(json.dumps({{
    "object_names": cmd.get_names("objects"),
    "segment_count": len(segments),
}}, indent=2), encoding="utf-8")
cmd.quit()
'''


def render_native(
    visualizer_data: Path,
    master: Path,
    pse: Path,
    pymol_python: str,
    plip_site: Path,
    timeout: int = 120,
    *,
    complex_pdb: Path,
) -> dict[str, object]:
    pse.parent.mkdir(parents=True, exist_ok=True)
    runner = pse.parent / "_plip_native_runner.py"
    result_json = pse.parent / "_plip_native_result.json"
    runner.write_text(RUNNER.format(
        plip_site=str(plip_site), visualizer_data=str(visualizer_data),
        master=str(master), pse=str(pse), result=str(result_json),
        scene_dir=str(Path(__file__).resolve().parent), complex_pdb=str(complex_pdb),
    ), encoding="utf-8")
    try:
        process = subprocess.run(
            [pymol_python, "-B", str(runner)], capture_output=True, text=True,
            timeout=timeout, cwd=str(pse.parent),
        )
        details = json.loads(result_json.read_text(encoding="utf-8")) if result_json.exists() else {}
        payload: dict[str, object] = {
            "source": "PLIP VisualizerData + official PyMOLVisualizer",
            "returncode": process.returncode,
            "stdout": process.stdout[-3000:],
            "stderr": process.stderr[-3000:],
            "pse": str(pse),
            **details,
        }
        payload["success"] = bool(details) and process.returncode == 0 and pse.is_file() and pse.stat().st_size > 0
        return payload
    finally:
        runner.unlink(missing_ok=True)
        result_json.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--visualizer-data", required=True)
    parser.add_argument("--master", required=True)
    parser.add_argument("--complex", required=True)
    parser.add_argument("--pse", required=True)
    parser.add_argument("--pymol-python", required=True)
    parser.add_argument("--plip-site", required=True)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    result = render_native(
        Path(args.visualizer_data).resolve(), Path(args.master).resolve(),
        Path(args.pse).resolve(), args.pymol_python,
        Path(args.plip_site).resolve(), args.timeout, complex_pdb=Path(args.complex).resolve(),
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
