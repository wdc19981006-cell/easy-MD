#!/usr/bin/env python
"""Rank every Vina pose and extract the global best PDBQT model."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path

ROW = re.compile(r"\s+(\d+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")


def parse_log(path: Path, run: int, seed: int | None = None) -> list[dict[str, object]]:
    poses: list[dict[str, object]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = ROW.fullmatch(line.rstrip())
        if not match:
            continue
        poses.append({
            "run": run,
            "seed": seed,
            "mode": int(match.group(1)),
            "affinity": float(match.group(2)),
            "rmsd_lb": float(match.group(3)),
            "rmsd_ub": float(match.group(4)),
        })
    return poses


def extract_model(source: Path, mode: int, target: Path) -> None:
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(True)
    selected: list[str] = []
    active = False
    saw_model = False
    for line in lines:
        if line.startswith("MODEL"):
            saw_model = True
            active = int(line.split()[1]) == mode
            continue
        if line.startswith("ENDMDL"):
            if active:
                break
            active = False
            continue
        if active or (not saw_model and mode == 1):
            selected.append(line)
    if not selected:
        raise ValueError(f"MODEL {mode} not found in {source}")
    target.write_text("".join(selected), encoding="utf-8")


def rank_poses(
    work_dir: Path, seeds: list[int], obabel: str | None = None, *,
    runs_dir: Path | None = None, pose_dir: Path | None = None,
) -> dict[str, object]:
    runs_dir = runs_dir or work_dir
    pose_dir = pose_dir or work_dir
    poses: list[dict[str, object]] = []
    for run, seed in enumerate(seeds, 1):
        log = runs_dir / f"run{run}" / "log.txt"
        poses.extend(parse_log(log, run, seed))
    if not poses:
        raise ValueError("No Vina score rows were found")
    poses.sort(key=lambda row: (row["affinity"], row["run"], row["mode"]))
    best = poses[0]
    best_pdbqt = pose_dir / "best_pose.pdbqt"
    extract_model(
        runs_dir / f"run{best['run']}" / "docked.pdbqt",
        int(best["mode"]),
        best_pdbqt,
    )
    best_pdb = pose_dir / "best_pose.pdb"
    if obabel:
        result = subprocess.run(
            [obabel, str(best_pdbqt), "-O", str(best_pdb)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"OpenBabel best-pose conversion failed: {result.stderr}")
    with (work_dir / "scores.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(poses[0]))
        writer.writeheader()
        writer.writerows(poses)
    result = {"best": best, "pose_count": len(poses), "poses": poses}
    (work_dir / "pose_ranking.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--seeds", default="42,2024,8888")
    parser.add_argument("--obabel")
    args = parser.parse_args()
    result = rank_poses(
        Path(args.work_dir).resolve(),
        [int(value) for value in args.seeds.split(",") if value.strip()],
        args.obabel,
    )
    print(json.dumps(result["best"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
