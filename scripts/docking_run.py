#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Vina docking run: CPU detection, parallel 3-seed execution, score parsing.

Usage:
    python docking_run.py --config config.txt --seeds 42,2024,8888
    python docking_run.py --vina PATH --cpu 5 --json

Output: run1/docked.pdbqt, run2/docked.pdbqt, run3/docked.pdbqt,
        run{1-3}/log.txt, docking_stats.json
"""

import os, sys, re, math, json, subprocess, argparse, time
from pathlib import Path

from common import configured_path, enter_work_dir, get_config, load_config


def detect_cpu(leave_free=1, num_runs=3):
    """Detect CPU cores and allocate threads."""
    ncpu = os.cpu_count()
    threads_per_run = max(1, (ncpu - leave_free) // num_runs)
    print(f"Detected {ncpu} logical CPU cores")
    print(f"Allocation: {threads_per_run} threads per Vina run ({num_runs} parallel)")
    print(f"Total Vina usage: {threads_per_run * num_runs} / {ncpu} cores")
    return ncpu, threads_per_run


def setup_run_dirs(config_file="config.txt", num_runs=3):
    """Create run directories and copy input files. Returns list of dirs."""
    dirs = []
    for i in range(1, num_runs + 1):
        d = f'run{i}'
        os.makedirs(d, exist_ok=True)
        for f in ['receptor.pdbqt', 'ligand.pdbqt']:
            if os.path.exists(f):
                import shutil
                shutil.copy2(f, os.path.join(d, f))
        if os.path.exists(config_file):
            import shutil
            shutil.copy2(config_file, os.path.join(d, "config.txt"))
        # Verify
        cfg = os.path.join(d, 'config.txt')
        if not os.path.exists(cfg):
            print(f"ERROR: Failed to create {cfg}")
            sys.exit(1)
        dirs.append(d)
    print("All run directories ready")
    return dirs


def parse_vina_log(log_path):
    """Parse Vina log file, return [(mode, affinity, rmsd_lb, rmsd_ub), ...]."""
    if not os.path.exists(log_path):
        return []
    with open(log_path, 'r') as f:
        text = f.read()
    table = []
    for line in text.split('\n'):
        m = re.match(r'\s+(\d+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)', line)
        if m:
            table.append((int(m.group(1)), float(m.group(2)),
                          float(m.group(3)), float(m.group(4))))
    return table


def compute_statistics(all_scores):
    """Compute mean, SD, and stability from per-run best scores.

    Args:
        all_scores: list of (run_num, seed, best_affinity)
    Returns:
        dict with best_affinity, mean_affinity, sd_affinity, stability, best_run, best_seed
    """
    if not all_scores:
        return {'best_affinity': 0, 'mean_affinity': 0, 'sd_affinity': 0,
                'stability': 'N/A', 'best_run': None, 'best_seed': None}

    best_affinity = min(s[2] for s in all_scores)
    scores_only = [s[2] for s in all_scores]
    mean_aff = sum(scores_only) / len(scores_only)

    if len(scores_only) >= 2:
        variance = sum((s - mean_aff) ** 2 for s in scores_only) / (len(scores_only) - 1)
        sd_aff = math.sqrt(variance)
    else:
        sd_aff = 0.0

    if sd_aff < 0.3:
        stability = "Excellent (SD < 0.3 kcal/mol)"
    elif sd_aff < 0.5:
        stability = "Good (SD < 0.5 kcal/mol)"
    elif sd_aff < 1.0:
        stability = "Moderate (SD < 1.0 kcal/mol)"
    else:
        stability = "High variability (SD >= 1.0 kcal/mol) - consider more runs"

    best_entry = min(all_scores, key=lambda x: x[2])
    return {
        'best_affinity': round(best_affinity, 3),
        'mean_affinity': round(mean_aff, 3),
        'sd_affinity': round(sd_aff, 3),
        'stability': stability,
        'best_run': best_entry[0],
        'best_seed': best_entry[1],
        'run_scores': scores_only,
    }


def report_results(all_scores, stats, all_tables):
    """Print formatted docking results report."""
    print("=" * 60)
    print("DOCKING RESULTS - Independent Vina Runs")
    print("=" * 60)
    for rn, sd, aff in all_scores:
        mark = " <-- BEST" if rn == stats['best_run'] else ""
        print(f"  Run {rn} (seed={sd}): Best = {aff:.3f} kcal/mol{mark}")
    print("")
    print(f"  Best affinity:    {stats['best_affinity']:.3f} kcal/mol")
    print(f"  Mean affinity:    {stats['mean_affinity']:.3f} kcal/mol")
    print(f"  Std deviation:    {stats['sd_affinity']:.3f} kcal/mol")
    print(f"  Stability:        {stats['stability']}")
    print("=" * 60)

    # Show full table for best run
    if stats['best_run'] and stats['best_run'] in all_tables:
        print(f"\n  Full {len(all_tables[stats['best_run']])}-mode scores (Best Run #{stats['best_run']}):")
        print(f"  {'Mode':<8} {'Affinity':<16} {'RMSD l.b.':<12} {'RMSD u.b.'}")
        print(f"  {'-'*8} {'-'*16} {'-'*12} {'-'*12}")
        for mode, aff, lb, ub in all_tables[stats['best_run']][:12]:
            star = " <--" if mode == 1 else ""
            print(f"  {mode:<8} {aff:<16.3f} {lb:<12.3f} {ub:.3f}{star}")


def main():
    parser = argparse.ArgumentParser(
        description="Run AutoDock Vina with multiple random seeds")
    parser.add_argument("--vina",
                        help="Path to Vina executable")
    parser.add_argument("--config", default="config.txt", help="Vina config file")
    parser.add_argument("--seeds",
                        help="Comma-separated random seeds")
    parser.add_argument("--cpu", type=int, default=None,
                        help="Threads per run (auto-detect if not set)")
    parser.add_argument("--num-runs", type=int)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--skill-config", help="Path to easy-md config.yaml")
    parser.add_argument("--structured", action="store_true",
                        help="Use docking/config.txt, shared inputs, and docking/runs/runN outputs")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--json", action="store_true",
                        help="Output stats as JSON")
    args = parser.parse_args()

    config = load_config(args.skill_config)
    enter_work_dir(args.work_dir, config)
    args.vina = args.vina or configured_path(config, "vina")
    if not args.vina:
        parser.error("Vina path is required via --vina or config.yaml")

    configured_seeds = get_config(config, "docking.seeds", [42, 2024, 8888])
    seed_text = args.seeds or ",".join(str(seed) for seed in configured_seeds)
    seeds = [int(s.strip()) for s in seed_text.split(',') if s.strip()]
    requested_runs = args.num_runs or int(
        get_config(config, "docking.num_runs", len(seeds))
    )
    if requested_runs < 1 or len(seeds) < requested_runs:
        parser.error("--num-runs must be positive and cannot exceed the seed count")
    seeds = seeds[:requested_runs]
    num_runs = requested_runs

    leave_free = int(get_config(config, "docking.cpu_leave_free", 1))
    ncpu, threads = detect_cpu(leave_free=leave_free, num_runs=num_runs)
    if args.cpu is not None:
        threads = args.cpu

    if args.structured:
        config_file = Path(args.config).resolve()
        if not config_file.is_file():
            parser.error(f"Vina config not found: {config_file}")
        docking = config_file.parent
        run_dirs = [docking / "runs" / f"run{i}" for i in range(1, num_runs + 1)]
        if not args.dry_run:
            for folder in run_dirs:
                folder.mkdir(parents=True, exist_ok=True)
    else:
        run_dirs = [Path(folder).resolve() for folder in setup_run_dirs(args.config, num_runs)]
        docking = Path.cwd()

    # Launch all runs in background (via subprocess)
    procs = []
    for i, (seed, folder) in enumerate(zip(seeds, run_dirs), 1):
        cwd = docking if args.structured else folder
        cfg = config_file.name if args.structured else "config.txt"
        cmd = [
            args.vina, '--config', cfg, '--out', str((folder / "docked.pdbqt").relative_to(cwd)),
            '--seed', str(seed), '--cpu', str(threads),
        ]
        if args.dry_run:
            print(f"cwd={cwd} command={json.dumps(cmd)}")
            continue
        print(f"Launching Run {i} (seed={seed})...")
        with (folder / "log.txt").open("w") as log_f:
            proc = subprocess.Popen(
                cmd, cwd=cwd,
                stdout=log_f, stderr=subprocess.STDOUT)
            procs.append((i, seed, proc))
    if args.dry_run:
        return

    # Wait for all to complete
    failed_runs = []
    for i, seed, proc in procs:
        ret = proc.wait()
        status = "OK" if ret == 0 else f"EXIT {ret}"
        print(f"  Run {i} (seed={seed}): {status}")
        if ret != 0:
            failed_runs.append((i, seed, ret))

    if failed_runs:
        print(f"ERROR: {len(failed_runs)} Vina run(s) failed: {failed_runs}")
        sys.exit(1)

    print("=== All Vina runs complete ===")

    # Parse results
    all_scores = []
    all_tables = {}
    for i, (seed, folder) in enumerate(zip(seeds, run_dirs), 1):
        table = parse_vina_log(folder / "log.txt")
        if table:
            all_tables[i] = table
            best = table[0][1]
            all_scores.append((i, seed, best))
        else:
            print(f"ERROR: No parsable score table for run {i}")

    if len(all_scores) != num_runs:
        print(
            f"ERROR: Parsed {len(all_scores)} of {num_runs} required Vina runs"
        )
        sys.exit(1)

    stats = compute_statistics(all_scores)
    stats['threads_per_run'] = threads
    stats['total_cpu'] = ncpu
    stats['seeds'] = seeds

    report_results(all_scores, stats, all_tables)

    # Save stats
    with (docking / "docking_stats.json").open("w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2)
    print("\nSaved docking_stats.json")

    if args.json:
        # Also output all tables
        tables_serializable = {
            str(k): [(m, round(a, 3), round(lb, 3), round(ub, 3))
                      for m, a, lb, ub in v]
            for k, v in all_tables.items()
        }
        print(json.dumps({'stats': stats, 'tables': tables_serializable},
                         indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
