#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Result packaging: split poses, merge multi-model PDB, reports, copy to target.

Usage:
    python package_results.py --target "D:/Desktop/test/ROCK2_Belumosudil"
        --protein "ROCK2" --pdb-id "8X92" --ligand "Belumosudil"
        --smiles "CC(C)..." --cid "11950170" --formula "C26H24N6O2" --mw 452.5

Output:
    {target}/{PROTEIN}_{LIGAND}/
    ├── original.pdb, protein_clean.pdb, receptor.pdbqt
    ├── ligand.sdf, ligand.pdbqt, config.txt
    ├── vina_logs/ (run1_log.txt, run2_log.txt, run3_log.txt)
    ├── all_poses.pdb, best_pose.pdb, best_pose.pdbqt
    ├── scores.csv, summary.txt, metadata.json, view.pml
    └── interaction_analysis/
        ├── complex_for_plip.pdb, plip_report.xml
        ├── interaction_summary.txt, interaction_summary.json
        ├── interaction_2D.png, interaction_2D.svg, interaction_2D.emf
        └── pymol_interaction.pml, pymol_interaction.pse
"""

import os, sys, re, json, argparse, shutil, time, subprocess
from pathlib import Path

from common import (
    configured_path,
    enter_work_dir,
    get_config,
    load_config,
    validate_cleanup_target,
)


def extract_best_pose_pdbqt(best_run):
    """Extract MODEL 1 from best run's docked.pdbqt → best_pose.pdbqt."""
    docked = f"run{best_run}/docked.pdbqt"
    if not os.path.exists(docked):
        print(f"WARNING: {docked} not found — cannot extract best_pose.pdbqt")
        return False
    inside_model1 = False
    lines = []
    with open(docked, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                model_num = int(line.split()[1])
                if model_num == 1:
                    inside_model1 = True
                    continue
                elif inside_model1:
                    break  # end of MODEL 1
            if inside_model1:
                lines.append(line)
    if lines:
        with open("best_pose.pdbqt", "w") as f:
            f.writelines(lines)
        print(f"Extracted best_pose.pdbqt ({len(lines)} atoms)")
        return True
    print("WARNING: Could not extract MODEL 1 from docked PDBQT")
    return False


def split_pdbqt(best_run, obabel_path="obabel"):
    """Split docked.pdbqt into individual PDBs using OpenBabel."""
    docked = f"run{best_run}/docked.pdbqt"
    if not os.path.exists(docked):
        print(f"ERROR: {docked} not found")
        return None
    result = subprocess.run(
        [obabel_path, docked, "-O", "best_pose.pdb", "-m"],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: obabel split failed: {result.stderr}")
        return None
    nposes = len([f for f in os.listdir(".")
                  if f.startswith("best_pose") and f.endswith(".pdb")])
    print(f"Split docked into {nposes} individual PDB files")
    return nposes


def parse_scores_from_log(best_run):
    """Parse Vina scores from best run log."""
    log_path = f"run{best_run}/log.txt"
    if not os.path.exists(log_path):
        return []
    with open(log_path, "r") as f:
        text = f.read()
    scores = []
    for line in text.split("\n"):
        m = re.match(r"\s+(\d+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)", line)
        if m:
            scores.append((int(m.group(1)), float(m.group(2)),
                           float(m.group(3)), float(m.group(4))))
    return scores


def write_multi_model_pdb(scores):
    """Merge individual PDBs into multi-model PDB with REMARK scores."""
    with open("all_poses.pdb", "w") as outf:
        for i in range(1, min(len(scores) + 1, 10)):
            pose_file = f"best_pose{i}.pdb"
            if not os.path.exists(pose_file):
                continue
            aff, lb, ub = scores[i-1][1], scores[i-1][2], scores[i-1][3]
            outf.write(f"MODEL     {i}\n")
            outf.write(f"REMARK Vina Score: {aff:.3f} kcal/mol  "
                       f"RMSD l.b.: {lb:.3f}  RMSD u.b.: {ub:.3f}\n")
            with open(pose_file, "r") as inf:
                for line in inf:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        outf.write(line)
            outf.write("ENDMDL\n")
    print(f"Wrote all_poses.pdb ({min(len(scores), 9)} models)")


def write_scores_csv(scores):
    """Write scores.csv."""
    with open("scores.csv", "w") as f:
        f.write("mode,affinity_kcal_mol,rmsd_lb,rmsd_ub\n")
        for mode, aff, lb, ub in scores:
            f.write(f"{mode},{aff:.3f},{lb:.3f},{ub:.3f}\n")
    print("Wrote scores.csv")


def _get_center_from_config(config_path="config.txt"):
    """Extract box center and size from Vina config."""
    center = [0.0, 0.0, 0.0]
    box_size = [0.0, 0.0, 0.0]
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            for line in f:
                for ai, ax in enumerate(["x", "y", "z"]):
                    if line.startswith(f"center_{ax}"):
                        try:
                            center[ai] = float(line.split("=")[1].strip())
                        except Exception:
                            pass
                for ai, ax in enumerate(["x", "y", "z"]):
                    if line.startswith(f"size_{ax}"):
                        try:
                            box_size[ai] = float(line.split("=")[1].strip())
                        except Exception:
                            pass
    return center, box_size


def write_metadata_json(args, dock_stats, redock_info, scores):
    """Write machine-readable metadata.json."""
    center, box_size = _get_center_from_config()
    metadata = {
        "workflow": {
            "name": "easy-MD", "version": "2.0",
            "toolchain": ["PDBFixer", "RDKit", "OpenBabel 3.1.1",
                          "AutoDock Vina 1.2.5"],
        },
        "protein": {
            "input_name": args.protein,
            "gene_symbol": args.gene_symbol or args.protein,
            "official_name": args.official_name or "",
            "pdb_id": args.pdb_id,
            "resolution": args.resolution or "?",
            "method": args.method or "X-RAY DIFFRACTION",
            "organism": args.organism or "Homo sapiens",
            "selected_chains": (args.chains.split(",")
                                if args.chains else ["A"]),
            "native_ligand": redock_info.get("native_ligand"),
        },
        "ligand": {
            "name": args.ligand, "cid": args.cid or "",
            "smiles": args.smiles or "",
            "formula": args.formula or "",
            "molecular_weight": args.mw or 0,
        },
        "docking": {
            "software": "AutoDock Vina", "version": "1.2.5",
            "mode": "semi-flexible", "receptor": "rigid",
            "ligand": "flexible",
            "exhaustiveness": args.exhaustiveness or 32,
            "num_modes": args.num_modes or 9,
            "energy_range": args.energy_range or 3,
            "num_runs": dock_stats.get("num_runs", 3),
            "seeds": dock_stats.get("seeds", [42, 2024, 8888]),
            "box_center": [round(c, 3) for c in center],
            "box_size": box_size,
            "cpu_total": dock_stats.get("total_cpu", "?"),
            "cpu_per_run": dock_stats.get("threads_per_run", "?"),
        },
        "results": {
            "best_affinity": dock_stats.get("best_affinity"),
            "mean_affinity": dock_stats.get("mean_affinity"),
            "sd_affinity": dock_stats.get("sd_affinity"),
            "stability": dock_stats.get("stability"),
            "run_scores": dock_stats.get("run_scores", []),
            "best_run": dock_stats.get("best_run"),
            "best_seed": dock_stats.get("best_seed"),
        },
        "validation": {
            "native_ligand": redock_info.get("native_ligand"),
            "redocking_rmsd": redock_info.get("redocking_rmsd"),
            "validation_status": redock_info.get(
                "validation_status", "Not available"),
        },
    }
    with open("metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print("Wrote metadata.json")


def write_view_pml(args, dock_stats, box_size):
    """Write PyMOL visualization script."""
    center, _ = _get_center_from_config()
    best_aff = dock_stats.get("best_affinity", 0)
    outdir = f"{args.protein}_{args.ligand}"
    pml = (
        f"# PyMOL session - {args.protein}_{args.ligand}\n"
        f"# easy-MD v2.0 | AutoDock Vina 1.2.5\n\n"
        f"load original.pdb, original\n"
        f"load protein_clean.pdb, receptor\n"
        f"load all_poses.pdb, poses\n"
        f"load best_pose.pdb, best\n\n"
        f"show cartoon, receptor\n"
        f"set cartoon_transparency, 0.3, receptor\n"
        f"color grey80, receptor\n"
        f"show sticks, best\n"
        f"color green, best\n"
        f"show lines, poses\n"
        f"color cyan, poses\n\n"
        f"pseudoatom dock_ctr, "
        f"pos=[{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}]\n"
        f"show spheres, dock_ctr\n"
        f"set sphere_scale, 0.5, dock_ctr\n"
        f"color red, dock_ctr\n"
        f"label dock_ctr, \"Box: {box_size[0]:.0f}x{box_size[1]:.0f}x{box_size[2]:.0f}A\"\n\n"
        f"cmd.center('best')\n"
        f"zoom center, 15.0\n"
        f"bg_color white\n"
        f"set antialias, 2\n"
        f"print \"easy-MD: {outdir} | Best: {best_aff:.2f} kcal/mol\"\n"
        f"print \"Use [ctrl]+[left]/[right] to switch poses (1-9)\"\n"
    )
    with open("view.pml", "w") as f:
        f.write(pml)
    print("Wrote view.pml")


def write_summary_txt(args, dock_stats, redock_info, scores):
    """Write human-readable summary.txt."""
    center, box_size = _get_center_from_config()
    redock_text = ""
    if os.path.exists("redocking_result.txt"):
        with open("redocking_result.txt", "r") as f:
            redock_text = f.read()

    protein_label = args.gene_symbol or args.protein
    native_lig = redock_info.get("native_ligand", "None (apo)")
    seeds = dock_stats.get("seeds", [42, 2024, 8888])
    run_scores = dock_stats.get("run_scores", [])
    best_run = dock_stats.get("best_run")

    L = []
    L.append("=" * 70)
    L.append(f"  Molecular Docking Report - "
             f"{args.protein}_{args.ligand}")
    L.append("  easy-MD v2.0 | AutoDock Vina 1.2.5 | "
             "Semi-Flexible Docking")
    L.append("=" * 70)
    L.append("")
    L.append("[PROTEIN]")
    L.append(f"  Input:      {args.protein}")
    L.append(f"  Gene:       {protein_label}")
    L.append(f"  Official:   {args.official_name or protein_label}")
    L.append(f"  PDB:        {args.pdb_id} | "
             f"{args.resolution or '?'} A | "
             f"{args.method or 'X-RAY'} | "
             f"{args.organism or 'Homo sapiens'}")
    L.append(f"  Chain(s):   {args.chains or 'A'}")
    L.append(f"  Native lig: {native_lig}")
    L.append("")
    L.append("[LIGAND]")
    L.append(f"  Name:       {args.ligand}")
    L.append(f"  CID:        {args.cid or '?'} | "
             f"Formula: {args.formula or '?'} | "
             f"MW: {args.mw or '?'} g/mol")
    L.append(f"  SMILES:     {args.smiles or '?'}")
    L.append("")
    L.append("[DOCKING]")
    L.append("  Software:   AutoDock Vina 1.2.5 | Mode: Semi-flexible")
    L.append(f"  Exhaust.:   {args.exhaustiveness or 32} | "
             f"Modes: {args.num_modes or 9} | "
             f"Energy range: {args.energy_range or 3} kcal/mol")
    L.append(f"  Runs:       {len(run_scores)} (seeds: "
             f"{', '.join(str(s) for s in seeds)})")
    L.append(f"  CPU:        "
             f"{dock_stats.get('threads_per_run', '?')} threads/run "
             f"({dock_stats.get('total_cpu', '?')} total)")
    L.append(f"  Box center: ({center[0]:.3f}, "
             f"{center[1]:.3f}, {center[2]:.3f})")
    L.append(f"  Box size:   {box_size[0]:.0f} x "
             f"{box_size[1]:.0f} x {box_size[2]:.0f} A")
    L.append("")
    L.append("[RESULTS - 3 Independent Runs]")
    L.append("  Run      Seed     Best Affinity")
    L.append("  ----------------------------------------")

    for i in range(min(3, len(run_scores))):
        rn = i + 1
        sd = seeds[i] if i < len(seeds) else 0
        aff = run_scores[i]
        mark = "  <-- BEST" if rn == best_run else ""
        L.append(f"  Run {rn}    {sd:<6}  {aff:.3f} kcal/mol{mark}")

    L.append("  ----------------------------------------")
    L.append(f"  Best:       "
             f"{dock_stats.get('best_affinity', 0):.3f} kcal/mol")
    L.append(f"  Mean:       "
             f"{dock_stats.get('mean_affinity', 0):.3f} kcal/mol")
    L.append(f"  SD:         "
             f"{dock_stats.get('sd_affinity', 0):.3f} kcal/mol")
    L.append(f"  Stability:  {dock_stats.get('stability', 'N/A')}")
    L.append("")
    L.append(f"[9-MODE SCORES - Best Run #{best_run or '?'}]")
    L.append("  Mode     Affinity         RMSD l.b.    RMSD u.b.")
    L.append("  --------------------------------------------------")

    for mode, aff, lb, ub in scores[:12]:
        star = " <--" if mode == 1 else ""
        L.append(f"  {mode:<8} {aff:<16.3f} {lb:<12.3f} {ub:.3f}{star}")

    L.append("")
    L.append("[VALIDATION]")
    if redock_text:
        for line in redock_text.strip().split("\n"):
            L.append(f"  {line}")
    else:
        L.append("  No redocking validation available.")

    L.append("")
    L.append("[OUTPUT]")
    L.append(f"  {args.protein}_{args.ligand}/")
    L.append("  |-- original.pdb       |-- best_pose.pdb")
    L.append("  |-- protein_clean.pdb  |-- best_pose.pdbqt")
    L.append("  |-- receptor.pdbqt     |-- scores.csv")
    L.append("  |-- ligand.sdf/.pdbqt  |-- summary.txt")
    L.append("  |-- config.txt         |-- metadata.json")
    L.append("  |-- all_poses.pdb      |-- view.pml")
    L.append("  |-- vina_logs/ (run1-3)")
    if os.path.isdir("interaction_analysis"):
        L.append("  |")
        L.append("  |-- interaction_analysis/")
        L.append("       |-- complex_for_plip.pdb")
        L.append("       |-- plip_report.xml")
        L.append("       |-- interaction_summary.txt")
        L.append("       |-- interaction_summary.json")
        L.append("       |-- interaction_2D.png")
        L.append("       |-- interaction_2D.svg")
        L.append("       |-- interaction_2D.emf")
        L.append("       |-- pymol_interaction.pml/.pse")
    L.append("")
    L.append("=" * 70)
    L.append("  easy-MD v2.0 | PDBFixer | RDKit | "
             "OpenBabel 3.1.1 | Vina 1.2.5")
    L.append("=" * 70)

    with open("summary.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("Wrote summary.txt")


def copy_to_target(target_dir, protein, ligand, pdb_id):
    """Copy result files to target directory."""
    outdir = os.path.join(target_dir, f"{protein}_{ligand}")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(outdir, "vina_logs"), exist_ok=True)

    files = [
        ("all_poses.pdb", "all_poses.pdb"),
        ("protein_clean.pdb", "protein_clean.pdb"),
        ("scores.csv", "scores.csv"),
        ("summary.txt", "summary.txt"),
        ("config.txt", "config.txt"),
        ("receptor.pdbqt", "receptor.pdbqt"),
        ("ligand.pdbqt", "ligand.pdbqt"),
        ("ligand.sdf", "ligand.sdf"),
        ("best_pose1.pdb", "best_pose.pdb"),
        ("best_pose.pdbqt", "best_pose.pdbqt"),
        ("metadata.json", "metadata.json"),
        ("view.pml", "view.pml"),
        (f"{pdb_id}.pdb", "original.pdb"),
        ("run1/log.txt", "vina_logs/run1_log.txt"),
        ("run2/log.txt", "vina_logs/run2_log.txt"),
        ("run3/log.txt", "vina_logs/run3_log.txt"),
    ]
    for src, dst in files:
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(outdir, dst))
        elif src == "best_pose.pdbqt":
            print("  Note: best_pose.pdbqt not found (pre-v2.1 run)")

    # Copy interaction_analysis/ directory (exists if PLIP ran)
    ia_dir = "interaction_analysis"
    if os.path.isdir(ia_dir):
        ia_dest = os.path.join(outdir, ia_dir)
        if os.path.exists(ia_dest):
            shutil.rmtree(ia_dest)
        shutil.copytree(ia_dir, ia_dest)
        print(f"  Copied interaction_analysis/ ({len(os.listdir(ia_dir))} files)")

    print(f"Copied results to: {outdir}")
    return outdir


def verify_target(outdir):
    """Verify all required files exist at target."""
    required = [
        "all_poses.pdb", "summary.txt", "scores.csv", "best_pose.pdb",
        "best_pose.pdbqt",
        "receptor.pdbqt", "ligand.pdbqt", "protein_clean.pdb",
        "metadata.json", "view.pml", "original.pdb",
    ]
    # interaction_analysis/ is optional — only verified if it exists at source
    ia_src = "interaction_analysis"
    ia_dst = os.path.join(outdir, "interaction_analysis")
    if os.path.isdir(ia_src):
        required.append("interaction_analysis/")
    all_ok = True
    report = []
    for f in required:
        if f == "interaction_analysis/":
            path = ia_dst
            exists = os.path.isdir(path)
            nf = len(os.listdir(path)) if exists else 0
            st = (f"OK ({nf} files)" if exists and nf > 0
                  else "MISSING or EMPTY!")
        else:
            path = os.path.join(outdir, f)
            exists = os.path.exists(path)
            size = os.path.getsize(path) if exists else 0
            st = (f"OK ({size} bytes)" if exists and size > 0
                  else "MISSING or EMPTY!")
        report.append(f"  {f}: {st}")
        if not exists or (f != "interaction_analysis/" and size == 0):
            all_ok = False
    return all_ok, "\n".join(report)


def cleanup_temp(tmp_dir):
    """Clean up temp directory with retry for Windows file locking."""
    try:
        tmp_dir = validate_cleanup_target(tmp_dir)
    except ValueError as exc:
        print(f"Refusing temp cleanup: {exc}")
        return False
    os.chdir("/")
    for root, dirs, files in os.walk(tmp_dir, topdown=False):
        for name in files:
            try:
                os.unlink(os.path.join(root, name))
            except Exception:
                pass
        for name in dirs:
            try:
                os.rmdir(os.path.join(root, name))
            except Exception:
                pass
    for i in range(10):
        try:
            shutil.rmtree(tmp_dir)
            print("Temp cleaned")
            return True
        except OSError:
            time.sleep(2)
    print("Note: temp not cleaned (files locked). "
          "Results are safe in target.")
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Package docking results")
    parser.add_argument("--target", required=True)
    parser.add_argument("--keep-temp", action="store_true")
    parser.add_argument("--protein", required=True)
    parser.add_argument("--gene-symbol")
    parser.add_argument("--official-name")
    parser.add_argument("--pdb-id", required=True)
    parser.add_argument("--resolution")
    parser.add_argument("--method", default="X-RAY DIFFRACTION")
    parser.add_argument("--organism", default="Homo sapiens")
    parser.add_argument("--chains", default="A")
    parser.add_argument("--ligand", required=True)
    parser.add_argument("--smiles")
    parser.add_argument("--cid")
    parser.add_argument("--formula")
    parser.add_argument("--mw", type=float)
    parser.add_argument("--exhaustiveness", type=int)
    parser.add_argument("--num-modes", type=int)
    parser.add_argument("--energy-range", type=int)
    parser.add_argument("--obabel")
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--config", help="Path to easy-md config.yaml")
    parser.add_argument("--docking-stats", default="docking_stats.json")
    parser.add_argument("--redocking-info",
                        default="redocking_info.json")
    args = parser.parse_args()

    config = load_config(args.config)
    args.target = str(Path(args.target).expanduser().resolve())
    tmp = enter_work_dir(args.work_dir, config)
    args.obabel = args.obabel or configured_path(config, "obabel") or "obabel"
    args.exhaustiveness = args.exhaustiveness or int(
        get_config(config, "docking.exhaustiveness", 32)
    )
    args.num_modes = args.num_modes or int(
        get_config(config, "docking.num_modes", 9)
    )
    args.energy_range = args.energy_range or int(
        get_config(config, "docking.energy_range", 3)
    )

    dock_stats = {}
    if os.path.exists(args.docking_stats):
        with open(args.docking_stats, "r") as f:
            dock_stats = json.load(f)

    redock_info = {}
    if os.path.exists(args.redocking_info):
        with open(args.redocking_info, "r") as f:
            redock_info = json.load(f)

    best_run = dock_stats.get("best_run", 1)

    nposes = split_pdbqt(best_run, args.obabel)
    if nposes is None:
        print("ERROR: Failed to split docked poses")
        sys.exit(1)

    # Extract single best_pose.pdbqt from best run
    extract_best_pose_pdbqt(best_run)

    scores = parse_scores_from_log(best_run)
    if not scores:
        print("WARNING: Could not parse scores from Vina log")

    write_multi_model_pdb(scores)
    write_scores_csv(scores)
    write_metadata_json(args, dock_stats, redock_info, scores)
    _, box_size = _get_center_from_config()
    write_view_pml(args, dock_stats, box_size)
    write_summary_txt(args, dock_stats, redock_info, scores)

    outdir = copy_to_target(args.target, args.protein,
                            args.ligand, args.pdb_id)

    all_ok, report = verify_target(outdir)
    print(f"\nVerification:\n{report}")
    if not all_ok:
        print("WARNING: Some files missing/empty. Retrying copy...")
        copy_to_target(args.target, args.protein,
                       args.ligand, args.pdb_id)
        all_ok, _ = verify_target(outdir)
        if not all_ok:
            print("ERROR: Verification still failing.")

    if not args.keep_temp and all_ok:
        cleanup_temp(tmp)
    elif not all_ok:
        print("Skipping temp cleanup - verification failed")

    print(f"\nDone! Results in: {outdir}")


if __name__ == "__main__":
    main()
