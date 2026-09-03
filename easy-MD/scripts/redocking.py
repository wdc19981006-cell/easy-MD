#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Redocking validation: extract native ligand, dock, compute RMSD.

Validates the docking protocol by removing the co-crystal ligand, docking it
back, and comparing the docked pose to the crystal pose.

Usage:
    python redocking.py --pdb original.pdb --native-info native_ligand_info.txt
    python redocking.py --skip  # Create placeholder when no native ligand

Output: redocking_result.txt, redocking_info.json
"""

import os, sys, argparse, json, subprocess
from pathlib import Path

from common import configured_path, enter_work_dir, get_config, load_config


def check_availability(native_info_file='native_ligand_info.txt'):
    """Check if redocking is possible."""
    if not os.path.exists(native_info_file):
        print("=" * 60)
        print("REDOCKING VALIDATION: NOT AVAILABLE")
        print("=" * 60)
        print("No co-crystal ligand found in this PDB structure.")
        print("=" * 60)
        with open('redocking_result.txt', 'w') as f:
            f.write("No redocking validation available - no co-crystal ligand in PDB.\n")
        return False, None, None, None
    with open(native_info_file) as f:
        lines = f.read().strip().split('\n')
    lig_name = lines[0]
    lig_chain = lines[1] if len(lines) > 1 else 'A'
    lig_residue = lines[2] if len(lines) > 2 else None
    print(f"Redocking target: {lig_name} (chain {lig_chain})")
    return True, lig_name, lig_chain, lig_residue


def extract_native_ligand(pdb_file, lig_name, lig_chain=None, lig_residue=None):
    """Extract native ligand HETATM records from original PDB."""
    native_pdb = f'native_ligand_{lig_name}.pdb'
    with open(pdb_file, 'r') as inf, open(native_pdb, 'w') as outf:
        for line in inf:
            if not line.startswith('HETATM'):
                continue
            if line[17:20].strip() != lig_name:
                continue
            if lig_chain and line[21:22].strip() != str(lig_chain).strip():
                continue
            if lig_residue and line[22:26].strip() != str(lig_residue).strip():
                continue
            outf.write(line)
        outf.write('END\n')
    # Count
    with open(native_pdb, 'r') as f:
        natoms = sum(1 for l in f if l.startswith('HETATM'))
    if natoms == 0:
        print(f"WARNING: Could not extract {lig_name} atoms from PDB")
        return None
    print(f"Extracted {natoms} atoms of native ligand {lig_name}")
    return native_pdb


def prepare_native_ligand(native_pdb, obabel_path='obabel'):
    """Convert native ligand PDB to PDBQT."""
    result = subprocess.run(
        [obabel_path, native_pdb, '-O', 'native_ligand.pdbqt',
         '--partialcharge', 'gasteiger'],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: obabel failed: {result.stderr}")
        return None
    print("Prepared native_ligand.pdbqt")
    return 'native_ligand.pdbqt'


def dock_native_ligand(vina_path, config='config.txt', cpu=1, seed=42):
    """Dock native ligand back into its binding site."""
    os.makedirs('redock', exist_ok=True)
    for f in ['receptor.pdbqt', 'config.txt']:
        if os.path.exists(f):
            import shutil
            shutil.copy2(f, os.path.join('redock', f))
    if os.path.exists('native_ligand.pdbqt'):
        import shutil
        shutil.copy2('native_ligand.pdbqt', os.path.join('redock', 'ligand.pdbqt'))
    log_file = os.path.join('redock', 'redock_log.txt')
    cmd = [
        vina_path, '--config', 'config.txt', '--out', 'docked_native.pdbqt',
        '--seed', str(seed), '--cpu', str(cpu),
    ]
    print(f"Running redocking Vina...")
    with open(log_file, 'w') as lf:
        proc = subprocess.Popen(cmd, cwd='redock', stdout=lf, stderr=subprocess.STDOUT)
        proc.wait()
    print(f"Vina redocking complete (exit {proc.returncode})")
    return 'redock/docked_native.pdbqt'


def split_poses(docked_pdbqt, obabel_path='obabel'):
    """Split docked PDBQT into individual PDB files."""
    result = subprocess.run(
        [obabel_path, docked_pdbqt, '-O', 'redock/docked_native.pdb', '-m'],
        capture_output=True, text=True, cwd='.')
    return 'redock/docked_native1.pdb'


def calculate_rmsd(crystal_pose, docked_pose):
    """Calculate RMSD between docked and crystal poses using RDKit."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    docked_mol = Chem.MolFromPDBFile(docked_pose, removeHs=False)
    crystal_mol = Chem.MolFromPDBFile(crystal_pose, removeHs=False)

    if docked_mol is None or crystal_mol is None:
        print("ERROR: Could not parse ligand structures for RMSD")
        return None

    try:
        docked_heavy = Chem.RemoveHs(docked_mol)
        crystal_heavy = Chem.RemoveHs(crystal_mol)
        rmsd = rdMolAlign.GetBestRMS(docked_heavy, crystal_heavy)
        print(f"Redocking symmetry-aware heavy-atom RMSD: {rmsd:.2f} A")
        return rmsd
    except Exception as e:
        print(f"Standard alignment failed: {e}")
        try:
            rmsd = rdMolAlign.CalcRMS(docked_mol, crystal_mol)
            print(f"Redocking RMSD (no alignment): {rmsd:.2f} A")
            return rmsd
        except Exception:
            print("RMSD calculation failed")
            return None


def assess_validation(rmsd, validated=2.0, borderline=3.0):
    """Assess protocol validation based on RMSD."""
    if rmsd is None:
        return "COULD NOT COMPUTE"
    if rmsd < validated:
        return f"VALIDATED (RMSD < {validated:.1f} A)"
    elif rmsd < borderline:
        return (
            f"BORDERLINE ({validated:.1f} <= RMSD < "
            f"{borderline:.1f} A)"
        )
    else:
        return f"NOT VALIDATED (RMSD >= {borderline:.1f} A)"


def main():
    parser = argparse.ArgumentParser(
        description="Redocking validation for docking protocol")
    parser.add_argument("--pdb", default="original.pdb",
                        help="Original PDB file (for extracting native ligand)")
    parser.add_argument("--native-info", default="native_ligand_info.txt",
                        help="Native ligand info file from protein_prepare.py")
    parser.add_argument("--vina")
    parser.add_argument("--obabel")
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--config", help="Path to easy-md config.yaml")
    parser.add_argument("--skip", action="store_true",
                        help="Skip redocking (no native ligand available)")
    args = parser.parse_args()

    config = load_config(args.config)
    enter_work_dir(args.work_dir, config)
    args.pdb = str(Path(args.pdb).expanduser().resolve())
    args.native_info = str(Path(args.native_info).expanduser().resolve())
    args.vina = args.vina or configured_path(config, "vina")
    args.obabel = args.obabel or configured_path(config, "obabel") or "obabel"
    if not args.vina and not args.skip:
        parser.error("Vina path is required via --vina or config.yaml")
    seed = int(get_config(config, "redocking.seed", 42))
    validated = float(get_config(config, "redocking.rmsd_validated", 2.0))
    borderline = float(get_config(config, "redocking.rmsd_borderline", 3.0))

    if args.skip:
        check_availability(args.native_info)
        with open('redocking_info.json', 'w') as f:
            json.dump({'native_ligand': None, 'redocking_rmsd': None,
                       'validation_status': 'Not available'}, f, indent=2)
        return

    available, lig_name, lig_chain, lig_residue = check_availability(args.native_info)
    if not available:
        with open('redocking_info.json', 'w') as f:
            json.dump({'native_ligand': None, 'redocking_rmsd': None,
                       'validation_status': 'Not available'}, f, indent=2)
        return

    # Extract native ligand
    native_pdb = extract_native_ligand(
        args.pdb, lig_name, lig_chain, lig_residue
    )
    if native_pdb is None:
        with open('redocking_info.json', 'w') as f:
            json.dump({'native_ligand': lig_name, 'redocking_rmsd': None,
                       'validation_status': 'Extraction failed'}, f, indent=2)
        return

    # Prepare
    prepare_native_ligand(native_pdb, args.obabel)

    # Dock
    dock_native_ligand(args.vina, cpu=args.cpu, seed=seed)

    # Split and get best pose
    docked_pose1 = split_poses('redock/docked_native.pdbqt', args.obabel)

    # RMSD
    rmsd = calculate_rmsd(native_pdb, docked_pose1)
    status = assess_validation(rmsd, validated, borderline)

    # Report
    print("\n" + "=" * 60)
    print("REDOCKING VALIDATION")
    print("=" * 60)
    print(f"  Native ligand:   {lig_name}")
    if rmsd is not None:
        print(f"  Redocking RMSD:  {rmsd:.2f} A")
    print(f"  Status:          {status}")
    if rmsd is not None and rmsd < validated:
        print("\n  The docking protocol reliably reproduces the crystal")
        print("  binding pose. Parameters are appropriate for this system.")
    print("=" * 60)

    # Save results
    with open('redocking_result.txt', 'w') as f:
        f.write(f"Native ligand: {lig_name}\n")
        f.write(f"Redocking RMSD: {rmsd:.2f} A\n" if rmsd is not None else "Redocking RMSD: N/A\n")
        f.write(f"Validation status: {status}\n")

    with open('redocking_info.json', 'w') as f:
        json.dump({
            'native_ligand': lig_name,
            'redocking_rmsd': round(rmsd, 2) if rmsd is not None else None,
            'validation_status': status,
        }, f, indent=2)

    print("\nSaved redocking_result.txt, redocking_info.json")


if __name__ == '__main__':
    main()
