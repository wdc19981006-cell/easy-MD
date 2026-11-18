#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Ligand preparation: SMILES -> multi-conformer 3D SDF -> PDBQT.

Enhanced v2.0 - structure validation, stereochemistry check, multi-conformer
MMFF optimization, drug-likeness check.

Usage:
    python ligand_prepare.py "<SMILES>" [--conformers 10] [--name "Belumosudil"]
    python ligand_prepare.py --sdf input.sdf  # Use existing SDF, skip RDKit

Output: ligand.sdf, ligand.pdbqt (in current working directory)
"""

import os, sys, argparse, subprocess
from pathlib import Path

from common import configured_path, enter_work_dir, get_config, load_config


def validate_structure(mol):
    """Validate SMILES: properties, stereochemistry, charge, drug-likeness."""
    from rdkit.Chem import Descriptors
    print("=" * 60)
    print("STRUCTURE VALIDATION")
    print("=" * 60)
    print(f"  Heavy atoms:      {mol.GetNumHeavyAtoms()}")
    print(f"  Rotatable bonds:  {Descriptors.NumRotatableBonds(mol)}")
    print(f"  HBD / HBA:        {Descriptors.NumHDonors(mol)} / {Descriptors.NumHAcceptors(mol)}")
    print(f"  TPSA:             {Descriptors.TPSA(mol):.1f} A^2")
    print(f"  LogP:             {Descriptors.MolLogP(mol):.2f}")

    # Stereochemistry
    undefined_chiral = []
    for atom in mol.GetAtoms():
        if atom.HasProp('_ChiralityPossible'):
            try:
                if atom.GetChiralTag() == 0:  # CHI_UNSPECIFIED
                    undefined_chiral.append(atom.GetIdx())
            except Exception:
                pass
    if undefined_chiral:
        print(f"  WARNING: {len(undefined_chiral)} undefined chiral centers at atoms {undefined_chiral}")
    else:
        print("  Stereochemistry:  All defined")

    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    print(f"  Formal charge:    {formal_charge}")
    if formal_charge != 0:
        print("  NOTE: Non-zero formal charge. Verify protonation state at pH 7.4.")

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    violations = []
    if mw > 500: violations.append(f'MW={mw:.0f}>500')
    if logp > 5: violations.append(f'LogP={logp:.1f}>5')
    if hbd > 5: violations.append(f'HBD={hbd}>5')
    if hba > 10: violations.append(f'HBA={hba}>10')
    if violations:
        print(f"  Lipinski:         {len(violations)} violation(s): {', '.join(violations)}")
    else:
        print("  Lipinski:         All rules satisfied")
    print("=" * 60)
    return mw, logp, hbd, hba


def generate_conformers(
    smiles,
    num_conformers=10,
    random_seed=42,
    prune_rms_thresh=0.5,
):
    """Generate multiple conformers, select lowest MMFF energy."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("ERROR: Invalid SMILES - cannot parse")
        sys.exit(1)

    mw, logp, hbd, hba = validate_structure(mol)
    mol = Chem.AddHs(mol)

    print(f"\nGenerating up to {num_conformers} conformers...")
    try:
        params = AllChem.ETKDGv3()
    except AttributeError:
        params = AllChem.ETKDG()
    params.randomSeed = int(random_seed)
    params.numThreads = 1
    params.pruneRmsThresh = float(prune_rms_thresh)

    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
    print(f"  Generated {len(conf_ids)} conformers")

    if len(conf_ids) == 0:
        print("  EmbedMultipleConfs failed - falling back to single conformer")
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        conf_ids = [0]

    # MMFF optimize all conformers and calculate the actual minimized energy.
    print("  MMFF optimization...")
    energies = []
    properties = AllChem.MMFFGetMoleculeProperties(mol)
    if properties is None:
        print("ERROR: MMFF parameters are unavailable for this ligand")
        sys.exit(1)
    for cid in conf_ids:
        result = AllChem.MMFFOptimizeMolecule(mol, confId=cid)
        force_field = AllChem.MMFFGetMoleculeForceField(
            mol, properties, confId=cid
        )
        energy = float(force_field.CalcEnergy())
        energies.append((cid, energy, result == 0))

    energies.sort(key=lambda x: x[1])
    best_cid, best_energy, best_converged = energies[0]

    print("\n  Conformer ranking (top 5):")
    for rank, (cid, energy, converged) in enumerate(energies[:5]):
        marker = " <-- SELECTED" if rank == 0 else ""
        print(f"    {rank+1}. Conf {cid}: {energy:.2f} kcal/mol "
              f"({'ok' if converged else 'not conv'}){marker}")
    if len(energies) > 5:
        print(f"    ... and {len(energies) - 5} more")

    print(f"\n  Selected: conformer {best_cid} (energy={best_energy:.2f} kcal/mol)")

    # Write SDF
    writer = Chem.SDWriter('ligand.sdf')
    writer.write(mol, confId=best_cid)
    writer.close()
    print("  Saved ligand.sdf")

    return mol, best_energy, mw, logp


def sdf_to_pdbqt(sdf_path='ligand.sdf', obabel_path='obabel', gen3d=False):
    """Convert ligand SDF to PDBQT via OpenBabel."""
    cmd = [obabel_path, sdf_path, '-O', 'ligand.pdbqt', '--partialcharge', 'gasteiger']
    if gen3d:
        cmd.insert(2, '--gen3d')
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: obabel failed: {result.stderr}")
        sys.exit(1)
    print("Converted to ligand.pdbqt")
    return 'ligand.pdbqt'


def main():
    parser = argparse.ArgumentParser(
        description="Ligand preparation: SMILES to multi-conformer 3D SDF + PDBQT")
    parser.add_argument("smiles", nargs='?', help="SMILES string of the ligand")
    parser.add_argument("--sdf", help="Use existing SDF file (skip RDKit generation)")
    parser.add_argument("--conformers", "-n", type=int)
    parser.add_argument("--name", default="ligand", help="Ligand name for logging")
    parser.add_argument("--obabel", help="OpenBabel executable path")
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--config", help="Path to easy-md config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    source_sdf = str(Path(args.sdf).expanduser().resolve()) if args.sdf else None
    enter_work_dir(args.work_dir, config)
    args.conformers = args.conformers or int(
        get_config(config, "ligand.num_conformers", 10)
    )
    args.obabel = args.obabel or configured_path(config, "obabel") or "obabel"

    if args.sdf:
        print(f"Using existing SDF: {source_sdf}")
        import shutil
        shutil.copy2(source_sdf, 'ligand.sdf')
        sdf_to_pdbqt('ligand.sdf', args.obabel)
    elif args.smiles:
        print(f"Ligand: {args.name}")
        print(f"SMILES: {args.smiles}")
        generate_conformers(
            args.smiles,
            args.conformers,
            get_config(config, "ligand.random_seed", 42),
            get_config(config, "ligand.prune_rms_thresh", 0.5),
        )
        sdf_to_pdbqt('ligand.sdf', args.obabel)
    else:
        parser.error("Either provide SMILES or --sdf")

    print("\nDone! Output: ligand.sdf, ligand.pdbqt")


if __name__ == '__main__':
    main()
