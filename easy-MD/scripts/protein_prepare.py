#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Protein preparation: HETATM classification, chain selection, pdbfixer, PDBQT.

Enhanced v2.0 - replaces blanket removeHeterogens() with intelligent classification.

Usage:
    python protein_prepare.py <pdb_file> [--chain A] [--ph 7.4] [--obabel PATH]
    python protein_prepare.py <pdb_file> --json  # Output classification as JSON

Output: protein_clean.pdb, receptor.pdbqt, native_ligand_info.txt
"""

import os, sys, argparse, json, subprocess, math
from pathlib import Path
import numpy as np

from common import (
    configured_path,
    enter_work_dir,
    get_config,
    load_config,
    normalize_name_list,
)

# --- HETATM Classification Lists ---
WATER_NAMES = {
    'HOH', 'WAT', 'H2O', 'DOD'
}
BUFFER_NAMES = {
    'PEG', 'GOL', 'EDO', 'SO4', 'PO4', 'ACT', 'EOH', 'MPD',
    'PG4', 'PGE', 'PG0', 'BOG', 'OCT', 'DMS', 'DMF', 'DMSO',
    'IPA', 'BU1', 'BU2', 'BU3', 'PX4', '1PE', '2PE', 'PE3',
    'TRS', 'BME', 'CIT', 'ACY', 'AZI', 'BCT', 'BEN',
    'CAC', 'CDL', 'CLR', 'CPS', 'DIO', 'DTD', 'ETF', 'FMT',
    'IMD', 'IOD', 'LDA', 'MES', 'MLA', 'MRD', 'MYR', 'N8E',
    'NA', 'NO3', 'P6G', 'PE8', 'PGO', 'PTY', 'SDS', 'SIN',
    'SPD', 'SPK', 'TAR', 'TAS', 'TCE', 'TMA', 'URE', 'XPE',
}
METAL_NAMES = {
    'MG', 'ZN', 'FE', 'MN', 'CA', 'NA', 'K', 'CO', 'NI', 'CU',
    'CD', 'HG', 'PT', 'AU', 'AG', 'SR', 'BA', 'CS', 'RB', 'LI',
}
COFACTOR_NAMES = {
    'FAD', 'NAD', 'NAP', 'NDP', 'HEM', 'HEME', 'FMN', 'COA',
    'PLP', 'TPP', 'SAM', 'SAH', 'FES', 'SF4', 'F3S', 'CLF',
    'BCL', 'BCB', 'BPB', 'CLA', 'PHO', 'MQ7', 'MQ8',
    'MQ9', 'U10', 'PLM', 'OLC', 'OLB', 'LHG', 'LMN', 'UND',
}


def classify_hetatm(pdb_file):
    """Scan PDB HETATM and classify every heterogen residue."""
    het_residues = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith('HETATM'):
                continue
            resn = line[17:20].strip()
            chain = line[21:22]
            try:
                resi = int(line[22:26].strip())
            except ValueError:
                resi = line[22:27].strip()
            key = (resn, chain, resi)
            het_residues.setdefault(key, []).append(line)

    native_ligands, waters, buffers, metals, cofactors = [], [], [], [], []
    for (resn, chain, resi), atoms in het_residues.items():
        rn = resn.upper().strip()
        if rn in WATER_NAMES:
            waters.append((resn, chain, resi, atoms))
        elif rn in BUFFER_NAMES:
            buffers.append((resn, chain, resi, atoms))
        elif rn in METAL_NAMES:
            metals.append((resn, chain, resi, atoms))
        elif rn in COFACTOR_NAMES:
            cofactors.append((resn, chain, resi, atoms))
        elif len(atoms) > 3:
            native_ligands.append((resn, chain, resi, atoms))
    return {
        'native_ligands': native_ligands, 'waters': waters,
        'buffers': buffers, 'metals': metals, 'cofactors': cofactors,
    }


def report_classification(cls):
    """Print HETATM classification report."""
    print("=" * 60)
    print("HETATM Classification Report")
    print("=" * 60)
    print(f"  Native ligand(s): {len(cls['native_ligands'])}")
    for resn, chain, resi, atoms in cls['native_ligands']:
        xs = [float(l[30:38]) for l in atoms]
        ys = [float(l[38:46]) for l in atoms]
        zs = [float(l[46:54]) for l in atoms]
        print(f"    {resn} (chain {chain}, res {resi}): "
              f"{len(atoms)} atoms, "
              f"center=({sum(xs)/len(xs):.1f}, {sum(ys)/len(ys):.1f}, {sum(zs)/len(zs):.1f})")
    print(f"  Water molecules:   {len(cls['waters'])}")
    print(f"  Buffer/additives:  {len(cls['buffers'])}")
    for resn, chain, resi, _ in cls['buffers'][:5]:
        print(f"    {resn} (chain {chain})")
    if len(cls['buffers']) > 5:
        print(f"    ... and {len(cls['buffers']) - 5} more")
    print(f"  Metal ions:        {len(cls['metals'])}")
    for resn, chain, resi, _ in cls['metals']:
        print(f"    {resn} (chain {chain})")
    print(f"  Cofactors:         {len(cls['cofactors'])}")
    for resn, chain, resi, _ in cls['cofactors']:
        print(f"    {resn} (chain {chain})")
    print("=" * 60)


def select_chains(pdb_file, classification):
    """Intelligently select protein chain(s) based on native ligand contacts."""
    chains = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21:22]
                chains[chain] = chains.get(chain, 0) + 1
    print(f"\nChain composition:")
    for cid, count in sorted(chains.items()):
        print(f"  Chain {cid}: {count} atoms")
    if len(chains) == 1:
        selected = list(chains.keys())
        print(f"Single chain -- using Chain {selected[0]}")
        return selected
    nl = classification['native_ligands']
    if nl:
        lig_chains = set(lg[1] for lg in nl)
        print(f"Native ligand(s) found in chain(s): {lig_chains}")
        selected = list(lig_chains)
        return selected
    selected = ['A']
    other = [c for c in chains if c != 'A']
    print("No native ligand found. Defaulting to Chain A.")
    if other:
        print(f"  Other chains available: {other}")
    return selected


def decide_metals(classification, selected_chains, proximity_cutoff=10.0):
    """Decide which metals to keep based on proximity to binding site."""
    metals_to_keep, metals_to_remove = [], []
    if not classification['metals']:
        return metals_to_keep, metals_to_remove
    nl = classification['native_ligands']
    if nl:
        all_coords = []
        for resn, chain, resi, atoms in nl:
            for line in atoms:
                all_coords.append([
                    float(line[30:38]), float(line[38:46]), float(line[46:54])
                ])
        lig_center = np.mean(all_coords, axis=0)
        for resn, chain, resi, atoms in classification['metals']:
            if chain not in set(selected_chains):
                metals_to_remove.append((resn, chain, resi, atoms))
                continue
            mc = np.array([
                float(atoms[0][30:38]), float(atoms[0][38:46]), float(atoms[0][46:54])
            ])
            dist = np.linalg.norm(mc - lig_center)
            if dist < proximity_cutoff:
                metals_to_keep.append((resn, chain, resi, atoms))
                print(f"  KEEP {resn} (chain {chain}): {dist:.1f} A from binding site")
            else:
                metals_to_remove.append((resn, chain, resi, atoms))
                print(f"  REMOVE {resn} (chain {chain}): {dist:.1f} A - too far")
    else:
        metals_to_remove = classification['metals']
        print("  No native ligand reference - removing all metal ions")
        for resn, chain, resi, _ in metals_to_remove:
            print(f"    REMOVE {resn} (chain {chain})")
    return metals_to_keep, metals_to_remove


def build_pre_clean_pdb(pdb_file, selected_chains, classification,
                        metals_to_keep, output_path):
    """Write pre-filtered PDB with protein + cofactors + kept metals only."""
    selected_set = set(selected_chains)
    keep_heterogens = set()
    for resn, chain, resi, atoms in classification['cofactors']:
        if chain in selected_set:
            keep_heterogens.add((resn, chain, str(resi)))
            print(f"  Keeping cofactor: {resn} (chain {chain})")
    for resn, chain, resi, atoms in metals_to_keep:
        if chain in selected_set:
            keep_heterogens.add((resn, chain, str(resi)))
            print(f"  Keeping metal: {resn} (chain {chain})")
    keep_atoms = []
    with open(pdb_file, 'r') as f:
        for raw in f:
            line = raw.rstrip('\n\r')
            if line.startswith('ATOM') and line[21:22] in selected_set:
                keep_atoms.append(line[:80])
            elif line.startswith('HETATM'):
                key = (line[17:20].strip(), line[21:22], line[22:26].strip())
                if key in keep_heterogens:
                    keep_atoms.append(line[:80])
    with open(output_path, 'w') as outf:
        for line in keep_atoms:
            outf.write(line + '\n')
        outf.write('END\n')
    print(f"\nPre-cleaned PDB written: {output_path}")
    print(f"  Kept: {len(keep_atoms)} atoms")
    return output_path


def run_pdbfixer(pre_clean_pdb, ph=7.4):
    """Run PDBFixer. Does NOT call removeHeterogens()."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    def prepare(path):
        instance = PDBFixer(filename=str(path))
        instance.findMissingResidues()
        instance.findNonstandardResidues()
        instance.replaceNonstandardResidues()
        instance.findMissingAtoms()
        instance.addMissingAtoms()
        instance.addMissingHydrogens(ph)
        return instance

    heterogen_lines = [
        line for line in Path(pre_clean_pdb).read_text(encoding='utf-8').splitlines()
        if line.startswith('HETATM')
    ]
    used_heterogen_fallback = False
    try:
        fixer = prepare(pre_clean_pdb)
    except (ValueError, KeyError, RuntimeError) as exc:
        if not heterogen_lines:
            raise
        # Some CCD components contain incomplete ideal-coordinate records that
        # PDBFixer cannot parse.  Repair the protein alone, then restore the
        # retained cofactor atoms at their exact crystallographic coordinates.
        protein_only = Path(pre_clean_pdb).with_suffix('.protein_only.pdb')
        source_lines = Path(pre_clean_pdb).read_text(encoding='utf-8').splitlines()
        protein_only.write_text(
            '\n'.join(line for line in source_lines if not line.startswith('HETATM')) + '\n',
            encoding='utf-8',
        )
        print(f"PDBFixer heterogen fallback: {exc}")
        fixer = prepare(protein_only)
        used_heterogen_fallback = True
    protein_clean = 'protein_clean.pdb'
    with open(protein_clean, 'w') as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)
    if used_heterogen_fallback:
        repaired = Path(protein_clean).read_text(encoding='utf-8').splitlines()
        while repaired and repaired[-1].strip() in {'END', ''}:
            repaired.pop()
        repaired.extend(heterogen_lines)
        repaired.append('END')
        Path(protein_clean).write_text('\n'.join(repaired) + '\n', encoding='utf-8')
        print(f"  Restored {len(heterogen_lines)} cofactor atom records from the crystal structure")
    print(f"Protein prepared: {protein_clean}")
    with open(protein_clean, 'r') as f:
        het = set(l[17:20].strip() for l in f if l.startswith('HETATM'))
    if het:
        print(f"  Heterogens in final receptor: {het}")
    return protein_clean


def convert_to_pdbqt(protein_clean_pdb, obabel_path):
    """Convert protein PDB to PDBQT via OpenBabel."""
    result = subprocess.run(
        [obabel_path, protein_clean_pdb, '-O', 'receptor.pdbqt',
         '-xr', '--partialcharge', 'gasteiger'],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: obabel failed: {result.stderr}")
        sys.exit(1)
    print("Converted to receptor.pdbqt")
    return 'receptor.pdbqt'


def save_native_ligand_info(classification, output_path='native_ligand_info.txt'):
    """Save native ligand identity for redocking workflow."""
    nl = classification['native_ligands']
    if nl:
        primary = max(nl, key=lambda x: len(x[3]))
        resn, chain, resi, atoms = primary
        with open(output_path, 'w') as f:
            f.write(f"{resn}\n{chain}\n{resi}\n")
        print(f"Saved native ligand info: {resn} (chain {chain})")
        return resn, chain, resi
    return None, None, None


def get_native_ligand_center(classification, selected_chains):
    """Get geometric center of primary native ligand for docking box."""
    nl = classification['native_ligands']
    chain_ligands = [n for n in nl if n[1] in set(selected_chains)]
    if not chain_ligands:
        chain_ligands = nl
    if not chain_ligands:
        return None, None, None, None
    primary = max(chain_ligands, key=lambda x: len(x[3]))
    resn, chain, resi, atoms = primary
    coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atoms])
    center = coords.mean(axis=0)
    dims = coords.max(axis=0) - coords.min(axis=0)
    dim_max = dims.max()
    print(f"\nNative ligand: {resn} (chain {chain})")
    print(f"  Center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    print(f"  Dimensions: {dims[0]:.1f} x {dims[1]:.1f} x {dims[2]:.1f} A")
    return center, dim_max, resn, chain


def main():
    parser = argparse.ArgumentParser(
        description="Protein preparation: HETATM classify, pdbfixer, PDBQT")
    parser.add_argument("pdb_file", help="Path to input PDB file")
    parser.add_argument("--chain", default=None,
                        help="Force specific chain (overrides auto-detection)")
    parser.add_argument("--site-resname",
                        help="Use this verified co-crystal/modified residue as the box marker")
    parser.add_argument("--pdb-id", default="PROT", help="PDB ID for output naming")
    parser.add_argument("--ph", type=float)
    parser.add_argument("--obabel", help="OpenBabel executable path")
    parser.add_argument("--metal-cutoff", type=float)
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--config", help="Path to easy-md config.yaml")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    pdb_file = str(Path(args.pdb_file).expanduser().resolve())
    enter_work_dir(args.work_dir, config)
    args.ph = args.ph if args.ph is not None else float(
        get_config(config, "protein.ph", 7.4)
    )
    args.metal_cutoff = (
        args.metal_cutoff
        if args.metal_cutoff is not None
        else float(get_config(config, "protein.metal_proximity_cutoff", 10.0))
    )
    args.obabel = args.obabel or configured_path(config, "obabel") or "obabel"

    global WATER_NAMES, BUFFER_NAMES, METAL_NAMES, COFACTOR_NAMES
    hetatm = config.get("hetatm", {})
    WATER_NAMES = normalize_name_list(hetatm.get("water_names")) or WATER_NAMES
    BUFFER_NAMES = normalize_name_list(hetatm.get("buffer_names")) or BUFFER_NAMES
    METAL_NAMES = normalize_name_list(hetatm.get("metal_names")) or METAL_NAMES
    COFACTOR_NAMES = normalize_name_list(hetatm.get("cofactor_names")) or COFACTOR_NAMES

    classification = classify_hetatm(pdb_file)
    if args.site_resname:
        preferred = [
            item for item in classification['native_ligands']
            if item[0].upper() == args.site_resname.upper()
        ]
        if not preferred:
            raise SystemExit(f"Verified site marker {args.site_resname} was not found in {pdb_file}")
        classification['native_ligands'] = preferred
        print(f"Using verified site marker: {args.site_resname.upper()}")
    report_classification(classification)

    if args.json:
        cls_out = {
            'native_ligands': [
                {'resn': r[0], 'chain': r[1], 'resi': str(r[2]), 'natoms': len(r[3])}
                for r in classification['native_ligands']
            ],
            'n_waters': len(classification['waters']),
            'n_buffers': len(classification['buffers']),
            'n_metals': len(classification['metals']),
            'n_cofactors': len(classification['cofactors']),
        }
        print(json.dumps(cls_out, indent=2))

    if args.chain:
        selected_chains = [args.chain]
        print(f"Using forced chain: {args.chain}")
    else:
        selected_chains = select_chains(pdb_file, classification)

    metals_to_keep, _ = decide_metals(classification, selected_chains, args.metal_cutoff)

    pre_clean = build_pre_clean_pdb(
        pdb_file, selected_chains, classification,
        metals_to_keep, f'{args.pdb_id}_pre_clean.pdb')

    protein_clean = run_pdbfixer(pre_clean, args.ph)
    convert_to_pdbqt(protein_clean, args.obabel)
    save_native_ligand_info(classification)

    print("\nDone! Output: protein_clean.pdb, receptor.pdbqt, native_ligand_info.txt")


if __name__ == '__main__':
    main()
