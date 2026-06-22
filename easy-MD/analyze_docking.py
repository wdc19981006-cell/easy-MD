"""Molecular docking result analysis script.

Usage:
    python analyze_docking.py <output_dir> [pdb_id] [ligand_name]

Example:
    python analyze_docking.py D:/桌面/测试/2O3S_ADP 2O3S ADP
    python analyze_docking.py /d/桌面/测试/4RWE_erlotinib 4RWE erlotinib
"""

import os
import sys
import argparse


def analyze(outdir, pdb_id="TARGET", ligand_name="LIGAND"):
    """Analyze docking results in outdir.

    Parameters
    ----------
    outdir : str
        Path to the docking output directory (contains protein_clean.pdb,
        best_pose.pdb or best_pose1.pdb, docked.pdbqt, log.txt, etc.)
    pdb_id : str
        PDB ID of the target protein (for report headers)
    ligand_name : str
        Name of the docked ligand (for report headers)
    """
    outdir = os.path.abspath(outdir)
    if not os.path.isdir(outdir):
        raise FileNotFoundError(f"Directory not found: {outdir}")

    print(f"Analyzing docking results in: {outdir}")
    print(f"  Protein: {pdb_id}")
    print(f"  Ligand:  {ligand_name}")

    import prolif
    from rdkit import Chem
    import matplotlib
    matplotlib.use('Agg')

    # ── Load protein via MDAnalysis ──────────────────────────────
    from MDAnalysis import Universe

    prot_pdb = os.path.join(outdir, "protein_clean.pdb")
    if not os.path.exists(prot_pdb):
        # Try chain A variant
        for f in os.listdir(outdir):
            if f.startswith("protein_clean") and f.endswith(".pdb"):
                prot_pdb = os.path.join(outdir, f)
                break
        else:
            raise FileNotFoundError(
                f"protein_clean.pdb not found in {outdir}. "
                "Make sure protein preparation (Workflow 1) completed."
            )

    prot_u = Universe(prot_pdb)
    prot = prolif.Molecule.from_mda(prot_u, res_name='UNK')
    print(f"Protein loaded: {prot_u.atoms.n_atoms} atoms from {os.path.basename(prot_pdb)}")

    # ── Load docked ligand ───────────────────────────────────────
    lig_path = os.path.join(outdir, "best_pose.pdb")
    if not os.path.exists(lig_path):
        lig_path = os.path.join(outdir, "best_pose1.pdb")
    if not os.path.exists(lig_path):
        raise FileNotFoundError(
            f"best_pose.pdb (or best_pose1.pdb) not found in {outdir}. "
            "Run Workflow 5a to split docked.pdbqt into individual PDB poses first."
        )

    lig_mol = Chem.MolFromPDBFile(lig_path, removeHs=False)
    if lig_mol is None:
        raise ValueError(f"RDKit could not parse ligand from {lig_path}")
    lig = prolif.Molecule.from_rdkit(lig_mol)
    print(f"Ligand loaded: {lig_mol.GetNumAtoms()} atoms (incl. H)")

    # ── Compute interaction fingerprint ──────────────────────────
    fp = prolif.Fingerprint([
        'HBDonor', 'HBAcceptor', 'PiStacking',
        'Hydrophobic', 'Cationic', 'Anionic', 'PiCation'
    ])
    fp.run_from_iterable([lig], prot)
    df = fp.to_dataframe()

    if df.empty:
        print("WARNING: No interaction dataframe returned — prolif may not have found contacts.")
        return None

    # ── Extract best binding energy from log ─────────────────────
    best_energy = None
    log_path = os.path.join(outdir, "log.txt")
    if os.path.exists(log_path):
        with open(log_path, 'r') as lf:
            for line in lf:
                if line.strip().startswith('1 '):
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            best_energy = float(parts[1])
                        except ValueError:
                            pass
                    break

    # ── Write interaction report ─────────────────────────────────
    int_file = os.path.join(outdir, "interactions.txt")
    with open(int_file, 'w', encoding='utf-8') as f:
        f.write(f"{pdb_id} + {ligand_name} — Molecular Docking Interaction Analysis\n")
        f.write("=" * 65 + "\n\n")

        if best_energy is not None:
            f.write(f"Best binding energy: {best_energy:.2f} kcal/mol\n\n")
            print(f"Best binding energy: {best_energy:.2f} kcal/mol")

        f.write("Interaction types detected:\n")
        f.write("-" * 35 + "\n")

        found_any = False
        type_count = 0
        for key in df.columns:
            interactions = df[key].iloc[0]
            if interactions and len(interactions) > 0:
                found_any = True
                type_count += 1
                res_list = []
                for inter in interactions:
                    try:
                        if hasattr(inter, 'resname') and hasattr(inter, 'resid'):
                            res_list.append(f"{inter.resname}{inter.resid}")
                    except Exception:
                        pass
                if res_list:
                    unique = sorted(set(res_list))
                    line = f"  {key}: {', '.join(unique[:20])}"
                    print(line)
                    f.write(line + "\n")
                else:
                    line = f"  {key}: {len(interactions)} contact(s)"
                    print(line)
                    f.write(line + "\n")

        if not found_any:
            msg = "  (No specific interactions detected at default cutoff distances)"
            print(msg)
            f.write(msg + "\n")

        f.write(f"\nTotal interaction types found: {type_count}\n")

    print(f"\nInteractions saved to: {int_file}")
    return int_file


def main():
    parser = argparse.ArgumentParser(
        description="Analyze molecular docking results (prolif interaction fingerprint)"
    )
    parser.add_argument(
        "output_dir",
        help="Path to the docking output directory"
    )
    parser.add_argument(
        "pdb_id", nargs="?", default="TARGET",
        help="PDB ID of the target protein (for report header)"
    )
    parser.add_argument(
        "ligand_name", nargs="?", default="LIGAND",
        help="Name of the docked ligand (for report header)"
    )
    args = parser.parse_args()

    analyze(args.output_dir, args.pdb_id, args.ligand_name)


if __name__ == '__main__':
    main()
