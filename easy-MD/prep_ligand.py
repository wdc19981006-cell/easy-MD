"""Ligand preparation: SMILES → 3D SDF → PDBQT.

Usage:  "D:/Python312/python.exe" prep_ligand.py "<SMILES>"
Output: ligand.sdf, ligand.pdbqt (in current working directory)
"""
import os, sys

# RDKit version-safe MMFFOptimizeMolecule wrapper
from rdkit import Chem
from rdkit.Chem import AllChem

if len(sys.argv) < 2:
    print("Usage: prep_ligand.py '<SMILES>'")
    sys.exit(1)

smiles = sys.argv[1]
print(f"SMILES: {smiles}")

print("Generating 3D structure...")
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("ERROR: Invalid SMILES")
    sys.exit(1)

mol = Chem.AddHs(mol)

status = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
print(f"  ETKDG embedding: status={status}")

# MMFFOptimizeMolecule API differs by RDKit version:
#   Old API → int (0=converged, 1=not converged)
#   New API → tuple (status, energy)
result = AllChem.MMFFOptimizeMolecule(mol)
if isinstance(result, tuple):
    print(f"  MMFF: converged={result[0]==0}, energy={result[1]:.2f}")
else:
    print(f"  MMFF: status={result} (0=converged)")

writer = Chem.SDWriter('ligand.sdf')
writer.write(mol)
writer.close()
print("Saved ligand.sdf")
print("Done!")
