"""OpenBabel preparation; Arpeggio remains the only chemical classifier."""
from __future__ import annotations
import math
from pathlib import Path

import gemmi
from .audit import Blocked, write_json
from .structures import write_structure


def heavy_atoms(st):
    result = {}
    for c in st[0]:
        for r in c:
            for a in r:
                if a.is_hydrogen():
                    continue
                key = (c.name, str(r.seqid), r.name, a.name)
                if key in result:
                    raise ValueError(f'Duplicate atom identity {key}')
                result[key] = (a.pos.x, a.pos.y, a.pos.z)
    return result


def protonate(pdb, output_dir, ph=7.4):
    from openbabel import openbabel as ob
    output_dir = Path(output_dir)
    ob.obErrorLog.ClearLog()
    conversion = ob.OBConversion()
    if not conversion.SetInAndOutFormats('pdb', 'pdb'):
        raise Blocked('OpenBabel native PDB reader/writer unavailable')
    mol = ob.OBMol()
    if not conversion.ReadFile(mol, str(pdb)):
        raise ValueError('OpenBabel failed to read prepared PDB')
    before = mol.NumAtoms()
    mol.AddHydrogens(False, True, ph)
    hydrogen_count = mol.NumAtoms() - before
    if hydrogen_count < before * 0.5:
        raise ValueError(f'Incomplete protonation: added {hydrogen_count} H for {before} heavy atoms')
    pdb_h = output_dir / 'complex_H.pdb'
    if not conversion.WriteFile(mol, str(pdb_h)):
        raise ValueError('Cannot write protonated PDB')
    conversion.CloseOutFile()
    st = gemmi.read_structure(str(pdb_h))
    # OpenBabel assigns duplicate generic H atom names; unique per-residue names
    # are required by Biopython/Arpeggio. Heavy-atom identifiers never change.
    names = []
    for c in st[0]:
        for r in c:
            hindex = 0
            for atom in r:
                if atom.is_hydrogen():
                    hindex += 1
                    previous = atom.name
                    atom.name = f'H{hindex}'
                    names.append({'chain': c.name, 'residue': str(r.seqid), 'old': previous, 'new': atom.name})
    original = heavy_atoms(gemmi.read_structure(str(pdb)))
    prepared = heavy_atoms(st)
    if original.keys() != prepared.keys():
        raise ValueError('Protonation changed heavy-atom identities')
    max_shift = max(math.dist(xyz, prepared[key]) for key, xyz in original.items())
    if max_shift > 0.002:
        raise ValueError(f'Protonation moved heavy atoms by {max_shift} A')
    write_structure(st, pdb_h)
    cif_h = output_dir / 'complex_H.cif'
    write_structure(st, cif_h)
    write_json(output_dir / 'protonation.json', {'engine': 'OpenBabel native PDB reader + AddHydrogens', 'openbabel_version': ob.OBReleaseVersion(),
               'ph': ph, 'heavy_atoms': before, 'added_hydrogens': hydrogen_count, 'heavy_atom_max_shift_A': max_shift,
               'hydrogen_name_mapping': names, 'minimization': False,
               'warnings': list(ob.obErrorLog.GetMessagesOfLevel(ob.obWarning)),
               'limitations': ['inferred protonation; no experimental hydrogen positions', 'missing heavy atoms are not reconstructed']})
    return cif_h
