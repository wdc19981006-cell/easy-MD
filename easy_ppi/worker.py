"""Process-isolated scientific adapters. No fallback chemistry or SASA implementation."""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from easy_ppi.audit import write_json


def sasa(folder, output):
    import freesasa
    classifier = freesasa.Classifier.getStandardClassifier('protor')
    parameters = freesasa.Parameters({'algorithm': freesasa.LeeRichards, 'probe-radius': 1.4, 'n-slices': 100, 'n-threads': 1})
    structures = {key: freesasa.Structure(str(folder / filename), classifier=classifier)
                  for key, filename in [('A', 'protein_A.pdb'), ('B', 'protein_B.pdb'), ('AB', 'complex.pdb')]}
    if structures['A'].nAtoms() + structures['B'].nAtoms() != structures['AB'].nAtoms():
        raise ValueError('SASA atom sets differ between complex and monomers')
    results = {key: freesasa.calc(st, parameters) for key, st in structures.items()}
    totals = {key: float(r.totalArea()) for key, r in results.items()}
    if not all(math.isfinite(x) and x > 0 for x in totals.values()):
        raise ValueError('Invalid SASA totals')
    per_residue = []
    bound = results['AB'].residueAreas()
    buried_sides = {}
    for chain in ('A', 'B'):
        buried_sides[chain] = 0.0
        for number, isolated in results[chain].residueAreas()[chain].items():
            complex_area = bound[chain][number].total
            delta = isolated.total - complex_area
            if delta < -0.5:
                raise ValueError('Negative residue burial: inconsistent SASA calculation')
            buried_sides[chain] += delta
            per_residue.append({'chain': chain, 'residue': number, 'residue_name': isolated.residueType,
                                'isolated_sasa_A2': isolated.total, 'complex_sasa_A2': complex_area, 'buried_A2': delta})
    total_buried = totals['A'] + totals['B'] - totals['AB']
    if total_buried < -0.5 or abs(sum(buried_sides.values()) - total_buried) > 0.01:
        raise ValueError('SASA conservation check failed')
    write_json(output, {'engine': 'FreeSASA', 'units': 'angstrom^2', 'parameters': {'algorithm': 'LeeRichards', 'probe_radius_A': 1.4,
                'n_slices': 100, 'classifier': 'protor', 'hydrogens': False, 'threads': 1}, 'sasa': totals,
                'buried_total_A2': total_buried, 'interface_area_half_A2': total_buried / 2,
                'buried_by_partner_A2': buried_sides, 'residues': per_residue})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['sasa'])
    parser.add_argument('folder', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    sasa(args.folder, args.output)
