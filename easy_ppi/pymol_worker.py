"""Runs under PyMOL's Python; consumes upstream interactions, never classifies."""
import json
import math
from pathlib import Path
import sys

def main(payload_path):
    import pymol
    pymol.finish_launching(['pymol', '-cq'])
    from pymol import cmd
    payload = json.loads(Path(payload_path).read_text(encoding='utf-8'))
    cmd.reinitialize()
    cmd.load(payload['pdb'], '_source')
    cmd.create('Protein_A', '_source and chain A')
    cmd.create('Protein_B', '_source and chain B')
    cmd.delete('_source')
    cmd.hide('everything')
    cmd.show('cartoon', 'Protein_A or Protein_B')
    cmd.color('cyan', 'Protein_A')
    cmd.color('salmon', 'Protein_B')
    selectors = [f'({"Protein_A" if c == "A" else "Protein_B"} and resi {r})' for c, r in payload['interface']]
    cmd.create('Interface_Residues', ' or '.join(selectors) if selectors else 'none')
    cmd.show('sticks', 'Interface_Residues')
    cmd.group('Interactions')
    endpoints = []
    draw_types = {'hbond', 'weak_hbond', 'ionic', 'hydrophobic', 'aromatic', 'polar'}
    for i, row in enumerate(payload['interactions']):
        types = draw_types.intersection(row['types'])
        if not types:
            continue
        sels = []
        for suffix in ('a', 'b'):
            obj = 'Protein_A' if row['chain_' + suffix] == 'A' else 'Protein_B'
            selection = f'{obj} and resi {row["residue_" + suffix]} and name {row["atom_" + suffix]}'
            if cmd.count_atoms(selection) != 1:
                raise ValueError('Interaction endpoint does not select one atom: ' + selection)
            sels.append(selection)
        xyz = [cmd.get_atom_coords(s) for s in sels]
        distance = math.sqrt(sum((x-y)**2 for x,y in zip(*xyz)))
        if abs(distance-row['distance_A']) > 0.05:
            raise ValueError('Arpeggio/PyMOL endpoint mismatch')
        name = 'contact_%05d' % i
        cmd.distance(name, sels[0], sels[1])
        cmd.hide('labels', name)
        cmd.set('dash_gap', 0.25, name)
        cmd.set('dash_length', 0.25, name)
        cmd.set('dash_width', 2.0, name)
        cmd.color('yellow' if 'hbond' in types else 'orange' if 'ionic' in types else 'gray60', name)
        cmd.group('Interactions', name)
        endpoints.append({'object': name, 'selections': sels, 'xyz': xyz, 'distance_A': distance, 'types': sorted(types)})
    cmd.bg_color('white')
    cmd.orient('Protein_A or Protein_B')
    cmd.zoom('Protein_A or Protein_B', 5)
    cmd.save(payload['pse'])
    cmd.png(payload['png'], width=1600, height=1200, dpi=200, ray=1, quiet=1)
    expected_atoms = {name: cmd.count_atoms(name) for name in ('Protein_A', 'Protein_B', 'Interface_Residues')}
    cmd.reinitialize()
    cmd.load(payload['pse'])
    names = set(cmd.get_names('all'))
    top = names - {x['object'] for x in endpoints}
    expected = {'Protein_A', 'Protein_B', 'Interface_Residues', 'Interactions'}
    if top != expected:
        raise ValueError('PSE top-level objects differ: ' + str(top))
    for name, count in expected_atoms.items():
        if cmd.count_atoms(name) != count:
            raise ValueError('PSE atom count changed on reload')
    raw = cmd.get_session()
    for row in endpoints:
        if cmd.get_type(row['object']) != 'object:measurement':
            raise ValueError('Saved interaction is not a measurement')
        # The session contains measurement endpoints even when get_raw_distances
        # is absent in a PyMOL edition. Check the saved object, not only inputs.
        saved = next(x for x in raw['names'] if x and x[0] == row['object'])
        actual = saved[5][2][0][1]
        wanted = [v for point in row['xyz'] for v in point]
        if len(actual) != 6 or max(abs(a-b) for a,b in zip(actual, wanted)) > 0.01:
            raise ValueError('Saved distance endpoints differ')
        if float(cmd.get('dash_gap', row['object'])) <= 0:
            raise ValueError('Interaction is not dashed')
    with open(payload['validation'], 'w', encoding='utf-8') as handle:
        json.dump({'status': 'SUCCESS', 'reopened_pse': True, 'top_level': sorted(top), 'atom_counts': expected_atoms,
                   'saved_endpoint_validation': True, 'pymol_version': list(cmd.get_version()), 'interactions': endpoints}, handle, indent=2)
    cmd.quit()

if __name__ == '__main__':
    main(sys.argv[1])
