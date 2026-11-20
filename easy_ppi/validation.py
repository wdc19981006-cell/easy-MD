"""Read-only artifact validation and geometric invariance controls."""
from pathlib import Path
import json
import math

import gemmi
import numpy as np
from scipy.spatial.transform import Rotation

from .audit import Blocked, sha256, write_json
from .backends import run_sasa
from .structures import make_structure, write_structure


def sasa_invariance(run, structure, baseline):
    if baseline is None:
        raise Blocked('Baseline FreeSASA calculation required')
    controls = {}
    for name in ('rotated_translated', 'chain_swap'):
        st = structure.clone()
        if name == 'rotated_translated':
            matrix = Rotation.from_euler('xyz', [37, 17, 91], degrees=True).as_matrix()
            for chain in st[0]:
                for residue in chain:
                    for atom in residue:
                        transformed = matrix @ np.array([atom.pos.x, atom.pos.y, atom.pos.z]) + np.array([22, -31, 17])
                        atom.pos = gemmi.Position(*transformed)
        else:
            for chain in st[0]:
                chain.name = 'B' if chain.name == 'A' else 'A'
                for residue in chain:
                    residue.subchain = chain.name
        folder = run.path / 'validation' / name
        folder.mkdir()
        write_structure(st, folder / 'complex.pdb')
        for chain in st[0]:
            write_structure(make_structure([chain]), folder / f'protein_{chain.name}.pdb')
        value = run_sasa(run, folder, label=name)
        delta = abs(value['buried_total_A2'] - baseline['buried_total_A2'])
        # Lee-Richards slice integration has a small orientation dependence;
        # PDB export also rounds coordinates to 0.001 A.
        tolerance = max(2.0, 0.01 * baseline['buried_total_A2']) if name == 'rotated_translated' else 0.01
        controls[name] = {'buried_total_A2': value['buried_total_A2'], 'difference_A2': delta, 'tolerance_A2': tolerance, 'pass': delta <= tolerance}
    write_json(run.path / 'validation/sasa_invariance.json', controls)
    if not all(x['pass'] for x in controls.values()):
        raise ValueError('FreeSASA invariance controls failed')


def chemistry_check(run, interactions):
    if not interactions:
        raise Blocked('Arpeggio output required')
    checks = {}
    # Preregistered literature-supported 1BRS interface contacts, by source
    # author numbering mapped A->A and D->B; counts are not hardcoded.
    for a, b in [('83', '39'), ('87', '39'), ('102', '39')]:
        found = []
        for row in interactions['atom_contacts']:
            sides = {row['chain_a']: row['residue_a'], row['chain_b']: row['residue_b']}
            if sides == {'A': a, 'B': b} and set(row['types']) & {'hbond', 'ionic', 'polar'}:
                found.append(row)
        checks[f'Barnase_{a}_Barstar_{b}'] = {'pass': bool(found), 'arpeggio_contacts': found}
    write_json(run.path / 'validation/1brs_hotspots.json', checks)
    if not all(x['pass'] for x in checks.values()):
        raise ValueError('Literature-supported 1BRS charged interface contacts not recovered')


def validate_artifacts(directory):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8'))
    failures = []
    for record in manifest:
        path = (directory / record['path']).resolve()
        if directory not in path.parents:
            raise ValueError('Manifest path escapes run directory')
        if not path.is_file() or path.stat().st_size != record['bytes'] or sha256(path) != record['sha256']:
            failures.append(record['path'])
    state = json.loads((directory / 'run.json').read_text(encoding='utf-8'))
    return {'artifact_integrity': 'PASS' if not failures else 'FAIL', 'failed_files': failures,
            'stages': {key: value['status'] for key, value in state['stages'].items()},
            'scientific_validation': state['scientific_validation'],
            'note': 'Integrity is separate from scientific validity; do not equate a valid manifest with successful docking'}


def reference_isolation(directory):
    directory = Path(directory)
    state = json.loads((directory / 'run.json').read_text(encoding='utf-8'))
    reference = str((directory / 'prepared/complex.pdb').resolve()).replace('\\', '/').casefold()
    checked = 0
    for command in state.get('commands', []):
        argv = command['argv']
        script_names = {Path(x).name for x in argv}
        if not script_names & {'lightdock3_setup.py', 'lightdock3.py'}:
            continue
        checked += 1
        if any(str(x).replace('\\', '/').casefold() == reference for x in argv):
            raise ValueError('Experimental reference path leaked into docking inputs')
    setups = list(directory.glob('docking/seed_*/setup.json'))
    for path in setups:
        setup = json.loads(path.read_text(encoding='utf-8'))
        if setup.get('restraints'):
            raise ValueError('Native/interface restraints supplied in supposedly unrestrained benchmark')
    if not setups or checked < 2 * len(setups):
        raise Blocked('Insufficient completed setup/run evidence for reference isolation audit')
    return {'status': 'PASS', 'docking_commands_checked': checked, 'setup_files_checked': len(setups),
            'scope': 'reference excluded from commands and restraints; independent centering/rotation recorded in input_transforms.json'}
