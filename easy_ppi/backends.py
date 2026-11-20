from __future__ import annotations

from collections import Counter
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys

from .audit import Blocked, write_csv, write_json

ROOT = Path(__file__).resolve().parent


def available(module):
    try:
        return importlib.util.find_spec(module) is not None
    except (ValueError, ModuleNotFoundError):
        return False


def executable(name):
    candidates = [Path(sys.executable).parent / name, Path(sys.executable).parent / (name + '.exe')]
    if sys.platform != 'win32':
        candidates.append(Path(sys.executable).parent / name)
    for p in candidates:
        if p.is_file():
            return [sys.executable, str(p)] if p.suffix == '.py' else [str(p)]
    found = shutil.which(name)
    if found:
        return [sys.executable, found] if found.endswith('.py') else [found]
    raise Blocked(f'Missing executable {name}; use the dedicated easy-PPI environment')


def arpeggio_rows(raw):
    if not isinstance(raw, list):
        raise ValueError('Unrecognized Arpeggio JSON schema; expected contact list')
    rows, other = [], []
    for contact in raw:
        bgn, end = contact.get('bgn', {}), contact.get('end', {})
        ca, cb = bgn.get('auth_asym_id'), end.get('auth_asym_id')
        if {ca, cb} != {'A', 'B'}:
            continue
        labels = contact.get('contact', [])
        if isinstance(labels, str):
            labels = [labels]
        if not isinstance(labels, list) or not labels:
            raise ValueError('Missing Arpeggio contact labels')
        if not bgn.get('auth_atom_id') or not end.get('auth_atom_id'):
            other.append(contact)
            continue
        def residue(atom):
            insertion = atom.get('pdbx_PDB_ins_code', '')
            return str(atom['auth_seq_id']) + ('' if insertion in ('?', '.', None, ' ') else insertion)
        distance = float(contact['distance'])
        if not math.isfinite(distance) or distance <= 0:
            raise ValueError('Invalid Arpeggio distance')
        rows.append({'chain_a': ca, 'residue_a': residue(bgn), 'atom_a': bgn['auth_atom_id'],
                     'chain_b': cb, 'residue_b': residue(end), 'atom_b': end['auth_atom_id'],
                     'distance_A': distance, 'types': labels, 'source': 'PDBe Arpeggio'})
    return rows, other


def run_arpeggio(run, folder=None, label='arpeggio', ph=7.4):
    if not available('arpeggio') or not available('openbabel'):
        raise Blocked('PDBe Arpeggio and OpenBabel imports are required')
    folder = Path(folder or run.path / 'prepared')
    out = run.path / 'analysis' / label
    out.mkdir()
    from .protonate import protonate
    hydrogenated = protonate(folder / 'complex.pdb', out, ph)
    # Fixed entrypoint, not arbitrary shell text; CLI performs the upstream chemistry.
    run.command([sys.executable, '-c', 'from arpeggio.scripts.process_protein_cli import main; main()',
                 hydrogenated, '-s', '/A//', '-o', out, '-ph', str(ph), '-wh'], timeout=1800)
    raw_file = out / 'complex_H.json'
    raw = json.loads(raw_file.read_text(encoding='utf-8'))
    rows, other = arpeggio_rows(raw)
    # An empty list is valid for a negative control; the pipeline tests geometry
    # separately before accepting an experimental interface.
    counts = dict(Counter(t for row in rows for t in row['types']))
    summary = {'engine': 'PDBe Arpeggio', 'ph': ph, 'atom_contacts': rows, 'non_atom_contacts': other,
               'counts_by_label': counts, 'ionic_interpretation': 'ionic contacts / candidate salt bridges, not experimentally measured salt bridges'}
    write_json(out / 'interactions.json', summary)
    write_csv(out / 'interactions.csv', [dict(r, types=';'.join(r['types'])) for r in rows],
              ['chain_a', 'residue_a', 'atom_a', 'chain_b', 'residue_b', 'atom_b', 'distance_A', 'types', 'source'])
    return summary


def run_sasa(run, folder=None, label='sasa'):
    if not available('freesasa'):
        raise Blocked('FreeSASA Python bindings are unavailable; no approximate replacement is used')
    output = run.path / 'analysis' / f'{label}.json'
    run.command([sys.executable, ROOT / 'worker.py', 'sasa', folder or run.path / 'prepared', output])
    data = json.loads(output.read_text(encoding='utf-8'))
    write_csv(output.with_suffix('.csv'), data['residues'], list(data['residues'][0]))
    return data


def parse_prodigy(output):
    dg = re.search(r'Predicted binding affinity[^:]*:\s*([-+\d.eE]+)', output)
    kd = re.search(r'Predicted dissociation constant[^:]*:\s*([-+\d.eE]+)', output)
    if not dg or not kd:
        raise ValueError('PRODIGY output lacks predicted DG/Kd')
    delta_g, dissociation = float(dg[1]), float(kd[1])
    if not math.isfinite(delta_g) or not math.isfinite(dissociation) or dissociation <= 0:
        raise ValueError('Invalid PRODIGY prediction')
    return {'delta_g_kcal_mol': delta_g, 'kd_M': dissociation, 'interpretation': 'empirical prediction; not experimental affinity'}


def run_prodigy(run, temperature=25):
    output = run.command(executable('prodigy') + [run.path / 'prepared/complex.pdb', '--selection', 'A', 'B', '--temperature', str(temperature)])
    result = dict(parse_prodigy(output), temperature_C=temperature, engine='PRODIGY')
    write_json(run.path / 'analysis/affinity_empirical.json', result)
    return result


def dockq_grade(score):
    return 'high' if score >= 0.8 else 'medium' if score >= 0.49 else 'acceptable' if score >= 0.23 else 'incorrect'


def run_dockq(run, model, reference, label='dockq'):
    if reference is None:
        result = {'status': 'SKIPPED', 'reason': 'no_experimental_reference', 'DockQ': None, 'fnat': None, 'iRMSD': None, 'LRMSD': None}
        write_json(run.path / 'validation' / f'{label}.json', result)
        run.state['stages'][label] = {'status': 'SKIPPED', 'reason': 'no_experimental_reference'}
        run.save()
        return result
    try:
        command = executable('DockQ')
    except Blocked:
        script = os.environ.get('EASY_PPI_DOCKQ_SCRIPT')
        if not script or not Path(script).is_file():
            raise
        command = [sys.executable, script]
    # v2.1.3's official script fallback can return numpy.float32 which its own
    # --json path cannot serialize. Preserve the native CLI text and derive a
    # plain-JSON record from the published short output, without altering DockQ.
    output = run.command(command + [model, reference, '--mapping', 'AB:AB', '--short', '--n_cpu', '1'])
    (run.path / 'validation' / f'{label}.raw.txt').write_text(output, encoding='utf-8')
    # Capture the four published metrics, not a number inferred from geometry.
    match = re.search(r'(?m)^DockQ\s+([\d.eE+-]+)\s+iRMSD\s+([\d.eE+-]+)\s+LRMSD\s+([\d.eE+-]+)\s+fnat\s+([\d.eE+-]+)', output)
    if not match:
        raise ValueError('Cannot parse DockQ v2 short output')
    score, irmsd, lrmsd, fnat = map(float, match.groups())
    if not all(map(math.isfinite, [score, irmsd, lrmsd, fnat])) or not 0 <= score <= 1 or not 0 <= fnat <= 1:
        raise ValueError('Invalid DockQ metrics')
    result = {'status': 'SUCCESS', 'DockQ': score, 'iRMSD': irmsd, 'LRMSD': lrmsd, 'fnat': fnat,
              'grade': dockq_grade(score), 'mapping': 'AB:AB', 'reference': str(reference), 'precision': 'upstream CLI printed precision'}
    write_json(run.path / 'validation' / f'{label}.json', result)
    return result
