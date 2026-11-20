"""LightDock adapter with score-only ranking and reference-only evaluation."""
from __future__ import annotations

import importlib.metadata
import gzip
import json
import os
from pathlib import Path
import shutil
import sys

import gemmi
import numpy as np

from .audit import Blocked, sha256, write_csv, write_json
from .backends import ROOT, available, dockq_grade, executable, run_dockq
from .structures import geometry, make_structure, scramble_partners, write_structure, coordinates


def lightdock_commands():
    # Pin v0.9.4: its CLI names differ from unreleased upstream master.
    if not available('lightdock'):
        raise Blocked('LightDock 0.9.4 unavailable. Windows native C-extension build or a Linux/WSL2 environment is required.')
    version = importlib.metadata.version('lightdock')
    if version != '0.9.4':
        raise Blocked(f'LightDock adapter targets 0.9.4; installed {version}; verify CLI before using another release')
    commands = {key: executable(name) for key, name in [('setup', 'lightdock3_setup.py'), ('run', 'lightdock3.py'),
                ('generate', 'lgd_generate_conformations.py'), ('cluster', 'lgd_cluster_bsas.py')]}
    if sys.version_info >= (3, 12):
        commands = {key: [sys.executable, str(ROOT / 'lightdock_compat.py'), args[-1]] for key, args in commands.items()}
    return commands


def read_gso(path):
    """The last GSO column is scoring; luciferin is not binding energy."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith('('):
            continue
        values = line[line.index(')')+1:].split()
        if len(values) not in (4, 6):
            raise ValueError('Unexpected LightDock GSO output schema')
        score = float(values[-1])
        if not np.isfinite(score):
            raise ValueError('Non-finite LightDock score')
        rows.append({'glowworm': len(rows), 'score': score})
    if not rows:
        raise ValueError('LightDock produced no final candidates')
    return rows


def normalize_windows_outputs(directory):
    """LightDock 0.9.4 writes os.linesep through text newline translation.

    On Windows this yields CR-CR-LF and false empty coordinate records. Keep
    exact upstream bytes in gzip, normalize line endings only, and audit it.
    """
    if os.name != 'nt':
        return
    records_path = directory / 'windows_newlines.json'
    records = json.loads(records_path.read_text()) if records_path.exists() else []
    for path in list(directory.rglob('*.dat')) + list(directory.rglob('*.out')):
        content = path.read_bytes()
        if b'\r\r\n' not in content:
            continue
        original = path.with_name(path.name + '.original.gz')
        if original.exists():
            raise ValueError(f'Already preserved original exists for {path}')
        with gzip.open(original, 'wb') as handle:
            handle.write(content)
        before = sha256(path)
        path.write_bytes(content.replace(b'\r\r\n', b'\n').replace(b'\r\n', b'\n'))
        records.append({'file': str(path.relative_to(directory)), 'before_sha256': before, 'after_sha256': sha256(path),
                        'original_bytes_gzip': str(original.relative_to(directory)), 'change': 'CR-CR-LF to LF only'})
    write_json(records_path, records)


def cluster_score_ranked(rows, threshold=4.0):
    """Greedy global ligand-CA RMSD clusters in the fixed receptor frame.

    No experimental coordinates are used. Keep every pose, including clashes.
    """
    representatives = []
    for row in rows:
        st = gemmi.read_structure(row['path'])
        ca = np.array([[a.pos.x, a.pos.y, a.pos.z] for c in st[0] if c.name == 'B' for r in c for a in r if a.name == 'CA'])
        if len(ca) == 0:
            raise ValueError('Candidate lacks ligand CA atoms')
        cluster = None
        for index, representative in enumerate(representatives):
            if representative.shape != ca.shape:
                raise ValueError('Different ligand atom sets across docking candidates')
            if float(np.sqrt(np.mean(np.sum((ca - representative)**2, axis=1)))) < threshold:
                cluster = index
                break
        if cluster is None:
            cluster = len(representatives)
            representatives.append(ca)
        row['global_cluster'] = cluster
    return {'method': 'greedy ligand CA RMSD in fixed receptor frame, no ligand superposition', 'threshold_A': threshold,
            'clusters': len(representatives)}


def dock(run, structure, seeds=(42, 2024, 8888), swarms=400, glowworms=300, steps=100, cores=1, reference=None, scoring='fastdfire'):
    commands = lightdock_commands()
    if min(swarms, glowworms, steps, cores) < 1 or not seeds:
        raise ValueError('Positive sampling parameters and at least one seed required')
    all_rows = []
    for seed in seeds:
        directory = run.path / 'docking' / f'seed_{seed}'
        directory.mkdir()
        scrambled, transforms = scramble_partners(structure, seed)
        write_json(directory / 'input_transforms.json', {'seed': seed, 'transforms': transforms, 'native_restraints': False})
        for chain in scrambled[0]:
            write_structure(make_structure([chain]), directory / f'protein_{chain.name}.pdb')
        receptor, ligand = directory / 'protein_A.pdb', directory / 'protein_B.pdb'
        run.command(commands['setup'] + [receptor, ligand, '-s', str(swarms), '-g', str(glowworms), '--seed_points', str(seed), '--noh', '--noxt'], cwd=directory, timeout=7200)
        setup = json.loads((directory / 'setup.json').read_text())
        normalize_windows_outputs(directory)
        # Actual generated swarm count is determined by the surface-covering setup.
        console = run.command(commands['run'] + ['setup.json', str(steps), '-s', scoring, '-sg', str(seed), '-c', str(cores)], cwd=directory, timeout=86400)
        if 'Lightdock has failed' in console:
            raise RuntimeError('LightDock caught an internal exception despite exit code 0; inspect native stdout/stderr logs')
        normalize_windows_outputs(directory)
        swarm_dirs = sorted(directory.glob('swarm_*'), key=lambda p: int(p.name.split('_')[-1]))
        if not swarm_dirs:
            raise ValueError('LightDock generated no swarms')
        for swarm in swarm_dirs:
            gso = swarm / f'gso_{steps}.out'
            if not gso.is_file():
                raise RuntimeError(f'LightDock did not generate {gso}; exit code alone does not establish success')
            rows = read_gso(gso)
            if len(rows) != glowworms:
                raise ValueError(f'{swarm.name}: expected {glowworms} candidates, got {len(rows)}')
            run.command(commands['generate'] + [receptor, ligand, gso, str(len(rows))], cwd=directory, timeout=7200)
            run.command(commands['cluster'] + [gso], cwd=swarm, timeout=7200)
            for row in rows:
                path = swarm / f'lightdock_{row["glowworm"]}.pdb'
                candidate = gemmi.read_structure(str(path))
                chain_ids = {c.name for c in candidate[0]}
                if chain_ids != {'A', 'B'}:
                    raise ValueError(f'Unexpected candidate chains: {chain_ids}')
                # Check that the rigid receptor is consistently centered for
                # meaningful inter-candidate clustering.
                rec = coordinates(candidate, 'A')
                if all_rows:
                    first = coordinates(gemmi.read_structure(all_rows[0]['path']), 'A')
                    if rec.shape != first.shape or np.max(np.abs(rec-first)) > 0.02:
                        raise ValueError('Receptor coordinates differ between rigid-body candidates')
                row.update(seed=seed, swarm=int(swarm.name.split('_')[-1]), path=str(path), **geometry(candidate))
                all_rows.append(row)
    # LightDock's upstream write_ranking_to_file sorts scoring in descending
    # order. Keep its convention, never interpret this value as affinity.
    all_rows.sort(key=lambda x: (-x['score'], x['seed'], x['swarm'], x['glowworm']))
    for i, row in enumerate(all_rows):
        row['rank'] = i + 1
    cluster_info = cluster_score_ranked(all_rows)
    write_json(run.path / 'docking/ranking.json', {'score_function': scoring, 'score_direction': 'higher_is_better', 'reference_used_for_ranking': False,
               'sampling': {'seeds': list(seeds), 'requested_swarms': swarms, 'glowworms': glowworms, 'steps': steps, 'anm': False},
               'clustering': cluster_info, 'candidates': all_rows})
    write_csv(run.path / 'docking/ranking.csv', all_rows, list(all_rows[0]))
    shutil.copy2(all_rows[0]['path'], run.path / 'docking/score_top1.pdb')
    if reference is None:
        run_dockq(run, all_rows[0]['path'], None)
        return all_rows
    # Reference never enters any setup, optimization or ranking command.
    request = run.path / 'validation/dockq_batch_request.json'
    metrics_path = run.path / 'validation/dockq_all_candidates.jsonl'
    write_json(request, {'reference': str(reference), 'candidates': [{'path': row['path'], 'rank': row['rank']} for row in all_rows],
                        'output': str(metrics_path), 'dockq_script': os.environ.get('EASY_PPI_DOCKQ_SCRIPT')})
    run.command([sys.executable, ROOT / 'dockq_batch.py', request], timeout=86400)
    metrics = [json.loads(line) for line in metrics_path.read_text(encoding='utf-8').splitlines()]
    if len(metrics) != len(all_rows):
        raise ValueError('DockQ did not evaluate every candidate')
    evaluated = []
    for row, result in zip(all_rows, metrics):
        if row['rank'] != result['rank'] or not 0 <= result['DockQ'] <= 1:
            raise ValueError('DockQ candidate correspondence/score invalid')
        evaluated.append({**row, **result, 'grade': dockq_grade(result['DockQ'])})
    write_csv(run.path / 'validation/candidate_metrics.csv', evaluated, list(evaluated[0]))
    best = max(evaluated, key=lambda x: x['DockQ'])
    summary = {'candidate_count': len(evaluated), 'score_top1': evaluated[0], 'oracle_best_for_sampling_only': best,
               'top_k_hits': {str(k): any(row['DockQ'] >= 0.23 for row in evaluated[:k]) for k in (1, 5, 10)},
               'sampling_validation': 'PASS' if best['DockQ'] >= 0.23 else 'FAIL',
               'ranking_validation': 'PASS' if any(row['DockQ'] >= 0.23 for row in evaluated[:10]) else 'FAIL',
               'bound_redocking_limitation': 'does not establish unbound docking performance'}
    write_json(run.path / 'validation/redocking.json', summary)
    return all_rows
