from __future__ import annotations

import json
from pathlib import Path
import re
import shutil

import gemmi

from .audit import Blocked, Run, sha256, write_json
from .backends import ROOT, run_arpeggio, run_dockq, run_prodigy, run_sasa
from .docking import dock
from ._inspect_environment import inspect as inspect_environment
from .sources import RCSB, entity_for_chain
from .structures import (assembly_pair, geometry, inspect_structure, make_structure,
                         prepare_pair, protein_residues, validate_pdb, write_structure)
from .validation import chemistry_check, sasa_invariance


def success(run, stage):
    return run.state['stages'].get(stage, {}).get('status') == 'SUCCESS'


def visualize(run, interactions, sasa, pymol_python):
    if interactions is None:
        raise Blocked('Arpeggio interactions required before rendering chemical interaction lines')
    if not pymol_python:
        raise Blocked('Set --pymol-python to a Python interpreter with PyMOL')
    interface = {(row['chain_' + end], row['residue_' + end]) for row in interactions['atom_contacts'] for end in ('a', 'b')}
    if sasa:
        interface.update((r['chain'], r['residue']) for r in sasa['residues'] if r['buried_A2'] > 1.0)
    out = run.path / 'visualization'
    payload = {'pdb': str(run.path / 'prepared/complex.pdb'), 'interactions': interactions['atom_contacts'],
               'interface': sorted(interface), 'pse': str(out / 'complex.pse'), 'png': str(out / 'complex.png'),
               'validation': str(out / 'session_validation.json')}
    write_json(out / 'scene.json', payload)
    # A reusable PML accompanies the editable PSE.
    (out / 'reopen.pml').write_text('load ' + json.dumps(payload['pse']) + '\n', encoding='utf-8')
    run.command([pymol_python, ROOT / 'pymol_worker.py', out / 'scene.json'], timeout=600)
    validation = json.loads((out / 'session_validation.json').read_text(encoding='utf-8'))
    if validation.get('status') != 'SUCCESS':
        raise ValueError('PyMOL saved-session validation failed')
    from PIL import Image
    with Image.open(out / 'complex.png') as image:
        image.verify()


def analyze_prepared(run, pymol_python=None, ph=7.4, temperature=25):
    interactions = sasa = None
    with run.stage('structural_sanity'):
        check = geometry(gemmi.read_structure(str(run.path / 'prepared/complex.pdb')))
        write_json(run.path / 'analysis/geometry.json', check)
        if check['severe_clash_atom_count_A_1_5A']:
            raise ValueError('Severe inter-protein heavy-atom overlap (<1.5 A); inspect/refine before interpreting affinity')
    with run.stage('pdb_tools'):
        validate_pdb(run, run.path / 'prepared/complex.pdb')
    with run.stage('arpeggio'):
        interactions = run_arpeggio(run, ph=ph)
        if not interactions['atom_contacts']:
            raise ValueError('No inter-protein Arpeggio contacts in the proposed interface')
    with run.stage('freesasa'):
        sasa = run_sasa(run)
        if sasa['buried_total_A2'] <= 1:
            raise ValueError('No buried interface')
    with run.stage('prodigy'):
        if not success(run, 'structural_sanity'):
            raise Blocked('PRODIGY interpretation blocked by severe structural clashes; no meaningful affinity asserted')
        if not success(run, 'arpeggio'):
            raise Blocked('Cannot interpret affinity without a verified contact-bearing interface')
        run_prodigy(run, temperature)
    with run.stage('visualization'):
        visualize(run, interactions, sasa, pymol_python)
    return interactions, sasa


def negative_control(run, st):
    folder = run.path / 'validation/separated_control'
    folder.mkdir()
    separated = st.clone()
    for c in separated[0]:
        if c.name == 'B':
            for residue in c:
                for atom in residue:
                    atom.pos.x += 200
    for suffix in ('pdb', 'cif'):
        write_structure(separated, folder / f'complex.{suffix}')
    for c in separated[0]:
        write_structure(make_structure([c]), folder / f'protein_{c.name}.pdb')
    write_json(folder / 'geometry.json', geometry(separated))
    with run.stage('negative_arpeggio'):
        contacts = run_arpeggio(run, folder, label='separated_arpeggio')
        if contacts['atom_contacts'] or contacts['non_atom_contacts']:
            raise ValueError('Separated negative control has inter-protein interactions')
    with run.stage('negative_sasa'):
        result = run_sasa(run, folder, label='separated_sasa')
        if abs(result['buried_total_A2']) > 1:
            raise ValueError('Separated negative control has nonzero burial')
    with run.stage('negative_dockq'):
        result = run_dockq(run, folder / 'complex.pdb', run.path / 'prepared/complex.pdb', label='separated_control')
        if result['DockQ'] >= 0.23:
            raise ValueError('Separated negative control passed DockQ')


def benchmark(output=None, pymol_python=None, seeds=(42, 2024, 8888), swarms=400, glowworms=300, steps=100, cores=1, scoring='fastdfire'):
    parameters = {'benchmark': '1BRS', 'chain_a': 'A', 'chain_b': 'D', 'assembly': '1', 'seeds': list(seeds),
                  'swarms': swarms, 'glowworms': glowworms, 'steps': steps, 'cores': cores, 'scoring': scoring,
                  'scope': 'experimental interface analysis and bound redocking, no native restraints'}
    run = Run(output, parameters)
    write_json(run.path / 'environment.json', inspect_environment())
    client = RCSB(run.path / 'raw')
    st = None
    interactions = sasa = None
    try:
        with run.stage('acquisition'):
            path, meta = client.entry('1BRS')
            write_json(run.path / 'raw/metadata.json', meta)
        with run.stage('preparation'):
            if not success(run, 'acquisition'):
                raise Blocked('1BRS acquisition failed')
            report = inspect_structure(path, meta)
            write_json(run.path / 'prepared/inspection.json', report)
            a, b = entity_for_chain(meta, 'A'), entity_for_chain(meta, 'D')
            if 'BARNASE' not in a['rcsb_polymer_entity']['pdbx_description'].upper() or 'BARSTAR' not in b['rcsb_polymer_entity']['pdbx_description'].upper():
                raise ValueError('1BRS identities do not match benchmark')
            st = prepare_pair(run, path, path, 'A', 'D', '1')
            write_json(run.path / 'route.json', {'mode': 'experimental_interface', 'redocking': 'benchmark-only exception', 'assembly_id': '1'})
        if success(run, 'preparation'):
            interactions, sasa = analyze_prepared(run, pymol_python)
            with run.stage('chemistry_validation'):
                chemistry_check(run, interactions)
            with run.stage('sasa_invariance'):
                sasa_invariance(run, st, sasa)
            with run.stage('dockq_positive_control'):
                result = run_dockq(run, run.path / 'prepared/complex.pdb', run.path / 'prepared/complex.pdb', 'self_comparison')
                if result['DockQ'] < 0.99:
                    raise ValueError('DockQ identical-structure positive control failed')
            negative_control(run, st)
            with run.stage('redocking'):
                rows = dock(run, st, seeds, swarms, glowworms, steps, cores, reference=run.path / 'prepared/complex.pdb', scoring=scoring)
            with run.stage('redocked_top1_analysis'):
                if not success(run, 'redocking'):
                    raise Blocked('LightDock did not complete and validate candidates')
                # Preserve a separate complete analysis run; never overwrite the
                # experimental preparation with the predicted top-ranked model.
                child = Run(run.path / 'docking/top1_analysis', {'source': 'LightDock score top1'})
                try:
                    model = gemmi.read_structure(str(run.path / 'docking/score_top1.pdb'))
                    for suffix in ('pdb', 'cif'):
                        write_structure(model, child.path / f'prepared/complex.{suffix}')
                    for c in model[0]:
                        write_structure(make_structure([c]), child.path / f'prepared/protein_{c.name}.pdb')
                    analyze_prepared(child, pymol_python)
                    if not all(success(child, x) for x in ('structural_sanity', 'arpeggio', 'freesasa', 'prodigy', 'visualization')):
                        raise Blocked('Predicted top1 interface analysis incomplete; inspect its run.json')
                finally:
                    child.finish()
        required_a = ('preparation', 'structural_sanity', 'pdb_tools', 'arpeggio', 'freesasa', 'prodigy', 'visualization', 'chemistry_validation', 'sasa_invariance', 'dockq_positive_control',
                      'negative_arpeggio', 'negative_sasa', 'negative_dockq')
        required_b = ('redocking', 'redocked_top1_analysis')
        run.state['scientific_validation'] = {
            'status': 'PASS' if all(success(run, x) for x in required_a + required_b) else 'INCOMPLETE',
            'phase1': 'PASS' if all(success(run, x) for x in required_a) else 'INCOMPLETE',
            'phase2': 'NOT_VALIDATED',
            'limitations': ['bound redocking only', 'protein-only interface excludes crystal waters',
                            'PRODIGY is an empirical prediction', '1BRS paper describes Barstar C40A/C82A; archive mutation annotation alone is insufficient']}
        run.state['scientific_validation']['sampling_tier'] = 'benchmark' if swarms >= 50 and glowworms >= 50 and steps >= 100 and len(seeds) >= 3 else 'smoke_only'
        metrics_file = run.path / 'validation/redocking.json'
        if metrics_file.exists():
            metrics = json.loads(metrics_file.read_text())
            run.state['scientific_validation']['phase2'] = {'sampling': metrics['sampling_validation'], 'ranking': metrics['ranking_validation'],
                'top1_analysis': run.state['stages'].get('redocked_top1_analysis', {}).get('status', 'NOT_VALIDATED')}
            if metrics['sampling_validation'] != 'PASS' or metrics['ranking_validation'] != 'PASS':
                run.state['scientific_validation']['status'] = 'FAIL'
        if run.state['scientific_validation']['sampling_tier'] == 'smoke_only' and run.state['scientific_validation']['status'] == 'PASS':
            run.state['scientific_validation']['status'] = 'SMOKE_PASS_NOT_FULL_VALIDATION'
    finally:
        run.finish()
    return run


def acquire_input(client, value, label, species=None):
    path = Path(value)
    if path.is_file():
        target = client.directory / (label + '_' + path.name)
        shutil.copy2(path, target)
        write_json(client.directory / f'{label}_source.json', {'source': str(path.resolve()), 'sha256': sha256(path), 'kind': 'local experimental structure'})
        # mmCIF entity/sequence metadata are inspected separately; local files
        # are never assumed to be an RCSB entry solely from their file name.
        return target, {}
    if re.fullmatch(r'[0-9][A-Za-z0-9]{3}', value):
        return client.entry(value)
    return client.resolve_name(value, species)


def select_chain(report, metadata, supplied=None, species=None):
    selection = metadata.get('easy_ppi_name_selection')
    allowed = selection['chains'] if selection else [c['chain'] for c in report['chains']]
    if supplied:
        base = supplied.split('@')[0]
        if base not in allowed:
            raise Blocked(f'Chain {supplied!r} does not match the resolved protein identity; allowed {allowed}')
        candidates = [c for c in report['chains'] if c['chain'] == base]
    else:
        entity_ids = {e['rcsb_polymer_entity_container_identifiers']['entity_id'] for e in metadata.get('entities', [])
                      if set(e['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']) & set(allowed)}
        if not selection and (len(entity_ids) > 1 or not entity_ids and len(allowed) > 1):
            raise Blocked('Multiple protein identities: specify a chain after inspecting input_inspection.json')
        candidates = [c for c in report['chains'] if c['chain'] in allowed]
    if species:
        candidates = [c for c in candidates if any(str(s.get('ncbi_taxonomy_id')) == species or s.get('ncbi_scientific_name', '').casefold() == species.casefold() for s in c.get('species') or [])]
    if not candidates:
        raise Blocked('No chain matches the selected identity/species')
    candidates.sort(key=lambda c: (len(c['missing_backbone']), len(c['missing_sidechains']), -c['residues'], c['chain']))
    return supplied or candidates[0]['chain']


def local_entity(report, chain):
    choices = [x for x in report['chains'] if x['chain'] == chain]
    if len(choices) != 1:
        raise Blocked(f'Chain {chain!r} is ambiguous or absent')
    return {'entity_poly': {'pdbx_seq_one_letter_code_can': choices[0].get('full_sequence') or choices[0]['sequence']},
            'rcsb_entity_source_organism': choices[0].get('species') or [], 'rcsb_polymer_entity_container_identifiers': {}}


def same_entity(a, b, require_construct=True):
    def accessions(e):
        return {x['database_accession'] for x in e.get('rcsb_polymer_entity_container_identifiers', {}).get('reference_sequence_identifiers', []) if x.get('database_name') == 'UniProt'}
    aa, ab = accessions(a), accessions(b)
    sa = a.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '').replace('\n', '')
    sb = b.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '').replace('\n', '')
    if aa and ab:
        identity = aa == ab
        # Shared UniProt identity does not make differing mutants equivalent.
        # Exact terminal construct extensions are tolerated; substitutions are
        # kept for review rather than silently replacing the requested protein.
        if require_construct and sa and sb and sa not in sb and sb not in sa:
            identity = False
    else:
        identity = bool(sa and sb and sa == sb)
    taxa_a = {str(x['ncbi_taxonomy_id']) for x in a.get('rcsb_entity_source_organism', []) if x.get('ncbi_taxonomy_id') is not None}
    taxa_b = {str(x['ncbi_taxonomy_id']) for x in b.get('rcsb_entity_source_organism', []) if x.get('ncbi_taxonomy_id') is not None}
    species_match = (not taxa_a or not taxa_b or bool(taxa_a & taxa_b)) if require_construct else bool(taxa_a and taxa_a == taxa_b)
    return identity and species_match


def automatic_route(client, pa, pb, ma, mb, ra, rb, ca, cb, assembly_id=None):
    if Path(pa).resolve() == Path(pb).resolve():
        assembly_pair(pa, ca, cb, assembly_id)
        return pa, pb, ca, cb, assembly_id, {'mode': 'experimental_interface', 'reason': 'both proteins in one verified biological assembly'}
    ea = entity_for_chain(ma, ca) if ma.get('entities') else local_entity(ra, ca)
    eb = entity_for_chain(mb, cb) if mb.get('entities') else local_entity(rb, cb)
    common = client.common_entries(ea, eb)
    candidates, review = [], []
    name_a = bool(ma.get('easy_ppi_name_selection'))
    name_b = bool(mb.get('easy_ppi_name_selection'))
    for pdb_id in common:
        pc, mc = client.entry(pdb_id)
        # A name does not specify a sequence construct. Select the actual
        # experimental complex for that identity and expose construct changes.
        # Explicit PDB/file inputs retain the strict construct check.
        valid_a = [e for e in mc['entities'] if same_entity(ea, e, require_construct=not name_a)]
        valid_b = [e for e in mc['entities'] if same_entity(eb, e, require_construct=not name_b)]
        if not valid_a or not valid_b:
            review.append({'pdb_id': pdb_id, 'reason': 'identity/species or explicitly requested construct mismatch'})
        for a in valid_a:
            for b in valid_b:
                for ach in a['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']:
                    for bch in b['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']:
                        if ach == bch:
                            continue
                        try:
                            _, _, evidence = assembly_pair(pc, ach, bch)
                            resolution = mc['entry'].get('rcsb_entry_info', {}).get('resolution_combined', [999])[0]
                            candidates.append((resolution, pdb_id, ach, bch, pc, evidence['assembly_id']))
                        except Blocked as exc:
                            review.append({'pdb_id': pdb_id, 'chains': [ach, bch], 'reason': str(exc)})
    discovery = {'common_experimental_entries': common, 'rejected_or_ambiguous': review,
                 'identity_policy': 'exact UniProt identity and source species; explicit structure inputs additionally require compatible constructs',
                 'name_only_construct_selection': {'A': name_a, 'B': name_b}}
    write_json(client.directory.parent / 'route_candidates.json', discovery)
    if candidates:
        candidates.sort(key=lambda x: x[:4])
        _, pdb_id, ca, cb, pc, aid = candidates[0]
        discovery.update(mode='experimental_interface', selected_pdb_id=pdb_id, selected_chains=[ca, cb], assembly_id=aid)
        _, selected_metadata = client.entry(pdb_id)
        changes = {}
        for label, original, chain, name_only in [('A', ea, ca, name_a), ('B', eb, cb, name_b)]:
            selected = entity_for_chain(selected_metadata, chain)
            original_sequence = original.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
            selected_sequence = selected.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
            changes[label] = {'name_only_input': name_only, 'same_full_sequence_as_initial_structure_candidate': original_sequence == selected_sequence,
                              'selected_archive_mutation_count': selected.get('entity_poly', {}).get('rcsb_mutation_count'),
                              'note': 'The actual experimental construct is analyzed; name input did not specify wild type or a particular construct' if name_only else 'Explicit structure construct matching enforced'}
        discovery['construct_review'] = changes
        return pc, pc, ca, cb, aid, discovery
    if common:
        # Entries may contain transformed copies or incomplete relevant chains.
        # Do not silently claim no complex simply because preparation failed.
        raise Blocked('Common experimental entries exist but no unambiguous prepared biological interface; inspect source/assembly evidence')
    discovery.update(mode='lightdock', reason='no common exact-identity experimental entry found; absence of evidence is not evidence of no interaction')
    return pa, pb, ca, cb, None, discovery


def prepare_reference(run, value, chain_a, chain_b):
    if not chain_a or not chain_b:
        raise Blocked('Experimental reference requires explicit --reference-chain-a and --reference-chain-b')
    child = Run(run.path / 'reference', {'reference_input': value, 'chain_a': chain_a, 'chain_b': chain_b})
    result = None
    try:
        client = RCSB(child.path / 'raw')
        with child.stage('experimental_reference'):
            path, metadata = acquire_input(client, value, 'reference')
            report = inspect_structure(path, metadata)
            write_json(child.path / 'prepared/inspection.json', report)
            prepare_pair(child, path, path, chain_a, chain_b, require_assembly=True)
            result = child.path / 'prepared/complex.pdb'
        if not success(child, 'experimental_reference'):
            raise Blocked('Reference is not a verified biological experimental complex; inspect reference/run.json')
    finally:
        child.finish()
    return result


def workflow(a, b, chain_a, chain_b, output=None, assembly_id=None, pymol_python=None, seeds=(42, 2024, 8888), swarms=400, glowworms=300, steps=100, cores=1, species_a=None, species_b=None,
             reference=None, reference_chain_a=None, reference_chain_b=None, scoring='fastdfire'):
    run = Run(output, {'protein_a': a, 'protein_b': b, 'chain_a': chain_a, 'chain_b': chain_b, 'assembly_id': assembly_id,
                       'seeds': list(seeds), 'swarms': swarms, 'glowworms': glowworms, 'steps': steps, 'cores': cores,
                       'reference': reference, 'reference_chain_a': reference_chain_a, 'reference_chain_b': reference_chain_b, 'scoring': scoring})
    write_json(run.path / 'environment.json', inspect_environment())
    try:
        client = RCSB(run.path / 'raw')
        with run.stage('acquisition'):
            pa, ma = acquire_input(client, a, 'a', species_a)
            if a == b:
                pb, mb = pa, ma
            else:
                pb, mb = acquire_input(client, b, 'b', species_b)
        with run.stage('inspection'):
            if not success(run, 'acquisition'):
                raise Blocked('Input acquisition/identity unresolved')
            ra, rb = inspect_structure(pa, ma), inspect_structure(pb, mb)
            write_json(run.path / 'prepared/input_inspection.json', {'A': ra, 'B': rb})
            chain_a = select_chain(ra, ma, chain_a, species_a)
            chain_b = select_chain(rb, mb, chain_b, species_b)
            write_json(run.path / 'prepared/chain_selection.json', {'A': chain_a, 'B': chain_b, 'species_a': species_a, 'species_b': species_b})
        with run.stage('route'):
            if not success(run, 'inspection'):
                raise Blocked('Inputs not inspected')
            pa, pb, ca, cb, aid, route = automatic_route(client, pa, pb, ma, mb, ra, rb, chain_a, chain_b, assembly_id)
            write_json(run.path / 'route.json', route)
        with run.stage('preparation'):
            if not success(run, 'route'):
                raise Blocked('Experimental-complex search unresolved')
            st = prepare_pair(run, pa, pb, ca, cb, aid, require_assembly=route['mode'] == 'experimental_interface')
        if success(run, 'preparation'):
            reference_path = None
            if reference:
                with run.stage('reference'):
                    reference_path = prepare_reference(run, reference, reference_chain_a, reference_chain_b)
            if route['mode'] == 'lightdock':
                with run.stage('lightdock'):
                    if reference and not success(run, 'reference'):
                        raise Blocked('Supplied experimental reference validation failed')
                    rows = dock(run, st, seeds, swarms, glowworms, steps, cores, reference_path, scoring=scoring)
                    predicted = gemmi.read_structure(rows[0]['path'])
                    # Original prepared monomers remain in docking input copies;
                    # predicted complex is kept in a separate analysis run.
                if success(run, 'lightdock'):
                    child = Run(run.path / 'docking/top1_analysis', {'source': 'LightDock score top1'})
                    try:
                        for suffix in ('pdb', 'cif'):
                            write_structure(predicted, child.path / f'prepared/complex.{suffix}')
                        for c in predicted[0]:
                            write_structure(make_structure([c]), child.path / f'prepared/protein_{c.name}.pdb')
                        analyze_prepared(child, pymol_python)
                    finally:
                        child.finish()
            else:
                analyze_prepared(run, pymol_python)
                with run.stage('dockq'):
                    if reference and not success(run, 'reference'):
                        raise Blocked('Supplied reference invalid; no DockQ computed')
                    run_dockq(run, run.path / 'prepared/complex.pdb', reference_path)
        run.state['scientific_validation'] = {'status': 'NOT_BENCHMARKED', 'note': 'Inspect per-stage statuses and benchmark evidence; no experimental affinity asserted'}
    finally:
        run.finish()
    return run
