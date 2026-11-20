from pathlib import Path
import json
import math
import tempfile
from unittest.mock import patch

import gemmi
import numpy as np
import pytest

from easy_ppi.audit import Blocked, Run
from easy_ppi.backends import arpeggio_rows, dockq_grade, parse_prodigy, run_dockq
from easy_ppi.docking import read_gso
from easy_ppi.pipeline import same_entity
from easy_ppi.structures import clean_chain, geometry, inspect_structure, make_structure, scramble_partners, write_structure


def small_structure():
    chains = []
    for name, offset in [('X', 0), ('Y', 3)]:
        c = gemmi.Chain(name)
        r = gemmi.Residue()
        r.name, r.seqid = 'ALA', gemmi.SeqId(10, 'A')
        for index, name in enumerate(['N', 'CA', 'C', 'O', 'CB']):
            a = gemmi.Atom()
            a.name, a.element = name, gemmi.Element('N' if name == 'N' else 'O' if name == 'O' else 'C')
            a.pos, a.occ = gemmi.Position(index * 1.3, offset, 0), 1
            r.add_atom(a)
        c.add_residue(r)
        chains.append(c)
    return make_structure(chains)


def test_refuses_existing_output(tmp_path):
    existing = tmp_path / 'occupied'
    existing.mkdir()
    marker = existing / 'valuable.txt'
    marker.write_text('keep me')
    with pytest.raises(FileExistsError):
        Run(existing)
    assert marker.read_text() == 'keep me'


def test_stage_failure_preserves_success_and_manifest(tmp_path):
    run = Run(tmp_path / 'run')
    with run.stage('a'):
        pass
    with run.stage('b'):
        raise Blocked('dependency absent')
    run.finish()
    assert run.state['stages']['a']['status'] == 'SUCCESS'
    assert run.state['stages']['b']['status'] == 'BLOCKED'
    assert run.state['scientific_validation']['status'] == 'NOT_VALIDATED'
    assert json.loads((run.path / 'manifest.json').read_text())


def test_no_reference_never_invokes_dockq(tmp_path):
    run = Run(tmp_path / 'run')
    with patch.object(run, 'command', side_effect=AssertionError('must not run')):
        data = run_dockq(run, 'model.pdb', None)
    assert data['DockQ'] is None
    assert data['status'] == 'SKIPPED'
    assert run.state['stages']['dockq']['status'] == 'SKIPPED'


def test_arpeggio_filters_chains_and_preserves_labels():
    def atom(c, n):
        return {'auth_asym_id': c, 'auth_seq_id': n, 'auth_atom_id': 'O', 'pdbx_PDB_ins_code': 'A' if n == 10 else '?'}
    raw = [{'bgn': atom('A', 10), 'end': atom('B', 35), 'distance': 2.9, 'contact': ['hbond', 'polar']},
           {'bgn': atom('A', 10), 'end': atom('A', 35), 'distance': 2.7, 'contact': ['hbond']}]
    rows, other = arpeggio_rows(raw)
    assert len(rows) == 1 and not other
    assert rows[0]['residue_a'] == '10A'
    assert rows[0]['types'] == ['hbond', 'polar']


def test_bad_arpeggio_schema_and_distance_fail():
    with pytest.raises(ValueError):
        arpeggio_rows({'contacts': []})
    with pytest.raises(ValueError):
        arpeggio_rows([{'bgn': {'auth_asym_id': 'A', 'auth_atom_id': 'O', 'auth_seq_id': 1},
                       'end': {'auth_asym_id': 'B', 'auth_atom_id': 'N', 'auth_seq_id': 2}, 'contact': ['hbond'], 'distance': float('nan')}])


def test_dockq_published_thresholds():
    assert [dockq_grade(x) for x in [0, .23, .49, .8, 1]] == ['incorrect', 'acceptable', 'medium', 'high', 'high']


def test_prodigy_requires_real_prediction():
    data = parse_prodigy('[++] Predicted binding affinity (kcal.mol-1): -12.5\n[++] Predicted dissociation constant (M) at 25.0 C: 6.8e-10')
    assert data['delta_g_kcal_mol'] == -12.5
    assert 'empirical' in data['interpretation']
    with pytest.raises(ValueError):
        parse_prodigy('program completed without errors')


def test_residue_mapping_preserves_insertion_code():
    chain, mapping, _ = clean_chain(small_structure(), 'X', 'A')
    assert mapping[0]['prepared_residue'] == '10A'
    assert mapping[0]['source_ins_code'] == 'A'
    assert chain[0].seqid.num == 10 and chain[0].seqid.icode == 'A'


def test_altloc_chooses_one_residue_conformer():
    st = small_structure()
    r = st[0]['X'][0]
    del r[4]
    for alt, occupancy, x in [('A', .6, 1), ('B', .4, 20)]:
        atom = gemmi.Atom()
        atom.name, atom.element, atom.altloc, atom.occ = 'CB', gemmi.Element('C'), alt, occupancy
        atom.pos = gemmi.Position(x, 0, 0)
        r.add_atom(atom)
    chain, mapping, decisions = clean_chain(st, 'X', 'A')
    assert len(chain[0]) == 5
    assert chain[0][-1].pos.x == 1
    assert mapping[0]['selected_altloc'] == 'A'


def test_missing_backbone_and_wrong_chain_are_blocked():
    st = small_structure()
    del st[0]['X'][0][0]
    with pytest.raises(Blocked, match='backbone'):
        clean_chain(st, 'X', 'A')
    with pytest.raises(Blocked, match='Choose one'):
        clean_chain(st, 'Z', 'A')


def test_predicted_local_file_is_rejected(tmp_path):
    path = tmp_path / 'prediction.cif'
    write_structure(small_structure(), path)
    with pytest.raises(Blocked, match='Experimental provenance'):
        inspect_structure(path)


def test_local_pdb_experimental_method_uses_gemmi_map_api(tmp_path):
    path = tmp_path / 'experimental.pdb'
    write_structure(small_structure(), path)
    path.write_text('EXPDTA    X-RAY DIFFRACTION\n' + path.read_text(), encoding='utf-8')
    report = inspect_structure(path)
    assert 'X-RAY DIFFRACTION' in report['experimental_methods']


def test_scramble_removes_original_relative_transform():
    st = small_structure()
    a, _, _ = clean_chain(st, 'X', 'A')
    b, _, _ = clean_chain(st, 'Y', 'B')
    pair = make_structure([a, b])
    before = geometry(pair)
    scrambled, transforms = scramble_partners(pair, 42)
    assert transforms['B']['rotation_matrix'] != np.eye(3).tolist()
    assert transforms['A']['subtract_centroid'] != transforms['B']['subtract_centroid']
    assert geometry(scrambled) != before
    scrambled2, transforms2 = scramble_partners(pair, 42)
    assert transforms == transforms2


def test_gso_parser_does_not_use_luciferin_as_score(tmp_path):
    path = tmp_path / 'gso.out'
    path.write_text('# header\n(0,0,0,1,0,0,0) 0 0 999.0 4 2.1 -12.5\n')
    assert read_gso(path) == [{'glowworm': 0, 'score': -12.5}]


def test_species_mismatch_cannot_be_identity_match():
    a = {'entity_poly': {'pdbx_seq_one_letter_code_can': 'ALAT'}, 'rcsb_entity_source_organism': [{'ncbi_taxonomy_id': 9606}]}
    b = {'entity_poly': {'pdbx_seq_one_letter_code_can': 'ALAT'}, 'rcsb_entity_source_organism': [{'ncbi_taxonomy_id': 10090}]}
    assert not same_entity(a, b)


def test_prepare_1brs_assembly_when_fixture_present():
    from easy_ppi.structures import assembly_pair
    fixture = Path(__file__).resolve().parents[2] / 'ppi_results/environment_check_20261118/rcsb/1BRS.cif'
    if not fixture.exists():
        migration_index = Path(__file__).resolve().parents[2] / 'references/result-location-map.json'
        if migration_index.exists():
            records = json.loads(migration_index.read_text(encoding='utf-8'))['entries']
            for record in records:
                if Path(record['source']).name == 'environment_check_20261118':
                    fixture = Path(record['destination']) / 'rcsb/1BRS.cif'
                    break
    if not fixture.exists():
        pytest.skip('Downloaded experimental fixture not present; run benchmark')
    st, mapping, info = assembly_pair(fixture, 'A', 'D', '1')
    assert info['assembly_id'] == '1'
    assert geometry(st)['severe_clash_atom_count_A_1_5A'] == 0
    assert len(mapping) == 195  # 108 Barnase + 87 Barstar coordinate residues
    with pytest.raises(Blocked):
        assembly_pair(fixture, 'A', 'E')
