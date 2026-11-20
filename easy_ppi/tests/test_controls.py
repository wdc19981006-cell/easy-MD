import json
from pathlib import Path
from unittest.mock import patch

import gemmi
import pytest

from easy_ppi.audit import Blocked, Run
from easy_ppi.pipeline import same_entity, select_chain
from easy_ppi.validation import validate_artifacts


def test_requested_species_is_enforced():
    report = {'chains': [{'chain': 'A', 'species': [{'ncbi_taxonomy_id': 1390, 'ncbi_scientific_name': 'Bacillus amyloliquefaciens'}],
                          'missing_backbone': [], 'missing_sidechains': [], 'residues': 110}]}
    assert select_chain(report, {}, supplied='A', species='1390') == 'A'
    with pytest.raises(Blocked):
        select_chain(report, {}, supplied='A', species='9606')


def test_identity_ambiguity_requires_explicit_chain():
    report = {'chains': [{'chain': x, 'missing_backbone': [], 'missing_sidechains': [], 'residues': 50} for x in ['A', 'B']]}
    with pytest.raises(Blocked, match='Multiple protein identities'):
        select_chain(report, {})


def test_shared_uniprot_does_not_hide_a_mutant():
    identifiers = {'reference_sequence_identifiers': [{'database_name': 'UniProt', 'database_accession': 'TEST'}]}
    a = {'entity_poly': {'pdbx_seq_one_letter_code_can': 'AAACCC'}, 'rcsb_polymer_entity_container_identifiers': identifiers}
    b = {'entity_poly': {'pdbx_seq_one_letter_code_can': 'AAAGCC'}, 'rcsb_polymer_entity_container_identifiers': identifiers}
    assert not same_entity(a, b)


def test_name_identity_can_use_actual_experimental_construct_but_not_a_fusion():
    identifiers = {'reference_sequence_identifiers': [{'database_name': 'UniProt', 'database_accession': 'TEST'}]}
    taxa = [{'ncbi_taxonomy_id': 1390}]
    a = {'entity_poly': {'pdbx_seq_one_letter_code_can': 'AAACCC'}, 'rcsb_polymer_entity_container_identifiers': identifiers,
         'rcsb_entity_source_organism': taxa}
    b = {'entity_poly': {'pdbx_seq_one_letter_code_can': 'AAAGCC'}, 'rcsb_polymer_entity_container_identifiers': identifiers,
         'rcsb_entity_source_organism': taxa}
    assert same_entity(a, b, require_construct=False)
    fused = dict(b, rcsb_polymer_entity_container_identifiers={'reference_sequence_identifiers': identifiers['reference_sequence_identifiers'] + [{'database_name': 'UniProt', 'database_accession': 'OTHER'}]})
    assert not same_entity(a, fused, require_construct=False)


def test_manifest_detects_changed_artifact(tmp_path):
    run = Run(tmp_path / 'run')
    artifact = run.path / 'prepared/test.pdb'
    artifact.write_text('original')
    run.finish()
    assert validate_artifacts(run.path)['artifact_integrity'] == 'PASS'
    artifact.write_text('tampered')
    assert validate_artifacts(run.path)['artifact_integrity'] == 'FAIL'


def test_failed_docking_cannot_become_success_from_a_manifest(tmp_path):
    run = Run(tmp_path / 'run')
    with run.stage('docking'):
        raise RuntimeError('native tool returned no candidates')
    run.finish()
    result = validate_artifacts(run.path)
    assert result['artifact_integrity'] == 'PASS'
    assert result['stages']['docking'] == 'FAILED'
    assert result['scientific_validation']['status'] == 'NOT_VALIDATED'


