"""Audited official RCSB REST acquisition. No predicted structure services."""
from __future__ import annotations

import json
from pathlib import Path
import re
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .audit import Blocked, sha256, write_json, now

DATA = 'https://data.rcsb.org/rest/v1/core'
SEARCH = 'https://search.rcsb.org/rcsbsearch/v2/query'


class RCSB:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.mount('https://', HTTPAdapter(max_retries=Retry(total=3, backoff_factor=1, status_forcelist=[429, 502, 503, 504])))
        self.records = []

    def get(self, url, filename, payload=None):
        path = self.directory / filename
        if path.exists():
            return path
        response = self.session.get(url, timeout=(15, 90)) if payload is None else self.session.post(url, json=payload, timeout=(15, 90))
        response.raise_for_status()
        with path.open('xb') as handle:
            handle.write(response.content)
        self.records.append({'url': url, 'query': payload, 'file': filename, 'time': now(), 'sha256': sha256(path)})
        write_json(self.directory / 'acquisition.json', self.records)
        return path

    def json(self, url, filename, payload=None):
        return json.loads(self.get(url, filename, payload).read_text(encoding='utf-8'))

    def entry(self, pdb_id):
        pdb_id = pdb_id.upper()
        if not re.fullmatch(r'[0-9][A-Z0-9]{3}', pdb_id):
            raise ValueError('Expected a four-character experimental PDB ID')
        entry = self.json(f'{DATA}/entry/{pdb_id}', f'{pdb_id}.entry.json')
        if not entry.get('exptl'):
            raise Blocked('RCSB entry has no experimental method')
        entities = []
        for eid in entry.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', []):
            entities.append(self.json(f'{DATA}/polymer_entity/{pdb_id}/{eid}', f'{pdb_id}.entity{eid}.json'))
        assemblies = []
        for aid in entry.get('rcsb_entry_container_identifiers', {}).get('assembly_ids', []):
            assemblies.append(self.json(f'{DATA}/assembly/{pdb_id}/{aid}', f'{pdb_id}.assembly{aid}.json'))
        path = self.get(f'https://files.rcsb.org/download/{pdb_id}.cif', f'{pdb_id}.cif')
        return path, {'entry': entry, 'entities': entities, 'assemblies': assemblies}

    def query(self, query, name, return_type='entry', rows=100):
        payload = {'query': query, 'return_type': return_type,
                   'request_options': {'paginate': {'start': 0, 'rows': rows}, 'results_content_type': ['experimental']}}
        # RCSB returns HTTP 204 for an empty result set.
        response = self.session.post(SEARCH, json=payload, timeout=(15, 90))
        response.raise_for_status()
        value = {} if response.status_code == 204 else response.json()
        write_json(self.directory / f'{name}.json', {'url': SEARCH, 'query': payload, 'response': value, 'time': now()})
        if value.get('total_count', 0) > rows:
            raise Blocked(f'RCSB result truncated ({value["total_count"]} > {rows}); narrow identity/species query')
        return [x['identifier'] for x in value.get('result_set', [])]

    def names(self, name, species=None):
        # Resolve a gene symbol first; never assume human for bacterial benchmarks.
        mygene_url = 'https://mygene.info/v3/query'
        params = {'q': name, 'fields': 'symbol,name,taxid', 'size': 10}
        # MyGene accepts supported common names/taxids, but does not map every
        # bacterial scientific name. Query all species there, then constrain
        # RCSB by the authoritative requested source taxonomy/name below.
        params['species'] = species if species and str(species).isdigit() else 'all'
        response = self.session.get(mygene_url, params=params, timeout=45)
        response.raise_for_status()
        gene = response.json()
        write_json(self.directory / ('gene_' + re.sub(r'\W', '_', name) + '.json'), {'url': response.url, 'response': gene})
        exact = {x['symbol'] for x in gene.get('hits', []) if x.get('symbol', '').casefold() == name.casefold()}
        symbol = next(iter(exact)) if len(exact) == 1 else name
        # A candidate list is reviewable, but full-text results alone are never identity proof.
        terms = [{'type': 'terminal', 'service': 'full_text', 'parameters': {'value': symbol}}]
        if species:
            if str(species).isdigit():
                terms.append({'type': 'terminal', 'service': 'text', 'parameters': {'attribute': 'rcsb_entity_source_organism.ncbi_taxonomy_id', 'operator': 'equals', 'value': int(species)}})
            else:
                terms.append(text_query('rcsb_entity_source_organism.ncbi_scientific_name', species))
        query = terms[0] if len(terms) == 1 else {'type': 'group', 'logical_operator': 'and', 'nodes': terms}
        return {'input': name, 'query_symbol': symbol, 'candidates': self.query(query, 'name_' + re.sub(r'\W', '_', name), rows=1000),
                'requires_identity_selection': True}

    def resolve_name(self, name, species=None):
        search = self.names(name, species)
        ids = search['candidates']
        if not ids:
            raise Blocked(f'No experimental RCSB candidates for {name!r}')
        query = '''query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id rcsb_entry_info{resolution_combined}
            polymer_entities{rcsb_id rcsb_polymer_entity{pdbx_description}
            entity_poly{pdbx_seq_one_letter_code_can}
            rcsb_polymer_entity_container_identifiers{auth_asym_ids reference_sequence_identifiers{database_name database_accession}}
            rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id rcsb_gene_name{value}}}}}'''
        payload = {'query': query, 'variables': {'ids': ids}}
        data = self.json('https://data.rcsb.org/graphql', 'name_entities_' + re.sub(r'\W', '_', name) + '.json', payload)
        if data.get('errors'):
            raise Blocked(f'RCSB identity metadata query failed: {data["errors"]}; candidates preserved')
        matches = []
        rejected_fusions = []
        identities = set()
        normalize = lambda x: re.sub(r'[^a-z0-9]', '', x.lower())
        requested = {normalize(name), normalize(search['query_symbol'])}
        for entry in data.get('data', {}).get('entries', []):
            if not entry:
                continue
            for entity in entry.get('polymer_entities') or []:
                description = entity.get('rcsb_polymer_entity', {}).get('pdbx_description', '')
                organism = entity.get('rcsb_entity_source_organism') or []
                gene_names = {normalize(g['value']) for org in organism for g in org.get('rcsb_gene_name') or []}
                if normalize(description) not in requested and not requested.intersection(gene_names):
                    continue
                if species and not any(str(org.get('ncbi_taxonomy_id')) == species or org.get('ncbi_scientific_name', '').casefold() == species.casefold() for org in organism):
                    continue
                refs = entity.get('rcsb_polymer_entity_container_identifiers', {}).get('reference_sequence_identifiers') or []
                accessions = tuple(sorted({x['database_accession'] for x in refs if x['database_name'] == 'UniProt'}))
                taxa = tuple(sorted({str(x['ncbi_taxonomy_id']) for x in organism if x.get('ncbi_taxonomy_id') is not None}))
                if len(accessions) > 1 or len(taxa) > 1:
                    rejected_fusions.append({'entity': entity['rcsb_id'], 'description': description, 'accessions': accessions, 'taxa': taxa,
                                             'reason': 'multi-reference or multi-organism construct; requires explicit selection'})
                    continue
                if not taxa:
                    continue
                sequence = entity.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '')
                identity = (accessions or (sequence,), taxa)
                identities.add(identity)
                resolution = entry.get('rcsb_entry_info', {}).get('resolution_combined') or [999]
                matches.append({'pdb_id': entry['rcsb_id'], 'entity_id': entity['rcsb_id'].split('_')[-1],
                                'description': description, 'resolution_A': resolution[0], 'taxa': taxa, 'accessions': accessions,
                                'chains': entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']})
        matches.sort(key=lambda x: (x['resolution_A'], x['pdb_id'], x['entity_id']))
        write_json(self.directory / ('resolved_' + re.sub(r'\W', '_', name) + '.json'), {'search': search, 'exact_identity_candidates': matches,
                   'unique_identity_count': len(identities), 'rejected_fusions': rejected_fusions,
                   'policy': 'exact normalized entity description or exact gene; UniProt plus source taxonomy identity'})
        if len(identities) != 1 or not matches:
            raise Blocked(f'Name {name!r} has {len(identities)} exact identities; select species/PDB/entity from resolved candidate JSON')
        selected = matches[0]
        path, metadata = self.entry(selected['pdb_id'])
        metadata['easy_ppi_name_selection'] = selected
        return path, metadata

    def common_entries(self, entity_a, entity_b):
        nodes = []
        for entity in (entity_a, entity_b):
            ids = entity.get('rcsb_polymer_entity_container_identifiers', {})
            accessions = [x['database_accession'] for x in ids.get('reference_sequence_identifiers', []) if x.get('database_name') == 'UniProt']
            if not accessions:
                # Exact experimental sequence identity; partial constructs may need manual review.
                sequence = entity.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', '').replace('\n', '')
                if not sequence:
                    raise Blocked('Cannot establish protein identity for complex discovery')
                nodes.append({'type': 'terminal', 'service': 'sequence', 'parameters': {'evalue_cutoff': 0.1, 'identity_cutoff': 1.0, 'sequence_type': 'protein', 'value': sequence}})
            else:
                nodes.append({'type': 'terminal', 'service': 'text', 'parameters': {'attribute': 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession', 'operator': 'in', 'value': accessions}})
        # Entity-scoped predicates in a single AND can wrongly demand both IDs in
        # one entity. Intersect separately searched ENTRY sets instead.
        sets = [set(self.query(node, f'identity_{i}', rows=1000)) for i, node in enumerate(nodes)]
        return sorted(sets[0] & sets[1])


def text_query(attribute, value):
    return {'type': 'terminal', 'service': 'text', 'parameters': {'attribute': attribute, 'operator': 'exact_match', 'value': value}}


def entity_for_chain(metadata, chain):
    chain = chain.split('@')[0]
    matches = [e for e in metadata.get('entities', []) if chain in e.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids', [])]
    if len(matches) != 1:
        raise Blocked(f'Chain {chain!r} does not identify one polymer entity')
    return matches[0]
