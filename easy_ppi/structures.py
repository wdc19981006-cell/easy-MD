from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path
import shutil
import re
import sys

import gemmi
import numpy as np
from scipy.spatial import cKDTree

from .audit import Blocked, sha256, write_csv, write_json

STANDARD = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
SIDECHAIN_ATOMS = {k: set(v.split()) for k, v in {
    'ALA': 'CB', 'ARG': 'CB CG CD NE CZ NH1 NH2', 'ASN': 'CB CG OD1 ND2', 'ASP': 'CB CG OD1 OD2',
    'CYS': 'CB SG', 'GLN': 'CB CG CD OE1 NE2', 'GLU': 'CB CG CD OE1 OE2', 'GLY': '',
    'HIS': 'CB CG ND1 CD2 CE1 NE2', 'ILE': 'CB CG1 CG2 CD1', 'LEU': 'CB CG CD1 CD2',
    'LYS': 'CB CG CD CE NZ', 'MET': 'CB CG SD CE', 'PHE': 'CB CG CD1 CD2 CE1 CE2 CZ',
    'PRO': 'CB CG CD', 'SER': 'CB OG', 'THR': 'CB OG1 CG2', 'TRP': 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
    'TYR': 'CB CG CD1 CD2 CE1 CE2 CZ OH', 'VAL': 'CB CG1 CG2'}.items()}
EXPERIMENTAL = ('X-RAY DIFFRACTION', 'SOLUTION NMR', 'SOLID-STATE NMR', 'ELECTRON MICROSCOPY', 'NEUTRON DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'FIBER DIFFRACTION')


def protein_residues(chain):
    return [r for r in chain if gemmi.find_tabulated_residue(r.name).is_amino_acid()]


def inspect_structure(path, metadata=None):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    if not len(st):
        raise Blocked('No coordinate models')
    block = gemmi.cif.read(str(path)).sole_block() if str(path).lower().endswith(('.cif', '.mmcif')) else None
    pdb_header = None
    if block is None:
        from Bio.PDB import parse_pdb_header
        pdb_header = parse_pdb_header(str(path))
    methods = [x['method'] for x in (metadata or {}).get('entry', {}).get('exptl', [])]
    if not methods:
        if block:
            methods = [gemmi.cif.as_string(value) for value in block.find_values('_exptl.method')]
        else:
            try:
                method = st.info['_exptl.method']
            except KeyError:
                method = (pdb_header or {}).get('structure_method', '').upper()
            methods = [method]
    if not any(any(allowed in method.upper() for allowed in EXPERIMENTAL) for method in methods):
        raise Blocked('Experimental provenance missing/unsupported; predicted structures are excluded')
    chains = []
    for chain in st[0]:
        residues = protein_residues(chain)
        if not residues:
            continue
        missing_backbone = [{'residue': str(r.seqid), 'name': r.name, 'atoms': sorted({'N', 'CA', 'C', 'O'} - {a.name for a in r if a.occ > 0})} for r in residues]
        missing_backbone = [x for x in missing_backbone if x['atoms']]
        entities = [e for e in (metadata or {}).get('entities', []) if chain.name in e.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids', [])]
        species = [s for e in entities for s in e.get('rcsb_entity_source_organism', [])]
        expected = next((e.get('entity_poly', {}).get('rcsb_sample_sequence_length') for e in entities), None)
        if not species and block:
            # Source organism, never expression host. Retain unknown explicitly.
            entity_id = next((e.name for e in st.entities if residues[0].subchain in e.subchains), None)
            for prefix, field, taxfield in [('_entity_src_gen.', 'pdbx_gene_src_scientific_name', 'pdbx_gene_src_ncbi_taxonomy_id'),
                                            ('_entity_src_nat.', 'pdbx_organism_scientific', 'pdbx_ncbi_taxonomy_id')]:
                values = block.get_mmcif_category(prefix)
                for i, eid in enumerate(values.get('entity_id', [])):
                    if eid == entity_id and field in values:
                        record = {'ncbi_scientific_name': values[field][i], 'source': 'local mmCIF source category'}
                        if taxfield in values and values[taxfield][i]:
                            record['ncbi_taxonomy_id'] = values[taxfield][i]
                        species.append(record)
        if not species and pdb_header:
            for mol_id, comp in pdb_header.get('compound', {}).items():
                chain_ids = [x.strip().casefold() for x in comp.get('chain', '').split(',')]
                source = pdb_header.get('source', {}).get(mol_id, {})
                if chain.name.casefold() in chain_ids and source.get('organism_scientific'):
                    species.append({'ncbi_scientific_name': source['organism_scientific'], 'ncbi_taxonomy_id': source.get('organism_taxid'), 'source': 'local PDB SOURCE/COMPND'})
        full_sequence = None
        if expected is None:
            for entity in st.entities:
                if any(r.subchain in entity.subchains for r in residues) and entity.full_sequence:
                    expected = len(entity.full_sequence)
                    full_sequence = ''.join(gemmi.find_tabulated_residue(x).one_letter_code for x in entity.full_sequence)
                    break
        missing_sidechains = [{'residue': str(r.seqid), 'name': r.name, 'atoms': sorted(SIDECHAIN_ATOMS.get(r.name, set()) - {a.name for a in r if a.occ > 0})} for r in residues]
        missing_sidechains = [x for x in missing_sidechains if x['atoms']]
        chains.append({'chain': chain.name, 'residues': len(residues), 'sequence': ''.join(gemmi.find_tabulated_residue(r.name).one_letter_code for r in residues),
                       'full_sequence': full_sequence, 'species': species or None, 'species_status': 'known' if species else 'unknown',
                       'expected_sequence_length': expected, 'unobserved_count': max(0, expected-len(residues)) if expected is not None else None,
                       'missing_backbone': missing_backbone, 'missing_sidechains': missing_sidechains,
                       'observed_fraction': len(residues) / expected if expected else None,
                       'nonstandard_residues': [{'residue': str(r.seqid), 'name': r.name} for r in residues if r.name not in STANDARD],
                       'altloc_atoms': sum(a.altloc != '\x00' for r in residues for a in r),
                       'zero_occupancy_atoms': sum(a.occ <= 0 for r in residues for a in r),
                       'partial_occupancy_atoms': sum(0 < a.occ < 0.999 for r in residues for a in r),
                       'mean_B_iso_A2': float(np.mean([a.b_iso for r in residues for a in r])),
                       'archive_mutation_count': next((e.get('entity_poly', {}).get('rcsb_mutation_count') for e in entities), None)})
    missing = []
    if block:
        table = block.find('_pdbx_unobs_or_zero_occ_residues.', ['auth_asym_id', 'auth_seq_id', 'auth_comp_id', 'PDB_ins_code'])
        missing = [{'chain': row[0], 'residue': row[1], 'name': row[2], 'insertion_code': row[3]} for row in table]
    else:
        missing = [line.rstrip() for line in Path(path).read_text(errors='replace').splitlines() if line.startswith('REMARK 465')]
    quality = (metadata or {}).get('entry', {}).get('rcsb_entry_info', {})
    refinement = (metadata or {}).get('entry', {}).get('refine', [])
    if not refinement and block:
        ref = block.get_mmcif_category('_refine.')
        refinement = [{k: values[i] for k, values in ref.items()} for i in range(len(next(iter(ref.values()), [])))]
    warnings = ['Only coordinate model 1 is analyzed'] if len(st) > 1 else []
    if any(c['missing_sidechains'] for c in chains):
        warnings.append('Missing sidechain heavy atoms; do not interpret missing interactions as absent chemistry')
    if any(c['species_status'] == 'unknown' for c in chains):
        warnings.append('Source species metadata unavailable for some chains')
    if st.resolution and st.resolution > 3.5:
        warnings.append('Resolution worse than 3.5 A: sidechain-level chemical interpretation requires caution')
    if any(c['observed_fraction'] is not None and c['observed_fraction'] < .9 for c in chains):
        warnings.append('Less than 90% coordinate coverage in at least one chain')
    return {'path': str(Path(path).resolve()), 'sha256': sha256(path), 'experimental_methods': methods,
            'experimental_provenance_status': 'RCSB_metadata_backed' if (metadata or {}).get('entry') else 'local_file_declared_experimental',
            'model_count': len(st),
            'selected_model': 1, 'resolution_A': quality.get('resolution_combined') or ([st.resolution] if st.resolution else None),
            'refinement': refinement,
            'chains': chains, 'missing_residue_records': missing,
            'assemblies': [{'name': a.name, 'author_determined': a.author_determined, 'software_determined': a.software_determined} for a in st.assemblies],
            'warnings': warnings}


def clean_chain(st, chain_id, new_id):
    matches = [c for c in st[0] if c.name == chain_id and protein_residues(c)]
    if len(matches) != 1:
        raise Blocked(f'Choose one protein chain: {chain_id!r}; available {[c.name for c in st[0] if protein_residues(c)]}')
    source = matches[0]
    chain = gemmi.Chain(new_id)
    mapping, decisions = [], []
    missing_heavy_atoms = []
    seen = set()
    for residue in protein_residues(source):
        if residue.name not in STANDARD:
            raise Blocked(f'Nonstandard residue {residue.name} {chain_id}:{residue.seqid}; requires explicit preparation policy')
        key = str(residue.seqid)
        if key in seen:
            raise Blocked(f'Duplicate/microheterogeneous residue {chain_id}:{key}')
        seen.add(key)
        scores = defaultdict(float)
        for atom in residue:
            if atom.altloc != '\x00' and not atom.is_hydrogen():
                scores[atom.altloc] += atom.occ
        chosen = sorted(scores, key=lambda k: (-scores[k], k))[0] if scores else '\x00'
        r = gemmi.Residue()
        r.name, r.seqid, r.het_flag = residue.name, residue.seqid, 'A'
        # Label numbering is distinct from author numbering and retained in mapping.
        r.subchain = new_id
        r.label_seq = residue.label_seq
        names = set()
        for atom in residue:
            if atom.is_hydrogen() or atom.occ <= 0 or atom.altloc not in ('\x00', chosen):
                continue
            if atom.name in names:
                raise Blocked(f'Duplicate selected atom {chain_id}:{key}:{atom.name}')
            names.add(atom.name)
            a = atom.clone()
            a.altloc = '\x00'
            r.add_atom(a)
        missing = {'N', 'CA', 'C', 'O'} - names
        if missing:
            raise Blocked(f'Missing backbone atoms {chain_id}:{key}: {sorted(missing)}')
        sidechain_missing = SIDECHAIN_ATOMS[residue.name] - names
        if sidechain_missing:
            missing_heavy_atoms.append({'residue': key, 'name': residue.name, 'atoms': sorted(sidechain_missing), 'policy': 'observed atoms retained; not rebuilt'})
        if not -999 <= r.seqid.num <= 9999:
            raise Blocked('Author residue number cannot be represented in PDB; explicit renumbering policy required')
        chain.add_residue(r)
        mapping.append({'source_chain': chain_id, 'source_auth_seq_id': residue.seqid.num, 'source_ins_code': residue.seqid.icode.strip(),
                        'source_label_asym_id': residue.subchain, 'source_label_seq_id': residue.label_seq,
                        'prepared_chain': new_id, 'prepared_residue': key, 'residue_name': residue.name, 'selected_altloc': chosen if chosen != '\x00' else ''})
        if scores:
            decisions.append({'residue': key, 'altloc': chosen, 'occupancy_sums': dict(scores)})
    excluded = Counter(r.name for r in source if not gemmi.find_tabulated_residue(r.name).is_amino_acid())
    return chain, mapping, {'altloc_choices': decisions, 'missing_sidechain_heavy_atoms': missing_heavy_atoms,
                            'excluded_nonprotein': dict(excluded), 'hydrogens': 'removed; Arpeggio will add H at configured pH',
                            'missing_sidechains_policy': 'not rebuilt; chemistry results require review'}


def make_structure(chains):
    st = gemmi.Structure()
    st.name = 'easy_ppi'
    model = gemmi.Model('1')
    for c in chains:
        model.add_chain(c)
    st.add_model(model)
    st.setup_entities()
    st.assign_label_seq_id()
    return st


def write_structure(st, path):
    path = Path(path)
    if path.suffix == '.pdb':
        st.assign_serial_numbers()
        st.write_pdb(str(path))
    else:
        document = st.make_mmcif_document()
        block = document.sole_block()
        # Gemmi cannot infer chemical component types without a monomer library.
        # Prepared chains contain standard amino acids only, checked above.
        table = block.find('_chem_comp.', ['id', 'type'])
        for row in table:
            if row[0] in STANDARD:
                row[1] = "'peptide linking'" if row[0] == 'GLY' else "'L-peptide linking'"
        document.write_file(str(path))


def coordinates(st, chain):
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for c in st[0] if c.name == chain for r in c for a in r if not a.is_hydrogen()], dtype=float)


def geometry(st):
    a, b = coordinates(st, 'A'), coordinates(st, 'B')
    if not len(a) or not len(b):
        raise Blocked('Empty protein partner')
    distances, _ = cKDTree(b).query(a)
    return {'min_interchain_distance_A': float(distances.min()), 'atoms_A_near_B_5A': int((distances <= 5).sum()),
            'severe_clash_atom_count_A_1_5A': int((distances < 1.5).sum()), 'contact_definition': 'geometry only; not chemical interaction classification'}


def assembly_pair(path, chain_a, chain_b, assembly_id=None):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    def selector(value):
        match = re.fullmatch(r'(.+?)(?:@(\d+))?', value)
        return match[1], int(match[2]) if match[2] else None
    source_a, copy_a = selector(chain_a)
    source_b, copy_b = selector(chain_b)
    options = []
    for assembly in st.assemblies:
        if assembly_id and assembly.name != str(assembly_id):
            continue
        assembled = gemmi.make_assembly(assembly, st[0], gemmi.HowToNameCopiedChain.Dup)
        ac = [c for c in assembled if c.name == source_a]
        bc = [c for c in assembled if c.name == source_b]
        ia, ib = copy_a, copy_b
        if source_a == source_b and ia is None and ib is None and len(ac) == 2:
            ia, ib = 1, 2
        if ia is None and len(ac) == 1:
            ia = 1
        if ib is None and len(bc) == 1:
            ib = 1
        if ia and ib and 1 <= ia <= len(ac) and 1 <= ib <= len(bc) and (source_a, ia) != (source_b, ib):
            ca, ma, da = clean_chain(make_structure([ac[ia-1]]), source_a, 'A')
            cb, mb, db = clean_chain(make_structure([bc[ib-1]]), source_b, 'B')
            for mapping, index in ((ma, ia), (mb, ib)):
                for row in mapping:
                    row['source_assembly_copy'] = index
            pair = make_structure([ca, cb])
            check = geometry(pair)
            if check['atoms_A_near_B_5A']:
                options.append((pair, ma + mb, {'assembly_id': assembly.name, 'source_chain_selectors': [f'{source_a}@{ia}', f'{source_b}@{ib}'],
                                'geometry': check, 'preparation_A': da, 'preparation_B': db}))
    if len(options) != 1:
        raise Blocked(f'Expected one contact-bearing biological assembly for {chain_a}/{chain_b}; found {len(options)}. Specify --assembly-id or review assembly copies.')
    return options[0]


def prepare_pair(run, path_a, path_b, chain_a, chain_b, assembly_id=None, require_assembly=True):
    if Path(path_a).resolve() == Path(path_b).resolve() and require_assembly:
        st, mapping, info = assembly_pair(path_a, chain_a, chain_b, assembly_id)
    else:
        a, ma, da = clean_chain(gemmi.read_structure(str(path_a)), chain_a, 'A')
        b, mb, db = clean_chain(gemmi.read_structure(str(path_b)), chain_b, 'B')
        st, mapping, info = make_structure([a, b]), ma + mb, {'preparation_A': da, 'preparation_B': db}
    folder = run.path / 'prepared'
    for suffix in ('pdb', 'cif'):
        write_structure(st, folder / f'complex.{suffix}')
    for chain in st[0]:
        write_structure(make_structure([chain]), folder / f'protein_{chain.name}.pdb')
    write_csv(folder / 'residue_mapping.csv', mapping, list(mapping[0]))
    write_json(folder / 'preparation.json', info)
    return st


def validate_pdb(run, path):
    import importlib.util
    if importlib.util.find_spec('pdbtools') is None:
        raise Blocked('pdb-tools is required')
    output = run.command([sys.executable, '-m', 'pdbtools.pdb_validate', path])
    (run.path / 'prepared' / 'pdb_validate.txt').write_text(output, encoding='utf-8')
    if 'ERROR' in output.upper():
        raise ValueError(f'pdb-tools validation failed: {output[-1000:]}')


def scramble_partners(st, seed):
    from scipy.spatial.transform import Rotation
    out = st.clone()
    transform = {}
    for name in ('A', 'B'):
        xyz = coordinates(out, name)
        center = xyz.mean(axis=0)
        rotation = np.eye(3) if name == 'A' else Rotation.random(random_state=seed).as_matrix()
        for c in out[0]:
            if c.name == name:
                for r in c:
                    for a in r:
                        p = rotation @ (np.array([a.pos.x, a.pos.y, a.pos.z])-center)
                        a.pos = gemmi.Position(*p)
        transform[name] = {'subtract_centroid': center.tolist(), 'rotation_matrix': rotation.tolist()}
    return out, transform
