import json, sys, urllib.request, ssl, argparse

from common import load_config

WEIGHTS = {
    'xray': 3, 'resolution_le_2_0': 3, 'resolution_le_2_5': 2,
    'resolution_le_3_0': 1, 'co_crystal_ligand': 3,
    'relevant_co_crystal_ligand': 4, 'binding_site_defined': 2,
    'active_site_evidence': 2, 'human_protein': 2, 'recent_deposition': 1,
}
PENALTIES = {
    'missing_gt_5pct': -2, 'missing_gt_10pct': -3,
    'mutations': -2, 'no_ligand': -2, 'non_xray': -4,
}
IGNORE_RES = {'HOH', 'WAT', 'GOL', 'EDO', 'PEG', 'SO4', 'PO4'}

def fetch_pdb_metadata(pdb_id):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    url = f'https://data.rcsb.org/rest/v1/core/entry/{pdb_id}'
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, context=ctx, timeout=15) as r:
            return json.loads(r.read().decode('utf-8'))
    except Exception as e:
        print(f'WARNING: REST fetch failed for {pdb_id}: {e}', file=sys.stderr)
        return None

def batch_fetch_metadata(pdb_ids):
    results = {}
    for pid in pdb_ids:
        data = fetch_pdb_metadata(pid)
        if data:
            results[pid] = data
    return results

def extract_ligands(entry_data):
    direct = entry_data.get("ligands", [])
    if direct:
        return sorted({str(item).upper() for item in direct if item and str(item).upper() not in IGNORE_RES})
    try:
        entities = entry_data.get('rcsb_entry_container_identifiers', {})
        non_poly = entities.get('non_polymer_entity_ids', [])
        return [e.get('rcsb_chem_comp_info', {}).get('chem_comp_id', '')
                for e in non_poly
                if isinstance(e, dict)
                if e.get('rcsb_chem_comp_info', {}).get('chem_comp_id', '') not in IGNORE_RES]
    except Exception:
        return []

def score_pdb_entry(entry_data, scoring=None):
    scoring = scoring or {}
    weights = {
        key: int(scoring.get(key, value))
        for key, value in WEIGHTS.items()
    }
    penalties = {
        "missing_gt_5pct": int(scoring.get(
            "penalty_missing_gt_5pct", PENALTIES["missing_gt_5pct"])),
        "missing_gt_10pct": int(scoring.get(
            "penalty_missing_gt_10pct", PENALTIES["missing_gt_10pct"])),
        "mutations": int(scoring.get(
            "penalty_mutations", PENALTIES["mutations"])),
        "no_ligand": int(scoring.get(
            "penalty_no_ligand", PENALTIES["no_ligand"])),
        "non_xray": int(scoring.get(
            "penalty_non_xray", PENALTIES["non_xray"])),
    }
    score = 0; pos = []; neg = []
    method = entry_data.get('exptl', [{}])[0].get('method', '?')
    res_info = entry_data.get('rcsb_entry_info', {})
    resolution = res_info.get('resolution_combined', [99])[0] or 99
    quality = entry_data.get("easy_md", {})
    organism = entry_data.get("organism", quality.get("organism", ""))
    is_human = bool(entry_data.get("is_human")) or "homo sapiens" in str(organism).lower()
    try:
        dep_year = int(entry_data.get('rcsb_accession_info', {}).get('deposit_date', '2000-01-01')[:4])
    except Exception:
        dep_year = 2000
    if method == 'X-RAY DIFFRACTION':
        score += weights['xray']; pos.append('X-ray diffraction')
    else:
        neg.append(f'Method: {method} (not X-RAY)'); score += penalties['non_xray']
    if resolution <= 2.0:
        score += weights['resolution_le_2_0']; pos.append(f'Resolution {resolution} A (atomic)')
    elif resolution <= 2.5:
        score += weights['resolution_le_2_5']; pos.append(f'Resolution {resolution} A (good)')
    elif resolution <= 3.0:
        score += weights['resolution_le_3_0']; pos.append(f'Resolution {resolution} A (acceptable)')
    ligands = extract_ligands(entry_data)
    if ligands:
        score += weights['co_crystal_ligand']; pos.append(f'Co-crystal ligand(s): {chr(44).join(sorted(ligands))}')
    else:
        score += penalties['no_ligand']; neg.append('No co-crystal ligand metadata')
    if quality.get("binding_site_defined"):
        score += weights["binding_site_defined"]; pos.append("Binding site defined")
    if quality.get("relevant_ligand_evidence"):
        score += weights["relevant_co_crystal_ligand"]
        pos.append(f"Relevant ligand evidence: {quality['relevant_ligand_evidence']}")
    if quality.get("active_site_evidence"):
        score += weights["active_site_evidence"]
        pos.append(f"Active-site evidence: {quality['active_site_evidence']}")
    if is_human:
        score += weights['human_protein']; pos.append('Human protein')
    recent_year = int(scoring.get("recent_deposition_year", 2015))
    if dep_year >= recent_year:
        score += weights['recent_deposition']; pos.append(f'Deposited {dep_year}')
    missing_pct = quality.get("missing_residue_pct")
    if isinstance(missing_pct, (int, float)) and missing_pct > 10:
        score += penalties["missing_gt_10pct"]; neg.append(f"Missing residues {missing_pct:.1f}%")
    elif isinstance(missing_pct, (int, float)) and missing_pct > 5:
        score += penalties["missing_gt_5pct"]; neg.append(f"Missing residues {missing_pct:.1f}%")
    if quality.get("mutations"):
        score += penalties["mutations"]; neg.append("Engineered mutation(s)")
    title = entry_data.get('struct', {}).get('title', 'Unknown')[:80]
    info = {'method': method, 'resolution': resolution, 'dep_year': dep_year,
            'ligands': ligands, 'title': title}
    return score, pos, neg, info

def rank_pdb_entries(pdb_metadata, scoring=None):
    ranked = []
    for pdb_id, data in pdb_metadata.items():
        s, pos, neg, info = score_pdb_entry(data, scoring)
        info['pdb_id'] = pdb_id
        ranked.append({'pdb_id': pdb_id, 'score': s, 'positives': pos,
                       'negatives': neg, 'info': info})
    ranked.sort(key=lambda x: (-x['score'], x['info']['resolution']))
    return ranked

def format_recommendation(ranked):
    if not ranked:
        return 'No candidates to recommend.'
    top = ranked[0]; info = top['info']
    lines = ['=' * 64, '  RECOMMENDED STRUCTURE', '=' * 64,
             f'  PDB ID:  {top["pdb_id"]:<45s} Score: {top["score"]}', '=' * 64]
    for p in top['positives']: lines.append(f'  + {p}')
    for n in top['negatives']: lines.append(f'  - {n}')
    lines.append('=' * 64)
    lines.append(f'  Title: {info["title"]}')
    lines.append('=' * 64)
    if len(ranked) > 1:
        lines.append(''); lines.append('Backup candidates:')
        for i, r in enumerate(ranked[1:], 2):
            inf = r['info']; ligs = ', '.join(inf['ligands']) if inf['ligands'] else 'none'
            lines.append(f'  Rank {i}: {r["pdb_id"]} - Score {r["score"]} | '
                         f'{inf["method"][:5]} | {inf["resolution"]} A | Ligand: {ligs}')
    return chr(10).join(lines)

def main():
    parser = argparse.ArgumentParser(description='Score and rank PDB entries for docking')
    parser.add_argument('--pdb-ids', required=True, help='Comma-separated PDB IDs')
    parser.add_argument('--json', action='store_true')
    parser.add_argument('--metadata-file', help='Pre-fetched JSON metadata')
    parser.add_argument('--config', help='Path to easy-md config.yaml')
    args = parser.parse_args()
    config = load_config(args.config)
    pdb_ids = [p.strip() for p in args.pdb_ids.split(',') if p.strip()]
    if args.metadata_file:
        with open(args.metadata_file, 'r', encoding='utf-8') as f:
            pdb_metadata = json.load(f)
    else:
        print(f'Fetching metadata for {len(pdb_ids)} PDB entries...', file=sys.stderr)
        pdb_metadata = batch_fetch_metadata(pdb_ids)
    ranked = rank_pdb_entries(pdb_metadata, config.get("pdb_scoring", {}))
    if args.json:
        print(json.dumps(ranked, indent=2, ensure_ascii=False))
    else:
        print(format_recommendation(ranked))
    if ranked:
        print(f'\n# TOP_PDB={ranked[0]["pdb_id"]}')

if __name__ == '__main__':
    main()
