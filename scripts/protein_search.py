#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Protein name normalization via MyGene.info.

Resolves user-input protein names to standardized gene symbols and builds
search term lists for downstream PDB search.

Usage:
    python protein_search.py "ROCK2"
    python protein_search.py "EGFR" --species human
"""

import json, sys, urllib.request, urllib.parse, ssl, argparse

MYGENE_URL = "https://mygene.info/v3/query"


def query_mygene(protein_name, species="human"):
    """Query MyGene.info for gene-level information."""
    query = urllib.parse.urlencode({
        "q": protein_name,
        "species": species,
        "fields": "symbol,name,alias,other_names,entrezgene",
        "size": 10,
    })
    url = f"{MYGENE_URL}?{query}"
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, context=ctx, timeout=15) as r:
            data = json.loads(r.read().decode('utf-8'))
    except Exception as e:
        print(f"ERROR: MyGene.info query failed: {e}", file=sys.stderr)
        return []
    hits = data.get('hits', [])
    if not hits:
        return []
    results = []
    for h in hits:
        symbol = h.get('symbol', '?')
        name = h.get('name', '?')
        aliases = h.get('alias', [])
        other = h.get('other_names', [])
        if isinstance(aliases, str):
            aliases = [aliases]
        if isinstance(other, str):
            other = [other]
        all_aliases = list(set(str(a) for a in (aliases + other) if a))
        results.append({'symbol': symbol, 'name': name, 'aliases': all_aliases})
    return results


def detect_multi_gene(hits):
    """Check for multiple distinct gene symbols."""
    seen = {}
    for h in hits:
        sym = h['symbol']
        if sym not in seen:
            seen[sym] = h['name']
    return list(seen.items())


def build_search_terms(gene_info):
    """Build deduplicated search term list."""
    terms = [gene_info['symbol'], gene_info['name']]
    terms.extend(gene_info.get('aliases', []))
    return list(set(t.strip() for t in terms if t and t.strip()))


def resolve_protein(protein_name, species="human", interactive=True):
    """Resolve protein name to standardized identifiers.

    Returns dict with: resolved, multi_gene, genes, selected, all_hits.
    exit code 2 = multi-gene ambiguity, needs user choice.
    """
    result = {
        'resolved': False, 'multi_gene': False, 'genes': [],
        'selected': None, 'all_hits': [], 'input': protein_name,
        'species': species,
    }
    hits = query_mygene(protein_name, species)
    if not hits:
        return result
    result['all_hits'] = hits
    genes = detect_multi_gene(hits)
    result['genes'] = genes
    result['resolved'] = True
    if len(genes) > 1:
        result['multi_gene'] = True
        if interactive:
            return result
        selected_gene = hits[0]
    else:
        selected_gene = hits[0]
    terms = build_search_terms(selected_gene)
    result['selected'] = {
        'symbol': selected_gene['symbol'],
        'name': selected_gene['name'],
        'aliases': selected_gene.get('aliases', []),
        'search_terms': terms,
    }
    return result


def format_multi_gene_prompt(result):
    """User-friendly multi-gene selection prompt."""
    lines = [f"Multiple genes match '{result['input']}':"]
    for i, (sym, name) in enumerate(result['genes'], 1):
        lines.append(f"  {i}. {sym} — {name}")
    lines.append("Please select the target gene [1/2/3]:")
    return '\n'.join(lines)


def format_identity_report(result):
    """Format protein identity confirmation report."""
    sel = result['selected']
    lines = [
        "Protein identity confirmed:",
        f"  Gene symbol:   {sel['symbol']}",
        f"  Official name: {sel['name']}",
    ]
    if sel['aliases']:
        lines.append(f"  Aliases:       {', '.join(sel['aliases'][:6])}")
    lines.append(f"  Search terms:  {len(sel['search_terms'])} terms for PDB search")
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Resolve protein name to gene symbol via MyGene.info")
    parser.add_argument("protein", help="Protein name or keyword")
    parser.add_argument("--species", default="human")
    parser.add_argument("--json", action="store_true")
    parser.add_argument("--no-interactive", action="store_true")
    args = parser.parse_args()
    result = resolve_protein(args.protein, args.species,
                             interactive=not args.no_interactive)
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return
    if not result['resolved']:
        print(f"NO_RESULTS for '{args.protein}'")
        sys.exit(1)
    if result['multi_gene'] and not args.no_interactive:
        print(format_multi_gene_prompt(result))
        sys.exit(2)
    print(format_identity_report(result))


if __name__ == '__main__':
    main()
