#!/usr/bin/env python
"""Search, enrich, and download RCSB PDB structures through stable REST APIs."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import ssl
import urllib.request
from pathlib import Path
from typing import Any

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/rest/v1/core"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
DOWNLOAD_URL = "https://files.rcsb.org/download"


def _json_request(url: str, payload: dict[str, Any] | None = None) -> dict[str, Any]:
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    request = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json", "User-Agent": "easy-md/4.0"},
        method="POST" if data else "GET",
    )
    with urllib.request.urlopen(request, context=ssl.create_default_context(), timeout=15) as response:
        return json.loads(response.read().decode("utf-8"))


def search_pdb_ids(term: str, rows: int = 25) -> list[str]:
    # Gene-symbol search avoids the large set of unrelated entries returned by
    # full text (for example PTGS2 appearing only in a citation).  The fallback
    # keeps natural-language protein names usable.
    gene_query = {
        "type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
            "operator": "exact_match", "value": term.upper(),
            "case_sensitive": True,
        },
    }
    full_text_query = {
        "type": "terminal", "service": "full_text", "parameters": {"value": term},
    }

    def request(query: dict[str, Any], return_type: str) -> list[str]:
        payload = {
            "query": query,
            "return_type": return_type,
            "request_options": {
                "paginate": {"start": 0, "rows": rows},
                "sort": [{"sort_by": "score", "direction": "desc"}],
            },
        }
        result = _json_request(SEARCH_URL, payload)
        return list(dict.fromkeys(
            item["identifier"].split("_")[0].upper()
            for item in result.get("result_set", [])
        ))

    try:
        exact = request(gene_query, "polymer_entity")
    except Exception:
        exact = []
    return exact or request(full_text_query, "entry")


def enrich_entries_bulk(pdb_ids: list[str]) -> dict[str, dict[str, Any]]:
    """Fetch ranking fields for many entries in one RCSB GraphQL request."""
    if not pdb_ids:
        return {}
    query = """
    query EasyMdEntries($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        exptl { method }
        struct { title }
        rcsb_accession_info { deposit_date }
        rcsb_entry_info { resolution_combined }
        polymer_entities {
          rcsb_polymer_entity_container_identifiers { auth_asym_ids }
          rcsb_entity_source_organism {
            ncbi_scientific_name
            rcsb_gene_name { value }
          }
          rcsb_polymer_entity { pdbx_description }
        }
        nonpolymer_entities {
          nonpolymer_comp { chem_comp { id } }
        }
      }
    }
    """
    payload = {"query": query, "variables": {"ids": [value.upper() for value in pdb_ids]}}
    response = _json_request(GRAPHQL_URL, payload)
    if response.get("errors"):
        raise RuntimeError(f"RCSB GraphQL error: {response['errors']}")
    result: dict[str, dict[str, Any]] = {}
    for entry in response.get("data", {}).get("entries", []) or []:
        if not entry:
            continue
        pdb_id = str(entry.get("rcsb_id", "")).upper()
        organisms: set[str] = set()
        chains: set[str] = set()
        gene_names: set[str] = set()
        descriptions: set[str] = set()
        mutation_count = 0
        for entity in entry.get("polymer_entities", []) or []:
            container = entity.get("rcsb_polymer_entity_container_identifiers", {}) or {}
            chains.update(container.get("auth_asym_ids", []) or [])
            polymer = entity.get("rcsb_polymer_entity", {}) or {}
            if polymer.get("pdbx_description"):
                descriptions.add(str(polymer["pdbx_description"]))
            for source in entity.get("rcsb_entity_source_organism", []) or []:
                if source.get("ncbi_scientific_name"):
                    organisms.add(str(source["ncbi_scientific_name"]))
                for gene in source.get("rcsb_gene_name", []) or []:
                    if gene.get("value"):
                        gene_names.add(str(gene["value"]).upper())
        ligands = sorted({
            str((entity.get("nonpolymer_comp", {}) or {}).get("chem_comp", {}).get("id", "")).upper()
            for entity in entry.get("nonpolymer_entities", []) or []
            if (entity.get("nonpolymer_comp", {}) or {}).get("chem_comp", {}).get("id")
        })
        entry["ligands"] = ligands
        entry["easy_md"] = {
            "organisms": sorted(organisms),
            "organism": "; ".join(sorted(organisms)),
            "is_human": any(value.casefold() == "homo sapiens" for value in organisms),
            "chains": sorted(chains),
            "gene_names": sorted(gene_names),
            "descriptions": sorted(descriptions),
            "mutations": mutation_count > 0,
            "mutation_count": mutation_count,
            "binding_site_defined": bool(ligands),
        }
        entry["is_human"] = entry["easy_md"]["is_human"]
        entry["organism"] = entry["easy_md"]["organism"]
        result[pdb_id] = entry
    return result


def enrich_entry(pdb_id: str) -> dict[str, Any]:
    pdb_id = pdb_id.upper()
    entry = _json_request(f"{DATA_URL}/entry/{pdb_id}")
    identifiers = entry.get("rcsb_entry_container_identifiers", {})
    organisms: set[str] = set()
    chains: set[str] = set()
    gene_names: set[str] = set()
    descriptions: set[str] = set()
    mutation_count = 0
    for entity_id in identifiers.get("polymer_entity_ids", []) or []:
        try:
            entity = _json_request(f"{DATA_URL}/polymer_entity/{pdb_id}/{entity_id}")
        except Exception:
            continue
        for source in entity.get("rcsb_entity_source_organism", []) or []:
            name = source.get("ncbi_scientific_name")
            if name:
                organisms.add(str(name))
        container = entity.get("rcsb_polymer_entity_container_identifiers", {})
        chains.update(container.get("auth_asym_ids", []) or [])
        polymer = entity.get("rcsb_polymer_entity", {})
        mutation_count += int(polymer.get("mutation_count") or 0)
        description = polymer.get("pdbx_description")
        if description:
            descriptions.add(str(description))
        for gene in polymer.get("rcsb_gene_name", []) or entity.get("rcsb_gene_name", []) or []:
            if isinstance(gene, dict):
                value = gene.get("value") or gene.get("name")
            else:
                value = gene
            if value:
                gene_names.add(str(value).upper())

    ligands: list[str] = []
    for entity_id in identifiers.get("non_polymer_entity_ids", []) or []:
        try:
            entity = _json_request(f"{DATA_URL}/nonpolymer_entity/{pdb_id}/{entity_id}")
        except Exception:
            continue
        comp_id = entity.get("pdbx_entity_nonpoly", {}).get("comp_id")
        if comp_id:
            ligands.append(str(comp_id).upper())

    entry["easy_md"] = {
        "organisms": sorted(organisms),
        "organism": "; ".join(sorted(organisms)),
        "is_human": any(name.casefold() == "homo sapiens" for name in organisms),
        "chains": sorted(chains),
        "gene_names": sorted(gene_names),
        "descriptions": sorted(descriptions),
        "mutations": mutation_count > 0,
        "mutation_count": mutation_count,
        "binding_site_defined": bool(ligands),
    }
    entry["ligands"] = sorted(set(ligands))
    entry["is_human"] = entry["easy_md"]["is_human"]
    entry["organism"] = entry["easy_md"]["organism"]
    return entry


def search_structures(term: str, rows: int = 25) -> dict[str, dict[str, Any]]:
    pdb_ids = search_pdb_ids(term, rows)
    try:
        return enrich_entries_bulk(pdb_ids)
    except Exception:
        pass
    result: dict[str, dict[str, Any]] = {}
    # RCSB core metadata is split across entry/entity endpoints.  Entries are
    # independent, so bounded concurrency reduces normal preparation from
    # minutes to seconds without weakening any scoring evidence.
    with ThreadPoolExecutor(max_workers=min(6, max(1, len(pdb_ids)))) as pool:
        futures = {pool.submit(enrich_entry, pdb_id): pdb_id for pdb_id in pdb_ids}
        for future in as_completed(futures):
            pdb_id = futures[future]
            try:
                result[pdb_id] = future.result()
            except Exception as exc:
                result[pdb_id] = {"error": str(exc)}
    return result


def download_pdb(pdb_id: str, output: str | Path) -> Path:
    target = Path(output).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        f"{DOWNLOAD_URL}/{pdb_id.upper()}.pdb",
        headers={"User-Agent": "easy-md/4.0"},
    )
    with urllib.request.urlopen(request, context=ssl.create_default_context(), timeout=60) as response:
        payload = response.read()
    if not payload.startswith((b"HEADER", b"ATOM", b"REMARK")):
        raise ValueError(f"RCSB did not return a PDB file for {pdb_id}")
    target.write_bytes(payload)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("term", nargs="?")
    parser.add_argument("--pdb-id")
    parser.add_argument("--rows", type=int, default=25)
    parser.add_argument("--download")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    if args.pdb_id:
        result = {args.pdb_id.upper(): enrich_entry(args.pdb_id)}
        if args.download:
            download_pdb(args.pdb_id, args.download)
    elif args.term:
        result = search_structures(args.term, args.rows)
    else:
        parser.error("provide a search term or --pdb-id")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
