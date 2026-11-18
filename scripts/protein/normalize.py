#!/usr/bin/env python
"""Normalize a protein name to a gene without silently resolving ambiguity."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from protein_search import query_mygene  # noqa: E402


# High-confidence protein aliases whose punctuation-insensitive forms otherwise
# collide with unrelated gene names in full-text search.  The mapped gene is
# still verified through MyGene; this table only chooses the query term.
CURATED_GENE_ALIASES = {
    "cox2": "PTGS2",
    "cyclooxygenase2": "PTGS2",
    "prostaglandinendoperoxidesynthase2": "PTGS2",
}


def _key(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).casefold())


def normalize_protein(name: str, species: str = "human") -> dict[str, object]:
    curated_symbol = CURATED_GENE_ALIASES.get(_key(name))
    query_term = curated_symbol or name
    hits = query_mygene(query_term, species)
    unique: list[dict[str, object]] = []
    seen: set[str] = set()
    for hit in hits:
        symbol = str(hit.get("symbol", "")).upper()
        if not symbol or symbol in seen:
            continue
        seen.add(symbol)
        unique.append(hit)

    input_key = _key(query_term)
    exact = [
        hit for hit in unique
        if _key(hit.get("symbol", "")) == input_key
        or _key(hit.get("name", "")) == input_key
        or input_key in {
            _key(value) for value in hit.get("aliases", [])
        }
    ]
    candidates = exact or unique
    result: dict[str, object] = {
        "input": name,
        "species": species,
        "resolved": False,
        "ambiguous": False,
        "selected": None,
        "candidates": candidates,
        "source": "curated alias + MyGene.info REST" if curated_symbol else "MyGene.info REST",
        "query_term": query_term,
        "curated_alias": curated_symbol,
    }
    if not candidates:
        return result
    if len(candidates) > 1:
        result["ambiguous"] = True
        return result
    selected = candidates[0]
    result["resolved"] = True
    result["selected"] = {
        "symbol": str(selected["symbol"]).upper(),
        "name": selected.get("name", ""),
        "aliases": selected.get("aliases", []),
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("protein")
    parser.add_argument("--species", default="human")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = normalize_protein(args.protein, args.species)
    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        print(result)
    if result["ambiguous"]:
        return 2
    return 0 if result["resolved"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
