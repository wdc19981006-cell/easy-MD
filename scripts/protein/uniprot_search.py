#!/usr/bin/env python
"""Reusable UniProt REST lookup for a resolved gene symbol."""

from __future__ import annotations

import argparse
import json
import ssl
import urllib.parse
import urllib.request

URL = "https://rest.uniprot.org/uniprotkb/search"


def search_uniprot(symbol: str, organism_id: int = 9606, size: int = 10) -> list[dict[str, object]]:
    params = urllib.parse.urlencode({
        "query": f"gene_exact:{symbol} AND organism_id:{organism_id}",
        "format": "json",
        "size": size,
        "fields": "accession,id,protein_name,gene_names,organism_name,reviewed",
    })
    request = urllib.request.Request(
        f"{URL}?{params}", headers={"User-Agent": "easy-md/4.0"}
    )
    context = ssl.create_default_context()
    with urllib.request.urlopen(request, context=context, timeout=30) as response:
        payload = json.loads(response.read().decode("utf-8"))
    return payload.get("results", [])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("symbol")
    parser.add_argument("--organism-id", type=int, default=9606)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    result = search_uniprot(args.symbol, args.organism_id)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result else 1


if __name__ == "__main__":
    raise SystemExit(main())
