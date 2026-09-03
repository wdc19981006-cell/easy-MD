#!/usr/bin/env python
"""Reusable PubChem lookup and 3D SDF download adapter."""

from __future__ import annotations

import argparse
import json
import ssl
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
PROPERTIES = "CanonicalSMILES,IsomericSMILES,MolecularFormula,MolecularWeight,IUPACName"


def search_ligand(name: str) -> dict[str, object]:
    encoded = urllib.parse.quote(name, safe="")
    url = f"{BASE}/name/{encoded}/property/{PROPERTIES}/JSON"
    request = urllib.request.Request(url, headers={"User-Agent": "easy-md/4.0"})
    with urllib.request.urlopen(request, context=ssl.create_default_context(), timeout=30) as response:
        payload = json.loads(response.read().decode("utf-8"))
    row = payload["PropertyTable"]["Properties"][0]
    return {
        "input": name,
        "cid": row.get("CID"),
        "canonical_smiles": row.get("ConnectivitySMILES") or row.get("CanonicalSMILES"),
        "isomeric_smiles": row.get("SMILES") or row.get("IsomericSMILES"),
        "formula": row.get("MolecularFormula"),
        "molecular_weight": row.get("MolecularWeight"),
        "iupac_name": row.get("IUPACName"),
        "source": "PubChem PUG REST",
    }


def download_ligand_sdf(name_or_cid: str | int, output: str | Path) -> Path:
    namespace = "cid" if str(name_or_cid).isdigit() else "name"
    encoded = urllib.parse.quote(str(name_or_cid), safe="")
    url = f"{BASE}/{namespace}/{encoded}/record/SDF?record_type=3d"
    request = urllib.request.Request(url, headers={"User-Agent": "easy-md/4.0"})
    with urllib.request.urlopen(request, context=ssl.create_default_context(), timeout=60) as response:
        payload = response.read()
    if b"V2000" not in payload and b"V3000" not in payload:
        raise ValueError("PubChem did not return a valid SDF")
    target = Path(output).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(payload)
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ligand")
    parser.add_argument("--download")
    args = parser.parse_args()
    result = search_ligand(args.ligand)
    if args.download:
        result["sdf"] = str(download_ligand_sdf(result.get("cid") or args.ligand, args.download))
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
