#!/usr/bin/env python
"""Rank enriched PDB candidates with the configured easy-md evidence weights."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from common import load_config  # noqa: E402
from pdb_selection import rank_pdb_entries  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata")
    parser.add_argument("--config")
    args = parser.parse_args()
    metadata = json.loads(Path(args.metadata).read_text(encoding="utf-8"))
    metadata = {key: value for key, value in metadata.items() if "error" not in value}
    config = load_config(args.config)
    ranked = rank_pdb_entries(metadata, config.get("pdb_scoring", {}))
    print(json.dumps(ranked, indent=2, ensure_ascii=False))
    return 0 if ranked else 1


if __name__ == "__main__":
    raise SystemExit(main())
