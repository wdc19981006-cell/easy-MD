#!/usr/bin/env python
"""Compatibility entry point for native-ligand or explicit docking-box generation."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from box_prepare import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
