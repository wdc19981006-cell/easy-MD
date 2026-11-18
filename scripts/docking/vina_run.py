#!/usr/bin/env python
"""Compatibility entry point for repeated seeded AutoDock Vina execution."""

from __future__ import annotations

import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from docking_run import main  # noqa: E402

if __name__ == "__main__":
    main()
