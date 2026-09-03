#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Validate easy-md configuration and installed dependencies without docking."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path

from common import configured_path, load_config


MODULES = {
    "rdkit": "rdkit",
    "pdbfixer": "pdbfixer",
    "openmm": "openmm",
    "biopython": "Bio",
    "numpy": "numpy",
    "plip": "plip",
    "pillow": "PIL",
    "pyyaml": "yaml",
}
CORE_TOOLS = {"python", "vina", "obabel"}
CORE_MODULES = {
    "rdkit",
    "pdbfixer",
    "openmm",
    "biopython",
    "numpy",
    "pyyaml",
}


def check_path(path: str | None) -> dict[str, object]:
    if not path:
        return {"path": None, "exists": False, "executable": False}
    resolved = Path(path).expanduser()
    return {
        "path": str(resolved),
        "exists": resolved.is_file(),
        "executable": resolved.is_file() and os.access(resolved, os.X_OK),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check easy-md paths and Python packages without docking"
    )
    parser.add_argument("--config", help="Path to easy-md config.yaml")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    config = load_config(args.config)
    tools = {
        key: check_path(configured_path(config, key))
        for key in ("python", "vina", "obabel", "pymol", "pymol_python")
    }
    modules = {
        label: importlib.util.find_spec(module_name) is not None
        for label, module_name in MODULES.items()
    }
    core_ok = all(tools[name]["exists"] for name in CORE_TOOLS) and all(
        modules[name] for name in CORE_MODULES
    )
    result = {
        "ok": core_ok,
        "optional": {
            "interaction_analysis": modules["plip"] and modules["pillow"],
            "pymol_session": (
                tools["pymol"]["exists"] and tools["pymol_python"]["exists"]
            ),
            "plip_native_visualization": (
                modules["plip"] and tools["pymol_python"]["exists"]
            ),
        },
        "tools": tools,
        "python_modules": modules,
    }

    if args.json:
        print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        print("easy-md preflight (no docking)")
        for name, item in tools.items():
            state = "OK" if item["exists"] else "MISSING"
            print(f"  {name:10s} {state:8s} {item['path']}")
        for name, available in modules.items():
            print(f"  {name:10s} {'OK' if available else 'MISSING'}")
        print(f"Core docking: {'OK' if result['ok'] else 'NOT READY'}")
        print(
            "Optional interaction analysis: "
            f"{'OK' if result['optional']['interaction_analysis'] else 'UNAVAILABLE'}"
        )
        print(
            "Optional PyMOL session: "
            f"{'OK' if result['optional']['pymol_session'] else 'UNAVAILABLE'}"
        )
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
