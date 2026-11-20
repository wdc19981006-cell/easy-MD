"""Enable versioned Git hooks for this checkout without changing global config."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main():
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    current = subprocess.run(["git", "config", "--local", "--get", "core.hooksPath"], capture_output=True, text=True).stdout.strip()
    if current and current != ".githooks":
        raise SystemExit(f"Existing hooksPath {current!r}; preserve it and integrate the guard manually")
    subprocess.run(["git", "config", "--local", "core.hooksPath", ".githooks"], check=True)
    subprocess.run(["git", "config", "--local", "easyMD.guardPython", sys.executable.replace("\\", "/")], check=True)
    print(f"Enabled pre-commit/pre-push hooks in {root}; Python={sys.executable}")


if __name__ == "__main__":
    main()
