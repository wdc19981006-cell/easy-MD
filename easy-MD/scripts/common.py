#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Shared configuration and work-directory helpers for easy-md."""

from __future__ import annotations

import os
import json
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml


SKILL_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = SKILL_ROOT / "config.yaml"
WORK_MARKER = ".easy-md-run"

STAGES = (
    "acquisition",
    "preparation",
    "docking",
    "interaction",
    "visualization",
    "packaging",
)


def load_config(config_path: str | os.PathLike[str] | None = None) -> dict[str, Any]:
    """Load the skill configuration."""
    path = Path(config_path).expanduser() if config_path else DEFAULT_CONFIG
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"easy-md config not found: {path}")
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"easy-md config must be a YAML mapping: {path}")
    return data


def get_config(config: dict[str, Any], dotted_key: str, default: Any = None) -> Any:
    """Read a dotted key from a nested config mapping."""
    value: Any = config
    for key in dotted_key.split("."):
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]
    return value


def configured_path(
    config: dict[str, Any],
    key: str,
    fallback: str | None = None,
) -> str | None:
    """Return a normalized path value from config.paths."""
    value = get_config(config, f"paths.{key}", fallback)
    if value in (None, ""):
        return None
    return str(Path(str(value)).expanduser())


def normalize_name_list(values: Any) -> set[str]:
    """Normalize YAML lists and tolerate legacy comma-delimited entries."""
    if values is None:
        return set()
    if isinstance(values, str):
        values = [values]
    result: set[str] = set()
    for value in values:
        result.update(part.strip().upper() for part in str(value).split(",") if part.strip())
    return result


def _unsafe_work_dir(path: Path) -> bool:
    home = Path.home().resolve()
    skill_root = SKILL_ROOT.resolve()
    anchors = {Path(path.anchor).resolve()} if path.anchor else set()
    protected = {home, skill_root, skill_root.parent, *anchors}
    return path in protected or len(path.parts) < 3


def resolve_work_dir(
    work_dir: str | os.PathLike[str] | None,
    config: dict[str, Any],
    *,
    create: bool = True,
    mark: bool = True,
) -> Path:
    """Resolve and validate a dedicated easy-md run directory.

    Existing non-empty directories without an easy-md marker are rejected.
    This prevents a typo in EASYMD_WORK_DIR from turning a broad directory
    into a cleanup target.
    """
    candidate = (
        work_dir
        or os.environ.get("EASYMD_WORK_DIR")
        or get_config(config, "paths.work_root")
    )
    if not candidate:
        raise ValueError(
            "A dedicated work directory is required. Pass --work-dir or set "
            "EASYMD_WORK_DIR."
        )

    path = Path(candidate).expanduser().resolve()
    if _unsafe_work_dir(path):
        raise ValueError(f"Refusing unsafe easy-md work directory: {path}")

    marker = path / WORK_MARKER
    if path.exists():
        if not path.is_dir():
            raise ValueError(f"Work directory path is not a directory: {path}")
        has_content = any(item.name != WORK_MARKER for item in path.iterdir())
        if has_content and not marker.is_file():
            raise ValueError(
                "Existing work directory is not marked as an easy-md run: "
                f"{path}"
            )
    elif create:
        path.mkdir(parents=True, exist_ok=False)
    else:
        raise FileNotFoundError(f"Work directory does not exist: {path}")

    if mark:
        marker.write_text("easy-md managed work directory\n", encoding="utf-8")
    return path


def enter_work_dir(
    work_dir: str | os.PathLike[str] | None,
    config: dict[str, Any],
) -> Path:
    """Resolve, mark, and change into a dedicated run directory."""
    path = resolve_work_dir(work_dir, config)
    os.chdir(path)
    return path


def validate_cleanup_target(path: str | os.PathLike[str]) -> Path:
    """Require a safe, marked run directory before recursive deletion."""
    resolved = Path(path).expanduser().resolve()
    if _unsafe_work_dir(resolved):
        raise ValueError(f"Refusing unsafe cleanup target: {resolved}")
    if not (resolved / WORK_MARKER).is_file():
        raise ValueError(
            f"Refusing cleanup without {WORK_MARKER} marker: {resolved}"
        )
    return resolved


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_slug(value: str, fallback: str = "item") -> str:
    """Return an ASCII-safe, deterministic directory component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", str(value)).strip("-._")
    return cleaned or fallback


def atomic_write_json(path: str | os.PathLike[str], value: Any) -> Path:
    """Write JSON atomically so interrupted runs do not corrupt state."""
    target = Path(path).resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{target.name}.", suffix=".tmp", dir=str(target.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return target


def new_status() -> dict[str, Any]:
    now = utc_now()
    return {
        "schema_version": 1,
        "created_at": now,
        "updated_at": now,
        "stages": {
            stage: {"status": "PENDING", "error": None, "artifacts": []}
            for stage in STAGES
        },
    }


def update_stage(
    status_path: str | os.PathLike[str],
    stage: str,
    state: str,
    *,
    error: str | None = None,
    artifacts: list[str] | None = None,
) -> dict[str, Any]:
    if stage not in STAGES:
        raise ValueError(f"Unknown easy-md stage: {stage}")
    state = state.upper()
    valid = {"PENDING", "RUNNING", "AWAITING_CONFIRMATION", "SUCCESS", "FAILED", "SKIPPED"}
    if state not in valid:
        raise ValueError(f"Invalid easy-md stage state: {state}")
    path = Path(status_path).resolve()
    status = json.loads(path.read_text(encoding="utf-8")) if path.exists() else new_status()
    current = status["stages"][stage]
    current["status"] = state
    current["error"] = error
    if artifacts is not None:
        current["artifacts"] = list(artifacts)
    current["updated_at"] = utc_now()
    status["updated_at"] = current["updated_at"]
    atomic_write_json(path, status)
    return status
