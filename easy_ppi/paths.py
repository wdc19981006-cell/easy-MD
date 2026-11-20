"""Configured output locations; never silently create results under the CWD."""
from __future__ import annotations

import json
import os
from pathlib import Path

DEFAULT_CONFIG_FILE = Path(__file__).resolve().with_name('config.json')


def default_output_root():
    override = os.environ.get('EASY_PPI_OUTPUT_ROOT')
    if override:
        value = override
    else:
        config_file = Path(os.environ.get('EASY_PPI_CONFIG') or DEFAULT_CONFIG_FILE).expanduser()
        configuration = json.loads(config_file.read_text(encoding='utf-8-sig'))
        value = configuration.get('paths', {}).get('output_root')
        if not value:
            raise ValueError(f'Set paths.output_root in {config_file}; project-local output fallback is disabled')
    root = Path(value).expanduser()
    if not root.is_absolute():
        raise ValueError(f'easy-PPI default output_root must be an absolute path: {value}')
    return root.resolve()
