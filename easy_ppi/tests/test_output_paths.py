import json
from pathlib import Path

import pytest

from easy_ppi.audit import Run
from easy_ppi.paths import default_output_root


def test_default_run_uses_configured_root_instead_of_cwd(tmp_path, monkeypatch):
    destination = tmp_path / 'desktop-results'
    configuration = tmp_path / 'config.json'
    configuration.write_text(json.dumps({'paths': {'output_root': str(destination)}}), encoding='utf-8')
    workspace = tmp_path / 'project'
    workspace.mkdir()
    monkeypatch.chdir(workspace)
    monkeypatch.delenv('EASY_PPI_OUTPUT_ROOT', raising=False)
    monkeypatch.setenv('EASY_PPI_CONFIG', str(configuration))
    run = Run()
    assert run.path.parent == destination
    assert not (workspace / 'ppi_results').exists()
    assert run.state['output_directory'] == str(run.path)


def test_output_root_environment_override(tmp_path, monkeypatch):
    monkeypatch.setenv('EASY_PPI_OUTPUT_ROOT', str(tmp_path / 'override'))
    monkeypatch.setenv('EASY_PPI_CONFIG', str(tmp_path / 'missing.json'))
    assert default_output_root() == tmp_path / 'override'


def test_bad_default_cannot_silently_write_to_project(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('EASY_PPI_OUTPUT_ROOT', 'ppi_results')
    with pytest.raises(ValueError, match='absolute'):
        Run()
    assert not (tmp_path / 'ppi_results').exists()


def test_explicit_run_directory_still_overrides_default(tmp_path, monkeypatch):
    monkeypatch.setenv('EASY_PPI_CONFIG', str(tmp_path / 'missing.json'))
    monkeypatch.delenv('EASY_PPI_OUTPUT_ROOT', raising=False)
    target = tmp_path / 'explicit'
    run = Run(target)
    assert run.path == target
