from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import uuid
import zipfile

from .paths import default_output_root


class Blocked(RuntimeError):
    """A missing prerequisite, never a successful scientific calculation."""


def now():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # NaN/Infinity must never enter a scientific result.
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows, fields):
    with Path(path).open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


class Run:
    def __init__(self, output=None, parameters=None):
        candidate = Path(output) if output is not None else default_output_root() / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '_' + uuid.uuid4().hex[:8])
        self.path = candidate.resolve()
        self.path.mkdir(parents=True, exist_ok=False)
        for folder in ('raw', 'prepared', 'analysis', 'docking', 'validation', 'visualization', 'logs'):
            (self.path / folder).mkdir()
        self.state = {'schema_version': 1, 'created_at': now(), 'parameters': parameters or {}, 'stages': {},
                      'scientific_validation': {'status': 'NOT_VALIDATED'}, 'commands': []}
        from . import __version__
        self.state['easy_ppi_version'] = __version__
        self.state['output_directory'] = str(self.path)
        source_root = Path(__file__).resolve().parent
        self.state['source_sha256'] = {str(p.relative_to(source_root)): sha256(p) for p in sorted(source_root.glob('*.py'))}
        provenance = self.path / 'provenance'
        provenance.mkdir()
        with zipfile.ZipFile(provenance / 'easy_ppi_source.zip', 'x', compression=zipfile.ZIP_DEFLATED) as archive:
            for source in sorted(source_root.iterdir()):
                if source.is_file() and source.suffix in ('.py', '.toml', '.txt', '.yml', '.md', '.json'):
                    archive.write(source, 'easy_ppi/' + source.name)
        self.save()

    def save(self):
        write_json(self.path / 'run.json', self.state)

    @contextmanager
    def stage(self, name):
        self.state['stages'][name] = {'status': 'RUNNING', 'start': now()}
        self.save()
        try:
            yield
        except Exception as exc:
            self.state['stages'][name].update(status='BLOCKED' if isinstance(exc, (Blocked, ModuleNotFoundError, FileNotFoundError)) else 'FAILED', error=str(exc))
        else:
            if self.state['stages'][name]['status'] == 'RUNNING':
                self.state['stages'][name]['status'] = 'SUCCESS'
        self.state['stages'][name]['end'] = now()
        self.save()

    def command(self, argv, cwd=None, timeout=600):
        argv = [str(x) for x in argv]
        index = len(self.state['commands'])
        record = {'argv': argv, 'cwd': str(cwd or self.path), 'start': now()}
        self.state['commands'].append(record)
        self.save()
        stdout_path = self.path / 'logs' / f'{index:04d}.stdout.txt'
        stderr_path = self.path / 'logs' / f'{index:04d}.stderr.txt'
        try:
            with stdout_path.open('wb') as out, stderr_path.open('wb') as err:
                result = subprocess.run(argv, cwd=cwd or self.path, stdout=out, stderr=err, timeout=timeout, check=False)
            record['returncode'] = result.returncode
        except (OSError, subprocess.TimeoutExpired) as exc:
            record['error'] = str(exc)
            raise Blocked(f'Cannot execute {argv[0]}: {exc}') from exc
        finally:
            record.update(end=now(), stdout=str(stdout_path.relative_to(self.path)), stderr=str(stderr_path.relative_to(self.path)))
            self.save()
        if result.returncode:
            tail = stderr_path.read_text(encoding='utf-8', errors='replace')[-2000:]
            raise RuntimeError(f'Command returned {result.returncode}: {argv}; {tail}')
        return stdout_path.read_text(encoding='utf-8', errors='replace')

    def finish(self):
        self.state['finished_at'] = now()
        self.save()
        from .report import write_report
        write_report(self.path, self.state)
        files = [{'path': str(p.relative_to(self.path)), 'bytes': p.stat().st_size, 'sha256': sha256(p)}
                 for p in sorted(self.path.rglob('*')) if p.is_file() and p.name != 'manifest.json']
        write_json(self.path / 'manifest.json', files)
