"""Read-only environment inventory; does not import the Vina workflow."""
import importlib.metadata as md
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

def inspect():
    from .paths import default_output_root
    packages = {}
    for name in ['gemmi', 'pdb-tools', 'lightdock', 'pdbe-arpeggio', 'freesasa', 'prodigy-prot', 'DockQ', 'numpy', 'biopython', 'openbabel']:
        try:
            dist = md.distribution(name)
            packages[name] = {'version': dist.version, 'license': (dist.metadata.get('License') or dist.metadata.get('License-Expression') or '')[:120], 'location': str(dist.locate_file(''))}
        except md.PackageNotFoundError:
            packages[name] = None
    executables = {name: shutil.which(name) for name in ['wsl', 'docker', 'conda', 'micromamba', 'gcc', 'cl', 'freesasa', 'pymol']}
    known_gcc = Path('D:/R/rtools45/x86_64-w64-mingw32.static.posix/bin/gcc.exe')
    if known_gcc.is_file():
        executables['rtools_gcc'] = str(known_gcc)
    wsl = subprocess.run(['wsl.exe', '--list', '--verbose'], capture_output=True) if os.name == 'nt' else None
    return {'python': sys.executable, 'platform': platform.platform(), 'packages': packages,
            'default_output_root': str(default_output_root()),
            'executables': executables, 'wsl': None if wsl is None else {'returncode': wsl.returncode, 'stdout': wsl.stdout.decode('utf-16', errors='replace'), 'stderr': wsl.stderr.decode('utf-16', errors='replace')},
            'pymol_python_exists': Path('D:/Interaction/PyMol/python.exe').exists(),
            'dockq_script': os.environ.get('EASY_PPI_DOCKQ_SCRIPT'),
            'scientific_environment': {key: os.environ.get(key) for key in ['BABEL_DATADIR', 'BABEL_LIBDIR', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']},
            'environment_scope': 'dedicated PPI interpreter; system-site packages may be inherited; Vina environment unchanged'}

if __name__ == '__main__':
    print(json.dumps(inspect(), indent=2, ensure_ascii=False))
