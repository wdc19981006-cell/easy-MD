"""Build unmodified upstream sources with an existing Windows GCC toolchain.

GCC >=14 turns a legacy NumPy C pointer conversion warning into an error.
Demote that one diagnostic for upstream LightDock 0.9.4; no math is changed.
This is a developer build helper, not part of scientific processing.
"""
from pathlib import Path
import os
import runpy
import sys

import setuptools

original = setuptools.setup

def configure(**kwargs):
    extensions = kwargs.get('ext_modules', [])
    regenerate = False
    for extension in extensions:
        sources = []
        for item in extension.sources:
            pyx = Path(item).with_suffix('.pyx')
            if pyx.is_file():
                original_c = Path(item)
                backup = original_c.with_suffix('.c.pre_cython')
                if original_c.is_file() and not backup.exists():
                    backup.write_bytes(original_c.read_bytes())
                sources.append(str(pyx))
                regenerate = True
            else:
                sources.append(item)
        extension.sources = sources
        if all(str(s).endswith('.c') for s in extension.sources):
            extension.extra_compile_args = list(extension.extra_compile_args or []) + ['-Wno-error=incompatible-pointer-types']
    if regenerate:
        from Cython.Build import cythonize
        print('Regenerating bundled Cython C from upstream .pyx for CPython 3.12 compatibility', flush=True)
        kwargs['ext_modules'] = cythonize(extensions, compiler_directives={'language_level': 3}, force=True)
    return original(**kwargs)

setuptools.setup = configure
source = Path(sys.argv[1]).resolve()
if not (source / 'setup.py').is_file():
    raise FileNotFoundError(source / 'setup.py')
if sys.version_info >= (3, 12) and (source / 'lightdock').is_dir() and not list((source / 'lightdock').rglob('*.pyx')):
    raise RuntimeError('LightDock sdist lacks .pyx sources: run python -m easy_ppi.fetch_lightdock_pyx before the Python 3.12 build')
os.chdir(source)
sys.path.insert(0, str(source))
sys.argv = ['setup.py', 'build_ext', '--compiler=mingw32', 'bdist_wheel']
print('Build compatibility flag: -Wno-error=incompatible-pointer-types', flush=True)
runpy.run_path(str(source / 'setup.py'), run_name='__main__')
