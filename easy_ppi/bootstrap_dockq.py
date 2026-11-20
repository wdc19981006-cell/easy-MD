"""Acquire the unmodified upstream DockQ script fallback with verified hashes.

This is optional: normal installations use the DockQ console entry point.
The official source includes operations_nocy specifically for script execution.
"""
from pathlib import Path
import hashlib
import io
import json
import sys
import tarfile

import requests

def main():
    root = Path(__file__).resolve().parent.parent / '.ppi-tools' / 'dockq-2.1.3'
    root.mkdir(parents=True, exist_ok=False)
    metadata = requests.get('https://pypi.org/pypi/DockQ/2.1.3/json', timeout=45)
    metadata.raise_for_status()
    data = metadata.json()
    package = next(x for x in data['urls'] if x['packagetype'] == 'sdist')
    response = requests.get(package['url'], timeout=90)
    response.raise_for_status()
    digest = hashlib.sha256(response.content).hexdigest()
    if digest != package['digests']['sha256']:
        raise ValueError('DockQ sdist SHA256 mismatch')
    (root / 'source.tar.gz').write_bytes(response.content)
    with tarfile.open(fileobj=io.BytesIO(response.content)) as archive:
        for member in archive.getmembers():
            # Copy regular Python sources/LICENSE only; never extract symlinks
            # or arbitrary archive paths to the host filesystem.
            if not member.isfile():
                continue
            if '/src/DockQ/' in member.name and member.name.endswith('.py') or member.name.endswith('/LICENSE'):
                name = Path(member.name).name
                content = archive.extractfile(member).read()
                (root / name).write_bytes(content)
    (root / 'provenance.json').write_text(json.dumps({'distribution': 'DockQ', 'version': '2.1.3', 'url': package['url'],
                'sha256': digest, 'modifications': 'none; official script-mode Python fallback', 'requires_dist': data['info']['requires_dist']}, indent=2), encoding='utf-8')
    print(root / 'DockQ.py')

if __name__ == '__main__':
    main()
