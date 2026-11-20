"""Official PyPI sdists, verified before safe extraction into .ppi-tools."""
from pathlib import Path
import hashlib
import io
import json
import sys
import tarfile
import requests

name, version = sys.argv[1:]
root = Path(__file__).resolve().parent.parent / '.ppi-tools' / 'build'
root.mkdir(parents=True, exist_ok=True)
response = requests.get(f'https://pypi.org/pypi/{name}/{version}/json', timeout=45)
response.raise_for_status()
data = response.json()
package = next(x for x in data['urls'] if x['packagetype'] == 'sdist')
download = requests.get(package['url'], timeout=180)
download.raise_for_status()
digest = hashlib.sha256(download.content).hexdigest()
assert digest == package['digests']['sha256']
target = root / f'{name}-{version}'
target.mkdir(exist_ok=False)
(target / 'provenance.json').write_text(json.dumps({'url': package['url'], 'sha256': digest, 'requires_dist': data['info'].get('requires_dist'), 'license': data['info'].get('license')}, indent=2), encoding='utf-8')
with tarfile.open(fileobj=io.BytesIO(download.content)) as archive:
    archive.extractall(target, filter='data')
print(target)
