"""Single-process upstream DockQ API evaluation, no reference-driven ranking."""
from pathlib import Path
import json
import sys

def main(request_path):
    request = json.loads(Path(request_path).read_text(encoding='utf-8'))
    script = request.get('dockq_script')
    if script:
        # The unmodified upstream script provides its official operations_nocy
        # fallback on this path. Keep it a separate scientific process.
        sys.path.insert(0, str(Path(script).parent))
        import importlib.util
        spec = importlib.util.spec_from_file_location('upstream_dockq', script)
        api = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(api)
    else:
        from DockQ import DockQ as api
    native = api.load_PDB(request['reference'])
    results = []
    output = Path(request['output'])
    # Preserve each completed candidate even if a later one fails.
    with output.open('x', encoding='utf-8') as handle:
        for index, row in enumerate(request['candidates']):
            model = api.load_PDB(row['path'])
            result, total = api.run_on_all_native_interfaces(model, native, chain_map={'A': 'A', 'B': 'B'})
            if not result or len(result) != 1:
                raise ValueError('Expected exactly one DockQ reference interface')
            interface = next(iter(result.values()))
            # numpy scalars are converted explicitly; no scientific formulas
            # or upstream geometric calculations are replaced.
            metrics = {key: float(interface[key]) for key in ('DockQ', 'iRMSD', 'LRMSD', 'fnat')}
            record = {'rank': row['rank'], 'path': row['path'], **metrics}
            handle.write(json.dumps(record, allow_nan=False) + '\n')
            handle.flush()
            if (index + 1) % 100 == 0:
                print(f'DockQ {index+1}/{len(request["candidates"])}', flush=True)

if __name__ == '__main__':
    main(sys.argv[1])
