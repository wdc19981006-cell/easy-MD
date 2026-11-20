from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

from .audit import write_json
from ._inspect_environment import inspect


def main():
    parser = argparse.ArgumentParser(description='easy-PPI: experimental protein structures only; isolated from Vina')
    sub = parser.add_subparsers(dest='command', required=True)
    doctor = sub.add_parser('doctor')
    doctor.add_argument('--output')
    validator = sub.add_parser('validate')
    validator.add_argument('run_directory')
    for name in ('run', 'benchmark'):
        p = sub.add_parser(name)
        p.add_argument('--output', help='Explicit new run directory; default root comes from easy_ppi/config.json (desktop codex-result/easy_PPI)')
        p.add_argument('--pymol-python', help='Interpreter containing PyMOL; no system changes')
        p.add_argument('--seeds', type=int, nargs='+', default=[42, 2024, 8888])
        p.add_argument('--swarms', type=int, default=400)
        p.add_argument('--glowworms', type=int, default=300)
        p.add_argument('--steps', type=int, default=100)
        p.add_argument('--cores', type=int, default=1)
        p.add_argument('--scoring', choices=['dfire', 'fastdfire'], default='fastdfire', help='Official LightDock potential; scores from different functions are not mixed')
        if name == 'run':
            p.add_argument('--protein-a', required=True)
            p.add_argument('--protein-b', required=True)
            p.add_argument('--chain-a')
            p.add_argument('--chain-b')
            p.add_argument('--assembly-id')
            p.add_argument('--species-a', help='Source taxonomy ID or exact scientific name')
            p.add_argument('--species-b', help='Source taxonomy ID or exact scientific name')
            p.add_argument('--reference', help='Experimental PDB ID or local complex; optional, never fabricated')
            p.add_argument('--reference-chain-a')
            p.add_argument('--reference-chain-b')
    inspect_cmd = sub.add_parser('inspect')
    inspect_cmd.add_argument('file')
    inspect_cmd.add_argument('--output')
    args = parser.parse_args()
    if args.command == 'validate':
        from .validation import validate_artifacts
        report = validate_artifacts(args.run_directory)
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return 0 if report['artifact_integrity'] == 'PASS' else 2
    if args.command == 'doctor':
        data = inspect()
        if args.output:
            write_json(args.output, data)
        print(json.dumps(data, ensure_ascii=False, indent=2))
        return 0
    if args.command == 'inspect':
        from .structures import inspect_structure
        data = inspect_structure(args.file)
        if args.output:
            write_json(args.output, data)
        print(json.dumps(data, ensure_ascii=False, indent=2))
        return 0
    from .pipeline import benchmark, workflow
    kwargs = dict(output=args.output, pymol_python=args.pymol_python, seeds=args.seeds, swarms=args.swarms,
                  glowworms=args.glowworms, steps=args.steps, cores=args.cores, scoring=args.scoring)
    if args.command == 'benchmark':
        run = benchmark(**kwargs)
    else:
        run = workflow(args.protein_a, args.protein_b, args.chain_a, args.chain_b, assembly_id=args.assembly_id,
                       species_a=args.species_a, species_b=args.species_b, reference=args.reference,
                       reference_chain_a=args.reference_chain_a, reference_chain_b=args.reference_chain_b, **kwargs)
    print(json.dumps({'run': str(run.path), 'stages': run.state['stages'], 'scientific_validation': run.state['scientific_validation']}, ensure_ascii=False, indent=2))
    if any(x['status'] in ('FAILED', 'BLOCKED') for x in run.state['stages'].values()):
        return 2
    if args.command == 'benchmark' and run.state['scientific_validation']['status'] != 'PASS':
        return 3
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
