# Command reference

## Primary pipeline

Use the configured absolute Python path and Windows-native arguments.

```text
pipeline.py --config CONFIG prepare
  --protein NAME --ligand NAME [--output-root DIR]
  [--pdb-id ID | --protein-pdb FILE]
  [--ligand-sdf FILE | --smiles TEXT]
  [--chain ID]
  [--center X Y Z --size SX SY SZ]
  [--redock]

pipeline.py execute --run-dir DIR --confirm
pipeline.py status --run-dir DIR
pipeline.py validate --run-dir DIR
repackage.py OLD_RUN [OLD_RUN ...] [--output-root DIR]
```

`prepare` uses `paths.output_root` from the configuration unless `--output-root`
is supplied. It performs no Vina run. `execute` rejects calls without `--confirm`.

`repackage.py` creates new result directories from trusted local runs that
already finished docking and PLIP. It uses the configured output root unless
overridden, refuses to overwrite results, and does not rerun either calculation.

## Subsystem CLIs

```text
protein/normalize.py NAME --species human --json
protein/uniprot_search.py SYMBOL --json
protein/pdb_search.py SYMBOL --rows 15 --json
protein/pdb_rank.py METADATA.json --config CONFIG
ligand/ligand_search.py NAME [--download FILE]
docking/pose_rank.py --work-dir DIR --seeds 42,2024,8888 --obabel EXE
analysis/plip_analysis.py ...
visualization/plip_pymol_native.py ...
visualization/pymol_generator.py ...
visualization/visualization_validate.py ...
```

All subsystem scripts support `--help`.

## Compatibility entry points

The existing `protein_search.py`, `pdb_selection.py`, `protein_prepare.py`,
`ligand_prepare.py`, `box_prepare.py`, `docking_run.py`, `redocking.py`,
`interaction_analysis.py`, the standalone `interaction_2d.py`, and `package_results.py`
remain available. Prefer `pipeline.py` for new runs because it owns stage
status, the structured output contract, and visualization validation.
