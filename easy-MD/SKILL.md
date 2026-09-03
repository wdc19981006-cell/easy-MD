---
name: easy-md
description: "Run an automated Windows protein-ligand molecular docking workflow with protein/gene normalization, human PDB selection, PubChem ligand acquisition, receptor and ligand preparation, repeated AutoDock Vina runs, optional redocking, PLIP-only interaction analysis, PLIP-native PyMOL visualization with a coordinate-driven fallback, validation, and structured packaging. Use for protein-small-molecule docking requests such as 'COX-2 with Aspirin' or 'Dock ROCK2 with Belumosudil', and for inspecting, resuming, validating, or packaging this workflow. Do not use for protein-protein docking, target prediction alone, molecular dynamics, or generic non-molecular uses of docking."
---

# easy-md

Run a reproducible protein-small-molecule docking workflow on Windows. Resolve
all paths relative to this file and read `config.yaml` before invoking scripts.

## Non-negotiable rules

- Use PLIP as the sole interaction classifier.
- Never let PyMOL or a custom distance cutoff infer hydrogen bonds,
  hydrophobic contacts, pi interactions, salt bridges, or halogen bonds.
- Prefer PLIP `VisualizerData` plus the official `PyMOLVisualizer`; use the
  generated coordinate-driven PML only when native rendering fails.
- Do not run Vina or redocking until the user has seen and explicitly approved
  the resolved protein, PDB, ligand, binding-site source, box, seeds,
  exhaustiveness, output directory, and expected runtime.
- Do not silently run blind docking. If no relevant native ligand exists,
  require experimentally supported `--center` and `--size` values.
- Disclose covalent-mechanism limitations before confirmation. Standard Vina
  does not model bond formation, reaction paths, or acyl transfer; a verified
  covalent residue may locate a pocket but does not make the protocol covalent.
- Do not overwrite a non-empty result directory or install dependencies without
  explicit permission.

## Preflight

Run:

```text
CONFIG=<skill-root>/config.yaml
PYTHON=<config paths.python>
PYTHON <skill-root>/scripts/preflight.py --config CONFIG --json
```

Continue core docking only when Python, Vina, OpenBabel, and core Python
modules pass. Report PLIP/PyMOL failures as stage-specific capability failures.

## MCP-first acquisition

When the `bio-server` tools are callable, use them first for PDB and PubChem
search/download and pass the resulting files to `pipeline.py` with
`--protein-pdb`, `--pdb-id`, and `--ligand-sdf`. Normalize common protein names
to a gene symbol before PDB search. Do not use the known-broken `get_pdb_info`
tool.

When MCP is unavailable, state that the bundled persistent REST adapters will
be used. Do not claim MCP usage when it was not callable.

## Prepare and pause

For the default name-only request, run:

```text
PYTHON <skill-root>/scripts/pipeline.py --config CONFIG prepare
  --protein PROTEIN
  --ligand LIGAND
  --output-root OUTPUT_ROOT
```

The command normalizes the identities, selects and downloads structures,
prepares receptor and ligand, builds the box, writes `run_status.json`, and
prints a `DOCKING CONFIRMATION` block. If protein identity is ambiguous, present
the candidates and wait for a gene selection.

Use optional overrides only when evidence or the user supplies them:

- `--pdb-id ID` or `--protein-pdb FILE`
- `--ligand-sdf FILE` or `--smiles TEXT`
- `--chain CHAIN`
- `--center X Y Z --size SX SY SZ`
- `--redock`

Stop after preparation and wait for explicit approval.

## Execute after approval

After the user approves the exact confirmation block, run:

```text
PYTHON <skill-root>/scripts/pipeline.py execute
  --run-dir RUN_DIR
  --confirm
```

The command performs repeated Vina runs, global pose ranking, optional
redocking, PLIP analysis, 2D rendering, PLIP-native PyMOL rendering, strict
fallback rendering when required, structured packaging, and validation.

Each stage records one of `PENDING`, `RUNNING`, `AWAITING_CONFIRMATION`,
`SUCCESS`, `FAILED`, or `SKIPPED`. A downstream failure must not replace a
previous stage's successful status or artifacts.

## Inspect or validate

```text
PYTHON <skill-root>/scripts/pipeline.py status --run-dir RUN_DIR
PYTHON <skill-root>/scripts/pipeline.py validate --run-dir RUN_DIR
```

Require non-empty docking artifacts, `interaction_summary.txt`,
`interaction_master.json`, PML, and PSE. Require the PLIP interaction count to
match the hidden `_easy_md_interaction_*` PyMOL marker count.

## References

- Read [references/environment-windows.md](references/environment-windows.md)
  for Windows, PowerShell, Git Bash, and path rules.
- Read [references/selection-and-validation.md](references/selection-and-validation.md)
  for structure ranking, box evidence, redocking, and interpretation limits.
- Read [references/workflow-commands.md](references/workflow-commands.md) for
  CLI details and compatibility entry points.
- Read [references/outputs.md](references/outputs.md) for the result contract
  and interaction schema.

## Completion report

Report the resolved identities, PDB selection evidence, binding-site source,
Vina parameters and scores, optional redocking status, PLIP counts, native or
fallback visualization source, per-stage status, validation result, and final
output path. Do not infer experimental affinity, efficacy, selectivity, or
mechanism from a docking score.
