---
name: easy-md
description: "Run protein-small-molecule docking with AutoDock Vina/PLIP, or protein-protein interface analysis and LightDock docking from experimental structures with easy-PPI. Use for tasks such as ROCK2 with Belumosudil or NRF2 with KEAP1, and for inspecting, validating or packaging either workflow. Excludes sequence-to-structure prediction, molecular dynamics and target prediction alone."
---

# easy-md

One skill with two independent scientific workflows. Resolve `<skill-root>`
relative to this `SKILL.md`; choose the workflow before checking dependencies,
acquiring structures or running commands.

## Select the workflow

| Request | Program | Required workflow reference |
| --- | --- | --- |
| Protein plus small molecule, e.g. ROCK2–Belumosudil | AutoDock Vina, `scripts/` | [Vina workflow](references/vina-workflow.md) |
| Two proteins or a protein plus a protein-derived peptide, e.g. NRF2–KEAP1 | easy-PPI, `easy_ppi/` | [PPI workflow](references/ppi-workflow.md) |
| Inspect, resume, validate or package an existing run | Use the program that created the run | Read its workflow reference |

Determine partner identity and molecular type; two names alone do not establish
that both are proteins. Ask only when identity or task type remains ambiguous.
For an existing run, use its metadata: Vina uses `run_status.json`, while PPI
uses `run.json`. Read only the relevant workflow reference and its linked
details. Do not run a scientific calculation merely to answer a status request.

## Shared constraints

- Keep program code, dependency environments, preparation methods, interaction
  classifiers and configurations independent. Vina uses PLIP; PPI uses PDBe
  Arpeggio. Do not apply Vina ligand preparation or its confirmation CLI to PPI.
- Prefer callable bio-server tools for PDB/PubChem acquisition. Normalize protein
  names to gene identity first; do not use the known-broken `get_pdb_info` tool.
  When unavailable, disclose use of the selected program's REST adapters.
- Use real local executable paths. Install missing dependencies within the
  user's authorization and preserve existing environments and configurations.
- Save actual calculations outside the source repository; never overwrite an
  existing nonempty result directory. Preserve raw inputs, numbering maps,
  parameters, candidates, scores, logs and validation evidence.
- Report stage failures and scientific limits. Program exit status, file
  integrity and scientific validity are separate conclusions. Never fabricate
  reference-based validation or interpret docking scores as experimental proof.

## Configuration and result locations

| Program | Configuration | Original workstation default |
| --- | --- | --- |
| Vina | `config.yaml` → `paths.output_root` | `D:/Desktop/codex-result/easy_MD` |
| easy-PPI | `easy_ppi/config.json` → `paths.output_root` | `D:/Desktop/codex-result/easy_PPI` |

On a new machine configure actual absolute output paths before running. Explicit
user-specified paths take precedence. PPI also supports `EASY_PPI_OUTPUT_ROOT`
and `EASY_PPI_CONFIG`; an invalid default must fail without falling back to
project-local `ppi_results/`. Use named partner/date directories for final
outputs, `easy_PPI/benchmarks` for benchmarks and `easy_PPI/_archive` for trials;
these are organization conventions, selected via explicit output arguments.
See [result storage](references/result-storage.md).
