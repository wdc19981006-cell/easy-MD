# easy-MD repository guidance

## Development and release repositories

Public `wdc19981006-cell/easy-MD` is the only released/shared Skill.
Develop and test in Private `wdc19981006-cell/easy-md-test` (renamed from easy-md-skill).
This workstation's source root uses `origin` for the test repository and `release` for the public repository; default pushes go to origin.
Read CONTRIBUTING.md and RELEASE_PROCESS.md before source commits or publication.
Install versioned hooks with the configured Python and tools/install_hooks.py; inspect the staged tree and outgoing commits with tools/repository_guard.py.
Never commit structures, molecular renderings, sessions, trajectories, raw archives/logs, environments, credentials or personal configurations, including renamed/embedded structures or LFS pointers.
Keep redacted method/error records and CSV/JSON metrics under development/ only in the test repository. Real run data stays in the existing desktop result roots.
Promote only explicitly reviewed stable changes through source-linked commits and a public PR; preserve both scientific workflows and distinguish software tests from scientific validation.
Do not rewrite history, delete the historical result branch, force push, move results or rerun completed calculations during repository maintenance.

Read `README.md`, `NEW_USER_GUIDE.md`, and `SKILL.md` before configuring or running this workflow. The repository root is the skill root; resolve scripts relative to it.

For a new user, check the operating system and actual Python environment first. This workflow targets Windows. The executable and output paths in `config.yaml` are examples from the original workstation: configure real local paths before invoking the pipeline. Do not assume the original author's programs or bio-server MCP are installed. Use bundled REST acquisition adapters when the MCP tools are unavailable.

A request to download alone authorizes downloading and inspection. Install dependencies when the user has authorized installation; the README first-use prompt provides that authorization for an isolated workflow environment. Preserve existing configurations and result directories. Check `scripts/preflight.py`, actual PLIP/OpenBabel imports, and the configured PyMOL Python before promising the full analysis and visualization workflow.

For protein-small-molecule docking, follow the prepare / confirmation / execute sequence in `SKILL.md`. Dependency installation consent is not approval for Vina or redocking. Present the resolved inputs, site evidence, box, seeds, exhaustiveness, output directory, and expected runtime, then wait for approval of that concrete run. Report ambiguity or missing site evidence instead of silently choosing blind docking.

In the Vina workflow, classify interactions only with PLIP. Preserve the four-object PyMOL scene and its geometry validation. Keep output data, virtual environments, caches and credentials out of source commits. Report actual stage failures and validation results. Do not promise a fixed six-minute runtime or claim a fresh-machine installation test that has not been performed.

# easy-MD / easy-PPI project instructions

## Result locations

- This repository is for source code. Save all actual molecular docking and protein interaction results under `D:/Desktop/codex-result`.
- Protein–small-molecule (AutoDock Vina) results: `D:/Desktop/codex-result/easy_MD`; configured by root `config.yaml` → `paths.output_root`.
- Protein–protein (easy-PPI) results: `D:/Desktop/codex-result/easy_PPI`; configured by `easy_ppi/config.json` → `paths.output_root`.
- Do not create new calculation results in project-local `results/` or `ppi_results/`. These folders contain relocation notices only after the 2026-11-19 migration.
- Put development/compatibility trials in `easy_PPI/_archive`, and scientific benchmarks in `easy_PPI/benchmarks`. Group final outputs by protein pair and date; keep validation data separate from final reports and views.
- Explicitly requested output directories take precedence. Never overwrite a nonempty result directory.

## Historical results

The local-only `references/result-location-map.json` (not distributed on GitHub) and the desktop result directory's migration index map historical paths to current locations. Original logs, scored files and manifests retain their historical text and numerical contents. Use the new location map and new reopening scripts; do not rerun archived workstation investigation scripts with their old hard-coded paths.

## Workflow boundaries

Keep the Vina workflow and the experimental-structure easy-PPI workflow independent. Use the root `SKILL.md` to select the workflow and read only its relevant reference. Preserve scientific failure reports and reference-free DockQ skips. Unit-test temporary files may use the test runner's temporary directories.

For protein-protein tasks, read `easy_ppi/README.zh-CN.md`, `easy_ppi/ARCHITECTURE.zh-CN.md`, and `easy_ppi/VALIDATION.zh-CN.md`. Use the independent environment and Arpeggio backend. The root `SKILL.md` routes both task types: read `references/vina-workflow.md` for protein-small-molecule work and `references/ppi-workflow.md` for protein-protein work. Vina-specific preparation, PLIP classification and confirmation commands do not apply to PPI analysis. Do not claim the failed 1BRS redocking benchmark passed. On other machines configure both output roots to real absolute paths first.
