# Protein–protein workflow (easy-PPI)

Use this mode for two proteins or an experimentally structured protein-derived
peptide. Run from the skill root with the independent PPI interpreter. Read
[usage and environment setup](../easy_ppi/README.zh-CN.md),
[architecture and acceptance criteria](../easy_ppi/ARCHITECTURE.zh-CN.md), and
[actual validation results](../easy_ppi/VALIDATION.zh-CN.md) as needed for the run.
Do not invoke the Vina/PLIP program for a PPI task.

## Prepare and route

1. Check the actual environment with `PPI_PYTHON -m easy_ppi doctor`; verify the
   needed backends in the independent environment and the chosen PyMOL Python.
   Configure `easy_ppi/config.json` or its documented environment overrides.
2. Resolve protein/gene identities and requested species. Prefer callable
   bio-server acquisition tools; otherwise disclose use of bundled REST adapters.
   Inputs may be names, PDB IDs or local PDB/mmCIF files. Unknown provenance is
   not proof of an experimental structure; sequence-prediction models are excluded.
3. Inspect experimental method, chain identity, biological assembly, resolution,
   missing residues/atoms and construct coverage. Expose ambiguities rather than
   silently choosing an unrelated chain, species or fusion construct. Preserve
   original files and auth/label chain and residue mappings.
4. Prefer an available experimental biological complex and analyze its actual
   interface directly. If only independent experimental structures are available,
   the program uses LightDock. A failed complex search means unknown availability,
   not confirmed absence. Reuse session authorization for an already specified run.

## Commands

`PPI_PYTHON` means the actual configured PPI interpreter, not the Vina Python or a
Windows Store placeholder. Run module commands from `<skill-root>`. `NEW_RUN` is
a new absolute output directory. The PyMOL interpreter may be separate.

```text
PPI_PYTHON -m easy_ppi doctor
PPI_PYTHON -m easy_ppi inspect EXPERIMENTAL_FILE
PPI_PYTHON -m easy_ppi run --protein-a INPUT_A --protein-b INPUT_B
  --output NEW_RUN --pymol-python PYMOL_PYTHON
PPI_PYTHON -m easy_ppi validate RUN_DIR
```

Pass chain/assembly/species overrides when established by the structure review.
For example, the existing NRF2–KEAP1 complex is `--protein-a 2FLU --chain-a P
--protein-b 2FLU --chain-b X --assembly-id 1`. Explain that it contains an NRF2
short peptide and a KEAP1 Kelch domain, rather than full-length proteins. Optional
`--reference` and reference-chain arguments are for an actual experimental
reference, never an invented structure. See module usage for docking parameters.

For the two-stage 1BRS test use `PPI_PYTHON -m easy_ppi benchmark` with an explicit
new output under the configured PPI root's `benchmarks` directory. This command
performs experimental-interface analysis and bound redocking; it is not needed
for routine direct analysis of an available experimental complex.

## Analyze, verify and deliver

- Gemmi/pdb-tools prepare structures and retain numbering maps. Report missing
  coordinates; do not silently complete loops or sidechains with predictions.
- PDBe Arpeggio is the interaction classifier. Ionic contacts are salt-bridge
  candidates; overlapping contact labels cannot be summed as independent forces.
- FreeSASA reports SASA(A), SASA(B), SASA(AB), total buried SASA `A+B-AB`, and
  interface area `(A+B-AB)/2`, with units and matching atom/parameter selections.
- PRODIGY affinity is an empirical estimate for the actual modeled fragments.
  Severe clashes or incomplete analyses must not be presented as valid affinity.
- Use DockQ only with an experimental reference; without one record SKIPPED and
  null. Reference comparison never determines LightDock optimization or ranking.
- Preserve every final LightDock candidate, native files, parameters, seeds,
  scores, clustering and logs. Reference-best is a sampling diagnostic, not the
  predicted score-top1. Validate geometry and chemistry alongside scores.
- PyMOL saves Protein_A, Protein_B, Interface_Residues and Interactions as separate
  top-level items; Arpeggio-derived lines use dashes. Reopen the PSE and verify
  atom counts and saved line endpoints. Deliver PSE/PNG plus CSV/JSON/PDB and logs.
- Require scientific controls appropriate to the run, such as separated-chain
  negative controls and SASA invariance. `validate` checks artifact integrity;
  inspect the independent scientific-validation status as well.

The recorded 1BRS interface analysis passed, but its 7,500-candidate redocking
failed the DockQ acceptance threshold. Preserve that FAIL; do not claim general
docking validation. Report source, species, chains, assembly, construct limits,
contact labels, area convention, empirical-affinity caveat, reference status,
stage/validation results and the actual output directory. Raw data remain local;
source publishing is separate from scientific result delivery.
