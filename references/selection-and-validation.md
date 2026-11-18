# Structure selection and scientific validation

## Contents

- MCP-first lookup
- Protein identity
- PDB ranking
- Docking-box selection
- Redocking
- Interpretation limits

## MCP-first lookup

When `bio-server` is callable, prefer:

- `search_pdb`
- `download_pdb`
- `search_ligand`
- `get_smiles`
- `download_ligand_sdf`

Do not use `get_pdb_info` when it returns the known string/object error. Use
the RCSB core entry REST endpoint for that metadata. If MCP is unavailable,
state the fallback before querying REST.

## Protein identity

Normalize a common protein name to a gene symbol before PDB search. Present
multiple genes for user selection rather than choosing the first hit.

Do not treat:

- Protein family names as a unique gene.
- Protein-small-molecule docking as protein-protein docking.
- Target prediction as evidence that a ligand binds a target.

## PDB ranking

Prefer a structure that:

- Uses X-ray diffraction when appropriate for the configured workflow.
- Has good resolution.
- Contains a relevant co-crystal ligand.
- Represents the requested organism and construct.
- Has a usable binding site with limited missing residues or mutations.

`pdb_selection.py` reads scoring weights from `config.yaml`. Human origin,
missing-residue percentage, mutation status, and binding-site status require
enriched metadata under the `easy_md` field or explicit top-level organism
metadata. Do not award those points when the evidence is absent.

Do not assume every non-polymer entity is a ligand. Exclude waters, buffers,
common crystallization agents, and irrelevant ions.

## Docking-box selection

Use a relevant native ligand whenever possible. Report the residue name,
chain, residue number, center, and dimensions.

If no native ligand exists:

1. Ask for experimentally supported coordinates.
2. Offer a separately validated pocket-detection workflow.
3. Do not run blind docking in this Skill; require an explicit, evidence-backed
   center and per-axis size.

## Redocking

Redocking is a protocol check, not proof of biological activity.

- Below the configured validated threshold: report protocol geometry as
  reproduced.
- Between validated and borderline thresholds: report borderline.
- Above the borderline threshold: report not validated.
- If atom mapping or ligand chemistry prevents RMSD calculation, report that
  validation could not be computed.

## Interpretation limits

- Vina affinity is a scoring estimate, not an experimental binding free energy.
- Do not convert scores directly into Kd, Ki, or efficacy.
- A favorable pose does not establish selectivity, cellular activity, safety,
  or mechanism.
- Inspect protonation, tautomerism, stereochemistry, cofactors, metals, and
  construct mutations before publication.
- Standard Vina is noncovalent. A covalent target mechanism or modified-residue
  site marker does not model reaction chemistry or a covalent bond.
