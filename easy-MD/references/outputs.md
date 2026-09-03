# Output and interaction contracts

## Contents

- Result tree
- Stage status
- Interaction master schema
- Visualization metadata

## Result tree

```text
<Protein>_<Ligand>/
  run_status.json
  metadata.json
  docking/
    receptor.pdbqt
    ligand.pdbqt
    config.txt
    best_pose.pdbqt
    docking.log
    scores.csv
    runs/
  structure/
    original.pdb
    protein_clean.pdb
    ligand.sdf
    best_pose.pdb
    complex_for_visualization.pdb
  analysis/
    interaction_summary.txt
    interaction_summary.json
    interaction_master.json
    plip_report.xml
    interaction_2D.png
    interaction_2D.svg
    interaction_2D.emf
  visualization/
    interaction_visualization.pse
    interaction_visualization.pml
    visualization_metadata.json
```

The run root may retain preparation logs, the `.easy-md-run` safety marker,
and intermediate files needed to resume a failed stage.

## Stage status

`run_status.json` stores independent `acquisition`, `preparation`, `docking`,
`interaction`, `visualization`, and `packaging` states. Preserve successful
stage artifacts when a later stage fails.

## Interaction master schema

`interaction_master.json` uses `source: PLIP` and a flat `interactions` list.
Every record includes:

```text
interaction_id
interaction_type
protein_chain
protein_residue_name
protein_residue_number
protein_atom_name
ligand_atom_name
protein_atom_names
ligand_atom_names
protein_atom_serials
ligand_atom_serials
distance
angle
protein_coordinate
ligand_coordinate
ring_center: {protein, ligand}
charge_center: {protein, ligand}
plip_source
```

Atom-pair interactions store atom names and coordinates. Group interactions
store member atom-name lists and PLIP ring/charge centers. Original atom
serials are retained only for traceability and are never the sole locator.

## Visualization metadata

`visualization_metadata.json` records the native attempt, selected source,
PLIP interaction count, PyMOL marker count, PSE byte size, object inventory,
validation checks, and any native/fallback errors.
