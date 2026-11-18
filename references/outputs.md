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
    config.txt
    docking.log
    scores.csv
    runs/run1/{docked.pdbqt,log.txt}
    runs/run2/{docked.pdbqt,log.txt}
    runs/run3/{docked.pdbqt,log.txt}
  structure/
    original.pdb
    protein_clean.pdb
    receptor.pdbqt
    ligand.sdf
    ligand.pdbqt
    best_pose.pdb
    best_pose.pdbqt
    complex_for_visualization.pdb
  analysis/
    interaction_summary.txt
    interaction_summary.json
    interaction_master.json
    plip_report.xml
    visualizer_data.pkl
  visualization/
    interaction_visualization.pse
    interaction_visualization.pml
    visualization_metadata.json
```

Run the shared Vina configuration from `docking/`. It points to the prepared
inputs in `../structure/`; seed folders contain only the pose outputs and
logs. Structure preparation and PLIP use temporary working directories.
The run root also contains the `.easy-md-run` safety marker. The standard
workflow does not produce 2D interaction diagrams.

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
PLIP interaction count, actual saved distance segment count, PSE byte size, object inventory,
validation checks, and any native/fallback errors.

The session contains `Protein`, `Contact_Residues`, `Ligand`, and
`Interactions`. Each can be toggled independently. Empty contact and
interaction objects appear as empty groups when applicable. The PSE opens
directly; to rerender the fallback PML, set PyMOL's working directory to
`visualization/` and run `@interaction_visualization.pml`.
