#!/usr/bin/env python
"""Convert PLIP interaction objects into a coordinate-rich master schema."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


def _coord(value: Any) -> list[float] | None:
    if value is None:
        return None
    try:
        return [round(float(item), 4) for item in list(value)[:3]]
    except Exception:
        return None


def parse_pdb_atoms(path: str | Path) -> dict[int, dict[str, Any]]:
    atoms: dict[int, dict[str, Any]] = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
            coordinate = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        except ValueError:
            continue
        atoms[serial] = {
            "serial": serial,
            "atom_name": line[12:16].strip(),
            "residue_name": line[17:20].strip(),
            "chain": line[21:22].strip(),
            "residue_number": line[22:26].strip(),
            "coordinate": _coord(coordinate),
        }
    return atoms


def _ids(value: Any) -> list[int]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return [int(item) for item in value]
    return [int(value)]


def _atom_names(atom_map: dict[int, dict[str, Any]], values: Iterable[int]) -> list[str]:
    return [atom_map[value]["atom_name"] for value in values if value in atom_map]


def _base(
    kind: str,
    item: Any,
    atom_map: dict[int, dict[str, Any]],
    protein_ids: Iterable[int] = (),
    ligand_ids: Iterable[int] = (),
) -> dict[str, Any]:
    pids = [int(value) for value in protein_ids]
    lids = [int(value) for value in ligand_ids]
    pnames = _atom_names(atom_map, pids)
    lnames = _atom_names(atom_map, lids)
    return {
        "interaction_id": "",
        "interaction_type": kind,
        "protein_chain": str(getattr(item, "reschain", "") or ""),
        "protein_residue_name": str(getattr(item, "restype", "") or ""),
        "protein_residue_number": str(getattr(item, "resnr", "") or ""),
        "protein_atom_name": pnames[0] if len(pnames) == 1 else None,
        "ligand_atom_name": lnames[0] if len(lnames) == 1 else None,
        "protein_atom_names": pnames,
        "ligand_atom_names": lnames,
        "protein_atom_serials": pids,
        "ligand_atom_serials": lids,
        "distance": None,
        "angle": None,
        "protein_coordinate": None,
        "ligand_coordinate": None,
        "ring_center": {"protein": None, "ligand": None},
        "charge_center": {"protein": None, "ligand": None},
        "plip_source": True,
    }


def parse_interactions(interaction_set: Any, complex_pdb: str | Path) -> dict[str, Any]:
    atoms = parse_pdb_atoms(complex_pdb)
    records: list[dict[str, Any]] = []

    for item in interaction_set.hydrophobic_contacts:
        pids, lids = [item.bsatom_orig_idx], [item.ligatom_orig_idx]
        row = _base("hydrophobic", item, atoms, pids, lids)
        row.update({
            "distance": round(float(item.distance), 3),
            "protein_coordinate": _coord(item.bsatom.coords),
            "ligand_coordinate": _coord(item.ligatom.coords),
        })
        records.append(row)

    for item in interaction_set.hbonds_pdon + interaction_set.hbonds_ldon:
        if item.protisdon:
            protein_atom, ligand_atom = item.d, item.a
            pids, lids = [item.d_orig_idx], [item.a_orig_idx]
        else:
            protein_atom, ligand_atom = item.a, item.d
            pids, lids = [item.a_orig_idx], [item.d_orig_idx]
        row = _base("hydrogen_bond", item, atoms, pids, lids)
        row.update({
            "distance": round(float(item.distance_ad), 3),
            "angle": round(float(item.angle), 3),
            "protein_coordinate": _coord(protein_atom.coords),
            "ligand_coordinate": _coord(ligand_atom.coords),
            "donor_side": "protein" if item.protisdon else "ligand",
        })
        records.append(row)

    for item in interaction_set.saltbridge_lneg + interaction_set.saltbridge_pneg:
        if item.protispos:
            protein_group, ligand_group = item.positive, item.negative
        else:
            protein_group, ligand_group = item.negative, item.positive
        pids = _ids(protein_group.atoms_orig_idx)
        lids = _ids(ligand_group.atoms_orig_idx)
        row = _base("salt_bridge", item, atoms, pids, lids)
        pcenter, lcenter = _coord(protein_group.center), _coord(ligand_group.center)
        row.update({
            "distance": round(float(item.distance), 3),
            "protein_coordinate": pcenter,
            "ligand_coordinate": lcenter,
            "charge_center": {"protein": pcenter, "ligand": lcenter},
            "protein_charge": "positive" if item.protispos else "negative",
        })
        records.append(row)

    for item in interaction_set.pistacking:
        pids = _ids(item.proteinring.atoms_orig_idx)
        lids = _ids(item.ligandring.atoms_orig_idx)
        row = _base("pi_stacking", item, atoms, pids, lids)
        pcenter, lcenter = _coord(item.proteinring.center), _coord(item.ligandring.center)
        row.update({
            "distance": round(float(item.distance), 3),
            "angle": round(float(item.angle), 3),
            "protein_coordinate": pcenter,
            "ligand_coordinate": lcenter,
            "ring_center": {"protein": pcenter, "ligand": lcenter},
            "stacking_type": str(item.type),
        })
        records.append(row)

    for item in interaction_set.pication_paro + interaction_set.pication_laro:
        if item.protcharged:
            pids = _ids(item.charge.atoms_orig_idx)
            lids = _ids(item.ring.atoms_orig_idx)
            pcenter, lcenter = _coord(item.charge.center), _coord(item.ring.center)
            ring_center = {"protein": None, "ligand": lcenter}
            charge_center = {"protein": pcenter, "ligand": None}
        else:
            pids = _ids(item.ring.atoms_orig_idx)
            lids = _ids(item.charge.atoms_orig_idx)
            pcenter, lcenter = _coord(item.ring.center), _coord(item.charge.center)
            ring_center = {"protein": pcenter, "ligand": None}
            charge_center = {"protein": None, "ligand": lcenter}
        row = _base("pi_cation", item, atoms, pids, lids)
        row.update({
            "distance": round(float(item.distance), 3),
            "protein_coordinate": pcenter,
            "ligand_coordinate": lcenter,
            "ring_center": ring_center,
            "charge_center": charge_center,
            "protein_charged": bool(item.protcharged),
        })
        records.append(row)

    for item in interaction_set.halogen_bonds:
        pids, lids = [item.acc_orig_idx], [item.don_orig_idx]
        row = _base("halogen_bond", item, atoms, pids, lids)
        row.update({
            "distance": round(float(item.distance), 3),
            "angle": round(float(item.don_angle), 3),
            "acceptor_angle": round(float(item.acc_angle), 3),
            "protein_coordinate": _coord(item.acc.o.coords),
            "ligand_coordinate": _coord(item.don.x.coords),
        })
        records.append(row)

    for item in interaction_set.water_bridges:
        if item.protisdon:
            pids, lids = [item.d_orig_idx], [item.a_orig_idx]
            pcoord, lcoord = _coord(item.d.coords), _coord(item.a.coords)
        else:
            pids, lids = [item.a_orig_idx], [item.d_orig_idx]
            pcoord, lcoord = _coord(item.a.coords), _coord(item.d.coords)
        row = _base("water_bridge", item, atoms, pids, lids)
        row.update({
            "distance": round(float(min(item.distance_aw, item.distance_dw)), 3),
            "angle": round(float(item.w_angle), 3),
            "protein_coordinate": pcoord,
            "ligand_coordinate": lcoord,
            "water_coordinate": _coord(item.water.coords),
            "water_atom_serial": int(item.water_orig_idx),
        })
        records.append(row)

    for item in interaction_set.metal_complexes:
        pids = [item.target_orig_idx] if item.location == "protein" else []
        lids = [item.target_orig_idx] if item.location == "ligand" else []
        row = _base("metal_coordination", item, atoms, pids, lids)
        metal_coord = _coord(item.metal.coords)
        target_coord = _coord(item.target.atom.coords)
        row.update({
            "distance": round(float(item.distance), 3),
            "protein_coordinate": target_coord if item.location == "protein" else metal_coord,
            "ligand_coordinate": target_coord if item.location == "ligand" else metal_coord,
            "metal_type": str(item.metal_type),
            "metal_location": str(item.location),
        })
        records.append(row)

    for index, row in enumerate(records, 1):
        row["interaction_id"] = f"I{index:04d}"
        row["residue"] = (
            f"{row['protein_residue_name']}{row['protein_residue_number']}"
            f"{row['protein_chain']}"
        )
    counts = Counter(row["interaction_type"] for row in records)
    return {
        "schema_version": 1,
        "source": "PLIP",
        "interaction_count": len(records),
        "counts": dict(sorted(counts.items())),
        "interactions": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Pre-generated JSON master input")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    payload = json.loads(Path(args.input).read_text(encoding="utf-8"))
    Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
