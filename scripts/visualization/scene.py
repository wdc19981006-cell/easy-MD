"""Build the four-object PyMOL session from PLIP interactions."""

import math


SCENE_OBJECTS = ("Protein", "Contact_Residues", "Ligand", "Interactions")
DASH_SETTINGS = {"dash_gap": 0.3, "dash_length": 0.15, "dash_width": 2.5}


def interaction_segments(master):
    segments = []
    for row in master["interactions"]:
        endpoints = [
            row.get(f"{side}_coordinate")
            or (row.get("ring_center") or {}).get(side)
            or (row.get("charge_center") or {}).get(side)
            for side in ("protein", "ligand")
        ]
        if row["interaction_type"] == "water_bridge":
            endpoints.insert(1, row.get("water_coordinate"))
        if not all(point and len(point) == 3 and all(map(math.isfinite, point)) for point in endpoints):
            raise ValueError(f"Missing or invalid PLIP coordinates for {row['interaction_id']}")
        segments.extend(zip(endpoints, endpoints[1:]))
    return segments


def measurement_segments(cmd, name):
    """Read saved distance geometry (PyMOL's ObjectDist session representation)."""
    if cmd.get_type(name) != "object:measurement":
        return []
    entries = cmd.get_session(name, partial=1, compress=0)["names"]
    entry = next(entry for entry in entries if entry and entry[0] == name)
    segments = []
    for state in entry[5][2]:
        if state and state[0]:
            coords = state[1]
            if len(coords) != 3 * state[0] or len(coords) % 6:
                raise ValueError(f"Invalid saved PyMOL distance geometry: {name}")
            segments.extend((coords[i:i + 3], coords[i + 3:i + 6])
                            for i in range(0, len(coords), 6))
    return segments


def match_segments(expected, actual, tolerance=0.02):
    """Match without relying on ordering or direction in the PyMOL session."""
    remaining = list(actual)
    if len(expected) != len(remaining):
        return False
    for start, end in expected:
        for i, (a, b) in enumerate(remaining):
            if (max(math.dist(start, a), math.dist(end, b)) <= tolerance
                    or max(math.dist(start, b), math.dist(end, a)) <= tolerance):
                remaining.pop(i)
                break
        else:
            return False
    return True


def contact_residues(master):
    return {(str(row.get("protein_chain") or ""),
             str(row["protein_residue_number"]),
             str(row.get("protein_residue_name") or ""))
            for row in master["interactions"] if row.get("protein_residue_number")
            and row.get("metal_location") not in {"ligand", "water"}}


def scene_atom_indices(cmd, master, source):
    site = master["ligand_site"]
    ligand, contacts = [], []
    residues = contact_residues(master)
    for atom in cmd.get_model(source).atom:
        key = atom.chain, atom.resi, atom.resn
        if key == (site["chain"], str(site["residue_number"]), site["resname"]):
            ligand.append(atom.index)
        elif key in residues or (*key[:2], "") in residues:
            contacts.append(atom.index)
    return ligand, contacts


def create_scene(cmd, master, complex_path, segments):
    """Create only three molecular objects and one toggle for measurements."""
    cmd.reinitialize()
    cmd.load(str(complex_path), "_complex", state=1)
    ligand, contacts = scene_atom_indices(cmd, master, "_complex")
    if not ligand:
        raise ValueError("The PLIP ligand site is absent from the complex")
    cmd.select_list("_ligand", "_complex", ligand, mode="index")
    cmd.create("Protein", "_complex and not _ligand", 1, 1)
    if contacts:
        cmd.select_list("_contacts", "_complex", contacts, mode="index")
        cmd.create("Contact_Residues", "_contacts", 1, 1)
        cmd.delete("_contacts")
    elif contact_residues(master):
        raise ValueError("No PLIP contact residues could be found in the complex")
    else:
        cmd.group("Contact_Residues")
    cmd.create("Ligand", "_ligand", 1, 1)
    cmd.delete("_ligand")
    cmd.delete("_complex")
    cmd.hide("everything", "all")
    cmd.show("cartoon", "Protein")
    cmd.show("sticks", "Contact_Residues or Ligand")
    cmd.color("gray80", "Protein")
    cmd.color("yellow", "Contact_Residues")
    cmd.color("cyan", "Ligand")
    cmd.util.cnc("Contact_Residues or Ligand")
    for start, end in segments:
        cmd.pseudoatom("_endpoint_a", pos=start)
        cmd.pseudoatom("_endpoint_b", pos=end)
        cmd.distance("Interactions", "_endpoint_a", "_endpoint_b", mode=0,
                     cutoff=math.dist(start, end) + 0.1, label=0, reset=0, state=1)
        cmd.delete("_endpoint_a")
        cmd.delete("_endpoint_b")
    if segments:
        cmd.set("dash_color", "gray30", "Interactions")
        cmd.set("dash_radius", 0, "Interactions")
        for key, value in DASH_SETTINGS.items():
            cmd.set(key, value, "Interactions")
        cmd.show("dashes", "Interactions")
        cmd.hide("labels", "Interactions")
    else:
        cmd.group("Interactions")
    for key, value in DASH_SETTINGS.items():
        cmd.set(key, value)
    cmd.bg_color("white")
    cmd.orient("Ligand")
    cmd.zoom("Ligand or Contact_Residues", 5)
