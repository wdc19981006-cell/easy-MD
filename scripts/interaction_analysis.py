#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Protein-Ligand Interaction Analysis — PLIP + 2D diagram + PyMOL visualization.

Analyzes the best docking pose for non-covalent interactions using the
Protein-Ligand Interaction Profiler (PLIP). Generates human-readable and
machine-readable reports, 2D interaction diagram, and PyMOL visualization.

IMPORTANT:
- Original docking files are NEVER modified.
- PLIP is the sole source of interaction detection.
- The 2D diagram reads from PLIP JSON — it does NOT re-compute interactions.
- PyMOL is used for 3D visualization only.
- Publication figures must use best_pose.pdb from Vina.

Usage:
    python interaction_analysis.py \
        --protein protein_clean.pdb \
        --best-pose-pdbqt run1/docked.pdbqt \
        --ligand-sdf ligand.sdf \
        --protein-name "CD38" --ligand-name "ADP" \
        --output-dir interaction_analysis \
        --obabel "D:/Interaction/OpenBabel-3.1.1/obabel.exe" \
        --pymol "D:/Interaction/PyMol/PyMOLWin.exe"

Output: interaction_analysis/
        ├── complex_for_plip.pdb
        ├── plip_report.xml
        ├── interaction_summary.txt
        ├── interaction_summary.json
        ├── interaction_2D.png
        ├── interaction_2D.svg
        ├── interaction_2D.emf
        ├── pymol_interaction.pml
        └── pymol_interaction.pse
"""

import os, sys, re, json, argparse, subprocess, shutil, time
from pathlib import Path

from common import configured_path, enter_work_dir, load_config


# ---------------------------------------------------------------------------
# Residue name generation from ligand name
# ---------------------------------------------------------------------------

def generate_resname(ligand_name):
    """Generate a 3-letter PDB residue name from the ligand name.

    Avoids collision with standard amino acid 3-letter codes (e.g. Aspirin
    would produce ASP which conflicts with aspartic acid).

    Rules:
    - Strip non-alphabetic chars, uppercase
    - Try first 3 letters; if it matches a forbidden code, try shifting
    - Fall back to 'LIG' if nothing works

    Examples: Aspirin->ASR, Belumosudil->BEL, ADP->ADP
    """
    # Standard 3-letter codes that PLIP would treat as protein residues
    _FORBIDDEN = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HOH", "WAT", "DOD", "ACE", "NME", "UNK", "PCA",
        "HYP", "CSD", "CSO", "CSW", "CSS", "CME", "MLY", "MLZ", "MSE",
        "SEP", "TPO", "PTR", "NLE", "NVA", "ABA", "AIB", "ORN", "CIT",
        "SAR", "BMT", "DAL", "DAR", "DAS", "DCY", "DGL", "DGN", "DHI",
        "DIL", "DLE", "DLY", "DME", "DPN", "DPR", "DSG", "DSN", "DSP",
        "DTH", "DTR", "DTY", "DVA",
    }

    clean = re.sub(r"[^a-zA-Z]", "", ligand_name).upper()
    if len(clean) == 0:
        return "LIG"

    # Try sliding windows of 3 chars until we find a non-colliding code
    # Prefer earlier windows (start of name)
    for start in range(len(clean) - 2):
        candidate = clean[start:start + 3]
        if candidate not in _FORBIDDEN:
            return candidate

    # If all 3-char windows collide, try skipping letters
    # (take characters at positions 0, 2, 4 or 0, 3, 6, etc.)
    for skip in range(2, len(clean)):
        indices = list(range(0, min(len(clean), skip * 3), skip))[:3]
        if len(indices) == 3:
            candidate = "".join(clean[i] for i in indices)
            if candidate not in _FORBIDDEN:
                return candidate

    # Last resort: pad with digits
    if len(clean) >= 2:
        candidate = clean[:2] + "1"
        if candidate not in _FORBIDDEN:
            return candidate

    return "LIG"


# ---------------------------------------------------------------------------
# Step 1: Convert best docking pose from PDBQT to PDB
# ---------------------------------------------------------------------------

def convert_best_pose(best_pose_pdbqt, obabel_path, work_dir):
    """Extract mode 1 from Vina docked PDBQT and convert to PDB via OpenBabel.

    The input PDBQT contains all 9 poses as MODEL/ENDMDL.  OpenBabel -m
    splits them; we keep only the first (mode 1 / best pose).
    """
    base = os.path.join(work_dir, "_best_pose_plip")
    # Use absolute path for -O (cwd is tmp_docking, not work_dir)
    result = subprocess.run(
        [obabel_path, best_pose_pdbqt, "-O", base + ".pdb", "-m"],
        capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: obabel PDBQT->PDB conversion failed: {result.stderr}")
        return None
    first_pose = base + "1.pdb"
    if not os.path.exists(first_pose):
        print(f"ERROR: obabel -m did not produce expected output: {first_pose}")
        print(f"  stdout: {result.stdout}")
        print(f"  stderr: {result.stderr}")
        return None
    natoms = sum(1 for l in open(first_pose) if l.startswith(("ATOM","HETATM")))
    print(f"  Best pose extracted: {natoms} atoms -> {os.path.basename(first_pose)}")
    return first_pose


# ---------------------------------------------------------------------------
# Step 2: Build protein + ligand complex for PLIP
# ---------------------------------------------------------------------------

def build_complex(protein_pdb, ligand_pdb, output_path,
                  lig_resname, lig_chain="L"):
    """Merge protein_clean.pdb and docked best-pose ligand into one PDB.

    Protein atoms are kept exactly as-is (ATOM/HETATM + TER).
    Ligand atoms are extracted from the docked pose PDB, assigned a distinct
    residue name (from ligand name, e.g. Aspirin->ASP) and chain so PLIP
    can unambiguously identify the ligand.

    Returns the number of protein atoms and ligand atoms written.
    """
    # --- read protein ---
    protein_lines = []
    with open(protein_pdb, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TER")):
                protein_lines.append(line.rstrip())

    # --- read ligand (best pose from obabel conversion) ---
    ligand_atoms = []
    with open(ligand_pdb, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                ligand_atoms.append(line.rstrip())

    if not ligand_atoms:
        print("ERROR: No ligand atoms found in docked pose PDB")
        return 0, 0

    # --- re-assign residue name & chain so PLIP recognizes the ligand ---
    # PDB column layout (1-indexed):
    #   1-6: record   7-11: serial   13-16: atname   17: altloc
    #   18-20: resname   22: chain   23-26: resseq
    # In 0-indexed Python string:  line[17:20]=resname  line[21]=chain
    new_ligand = []
    for line in ligand_atoms:
        before = line[:17]
        after_chain = line[22:] if len(line) > 22 else ""
        new_line = before + f"{lig_resname:3s}" + line[20:21] + lig_chain + after_chain
        new_ligand.append(new_line)

    # --- write merged complex ---
    with open(output_path, "w") as f:
        for line in protein_lines:
            f.write(line + "\n")
        f.write("TER\n")
        for line in new_ligand:
            f.write(line + "\n")
        f.write("TER\n")
        f.write("END\n")

    print(f"  Complex written: {len(protein_lines)} protein atoms + "
          f"{len(new_ligand)} ligand atoms -> {output_path}")
    return len(protein_lines), len(new_ligand)


# ---------------------------------------------------------------------------
# Step 3: Run PLIP analysis
# ---------------------------------------------------------------------------

def _extract_reslabel(residue_obj):
    """Build a human-readable residue label like 'ASP123A' from a PLIP object."""
    chain = getattr(residue_obj, "reschain", "")
    return f"{residue_obj.restype}{residue_obj.resnr}{chain}"


def run_plip(complex_pdb):
    """Run PLIP on a protein-ligand complex PDB.

    Returns a dictionary of interaction lists keyed by interaction type:
        hbonds, hydrophobic, pistacking, saltbridges, pication, halogen
    """
    from plip.structure.preparation import PDBComplex

    mol = PDBComplex()
    mol.load_pdb(complex_pdb)

    # auto-detect ligands — our LIG, chain L should be picked up
    mol.analyze()

    ligand_keys = list(mol.interaction_sets.keys())
    if not ligand_keys:
        print("WARNING: PLIP found no ligands / no interaction sets")
        return {}

    # Use the interaction set for our docked ligand (last one added)
    lig_key = ligand_keys[-1]
    iset = mol.interaction_sets[lig_key]
    print(f"  PLIP analyzed ligand: {lig_key} ({len(ligand_keys)} ligands found)")

    interactions = {}

    # --- hydrogen bonds (PLIP v3: hbonds_ldon + hbonds_pdon) ---
    interactions["hbonds"] = []
    for hb in iset.hbonds_ldon + iset.hbonds_pdon:
        interactions["hbonds"].append({
            "residue": _extract_reslabel(hb),
            "residue_chain": hb.reschain,
            "atom_donor": hb.d_orig_idx,
            "atom_acceptor": hb.a_orig_idx,
            "distance_hb": round(hb.distance_ah, 2),
            "distance_da": round(hb.distance_ad, 2),
            "angle": round(hb.angle, 1),
            "type": hb.type,
        })

    # --- hydrophobic contacts ---
    interactions["hydrophobic"] = []
    for hc in iset.hydrophobic_contacts:
        interactions["hydrophobic"].append({
            "residue": _extract_reslabel(hc),
            "residue_chain": hc.reschain,
            "distance": round(hc.distance, 2),
        })

    # --- pi-pi stacking ---
    interactions["pistacking"] = []
    for ps in iset.pistacking:
        interactions["pistacking"].append({
            "residue": _extract_reslabel(ps),
            "residue_chain": ps.reschain,
            "distance": round(ps.distance, 2),
            "angle": round(ps.angle, 1),
            "type": ps.type,
        })

    # --- salt bridges (PLIP v3: saltbridge_lneg + saltbridge_pneg) ---
    interactions["saltbridges"] = []
    for sb in iset.saltbridge_lneg + iset.saltbridge_pneg:
        interactions["saltbridges"].append({
            "residue": _extract_reslabel(sb),
            "residue_chain": sb.reschain,
            "distance": round(sb.distance, 2),
            "protein_charge": "positive" if getattr(sb, "protispos",
                                       not sb.__class__.__name__.endswith("lneg")) else "negative",
        })

    # --- pi-cation (PLIP v3: pication_laro + pication_paro) ---
    interactions["pication"] = []
    for pc in iset.pication_laro + iset.pication_paro:
        interactions["pication"].append({
            "residue": _extract_reslabel(pc),
            "residue_chain": pc.reschain,
            "distance": round(pc.distance, 2),
        })

    # --- halogen bonds ---
    interactions["halogen"] = []
    for hx in iset.halogen_bonds:
        interactions["halogen"].append({
            "residue": _extract_reslabel(hx),
            "residue_chain": hx.reschain,
            "distance": round(hx.distance, 2),
            "angle": round(hx.angle, 1),
        })

    # print counts
    for k, v in interactions.items():
        if v:
            print(f"    {k}: {len(v)} found")
    return interactions


def generate_plip_xml(complex_pdb, output_dir):
    """Generate PLIP XML report via command-line tool (supplementary)."""
    xml_path = os.path.join(output_dir, "plip_report.xml")
    complex_abs = os.path.abspath(complex_pdb)
    outdir_abs = os.path.abspath(output_dir)
    try:
        # PLIP CLI: -f input.pdb -x (XML mode) -o output_prefix
        # When -x is used, PLIP writes 'report.xml' in CWD.
        result = subprocess.run(
            ["plip", "-f", complex_abs, "-x"],
            capture_output=True, text=True, timeout=120,
            cwd=outdir_abs)
        expected = os.path.join(outdir_abs, "report.xml")
        if os.path.exists(expected):
            if os.path.exists(xml_path):
                os.unlink(xml_path)
            os.rename(expected, xml_path)
            print(f"  PLIP XML report: {xml_path}")
            return xml_path
        else:
            # maybe PLIP wrote to original CWD
            if os.path.exists("report.xml"):
                shutil.move("report.xml", xml_path)
                print(f"  PLIP XML report: {xml_path}")
                return xml_path
    except FileNotFoundError:
        print("  PLIP CLI not available (XML report skipped)")
    except Exception as e:
        print(f"  PLIP XML generation skipped: {e}")
    return None


# ---------------------------------------------------------------------------
# Step 4: Write interaction reports
# ---------------------------------------------------------------------------

def write_interaction_summary(interactions, protein_name, ligand_name,
                               best_affinity, output_txt, output_json):
    """Write human-readable interaction_summary.txt and interaction_summary.json."""

    hbonds      = interactions.get("hbonds", [])
    hydrophobic = interactions.get("hydrophobic", [])
    pistacking  = interactions.get("pistacking", [])
    saltbridges = interactions.get("saltbridges", [])
    pication    = interactions.get("pication", [])
    halogen     = interactions.get("halogen", [])

    # ---- text report ----
    L = []
    L.append("=" * 60)
    L.append("  Protein-Ligand Interaction Analysis")
    L.append("  easy-MD v2.0 | PLIP (Protein-Ligand Interaction Profiler)")
    L.append("=" * 60)
    L.append("")
    L.append(f"Protein:")
    L.append(f"  {protein_name}")
    L.append("")
    L.append(f"Ligand:")
    L.append(f"  {ligand_name}")
    L.append("")
    L.append(f"Docking score:")
    L.append(f"  {best_affinity:.3f} kcal/mol")
    L.append("")

    # Hydrogen Bonds
    L.append("-" * 60)
    L.append("Hydrogen Bonds:")
    if hbonds:
        for hb in hbonds:
            L.append(f"  {hb['residue']} ({hb['residue_chain']}) - "
                     f"Distance: {hb['distance_hb']} A, "
                     f"Angle: {hb['angle']} deg, "
                     f"Type: {hb['type']}")
    else:
        L.append("  None detected")
    L.append("")

    # Hydrophobic Interactions
    L.append("-" * 60)
    L.append("Hydrophobic Interactions:")
    if hydrophobic:
        residues = sorted(set(h["residue"] for h in hydrophobic))
        L.append(f"  Residues: {', '.join(residues)}")
        L.append(f"  Total contacts: {len(hydrophobic)}")
    else:
        L.append("  None detected")
    L.append("")

    # Pi-Pi Interactions
    L.append("-" * 60)
    L.append("Pi-Pi Interactions:")
    if pistacking:
        for ps in pistacking:
            L.append(f"  {ps['residue']} ({ps['residue_chain']}) - "
                     f"Distance: {ps['distance']} A, "
                     f"Angle: {ps['angle']} deg, "
                     f"Type: {ps['type']}")
    else:
        L.append("  None detected")
    L.append("")

    # Salt Bridges
    L.append("-" * 60)
    L.append("Salt Bridges:")
    if saltbridges:
        for sb in saltbridges:
            L.append(f"  {sb['residue']} ({sb['residue_chain']}) - "
                     f"Distance: {sb['distance']} A, "
                     f"Protein charge: {sb['protein_charge']}")
    else:
        L.append("  None detected")
    L.append("")

    # Pi-Cation (if any)
    if pication:
        L.append("-" * 60)
        L.append("Pi-Cation Interactions:")
        for pc in pication:
            L.append(f"  {pc['residue']} ({pc['residue_chain']}) - "
                     f"Distance: {pc['distance']} A")
        L.append("")

    # Halogen Bonds (if any)
    if halogen:
        L.append("-" * 60)
        L.append("Halogen Bonds:")
        for hx in halogen:
            L.append(f"  {hx['residue']} ({hx['residue_chain']}) - "
                     f"Distance: {hx['distance']} A, "
                     f"Angle: {hx['angle']} deg")
        L.append("")

    # Summary
    L.append("=" * 60)
    L.append("Summary:")
    L.append(f"  Total hydrogen bonds:        {len(hbonds)}")
    L.append(f"  Total hydrophobic contacts:   {len(hydrophobic)}")
    L.append(f"  Total pi interactions:        {len(pistacking)}")
    L.append(f"  Total salt bridges:           {len(saltbridges)}")
    if pication:
        L.append(f"  Total pi-cation interactions: {len(pication)}")
    if halogen:
        L.append(f"  Total halogen bonds:          {len(halogen)}")
    L.append("")
    L.append("NOTE: This analysis is based on the best docking pose.")
    L.append("Publication figures MUST use original best_pose.pdb from")
    L.append("AutoDock Vina, NOT PLIP-processed coordinates.")
    L.append("=" * 60)

    with open(output_txt, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"  Wrote {output_txt}")

    # ---- JSON report ----
    json_data = {
        "protein": protein_name,
        "ligand": ligand_name,
        "best_affinity_kcal_mol": best_affinity,
        "interactions": {
            "hydrogen_bonds":      {"count": len(hbonds),      "details": hbonds},
            "hydrophobic_contacts": {"count": len(hydrophobic),
                                     "residues": sorted(set(h["residue"] for h in hydrophobic)),
                                     "total_contacts": len(hydrophobic),
                                     "details": hydrophobic},
            "pi_stacking":         {"count": len(pistacking),  "details": pistacking},
            "salt_bridges":        {"count": len(saltbridges), "details": saltbridges},
            "pi_cation":           {"count": len(pication),    "details": pication},
            "halogen_bonds":       {"count": len(halogen),     "details": halogen},
        },
        "summary": {
            "total_hydrogen_bonds":  len(hbonds),
            "total_hydrophobic":     len(hydrophobic),
            "total_pi_interactions": len(pistacking),
            "total_salt_bridges":    len(saltbridges),
            "total_pi_cation":       len(pication),
            "total_halogen_bonds":   len(halogen),
        },
    }
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(json_data, f, indent=2, ensure_ascii=False)
    print(f"  Wrote {output_json}")

    return json_data


# ---------------------------------------------------------------------------
# Step 5: PyMOL visualization
# ---------------------------------------------------------------------------

# regex for extracting residue number from label like "ASP123" or "ASP123A"
_RESNUM_RE = re.compile(r"^[A-Z]{1,3}(\d+)")


def _extract_resnums(residue_labels):
    """Return sorted unique residue numbers from labels like ['ASP123','TYR45A']."""
    nums = set()
    for r in residue_labels:
        m = _RESNUM_RE.match(r)
        if m:
            nums.add(int(m.group(1)))
    return sorted(nums)


def write_pymol_pml(protein_pdb, interactions, output_pml,
                     protein_name, ligand_name, lig_resname="LIG"):
    """Write pymol_interaction.pml with pre-configured visualization.

    The user opens this script in PyMOL and can freely adjust colors,
    angles, representations, and export publication images.
    """
    hbonds      = interactions.get("hbonds", [])
    hydrophobic = interactions.get("hydrophobic", [])
    pistacking  = interactions.get("pistacking", [])
    saltbridges = interactions.get("saltbridges", [])
    pication    = interactions.get("pication", [])
    halogen     = interactions.get("halogen", [])

    # Collect all interacting residues
    all_interacting = set()
    for lst in [hbonds, hydrophobic, pistacking, saltbridges, pication, halogen]:
        for item in lst:
            all_interacting.add(item["residue"])
    inter_resnums = _extract_resnums(all_interacting)

    pml = []
    pml.append(f"# PyMOL Interaction Visualization")
    pml.append(f"# {protein_name} - {ligand_name}")
    pml.append(f"# easy-MD v2.0 | PLIP Interaction Analysis")
    pml.append(f"#")
    pml.append(f"# Open this file in PyMOL, then freely adjust:")
    pml.append(f"#   - colors (color, set_color)")
    pml.append(f"#   - angles (turn, rotate)")
    pml.append(f"#   - display (show, hide)")
    pml.append(f"#   - export  (ray, png)")
    pml.append("")

    # Load complex (protein + docked ligand) from same directory
    pml.append(f"load complex_for_plip.pdb, complex")
    pml.append("")

    # Protein: cartoon, semi-transparent grey
    pml.append("# --- Protein ---")
    pml.append("show cartoon, complex")
    pml.append("set cartoon_transparency, 0.3, complex")
    pml.append("color grey80, complex")
    pml.append("")

    # Best pose ligand: sticks, cyan
    pml.append("# --- Ligand (best docking pose) ---")
    pml.append(f"show sticks, complex and resn {lig_resname}")
    pml.append(f"color cyan, complex and resn {lig_resname}")
    pml.append("")

    # Interaction residues: sticks, yellow
    if inter_resnums:
        pml.append("# --- Interacting residues (yellow sticks) ---")
        sel_parts = []
        for n in inter_resnums:
            sel_parts.append(f"resi {n}")
        sel_str = " or ".join(sel_parts)
        pml.append(f"select inter_res, complex and ({sel_str})")
        pml.append("show sticks, inter_res")
        pml.append("color yellow, inter_res")
        pml.append("")

    # Hydrogen bonds: dashed lines between donor and acceptor
    if hbonds:
        pml.append("# --- Hydrogen bonds (dashed magenta lines) ---")
        for i, hb in enumerate(hbonds):
            resn_match = _RESNUM_RE.match(hb["residue"])
            if resn_match:
                rnum = resn_match.group(1)
                pml.append(f"distance hbond_{i+1}, "
                           f"complex and resi {rnum}, "
                           f"complex and resn {lig_resname}")
        pml.append("color magenta, hbond_*")
        pml.append("hide labels, hbond_*")
        pml.append("")

    # Salt bridges: also shown as dashed lines
    if saltbridges:
        pml.append("# --- Salt bridges (dashed orange lines) ---")
        for i, sb in enumerate(saltbridges):
            resn_match = _RESNUM_RE.match(sb["residue"])
            if resn_match:
                rnum = resn_match.group(1)
                pml.append(f"distance saltbr_{i+1}, "
                           f"complex and resi {rnum}, "
                           f"complex and resn {lig_resname}")
        pml.append("color orange, saltbr_*")
        pml.append("hide labels, saltbr_*")
        pml.append("")

    # Labels: key interacting residues (max ~15 to avoid clutter)
    if inter_resnums:
        pml.append("# --- Labels on CA atoms of interacting residues ---")
        for n in inter_resnums[:15]:
            pml.append(f"label complex and resi {n} and name CA, \"%s-%s\" % "
                       f"(one_letter[resn], resi)")
        pml.append("")

    # Final view settings
    pml.append("# --- Final settings ---")
    pml.append(f"cmd.center('complex and resn {lig_resname}')")
    pml.append("zoom center, 15.0")
    pml.append("bg_color white")
    pml.append("set antialias, 2")
    pml.append("set dash_gap, 0.4")
    pml.append("set dash_length, 0.2")
    pml.append("set dash_radius, 0.1")
    pml.append("set label_size, 14")
    pml.append("set stick_radius, 0.2")
    pml.append("set valence, 1")
    pml.append("set depth_cue, 0")
    pml.append("")

    pml.append(f"print 'easy-MD: {protein_name} - {ligand_name} interaction analysis'")
    pml.append("print '  Protein: cartoon, grey'")
    pml.append("print '  Ligand: sticks, cyan'")
    pml.append("print '  Interacting residues: sticks, yellow'")
    pml.append("print '  H-bonds: magenta dashes'")
    if saltbridges:
        pml.append("print '  Salt bridges: orange dashes'")

    with open(output_pml, "w") as f:
        f.write("\n".join(pml) + "\n")
    print(f"  Wrote {output_pml}")
    return True


def save_pymol_session(pml_file, pse_file, pymol_path):
    """Attempt to convert .pml to .pse via PyMOL batch mode.

    This is non-critical — if PyMOLWin.exe doesn't support batch mode on
    this Windows installation the .pml file is still fully usable.
    """
    wrapper_pml = pml_file.replace(".pml", "_save.pml")
    pml_dir = os.path.dirname(os.path.abspath(pml_file))
    pse_abs = os.path.abspath(pse_file)

    with open(pml_file, "r") as src, open(wrapper_pml, "w") as dst:
        dst.write(src.read())
        dst.write(f"\nsave {pse_abs}\n")

    try:
        result = subprocess.run(
            [pymol_path, "-cq", wrapper_pml],
            capture_output=True, text=True, timeout=120,
            cwd=pml_dir)
        if os.path.exists(pse_file):
            print(f"  Saved PyMOL session: {pse_file}")
            return True
        else:
            print(f"  PyMOL batch did not produce .pse (non-critical)")
            print(f"    stdout: {result.stdout[:300]}")
            print(f"    stderr: {result.stderr[:300]}")
            return False
    except FileNotFoundError:
        print(f"  PyMOL not found at {pymol_path} — .pse skipped")
        return False
    except Exception as e:
        print(f"  PyMOL batch conversion failed (non-critical): {e}")
        return False
    finally:
        if os.path.exists(wrapper_pml):
            try:
                os.unlink(wrapper_pml)
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Cleanup helpers
# ---------------------------------------------------------------------------

def _cleanup_temp_files(output_dir):
    """Remove internal scratch files from the output directory."""
    for f in os.listdir(output_dir):
        if f.startswith("_best_pose_plip") and f.endswith(".pdb"):
            try:
                os.unlink(os.path.join(output_dir, f))
            except Exception:
                pass


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Protein-ligand interaction analysis (PLIP + PyMOL)")
    parser.add_argument("--protein", required=True,
                        help="Path to protein_clean.pdb")
    parser.add_argument("--best-pose-pdbqt", required=True,
                        help="Path to docked PDBQT from best Vina run "
                             "(e.g. run2/docked.pdbqt — mode 1 is extracted)")
    parser.add_argument("--ligand-sdf", required=True,
                        help="Path to ligand SDF (chemical topology reference)")
    parser.add_argument("--docking-stats", default="docking_stats.json",
                        help="Path to docking_stats.json for affinity score")
    parser.add_argument("--output-dir", default="interaction_analysis",
                        help="Output directory name (created in CWD)")
    parser.add_argument("--protein-name", default="Protein",
                        help="Protein name for reports")
    parser.add_argument("--ligand-name", default="Ligand",
                        help="Ligand name for reports")
    parser.add_argument("--obabel")
    parser.add_argument("--pymol")
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--config", help="Path to easy-md config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    args.protein = str(Path(args.protein).expanduser().resolve())
    args.best_pose_pdbqt = str(
        Path(args.best_pose_pdbqt).expanduser().resolve()
    )
    args.ligand_sdf = str(Path(args.ligand_sdf).expanduser().resolve())
    enter_work_dir(args.work_dir, config)
    args.obabel = args.obabel or configured_path(config, "obabel") or "obabel"
    args.pymol = args.pymol or configured_path(config, "pymol")

    # --- load affinity from docking stats ---
    best_affinity = 0.0
    if os.path.exists(args.docking_stats):
        with open(args.docking_stats, "r") as f:
            dock_stats = json.load(f)
        best_affinity = dock_stats.get("best_affinity", 0.0)

    # --- create output directory ---
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("INTERACTION ANALYSIS")
    print(f"  Protein:  {args.protein_name}")
    print(f"  Ligand:   {args.ligand_name}")
    print(f"  Affinity: {best_affinity:.3f} kcal/mol")
    print("=" * 60)

    # Generate a meaningful PDB residue name from the ligand name
    # e.g. Aspirin->ASP, Belumosudil->BEL, ADP->ADP
    lig_resname = generate_resname(args.ligand_name)
    print(f"\n  Ligand resname: {lig_resname} (chain L)")

    # ---- [1/6] Convert best pose PDBQT -> PDB ----
    print("\n[1/6] Converting best pose PDBQT to PDB...")
    best_pose_pdb = convert_best_pose(args.best_pose_pdbqt,
                                      args.obabel, args.output_dir)
    if best_pose_pdb is None:
        print("ERROR: Failed to convert best docking pose")
        sys.exit(1)

    # ---- [2/6] Build complex for PLIP ----
    print("\n[2/6] Building protein-ligand complex for PLIP...")
    complex_pdb = os.path.join(args.output_dir, "complex_for_plip.pdb")
    nprot, nlig = build_complex(args.protein, best_pose_pdb, complex_pdb,
                                  lig_resname=lig_resname)
    if nlig == 0:
        print("ERROR: Failed to build complex (no ligand atoms)")
        sys.exit(1)

    # ---- [3/6] Run PLIP ----
    print("\n[3/6] Running PLIP interaction detection...")
    try:
        interactions = run_plip(complex_pdb)
    except Exception as e:
        print(f"ERROR: PLIP analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Generate XML report via CLI (supplementary)
    generate_plip_xml(complex_pdb, args.output_dir)

    # ---- [4/6] Write reports ----
    print("\n[4/6] Writing interaction reports...")
    out_txt  = os.path.join(args.output_dir, "interaction_summary.txt")
    out_json = os.path.join(args.output_dir, "interaction_summary.json")
    json_data = write_interaction_summary(
        interactions, args.protein_name, args.ligand_name,
        best_affinity, out_txt, out_json)

    # ---- [5/6] 2D interaction diagram ----
    print("\n[5/6] Generating 2D interaction diagram...")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ia2d_py = os.path.join(script_dir, "interaction_2d.py")
    if os.path.exists(ia2d_py):
        ia2d_cmd = [
            sys.executable, ia2d_py,
            "--json", os.path.abspath(out_json),
            "--sdf", os.path.abspath(args.ligand_sdf),
            "--complex", os.path.abspath(complex_pdb),
            "--output-dir", os.path.abspath(args.output_dir),
            "--protein-name", args.protein_name,
            "--ligand-name", args.ligand_name,
            "--ligand-resname", lig_resname,
            "--best-affinity", str(best_affinity),
        ]
        result = subprocess.run(ia2d_cmd, capture_output=True, text=True,
                                timeout=300)
        print(result.stdout)
        if result.returncode != 0:
            print(f"  WARNING: 2D diagram generation had issues:")
            print(f"    {result.stderr[:500]}")
    else:
        print(f"  WARNING: interaction_2d.py not found at {ia2d_py}")
        print("    Skipping 2D diagram generation.")

    # ---- [6/6] PyMOL visualization ----
    print("\n[6/6] Generating PyMOL visualization files...")
    pml_path = os.path.join(args.output_dir, "pymol_interaction.pml")

    # Use relative paths so the .pml / .pse are portable if the directory moves
    protein_rel = os.path.basename(args.protein)
    write_pymol_pml(protein_rel, interactions, pml_path,
                     args.protein_name, args.ligand_name,
                     lig_resname=lig_resname)

    pse_path = os.path.join(args.output_dir, "pymol_interaction.pse")
    save_pymol_session(pml_path, pse_path, args.pymol)

    # ---- cleanup scratch files ----
    _cleanup_temp_files(args.output_dir)

    # ---- final summary ----
    s = json_data["summary"]
    print("\n" + "=" * 60)
    print("INTERACTION ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"  Output directory: {os.path.abspath(args.output_dir)}")
    print(f"  Hydrogen bonds:     {s['total_hydrogen_bonds']}")
    print(f"  Hydrophobic:        {s['total_hydrophobic']}")
    print(f"  Pi interactions:    {s['total_pi_interactions']}")
    print(f"  Salt bridges:       {s['total_salt_bridges']}")
    print("=" * 60)
    print("  Output files:")
    for f in sorted(os.listdir(args.output_dir)):
        if not f.startswith("_"):
            fpath = os.path.join(args.output_dir, f)
            size = os.path.getsize(fpath) if os.path.isfile(fpath) else 0
            print(f"    {f}  ({size:,} bytes)")
    print("=" * 60)


if __name__ == "__main__":
    main()
