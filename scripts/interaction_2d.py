#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""2D Protein-Ligand Interaction Diagram Generator.

Reads PLIP interaction_summary.json + ligand SDF and produces
publication-quality 2D interaction diagrams.

IMPORTANT: All interaction data comes from PLIP via interaction_summary.json.
This script does NOT re-compute interactions, does NOT call LigPlot+,
and does NOT run any independent force-field analysis.

Output: interaction_2D.png, interaction_2D.svg, interaction_2D.emf

Usage:
    python interaction_2d.py \
        --json interaction_analysis/interaction_summary.json \
        --sdf ligand.sdf \
        --complex interaction_analysis/complex_for_plip.pdb \
        --output-dir interaction_analysis \
        --ligand-name "Aspirin" --ligand-resname "ASP" \
        --protein-name "PTGS2"
"""

import os, sys, json, argparse, math, io, struct, re
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from PIL import Image, ImageDraw, ImageFont

# ---------------------------------------------------------------------------
# Color & style constants
# ---------------------------------------------------------------------------
COLORS = {
    "hbond":       (34,  139, 34),
    "hydrophobic": (220, 20,  60),
    "pistacking":  (65,  105, 225),
    "saltbridge":  (255, 140, 0),
    "pication":    (147, 112, 219),
    "halogen":     (0,   206, 209),
    "residue_bg":  (248, 248, 248),
    "residue_text": (30, 30, 30),
    "legend_text":  (60, 60, 60),
    "bg":           (255, 255, 255),
}

CANVAS_W = 2800
CANVAS_H = 2400
LIGAND_BOX_W = 1000
LIGAND_BOX_H = 1000
LIGAND_CENTER = (CANVAS_W // 2, 1050)
DPI = 300

RESIDUE_FONT_SIZE = 28
LEGEND_FONT_SIZE = 24
TITLE_FONT_SIZE = 36
SUBTITLE_FONT_SIZE = 22

LABEL_RADIUS_MIN = 480
LABEL_RADIUS_MAX = 750


# ---------------------------------------------------------------------------
# Font loading
# ---------------------------------------------------------------------------

def _get_font(size, bold=False):
    """Load a TrueType font, falling back to default."""
    candidates = []
    if bold:
        candidates = [
            "C:/Windows/Fonts/arialbd.ttf",
            "C:/Windows/Fonts/calibrib.ttf",
            "C:/Windows/Fonts/segoeuib.ttf",
        ]
    candidates += [
        "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibri.ttf",
        "C:/Windows/Fonts/segoeui.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]
    for fp in candidates:
        if os.path.exists(fp):
            try:
                return ImageFont.truetype(fp, size)
            except Exception:
                continue
    return ImageFont.load_default()


# ---------------------------------------------------------------------------
# Step 1: Build atom index mapping from complex PDB
# ---------------------------------------------------------------------------

def build_atom_map(complex_pdb, lig_resname, lig_chain):
    """Map PDB atom serial numbers -> 0-based ligand atom indices.

    Parses complex_for_plip.pdb to find all ligand HETATM records.
    Returns: {pdb_serial_num: ligand_atom_index}, num_ligand_atoms
    """
    serial_to_idx = {}
    lig_count = 0
    with open(complex_pdb, "r") as f:
        for line in f:
            if not (line.startswith("HETATM") or line.startswith("ATOM")):
                continue
            if len(line) < 22:
                continue
            try:
                serial = int(line[6:11].strip())
            except ValueError:
                continue
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            if resname == lig_resname and chain == lig_chain:
                serial_to_idx[serial] = lig_count
                lig_count += 1
    print(f"  Atom map: {lig_count} ligand atoms indexed from complex PDB")
    return serial_to_idx, lig_count


# ---------------------------------------------------------------------------
# Step 2: Draw ligand with RDKit -> PIL Image
# ---------------------------------------------------------------------------

def draw_ligand_image(mol, size=LIGAND_BOX_W):
    """Render ligand as a transparent-background PIL Image."""
    mol = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(size, size)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.bondLineWidth = 3
    opts.addStereoAnnotation = True
    opts.baseFontSize = 0.8
    opts.atomLabelDeuteriumTritium = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    png_bytes = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png_bytes)).convert("RGBA")
    return img, mol


def get_atom_2d_positions(mol):
    """Return dict: atom_idx -> (x, y) in RDKit normalized coordinates."""
    conf = mol.GetConformer()
    positions = {}
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        positions[i] = (p.x, p.y)
    return positions


# ---------------------------------------------------------------------------
# Step 3: Residue label placement
# ---------------------------------------------------------------------------

def place_residue_labels(interactions, atom_map, ligand_mol,
                         ligand_img_w, ligand_img_h,
                         lig_center_x, lig_center_y):
    """Compute label positions and interaction line endpoints.

    Returns list of dicts with keys: residue, interaction, interaction_kind,
    anchor_x, anchor_y, label_x, label_y, distance, detail
    """
    atom_pos = get_atom_2d_positions(ligand_mol)

    # Normalize RDKit coords to [0,1]
    xs = [p[0] for p in atom_pos.values()]
    ys = [p[1] for p in atom_pos.values()]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    xspan = (xmax - xmin) or 1.0
    yspan = (ymax - ymin) or 1.0

    def norm_pos(idx):
        x, y = atom_pos.get(idx, (0.5, 0.5))
        return ((x - xmin) / xspan, (y - ymin) / yspan)

    def to_canvas(nx, ny):
        px = lig_center_x - ligand_img_w // 2 + nx * ligand_img_w
        py = lig_center_y - ligand_img_h // 2 + (1.0 - ny) * ligand_img_h
        return px, py

    labels = []

    # --- Hydrogen bonds ---
    hbonds = interactions.get("hydrogen_bonds", {}).get("details", [])
    for hb in hbonds:
        residue = hb["residue"]
        lig_atom_idx = None
        ligand_serials = hb.get("ligand_atom_serials", [])
        don = hb.get("atom_donor") or (ligand_serials[0] if ligand_serials else None)
        acc = hb.get("atom_acceptor")
        if don is not None and don in atom_map:
            lig_atom_idx = atom_map[don]
        elif acc is not None and acc in atom_map:
            lig_atom_idx = atom_map[acc]

        if lig_atom_idx is not None and lig_atom_idx in atom_pos:
            nx, ny = norm_pos(lig_atom_idx)
        else:
            angle = (hash(residue) % 360) * math.pi / 180.0
            nx, ny = 0.5 + 0.35 * math.cos(angle), 0.5 + 0.35 * math.sin(angle)

        ax, ay = to_canvas(nx, ny)
        angle = math.atan2(ay - lig_center_y, ax - lig_center_x)
        dist = math.hypot(ax - lig_center_x, ay - lig_center_y)
        label_dist = max(dist * 1.35, LABEL_RADIUS_MIN)
        label_dist = min(label_dist, LABEL_RADIUS_MAX)
        lx = lig_center_x + label_dist * math.cos(angle)
        ly = lig_center_y + label_dist * math.sin(angle)

        hb_distance = float(hb.get("distance_hb", hb.get("distance", 0.0)))
        labels.append({
            "residue": residue,
            "interaction": "H-bond",
            "interaction_kind": "hbond",
            "anchor_x": ax, "anchor_y": ay,
            "label_x": lx, "label_y": ly,
            "distance": f"{hb_distance:.1f} A",
            "detail": f"{hb_distance:.1f} A, {float(hb.get('angle') or 0):.0f}\u00b0",
        })

    # --- Salt bridges ---
    saltbridges = interactions.get("salt_bridges", {}).get("details", [])
    for sb in saltbridges:
        residue = sb["residue"]
        angle = (hash(residue + "sb") % 360) * math.pi / 180.0
        nx, ny = 0.5 + 0.38 * math.cos(angle), 0.5 + 0.38 * math.sin(angle)
        ax, ay = to_canvas(nx, ny)
        dist = math.hypot(ax - lig_center_x, ay - lig_center_y)
        label_dist = max(dist * 1.4, LABEL_RADIUS_MIN)
        label_dist = min(label_dist, LABEL_RADIUS_MAX)
        lx = lig_center_x + label_dist * math.cos(angle)
        ly = lig_center_y + label_dist * math.sin(angle)
        labels.append({
            "residue": residue,
            "interaction": "Salt bridge",
            "interaction_kind": "saltbridge",
            "anchor_x": ax, "anchor_y": ay,
            "label_x": lx, "label_y": ly,
            "distance": f"{sb['distance']:.1f} A",
            "detail": f"{sb['distance']:.1f} A",
        })

    # --- Pi-stacking ---
    pistacking = interactions.get("pi_stacking", {}).get("details", [])
    for ps in pistacking:
        residue = ps["residue"]
        angle = (hash(residue + "pi") % 360) * math.pi / 180.0
        nx, ny = 0.5 + 0.36 * math.cos(angle), 0.5 + 0.36 * math.sin(angle)
        ax, ay = to_canvas(nx, ny)
        dist = math.hypot(ax - lig_center_x, ay - lig_center_y)
        label_dist = max(dist * 1.38, LABEL_RADIUS_MIN)
        label_dist = min(label_dist, LABEL_RADIUS_MAX)
        lx = lig_center_x + label_dist * math.cos(angle)
        ly = lig_center_y + label_dist * math.sin(angle)
        labels.append({
            "residue": residue,
            "interaction": "Pi-stacking",
            "interaction_kind": "pistacking",
            "anchor_x": ax, "anchor_y": ay,
            "label_x": lx, "label_y": ly,
            "distance": f"{ps['distance']:.1f} A",
            "detail": f"{ps['distance']:.1f} A, {ps.get('type', '')}",
        })

    # --- Pi-cation ---
    pication = interactions.get("pi_cation", {}).get("details", [])
    for pc in pication:
        residue = pc["residue"]
        angle = (hash(residue + "pc") % 360) * math.pi / 180.0
        nx, ny = 0.5 + 0.37 * math.cos(angle), 0.5 + 0.37 * math.sin(angle)
        ax, ay = to_canvas(nx, ny)
        dist = math.hypot(ax - lig_center_x, ay - lig_center_y)
        label_dist = max(dist * 1.36, LABEL_RADIUS_MIN)
        label_dist = min(label_dist, LABEL_RADIUS_MAX)
        lx = lig_center_x + label_dist * math.cos(angle)
        ly = lig_center_y + label_dist * math.sin(angle)
        labels.append({
            "residue": residue,
            "interaction": "Pi-cation",
            "interaction_kind": "pication",
            "anchor_x": ax, "anchor_y": ay,
            "label_x": lx, "label_y": ly,
            "distance": f"{pc['distance']:.1f} A",
            "detail": f"{pc['distance']:.1f} A",
        })

    # --- Hydrophobic contacts (grouped by residue) ---
    hydrophobic = interactions.get("hydrophobic_contacts", {}).get("details", [])
    if hydrophobic:
        hp_by_res = defaultdict(list)
        for hc in hydrophobic:
            hp_by_res[hc["residue"]].append(hc.get("distance", 0))

        n_hp = len(hp_by_res)
        for i, (residue, dists) in enumerate(sorted(hp_by_res.items())):
            avg_dist = sum(dists) / len(dists) if dists else 0
            angle = (360.0 * i / n_hp + 180.0 / n_hp) * math.pi / 180.0
            nx, ny = 0.5 + 0.4 * math.cos(angle), 0.5 + 0.4 * math.sin(angle)
            ax, ay = to_canvas(nx, ny)
            dist = math.hypot(ax - lig_center_x, ay - lig_center_y)
            label_dist = max(dist * 1.42, LABEL_RADIUS_MIN + 40)
            label_dist = min(label_dist, LABEL_RADIUS_MAX)
            lx = lig_center_x + label_dist * math.cos(angle)
            ly = lig_center_y + label_dist * math.sin(angle)
            labels.append({
                "residue": residue,
                "interaction": "Hydrophobic",
                "interaction_kind": "hydrophobic",
                "anchor_x": ax, "anchor_y": ay,
                "label_x": lx, "label_y": ly,
                "distance": f"{avg_dist:.1f} A" if avg_dist > 0 else "",
                "detail": (f"{len(dists)} contacts, ~{avg_dist:.1f} A"
                           if avg_dist > 0 else f"{len(dists)} contacts"),
            })

    # --- Halogen bonds ---
    halogen = interactions.get("halogen_bonds", {}).get("details", [])
    for hx in halogen:
        residue = hx["residue"]
        angle = (hash(residue + "hx") % 360) * math.pi / 180.0
        nx, ny = 0.5 + 0.36 * math.cos(angle), 0.5 + 0.36 * math.sin(angle)
        ax, ay = to_canvas(nx, ny)
        dist = math.hypot(ax - lig_center_x, ay - lig_center_y)
        label_dist = max(dist * 1.38, LABEL_RADIUS_MIN)
        label_dist = min(label_dist, LABEL_RADIUS_MAX)
        lx = lig_center_x + label_dist * math.cos(angle)
        ly = lig_center_y + label_dist * math.sin(angle)
        labels.append({
            "residue": residue,
            "interaction": "Halogen bond",
            "interaction_kind": "halogen",
            "anchor_x": ax, "anchor_y": ay,
            "label_x": lx, "label_y": ly,
            "distance": f"{hx['distance']:.1f} A",
            "detail": f"{hx['distance']:.1f} A, {hx['angle']:.0f}\u00b0",
        })

    # --- Deduplicate: one label per residue (prefer more specific interactions) ---
    priority = {"hbond": 0, "saltbridge": 1, "pistacking": 2,
                "pication": 3, "halogen": 4, "hydrophobic": 5}
    seen = {}
    for lbl in labels:
        r = lbl["residue"]
        if r not in seen or priority.get(lbl["interaction_kind"], 9) < \
           priority.get(seen[r]["interaction_kind"], 9):
            seen[r] = lbl

    result = list(seen.values())
    print(f"  Placed {len(result)} residue labels")
    return result


# ---------------------------------------------------------------------------
# Step 4: Render the composite diagram
# ---------------------------------------------------------------------------

def draw_interaction_diagram(ligand_img, ligand_mol, residue_labels,
                              protein_name, ligand_name, best_affinity,
                              output_png, output_svg, output_emf):
    """Composite: ligand in center, residues around, interaction lines."""

    canvas = Image.new("RGBA", (CANVAS_W, CANVAS_H), COLORS["bg"] + (255,))
    draw = ImageDraw.Draw(canvas)

    # Paste ligand
    lig_w, lig_h = ligand_img.size
    lig_x = LIGAND_CENTER[0] - lig_w // 2
    lig_y = LIGAND_CENTER[1] - lig_h // 2
    canvas.paste(ligand_img, (lig_x, lig_y), ligand_img)

    font_res = _get_font(RESIDUE_FONT_SIZE, bold=True)
    font_res_sm = _get_font(RESIDUE_FONT_SIZE - 4)
    font_legend = _get_font(LEGEND_FONT_SIZE)
    font_title = _get_font(TITLE_FONT_SIZE, bold=True)
    font_sub = _get_font(SUBTITLE_FONT_SIZE)
    font_detail = _get_font(RESIDUE_FONT_SIZE - 8)

    # Draw interaction lines + residue labels
    for lbl in residue_labels:
        kind = lbl["interaction_kind"]
        color = COLORS.get(kind, (100, 100, 100))
        ax, ay = lbl["anchor_x"], lbl["anchor_y"]
        lx, ly = lbl["label_x"], lbl["label_y"]

        if kind == "hbond":
            _draw_dashed_line(draw, (ax, ay), (lx, ly), color,
                              dash_len=10, gap_len=6, width=3)
        elif kind == "saltbridge":
            _draw_dashed_line(draw, (ax, ay), (lx, ly), color,
                              dash_len=14, gap_len=6, width=3)
        elif kind in ("pistacking", "pication", "halogen"):
            _draw_dashed_line(draw, (ax, ay), (lx, ly), color,
                              dash_len=8, gap_len=4, width=2)
        elif kind == "hydrophobic":
            _draw_hydrophobic_arc(draw, (ax, ay), (lx, ly), color)

        _draw_residue_label(draw, lx, ly, lbl, color, font_res, font_res_sm)

    # Legend
    legend_y = CANVAS_H - 180
    legend_items = [
        ("H-bond",      COLORS["hbond"],       "dashed"),
        ("Salt bridge", COLORS["saltbridge"],  "dashed"),
        ("Hydrophobic", COLORS["hydrophobic"], "arc"),
        ("Pi-stacking", COLORS["pistacking"],  "dashed"),
        ("Pi-cation",   COLORS["pication"],    "dashed"),
    ]
    if any(l["interaction_kind"] == "halogen" for l in residue_labels):
        legend_items.append(("Halogen", COLORS["halogen"], "dashed"))

    legend_x_start = CANVAS_W // 2 - 520
    legend_x = legend_x_start
    for name, clr, style in legend_items:
        x0, y0 = legend_x, legend_y
        if style == "dashed":
            _draw_dashed_line(draw, (x0, y0 + 8), (x0 + 50, y0 + 8),
                              clr, dash_len=8, gap_len=4, width=3)
        elif style == "arc":
            _draw_hydrophobic_arc_small(draw, x0 + 25, y0 + 8, clr)
        draw.text((x0 + 58, y0 - 4), name, fill=COLORS["legend_text"],
                  font=font_legend)
        legend_x += 200

    # Title
    title = f"{protein_name}  -  {ligand_name}"
    draw.text((CANVAS_W // 2, 35), title, fill=(20, 20, 20),
              font=font_title, anchor="mt")

    # Subtitle
    aff_text = f"Best affinity: {best_affinity:.2f} kcal/mol" if best_affinity else ""
    subtitle = f"Protein-Ligand Interaction Diagram  |  {aff_text}"
    draw.text((CANVAS_W // 2, 80), subtitle, fill=(100, 100, 100),
              font=font_sub, anchor="mt")

    # Save PNG
    canvas = canvas.convert("RGB")
    canvas.save(output_png, "PNG", dpi=(DPI, DPI))
    print(f"  Saved {output_png}")

    # Generate SVG
    _generate_svg(ligand_img, ligand_mol, residue_labels,
                  protein_name, ligand_name, best_affinity,
                  output_svg)
    print(f"  Saved {output_svg}")

    # Generate EMF
    _generate_emf(output_png, output_emf, CANVAS_W, CANVAS_H)
    print(f"  Saved {output_emf}")

    return True


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def _draw_dashed_line(draw, p1, p2, color, dash_len=10, gap_len=6, width=2):
    """Draw a dashed line from p1 to p2."""
    x1, y1 = p1
    x2, y2 = p2
    total = math.hypot(x2 - x1, y2 - y1)
    if total < 1:
        return
    dx = (x2 - x1) / total
    dy = (y2 - y1) / total
    pos = 0.0
    on = True
    while pos < total:
        seg_len = dash_len if on else gap_len
        seg_len = min(seg_len, total - pos)
        if on:
            sx, sy = x1 + pos * dx, y1 + pos * dy
            ex, ey = x1 + (pos + seg_len) * dx, y1 + (pos + seg_len) * dy
            draw.line([(sx, sy), (ex, ey)], fill=color, width=width)
        pos += seg_len
        on = not on


def _draw_hydrophobic_arc(draw, anchor, label_pos, color):
    """Draw a hydrophobic 'eyelash' indicator: arc + radiating spokes."""
    ax, ay = anchor
    lx, ly = label_pos
    cx = (ax + lx) / 2
    cy = (ay + ly) / 2
    angle = math.atan2(ay - ly, ax - lx)
    radius = 18

    steps = 16
    arc_start = angle - math.pi * 0.4
    arc_end = angle + math.pi * 0.4
    pts = []
    for i in range(steps + 1):
        a = arc_start + (arc_end - arc_start) * i / steps
        pts.append((cx + radius * math.cos(a), cy + radius * math.sin(a)))
    for i in range(len(pts) - 1):
        draw.line([pts[i], pts[i + 1]], fill=color, width=2)

    for i in range(3):
        spoke_a = angle - 0.35 + i * 0.35
        sx = cx + (radius * 0.5) * math.cos(spoke_a)
        sy = cy + (radius * 0.5) * math.sin(spoke_a)
        ex = cx + (radius * 1.6) * math.cos(spoke_a)
        ey = cy + (radius * 1.6) * math.sin(spoke_a)
        draw.line([(sx, sy), (ex, ey)], fill=color, width=1)


def _draw_hydrophobic_arc_small(draw, cx, cy, color):
    """Small hydrophobic indicator for legend."""
    radius = 10
    steps = 12
    angle_up = -math.pi / 2
    arc_start = angle_up - math.pi * 0.4
    arc_end = angle_up + math.pi * 0.4
    pts = []
    for i in range(steps + 1):
        a = arc_start + (arc_end - arc_start) * i / steps
        pts.append((cx + radius * math.cos(a), cy + radius * math.sin(a)))
    for i in range(len(pts) - 1):
        draw.line([pts[i], pts[i + 1]], fill=color, width=2)
    for i in range(3):
        sa = angle_up - 0.35 + i * 0.35
        sx = cx + 5 * math.cos(sa)
        sy = cy + 5 * math.sin(sa)
        ex = cx + 16 * math.cos(sa)
        ey = cy + 16 * math.sin(sa)
        draw.line([(sx, sy), (ex, ey)], fill=color, width=1)


def _draw_residue_label(draw, lx, ly, lbl, color, font, font_sm):
    """Draw a residue label box with interaction type indicator."""
    residue = lbl["residue"]
    detail = lbl.get("detail", "")

    bbox = draw.textbbox((0, 0), residue, font=font)
    tw = bbox[2] - bbox[0]
    th = bbox[3] - bbox[1]
    pad_x, pad_y = 14, 8

    rx1 = lx - tw // 2 - pad_x
    ry1 = ly - th // 2 - pad_y
    rx2 = lx + tw // 2 + pad_x
    ry2 = ly + th // 2 + pad_y

    draw.rounded_rectangle(
        [rx1, ry1, rx2, ry2], radius=8,
        fill=COLORS["residue_bg"],
        outline=color, width=2)

    draw.text((lx, ly), residue, fill=COLORS["residue_text"],
              font=font, anchor="mm")

    if detail:
        draw.text((lx, ry2 + 4), detail, fill=COLORS["legend_text"],
                  font=font_sm, anchor="mt")


# ---------------------------------------------------------------------------
# Step 5: SVG generation
# ---------------------------------------------------------------------------

def _generate_svg(ligand_img, ligand_mol, residue_labels,
                  protein_name, ligand_name, best_affinity,
                  output_path):
    """Generate SVG with RDKit ligand + manual annotations."""
    mol = Chem.Mol(ligand_mol)
    AllChem.Compute2DCoords(mol)

    lig_size = LIGAND_BOX_W
    drawer = rdMolDraw2D.MolDraw2DSVG(lig_size, lig_size)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    opts.bondLineWidth = 3
    opts.baseFontSize = 0.8
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    lig_svg = drawer.GetDrawingText()

    inner_match = re.search(r"<svg[^>]*>(.*)</svg>", lig_svg, re.DOTALL)
    lig_inner = inner_match.group(1) if inner_match else lig_svg

    color_css = {
        "hbond": "#228B22", "saltbridge": "#FF8C00",
        "pistacking": "#4169E1", "pication": "#9370DB",
        "halogen": "#00CED1", "hydrophobic": "#DC143C",
    }

    lig_left = LIGAND_CENTER[0] - lig_size // 2
    lig_top = LIGAND_CENTER[1] - lig_size // 2

    anno_parts = []
    for lbl in residue_labels:
        kind = lbl["interaction_kind"]
        clr = color_css.get(kind, "#666666")
        ax, ay = lbl["anchor_x"], lbl["anchor_y"]
        lx, ly = lbl["label_x"], lbl["label_y"]
        residue = lbl["residue"]
        detail = lbl.get("detail", "")

        if kind == "hbond":
            dash = 'stroke-dasharray="8,5"'
        elif kind == "saltbridge":
            dash = 'stroke-dasharray="12,6"'
        elif kind in ("pistacking", "pication", "halogen"):
            dash = 'stroke-dasharray="6,4"'
        else:
            dash = ""

        anno_parts.append(
            f'<line x1="{ax:.1f}" y1="{ay:.1f}" x2="{lx:.1f}" y2="{ly:.1f}" '
            f'stroke="{clr}" stroke-width="2.5" {dash}/>')

        if kind == "hydrophobic":
            cx = (ax + lx) / 2
            cy = (ay + ly) / 2
            angle = math.atan2(ay - ly, ax - lx)
            r = 18
            a_start = angle - math.pi * 0.4
            a_end = angle + math.pi * 0.4
            sx = cx + r * math.cos(a_start)
            sy = cy + r * math.sin(a_start)
            ex = cx + r * math.cos(a_end)
            ey = cy + r * math.sin(a_end)
            large = 1 if (a_end - a_start) > math.pi else 0
            anno_parts.append(
                f'<path d="M {sx:.1f},{sy:.1f} A {r:.1f},{r:.1f} 0 {large},1 '
                f'{ex:.1f},{ey:.1f}" stroke="{clr}" stroke-width="2" fill="none"/>')
            for i in range(3):
                sa = angle - 0.35 + i * 0.35
                sx2 = cx + r * 0.5 * math.cos(sa)
                sy2 = cy + r * 0.5 * math.sin(sa)
                ex2 = cx + r * 1.6 * math.cos(sa)
                ey2 = cy + r * 1.6 * math.sin(sa)
                anno_parts.append(
                    f'<line x1="{sx2:.1f}" y1="{sy2:.1f}" x2="{ex2:.1f}" '
                    f'y2="{ey2:.1f}" stroke="{clr}" stroke-width="1"/>')

        tw = len(residue) * (RESIDUE_FONT_SIZE * 0.65)
        th = RESIDUE_FONT_SIZE * 1.4
        pad_x, pad_y = 12, 6
        rx1 = lx - tw / 2 - pad_x
        ry1 = ly - th / 2 - pad_y
        rx2 = lx + tw / 2 + pad_x
        ry2 = ly + th / 2 + pad_y

        anno_parts.append(
            f'<rect x="{rx1:.1f}" y="{ry1:.1f}" width="{rx2 - rx1:.1f}" '
            f'height="{ry2 - ry1:.1f}" rx="8" ry="8" '
            f'fill="#f8f8f8" stroke="{clr}" stroke-width="1.5"/>')
        anno_parts.append(
            f'<text x="{lx:.1f}" y="{ly + th * 0.15:.1f}" '
            f'text-anchor="middle" font-family="Arial,sans-serif" '
            f'font-size="{RESIDUE_FONT_SIZE}" font-weight="bold" '
            f'fill="#1e1e1e">{residue}</text>')

        if detail:
            anno_parts.append(
                f'<text x="{lx:.1f}" y="{ry2 + RESIDUE_FONT_SIZE - 4:.1f}" '
                f'text-anchor="middle" font-family="Arial,sans-serif" '
                f'font-size="{RESIDUE_FONT_SIZE - 8}" fill="#646464">{detail}</text>')

    # Legend
    leg_x = CANVAS_W // 2 - 500
    leg_y = CANVAS_H - 160
    items = [
        ("H-bond", "#228B22"), ("Salt bridge", "#FF8C00"),
        ("Hydrophobic", "#DC143C"), ("Pi-stacking", "#4169E1"),
        ("Pi-cation", "#9370DB"),
    ]
    for i, (name, clr) in enumerate(items):
        lx2 = leg_x + i * 195
        anno_parts.append(
            f'<line x1="{lx2}" y1="{leg_y + 6}" x2="{lx2 + 45}" y2="{leg_y + 6}" '
            f'stroke="{clr}" stroke-width="3" stroke-dasharray="7,4"/>')
        anno_parts.append(
            f'<text x="{lx2 + 52}" y="{leg_y + 12}" '
            f'font-family="Arial,sans-serif" font-size="{LEGEND_FONT_SIZE}" '
            f'fill="#646464">{name}</text>')

    # Title
    aff_str = f"Best affinity: {best_affinity:.2f} kcal/mol" if best_affinity else ""
    anno_parts.append(
        f'<text x="{CANVAS_W // 2}" y="38" text-anchor="middle" '
        f'font-family="Arial,sans-serif" font-size="{TITLE_FONT_SIZE}" '
        f'font-weight="bold" fill="#141414">{protein_name} - {ligand_name}</text>')
    if aff_str:
        anno_parts.append(
            f'<text x="{CANVAS_W // 2}" y="78" text-anchor="middle" '
            f'font-family="Arial,sans-serif" font-size="{SUBTITLE_FONT_SIZE}" '
            f'fill="#888888">Protein-Ligand Interaction Diagram  |  {aff_str}</text>')

    full_svg = f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{CANVAS_W}" height="{CANVAS_H}"
     viewBox="0 0 {CANVAS_W} {CANVAS_H}">
<rect width="100%" height="100%" fill="white"/>
<g transform="translate({lig_left},{lig_top})">
{lig_inner}
</g>
<g id="interactions">
{chr(10).join(f"  {p}" for p in anno_parts)}
</g>
</svg>"""

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(full_svg)


# ---------------------------------------------------------------------------
# Step 6: EMF generation (Windows Enhanced Metafile)
# ---------------------------------------------------------------------------

def _generate_emf(png_path, emf_path, width_px, height_px):
    """Create an EMF file embedding the high-res PNG diagram.

    Uses raw EMF binary format with EMR_STRETCHDIBITS to embed a
    device-independent bitmap. Compatible with Microsoft Office.
    Reference: [MS-EMF] Enhanced Metafile Format Specification.
    """
    png_img = Image.open(png_path).convert("RGBA")
    rgb_img = Image.new("RGB", png_img.size, (255, 255, 255))
    rgb_img.paste(png_img, mask=png_img.split()[3])
    w, h = rgb_img.size

    # Build bottom-up DIB data (BGR, DWORD-aligned rows)
    row_size = ((w * 24 + 31) // 32) * 4
    pixel_data = bytearray()
    for y in range(h - 1, -1, -1):
        row = bytearray(row_size)
        for x in range(w):
            r, g, b = rgb_img.getpixel((x, y))
            row[x * 3] = b
            row[x * 3 + 1] = g
            row[x * 3 + 2] = r
        pixel_data.extend(row)

    # BITMAPINFOHEADER (40 bytes)
    bih = struct.pack("<IiiHHiiiiII",
                       40, w, h, 1, 24, 0, len(pixel_data),
                       2835, 2835, 0, 0)

    dib_data = bih + bytes(pixel_data)

    # EMR_STRETCHDIBITS fixed part (84 bytes after record header)
    fixed_part = struct.pack("<iiii", 0, 0, w, h)          # rclBounds
    fixed_part += struct.pack("<iiii", 0, 0, w - 1, h - 1)
    fixed_part += struct.pack("<iiii", 0, 0, w, h)          # xDest,yDest,cxDest,cyDest
    fixed_part += struct.pack("<iiii", 0, 0, w, h)          # xSrc,ySrc,cxSrc,cySrc
    fixed_part += struct.pack("<IIII", 0x00CC0020, 0, 0, 0) # dwRop=SRCCOPY, etc.

    record_type = 81  # EMR_STRETCHDIBITS
    record_size = 8 + len(fixed_part) + len(dib_data)
    sdb_record = struct.pack("<II", record_type, record_size)
    sdb_record += fixed_part + dib_data

    # EMR_EOF record
    eof_record = struct.pack("<IIiiII", 14, 20, 0, 0, 0, 0)

    # EMF Header
    himetric_scale = 2540.0 / 72.0
    frame_right = max(1, int(w * himetric_scale))
    frame_bottom = max(1, int(h * himetric_scale))

    header_record = struct.pack("<II", 1, 108)  # EMR_HEADER, size=108
    header_record += struct.pack("<iiii", 0, 0, w - 1, h - 1)  # bounds (inclusive)
    header_record += struct.pack("<iiii", 0, 0, frame_right, frame_bottom)  # frame
    header_record += b" EMF"                                         # signature
    header_record += struct.pack("<I", 0x00010000)                   # version
    total_size = 108 + record_size + 20
    header_record += struct.pack("<I", total_size)                   # total bytes
    header_record += struct.pack("<I", 3)                            # 3 records
    header_record += struct.pack("<H", 1)                            # handles
    header_record += struct.pack("<H", 0)                            # reserved
    header_record += struct.pack("<II", 0, 0)                        # description
    header_record += struct.pack("<I", 0)                            # nPalEntries
    header_record += struct.pack("<ii", w, h)                        # szlDevice (pixels)
    mm_w = int(w * 25.4 / 72.0)
    mm_h = int(h * 25.4 / 72.0)
    header_record += struct.pack("<ii", mm_w, mm_h)                  # szlMillimeters
    header_record += struct.pack("<II", 0, 0)                        # cbPixelFormat, offPixelFormat
    header_record += struct.pack("<I", 0)                            # bOpenGL
    um_w = int(w * 25400.0 / 72.0)
    um_h = int(h * 25400.0 / 72.0)
    header_record += struct.pack("<ii", um_w, um_h)                  # szlMicrometers
    assert len(header_record) == 108, \
        f"Header is {len(header_record)} bytes, expected 108"

    with open(emf_path, "wb") as f:
        f.write(header_record)
        f.write(sdb_record)
        f.write(eof_record)

    emf_size = os.path.getsize(emf_path)
    print(f"  EMF file: {emf_size:,} bytes ({w}x{h})")


# ---------------------------------------------------------------------------
# Residue name generation
# ---------------------------------------------------------------------------

def generate_resname(ligand_name):
    """Generate a 3-letter PDB residue name from the ligand name.

    Avoids collision with standard amino acid 3-letter codes.

    Examples: Aspirin->ASR, Belumosudil->BEL, Ibu->IBU, ADP->ADP
    """
    _FORBIDDEN = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "HOH", "WAT", "DOD", "ACE", "NME", "UNK", "PCA",
        "HYP", "CSD", "CSO", "CSW", "CSS", "CME", "MLY", "MLZ", "MSE",
        "SEP", "TPO", "PTR", "NLE", "NVA", "ABA", "AIB", "ORN", "CIT",
        "SAR", "BMT",
    }
    clean = re.sub(r"[^a-zA-Z]", "", ligand_name).upper()
    if len(clean) == 0:
        return "LIG"
    for start in range(len(clean) - 2):
        candidate = clean[start:start + 3]
        if candidate not in _FORBIDDEN:
            return candidate
    for skip in range(2, len(clean)):
        indices = list(range(0, min(len(clean), skip * 3), skip))[:3]
        if len(indices) == 3:
            candidate = "".join(clean[i] for i in indices)
            if candidate not in _FORBIDDEN:
                return candidate
    if len(clean) >= 2:
        candidate = clean[:2] + "1"
        if candidate not in _FORBIDDEN:
            return candidate
    return "LIG"


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="2D Protein-Ligand Interaction Diagram (from PLIP data)")
    parser.add_argument("--json", required=True,
                        help="Path to interaction_summary.json from PLIP")
    parser.add_argument("--sdf", required=True,
                        help="Path to ligand SDF file")
    parser.add_argument("--complex", required=True,
                        help="Path to complex_for_plip.pdb (for atom mapping)")
    parser.add_argument("--output-dir", default="interaction_analysis",
                        help="Output directory")
    parser.add_argument("--protein-name", default="Protein",
                        help="Protein name for display")
    parser.add_argument("--ligand-name", default="Ligand",
                        help="Ligand name for display")
    parser.add_argument("--ligand-resname", default=None,
                        help="PDB residue name in complex (auto-generated if omitted)")
    parser.add_argument("--ligand-chain", default="L",
                        help="PDB chain for ligand in complex")
    parser.add_argument("--best-affinity", type=float, default=0.0,
                        help="Best docking affinity (kcal/mol)")
    args = parser.parse_args()

    resname = args.ligand_resname or generate_resname(args.ligand_name)

    print("=" * 60)
    print("2D INTERACTION DIAGRAM")
    print(f"  Protein:   {args.protein_name}")
    print(f"  Ligand:    {args.ligand_name}  (resname: {resname})")
    print(f"  Affinity:  {args.best_affinity:.3f} kcal/mol")
    print("=" * 60)

    # Load PLIP data
    if not os.path.exists(args.json):
        print(f"ERROR: interaction_summary.json not found: {args.json}")
        sys.exit(1)
    with open(args.json, "r") as f:
        plip_data = json.load(f)
    interactions = plip_data.get("interactions", {})

    # Load ligand
    if not os.path.exists(args.sdf):
        print(f"ERROR: ligand SDF not found: {args.sdf}")
        sys.exit(1)
    # RDKit's filename-based suppliers can fail on non-ASCII Windows paths.
    # Let Python open the Unicode path and pass a binary stream to RDKit.
    with open(args.sdf, "rb") as sdf_handle:
        supplier = Chem.ForwardSDMolSupplier(sdf_handle, removeHs=False)
        mol = next((candidate for candidate in supplier if candidate is not None), None)
    if mol is None:
        print("ERROR: Cannot parse ligand SDF")
        sys.exit(1)
    print(f"  Ligand: {mol.GetNumAtoms()} heavy atoms")

    # Build atom index mapping
    atom_map, lig_count = build_atom_map(args.complex, resname, args.ligand_chain)
    if lig_count == 0:
        print("WARNING: No ligand atoms found in complex PDB")

    # Draw ligand
    print("\n[1/4] Drawing ligand structure...")
    lig_img, lig_mol = draw_ligand_image(mol)

    # Compute label positions
    print("\n[2/4] Computing residue label positions...")
    residues = place_residue_labels(
        interactions, atom_map, lig_mol,
        LIGAND_BOX_W, LIGAND_BOX_H,
        LIGAND_CENTER[0], LIGAND_CENTER[1])

    # Render final diagram
    print("\n[3/4] Rendering 2D diagram...")
    os.makedirs(args.output_dir, exist_ok=True)
    out_png = os.path.join(args.output_dir, "interaction_2D.png")
    out_svg = os.path.join(args.output_dir, "interaction_2D.svg")
    out_emf = os.path.join(args.output_dir, "interaction_2D.emf")

    draw_interaction_diagram(
        lig_img, lig_mol, residues,
        args.protein_name, args.ligand_name,
        args.best_affinity,
        out_png, out_svg, out_emf)

    # Summary
    print(f"\n[4/4] Done! Output files:")
    for p in [out_png, out_svg, out_emf]:
        if os.path.exists(p):
            sz = os.path.getsize(p)
            print(f"  {os.path.basename(p):30s} {sz:>10,} bytes")
        else:
            print(f"  {os.path.basename(p):30s} MISSING")

    kinds = {}
    for r in residues:
        k = r["interaction_kind"]
        kinds[k] = kinds.get(k, 0) + 1
    print(f"\n  Residue labels: {len(residues)}")
    for k, v in sorted(kinds.items()):
        print(f"    {k}: {v}")

    print("=" * 60)


if __name__ == "__main__":
    main()
