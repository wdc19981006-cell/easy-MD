from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analysis.interaction_parser import parse_interactions
from common import atomic_write_json, new_status, update_stage
from docking.pose_rank import rank_poses
from pdb_selection import rank_pdb_entries
from pipeline import _verified_modified_site
from protein.normalize import normalize_protein
from visualization.scene import interaction_segments


class NormalizeTests(unittest.TestCase):
    def test_cox2_alias_resolves_to_ptgs2(self):
        hits = [
            {"symbol": "PTGS2", "name": "prostaglandin-endoperoxide synthase 2", "aliases": ["COX2", "COX-2"]},
            {"symbol": "JMJD2D", "name": "lysine demethylase 4D", "aliases": []},
        ]
        with patch("protein.normalize.query_mygene", return_value=hits):
            result = normalize_protein("COX-2")
        self.assertTrue(result["resolved"])
        self.assertFalse(result["ambiguous"])
        self.assertEqual(result["selected"]["symbol"], "PTGS2")


class RankingTests(unittest.TestCase):
    def test_human_xray_liganded_structure_wins(self):
        good = {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.0]},
            "rcsb_accession_info": {"deposit_date": "2022-01-01"},
            "struct": {"title": "human PTGS2"},
            "ligands": ["AAS"], "is_human": True, "organism": "Homo sapiens",
            "easy_md": {"binding_site_defined": True, "mutations": False},
        }
        weak = {
            "exptl": [{"method": "ELECTRON MICROSCOPY"}],
            "rcsb_entry_info": {"resolution_combined": [3.5]},
            "rcsb_accession_info": {"deposit_date": "2010-01-01"},
            "struct": {"title": "other"}, "ligands": [], "is_human": False,
            "easy_md": {},
        }
        ranked = rank_pdb_entries({"GOOD": good, "WEAK": weak}, {})
        self.assertEqual(ranked[0]["pdb_id"], "GOOD")

    def test_relevant_ligand_and_active_site_evidence_are_scored(self):
        entry = {
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.4]},
            "rcsb_accession_info": {"deposit_date": "2015-01-01"},
            "struct": {"title": "Aspirin acetylated enzyme"},
            "ligands": ["OAS"], "is_human": True,
            "easy_md": {
                "binding_site_defined": True,
                "relevant_ligand_evidence": "entry title",
                "active_site_evidence": "verified OAS",
            },
        }
        ranked = rank_pdb_entries({"SITE": entry}, {})
        self.assertIn("Relevant ligand evidence: entry title", ranked[0]["positives"])
        self.assertIn("Active-site evidence: verified OAS", ranked[0]["positives"])


class SiteEvidenceTests(unittest.TestCase):
    def test_aspirin_site_requires_verified_oas_modres(self):
        with tempfile.TemporaryDirectory() as temp:
            pdb = Path(temp) / "site.pdb"
            pdb.write_text(
                "MODRES 5F19 OAS A  530  SER  MODIFIED RESIDUE\n",
                encoding="utf-8",
            )
            marker = _verified_modified_site(
                pdb, "Aspirin", "The Crystal Structure of Aspirin Acetylated Human COX-2"
            )
        self.assertEqual(marker, ("OAS", "A"))


class PoseRankingTests(unittest.TestCase):
    def test_global_best_pose_is_extracted(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for run, score in ((1, -7.0), (2, -8.1), (3, -7.4)):
                folder = root / f"run{run}"
                folder.mkdir()
                folder.joinpath("log.txt").write_text(
                    f"    1    {score:.1f}    0.000    0.000\n    2    -6.0    1.0    2.0\n",
                    encoding="utf-8",
                )
                folder.joinpath("docked.pdbqt").write_text(
                    "MODEL 1\nATOM      1  C   LIG L   1       0.000   0.000   0.000\nENDMDL\n"
                    "MODEL 2\nATOM      1  C   LIG L   1       1.000   0.000   0.000\nENDMDL\n",
                    encoding="utf-8",
                )
            result = rank_poses(root, [42, 2024, 8888])
            self.assertEqual(result["best"]["run"], 2)
            self.assertTrue((root / "best_pose.pdbqt").is_file())


class InteractionSchemaTests(unittest.TestCase):
    def test_hydrophobic_record_contains_names_and_coordinates(self):
        atom = lambda coords: SimpleNamespace(coords=coords)
        contact = SimpleNamespace(
            reschain="A", restype="LEU", resnr=10, distance=3.5,
            bsatom_orig_idx=1, ligatom_orig_idx=2,
            bsatom=atom([0, 0, 0]), ligatom=atom([1, 1, 1]),
        )
        empty = SimpleNamespace(
            hydrophobic_contacts=[contact], hbonds_pdon=[], hbonds_ldon=[],
            saltbridge_lneg=[], saltbridge_pneg=[], pistacking=[],
            pication_paro=[], pication_laro=[], halogen_bonds=[],
            water_bridges=[], metal_complexes=[],
        )
        with tempfile.TemporaryDirectory() as temp:
            pdb = Path(temp) / "complex.pdb"
            pdb.write_text(
                "ATOM      1  CD1 LEU A  10       0.000   0.000   0.000  1.00 20.00           C\n"
                "HETATM    2  C1  LIG L   1       1.000   1.000   1.000  1.00 20.00           C\n",
                encoding="utf-8",
            )
            master = parse_interactions(empty, pdb)
        row = master["interactions"][0]
        self.assertEqual(row["protein_atom_name"], "CD1")
        self.assertEqual(row["ligand_atom_name"], "C1")
        self.assertEqual(row["protein_coordinate"], [0.0, 0.0, 0.0])
        self.assertTrue(row["plip_source"])


class VisualizationTests(unittest.TestCase):
    def test_water_bridge_draws_through_reported_water(self):
        protein, water, ligand = [0, 0, 0], [1, 1, 1], [2, 0, 0]
        row = {
            "interaction_id": "I0001", "interaction_type": "water_bridge",
            "protein_coordinate": protein, "ligand_coordinate": ligand, "water_coordinate": water,
        }
        self.assertEqual(interaction_segments({"interactions": [row]}), [(protein, water), (water, ligand)])
        del row["water_coordinate"]
        with self.assertRaisesRegex(ValueError, "I0001"):
            interaction_segments({"interactions": [row]})


class StatusTests(unittest.TestCase):
    def test_downstream_failure_preserves_success(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "run_status.json"
            atomic_write_json(path, new_status())
            update_stage(path, "docking", "SUCCESS", artifacts=["best_pose.pdbqt"])
            update_stage(path, "visualization", "FAILED", error="PyMOL failed")
            status = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(status["stages"]["docking"]["status"], "SUCCESS")
        self.assertEqual(status["stages"]["visualization"]["status"], "FAILED")


if __name__ == "__main__":
    unittest.main()
