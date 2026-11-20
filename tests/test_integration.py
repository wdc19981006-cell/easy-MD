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

import pipeline
from box_prepare import calculate_box
from common import atomic_write_json, new_status
from visualization.pymol_generator import generate_pml, run_pymol
from visualization.visualization_validate import validate


class BoxIntegrationTests(unittest.TestCase):
    def test_box_sizes_are_computed_per_axis(self):
        import numpy as np

        coordinates = np.asarray([[0.0, 0.0, 0.0], [2.0, 4.0, 6.0]])
        _, sizes = calculate_box(coordinates, padding=1.0, minimum=0.0, maximum=99.0)
        self.assertEqual(sizes.tolist(), [4.0, 6.0, 8.0])


class PipelineIntegrationTests(unittest.TestCase):
    def test_execute_requires_confirmation_before_reading_run(self):
        result = pipeline.execute(SimpleNamespace(run_dir="missing", confirm=False))
        self.assertEqual(result, 2)

    def test_mocked_execute_completes_every_stage(self):
        with tempfile.TemporaryDirectory() as temp:
            run_dir = Path(temp) / "PTGS2_Aspirin"
            for name in ("docking", "structure", "analysis", "visualization"):
                (run_dir / name).mkdir(parents=True, exist_ok=True)
            config = Path(temp) / "config.yaml"
            config.write_text(
                "paths:\n  python: null\n  obabel: obabel\n  pymol: pymol\n  pymol_python: pymol-python\n"
                "docking:\n  seeds: [42, 2024, 8888]\ninteraction:\n  pymol_timeout: 5\n"
                "validation:\n  min_pse_bytes: 16\n",
                encoding="utf-8",
            )
            metadata = {
                "config": str(config), "redock": False,
                "protein": {"gene_symbol": "PTGS2"}, "docking": {"seeds": [42, 2024, 8888]},
                "ligand": {"name": "Aspirin"},
            }
            atomic_write_json(run_dir / "metadata.json", metadata)
            atomic_write_json(run_dir / "run_status.json", new_status())
            for path in (
                run_dir / "structure" / "receptor.pdbqt",
                run_dir / "structure" / "ligand.pdbqt",
                run_dir / "docking" / "config.txt",
                run_dir / "structure" / "original.pdb",
                run_dir / "structure" / "protein_clean.pdb",
                run_dir / "structure" / "ligand.sdf",
            ):
                path.write_text("fixture\n", encoding="utf-8")

            def fake_run(command, cwd, log=None):
                self.assertIn("--structured", command)
                if log:
                    log.parent.mkdir(parents=True, exist_ok=True)
                    log.write_text("mocked\n", encoding="utf-8")
                for number in (1, 2, 3):
                    folder = run_dir / "docking" / "runs" / f"run{number}"
                    folder.mkdir(parents=True)
                    (folder / "docked.pdbqt").write_text("MODEL 1\nREMARK VINA RESULT: -7.0 0 0\nENDMDL\n", encoding="utf-8")
                    (folder / "log.txt").write_text("   1 -7.0 0 0\n", encoding="utf-8")
                return SimpleNamespace(returncode=0, stdout="", stderr="")

            def fake_analyze(*_args, **_kwargs):
                master = run_dir / "analysis" / "interaction_master.json"
                summary = run_dir / "analysis" / "interaction_summary.txt"
                summary_json = run_dir / "analysis" / "interaction_summary.json"
                xml = run_dir / "analysis" / "plip_report.xml"
                visualizer = run_dir / "analysis" / "visualizer_data.pkl"
                complex_pdb = run_dir / "structure" / "complex_for_visualization.pdb"
                atomic_write_json(master, {"interaction_count": 0, "interactions": []})
                summary.write_text("PLIP interactions: 0\n", encoding="utf-8")
                atomic_write_json(summary_json, {"interaction_count": 0})
                xml.write_text("<report/>\n", encoding="utf-8")
                visualizer.write_bytes(b"fixture")
                complex_pdb.write_text("END\n", encoding="utf-8")
                (complex_pdb.parent / "best_pose.pdb").write_text("END\n", encoding="utf-8")
                return {
                    "master": str(master), "summary_txt": str(summary),
                    "summary_json": str(summary_json), "xml": str(xml),
                    "visualizer_data": str(visualizer), "complex": str(complex_pdb),
                    "interaction_count": 0, "ligand_resname": "SPI",
                }

            def fake_pml(_master, _complex, pml, pse):
                pml.write_text("quit\n", encoding="utf-8")
                return {"rendered_interaction_count": 0, "pse": str(pse)}

            def fake_native(_data, _master, pse, *_args, **_kwargs):
                pse.write_bytes(b"PSE" * 100)
                return {"success": True, "segment_count": 0}

            with (
                # This is a fully mocked orchestration test. Importing a real
                # PLIP installation would make it depend on the workstation.
                patch.dict(sys.modules, {"plip": SimpleNamespace(__file__=str(ROOT / "tests" / "mock_plip.py"))}),
                patch.object(pipeline, "_run", side_effect=fake_run),
                patch.object(pipeline, "analyze", side_effect=fake_analyze),
                patch.object(pipeline, "generate_pml", side_effect=fake_pml),
                patch.object(pipeline, "render_native", side_effect=fake_native),
                patch.object(pipeline, "validate_visualization", return_value={"success": True, "rendered_segment_count": 0}),
            ):
                result = pipeline.execute(SimpleNamespace(run_dir=str(run_dir), confirm=True))

            self.assertEqual(result, 0)
            self.assertEqual(sorted(p.name for p in run_dir.iterdir() if p.is_dir()),
                             ["analysis", "docking", "structure", "visualization"])
            self.assertFalse(list(run_dir.rglob("interaction_2*")))
            self.assertEqual(list(run_dir.rglob("best_pose.pdbqt")), [run_dir / "structure" / "best_pose.pdbqt"])
            status = json.loads((run_dir / "run_status.json").read_text(encoding="utf-8"))
            for stage in ("docking", "interaction", "visualization", "packaging"):
                self.assertEqual(status["stages"][stage]["status"], "SUCCESS")


class PyMOLRuntimeIntegrationTests(unittest.TestCase):
    def test_fallback_preserves_receptor_heterogen_and_checks_distance_endpoints(self):
        pymol_python = Path("D:/Interaction/PyMol/python.exe")
        if not pymol_python.is_file():
            self.skipTest("PyMOL runtime is not installed")
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            complex_pdb = root / "complex.pdb"
            complex_pdb.write_text(
                "HETATM    1  O2  FUC A  10       0.000   0.000   0.000  1.00 20.00           O\n"
                "HETATM    2  C1  SPI L   1       2.000   0.000   0.000  1.00 20.00           C\nEND\n",
                encoding="utf-8",
            )
            master = root / "master.json"
            report = {
                "interaction_count": 1,
                "ligand_site": {"resname": "SPI", "chain": "L", "residue_number": 1},
                "interactions": [{
                    "interaction_id": "I0001", "interaction_type": "hydrogen_bond",
                    "protein_chain": "A", "protein_residue_number": "10", "protein_residue_name": "FUC",
                    "protein_coordinate": [0.0, 0.0, 0.0],
                    "ligand_coordinate": [2.0, 0.0, 0.0],
                    "ring_center": {}, "charge_center": {},
                }],
            }
            atomic_write_json(master, report)
            pml, pse = root / "view.pml", root / "view.pse"
            generate_pml(master, complex_pdb, pml, pse)
            process = run_pymol(str(pymol_python), pml)
            checked = validate(master, pse, pml, str(pymol_python), min_pse_bytes=100, complex_pdb=complex_pdb)
            report["interactions"][0]["ligand_coordinate"] = [3.0, 0.0, 0.0]
            atomic_write_json(master, report)
            mismatched = validate(master, pse, pml, str(pymol_python), min_pse_bytes=100)
        self.assertEqual(process.returncode, 0)
        self.assertTrue(checked["success"], checked)
        self.assertEqual(checked["contact_residues"], [["A", "10", "FUC"]])
        self.assertEqual(checked["rendered_segment_count"], 1)
        self.assertFalse(mismatched["checks"]["interactions_match_plip_geometry"])


if __name__ == "__main__":
    unittest.main()
