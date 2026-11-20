"""Exercise staged Git blobs and outgoing commits in isolated repositories."""
from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
GUARD = ROOT / "tools" / "repository_guard.py"


class RepositoryGuardTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name)
        self.git("init", "-b", "main")
        self.git("config", "user.name", "Guard Test")
        self.git("config", "user.email", "guard@example.invalid")
        self.write("README.md", "Synthetic guard test; no scientific computation.\n")
        self.git("add", "README.md")
        self.git("commit", "-m", "base")

    def git(self, *args):
        result = subprocess.run(["git", *args], cwd=self.repo,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
        return result.stdout.decode().strip()

    def write(self, name, text):
        path = self.repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def check_guard(self, *args, stdin=None):
        return subprocess.run([sys.executable, str(GUARD), *args], cwd=self.repo,
                              input=stdin, capture_output=True, text=True)

    def test_force_added_uppercase_structure_is_blocked_and_deletion_is_allowed(self):
        self.write(".gitignore", "*.PDB\n")
        self.write("PROTEIN.PDB", "not even coordinates\n")
        self.git("add", "-f", "PROTEIN.PDB")
        self.assertEqual(self.check_guard().returncode, 1)
        self.git("rm", "--cached", "PROTEIN.PDB")
        self.assertEqual(self.check_guard().returncode, 0)

    def test_reads_staged_secret_even_if_worktree_has_been_cleaned(self):
        token = "ghp_" + "X" * 40
        self.write("development/records/redacted.txt", token)
        self.git("add", "development")
        self.write("development/records/redacted.txt", "REDACTED\n")
        result = self.check_guard()
        self.assertEqual(result.returncode, 1)

        self.assertNotIn(token, result.stdout + result.stderr)

    def test_renamed_coordinate_file_is_blocked(self):
        self.write("coordinates.txt", "".join(f"ATOM  {i:5d}  CA  ALA A   1       0.000   0.000   0.000\n" for i in range(1, 4)))
        self.git("add", "coordinates.txt")
        self.assertEqual(self.check_guard().returncode, 1)

    def test_source_visualization_and_safe_score_summary_are_allowed(self):
        self.write("scripts/visualization/example.py", 'filename = "model.pdb"\n')
        self.write("development/records/summary.csv", "pdb_id,seed,DockQ,status\n1BRS,42,0.05,FAIL\n")
        self.git("add", ".")
        self.assertEqual(self.check_guard().returncode, 0)

    def test_nested_archive_and_environment_are_blocked(self):
        self.write("structure.CIF.GZ", "compressed data disguised as text\n")
        self.write(".env.local", "value=example\n")
        self.git("add", ".")
        self.assertEqual(self.check_guard().returncode, 1)

    def test_binary_under_a_text_suffix_is_blocked(self):
        (self.repo / "report.json").write_bytes(b"binary\x00payload")
        self.git("add", ".")
        self.assertEqual(self.check_guard().returncode, 1)

    def test_pre_push_rejects_bad_intermediate_commit_even_if_tip_removed_it(self):
        base = self.git("rev-parse", "HEAD")
        self.git("update-ref", "refs/remotes/origin/main", base)
        self.write("pose.sdf", "not permitted\n")
        self.git("add", ".")
        self.git("commit", "-m", "unsafe intermediate")
        self.git("rm", "pose.sdf")
        self.git("commit", "-m", "remove unsafe file")
        tip = self.git("rev-parse", "HEAD")
        self.assertEqual(self.check_guard("--ref", tip).returncode, 0)
        result = self.check_guard("--pre-push", "origin", stdin=f"refs/heads/main {tip} refs/heads/main {base}\n")
        self.assertEqual(result.returncode, 1)
        self.assertEqual(self.check_guard("--range", f"{base}..{tip}").returncode, 1)

    def test_personal_path_in_raw_transcript_requires_redaction(self):
        self.write("transcript.txt", "C:" + "/Users/" + "ResearcherPrivate/run\n")
        self.git("add", ".")
        self.assertEqual(self.check_guard().returncode, 1)

    def test_installed_pre_commit_hook_blocks_real_commit(self):
        shutil.copytree(ROOT / "tools", self.repo / "tools", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copytree(ROOT / ".githooks", self.repo / ".githooks")
        subprocess.run([sys.executable, str(self.repo / "tools/install_hooks.py")],
                       cwd=self.repo, capture_output=True, check=True)
        self.write("pose.PDBQT", "unsafe\n")
        self.git("add", ".")
        result = subprocess.run(["git", "commit", "-m", "must be blocked"],
                                cwd=self.repo, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("BLOCK pose.PDBQT", result.stderr)
        self.git("rm", "--cached", "pose.PDBQT")
        self.git("commit", "-m", "guarded safe source")


if __name__ == "__main__":
    unittest.main()
