"""Read every reachable commit and ref; write findings without changing history."""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

from repository_guard import entries, git, scan


def audit(repo: Path) -> dict:
    previous = Path.cwd()
    os.chdir(repo)
    try:
        refs = {}
        for line in git("for-each-ref", "--format=%(refname) %(objectname)").decode().splitlines():
            ref, sha = line.split()
            refs[ref] = {"sha": sha, "commits": set(git("rev-list", ref).decode().split())}
        commits = git("rev-list", "--all").decode().split()
        by_blob_path = {}
        tree_cache = {}
        for sha in commits:
            tree = git("rev-parse", f"{sha}^{{tree}}").decode().strip()
            if tree not in tree_cache:
                tree_cache[tree] = scan(entries(sha))
            for item in tree_cache[tree]:
                key = (item["path"], item["blob"])
                found = by_blob_path.setdefault(key, {**item, "commits": []})
                found["commits"].append(sha)
        findings = []
        for item in by_blob_path.values():
            item["reachable_from"] = [ref for ref, info in refs.items() if set(item["commits"]) & info["commits"]]
            findings.append(item)
        return {
            "refs": {ref: info["sha"] for ref, info in refs.items()},
            "commit_count": len(commits), "finding_count": len(findings),
            "findings": sorted(findings, key=lambda x: (x["path"], x["blob"])),
            "scope": "all refs returned by the complete Git mirror; unreachable/deleted remote objects are not observable",
            "action": "read-only; no deletions, history rewriting or force pushes",
        }
    finally:
        os.chdir(previous)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    report = audit(args.repo.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "findings"}, ensure_ascii=False))


if __name__ == "__main__":
    main()
