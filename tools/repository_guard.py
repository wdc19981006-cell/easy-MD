"""Check Git index/trees, including force-added files; never print secret values."""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import PurePosixPath

BLOCKED_SUFFIXES = set("""
.pdb .ent .cif .mmcif .bcif .mmtf .pdbqt .sdf .mol .mol2 .xyz .gro .crd
.mae .maegz .pse .psw .pmo .dcd .xtc .trr .tpr .nc .h5 .hdf5 .npy .npz
.png .jpg .jpeg .gif .webp .bmp .tif .tiff .svg .mp4 .avi .mov
.zip .gz .bz2 .xz .7z .rar .tar .tgz .zst
.exe .dll .so .dylib .pyd .pyc .pyo .whl .egg .pkl .pickle .bin
""".split())
BLOCKED_DIRS = {
    ".git", ".venv", "venv", ".ppi-venv", ".ppi-tools", "__pycache__",
    ".pytest_cache", ".mypy_cache", ".repository-admin", "node_modules",
    "results", "ppi_results", "runs", "work", "raw", "prepared",
    "credentials", ".aws", ".ssh", "site-packages",
}
BLOCKED_NAMES = {"settings.local.json", "credentials.json", "credentials",
                 "id_rsa", "id_ed25519", "result-location-map.json", ".netrc"}
SECRET_PATTERNS = (
    re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9_-]{20,})"),
    re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    re.compile(rb"AKIA[0-9A-Z]{16}"),
    re.compile(rb"https?://[^\s/:]+:[^\s/@]+@", re.I),
)
ASSIGNMENT = re.compile(
    r'''(?im)^\s*["']?(?:api[_-]?key|access[_-]?token|auth[_-]?token|password|secret[_-]?key)["']?\s*[:=]\s*["']?([^\s"',;#}]+)'''
)
PERSONAL_PATH = re.compile(r'''(?i)(?:[a-z]:[/\\]Users[/\\]|/home/)([^/\\\s<>"'()\[\]{}|;]+)''')


def git(*args: str, data: bytes | None = None) -> bytes:
    return subprocess.run(["git", *args], input=data, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True).stdout


def path_issues(path: str) -> list[str]:
    p = PurePosixPath(path.replace("\\", "/").lower())
    issues = []
    if any(part in BLOCKED_DIRS or part.startswith((".public-release-", ".github-upload-"))
           for part in p.parts[:-1]):
        issues.append("local environment or raw-result directory")
    if any(s in BLOCKED_SUFFIXES for s in p.suffixes):
        issues.append("structure, image, archive, trajectory or binary file")
    name = p.name
    if name in BLOCKED_NAMES or name.startswith(("config.local.", "config.personal.")):
        issues.append("personal configuration or credentials")
    if name == ".env" or (name.startswith(".env.") and name != ".env.example"):
        issues.append("environment secrets file")
    if name.endswith((".log", ".out", ".err")):
        issues.append("raw log; put reviewed, redacted text in development/records instead")
    if name.endswith((".pem", ".key", ".p12", ".pfx")):
        issues.append("credential or key file")
    return issues


def content_issues(data: bytes) -> list[str]:
    issues = []
    if len(data) > 2_000_000:
        issues.append("file exceeds 2 MB; store raw outputs locally")
    if b"\0" in data:
        issues.append("binary content")
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        return issues + ["non-UTF-8 content"]
    if any(pattern.search(data) for pattern in SECRET_PATTERNS):
        issues.append("credential-like content (value redacted)")
    for match in ASSIGNMENT.finditer(text):
        value = match.group(1).lower()
        if value not in {"null", "none", "false", "true", "redacted", "changeme", "example"} and not value.startswith(("<", "${", "os.", "environ", "your_")):
            issues.append("non-placeholder credential assignment (value redacted)")
            break
    for match in PERSONAL_PATH.finditer(text):
        if match.group(1).lower() not in {"public", "default", "example", "username", "user", "yourname"}:
            issues.append("personal user-directory path; redact before committing")
            break
    # Detect actual coordinate files renamed to .txt/.csv/.json.
    if len(re.findall(r"(?m)^(?:ATOM  |HETATM)\s*\d+\s+", text)) >= 3 or re.search(r"(?m)^_atom_site\.(?:Cartn_[xyz]|label_atom_id)\b", text):
        issues.append("coordinate content disguised as text")
    if text.startswith("version https://git-lfs.github.com/spec/v1\n"):
        issues.append("LFS pointer; structures/binaries are not permitted through LFS")
    return issues


def entries(ref: str | None = None) -> list[tuple[str, str, str]]:
    if ref is None:
        rows = git("ls-files", "--stage", "-z").split(b"\0")
        result = []
        for row in filter(None, rows):
            meta, path = row.split(b"\t", 1)
            mode, sha, stage = meta.decode().split()
            if stage != "0":
                raise ValueError("Unresolved merge entries")
            result.append((path.decode("utf-8"), sha, mode))
        return result
    result = []
    for row in filter(None, git("ls-tree", "-r", "-z", ref).split(b"\0")):
        meta, path = row.split(b"\t", 1)
        mode, kind, sha = meta.decode().split()
        result.append((path.decode("utf-8"), sha, mode))
    return result


def scan(items: list[tuple[str, str, str]]) -> list[dict]:
    findings = []
    cache = {}
    for path, sha, mode in items:
        reasons = path_issues(path)
        if mode not in {"100644", "100755"}:
            reasons.append("symlink, submodule or non-regular file")
        else:
            if sha not in cache:
                cache[sha] = content_issues(git("cat-file", "blob", sha))
            reasons.extend(cache[sha])
        if reasons:
            findings.append({"path": path, "blob": sha, "reasons": sorted(set(reasons))})
    return findings


def pre_push(remote: str) -> int:
    commits = set()
    for line in sys.stdin:
        local_ref, local_sha, remote_ref, remote_sha = line.split()
        if set(local_sha) == {"0"}:
            continue
        if set(remote_sha) == {"0"}:
            args = (local_sha, "--not", f"--remotes={remote}")
        else:
            args = (f"{remote_sha}..{local_sha}",)
        commits.update(git("rev-list", *args).decode().split())
    failed = False
    for sha in sorted(commits):
        findings = scan(entries(sha))
        for item in findings:
            print(f"BLOCK {sha[:12]} {item['path']}: {', '.join(item['reasons'])}", file=sys.stderr)
        failed |= bool(findings)
    if not failed:
        print(f"Repository guard: {len(commits)} outgoing commit(s) checked")
    return int(failed)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--ref", help="Check the full tree of a commit (CI)")
    group.add_argument("--range", dest="commit_range", help="Check every tree in BASE..HEAD")
    group.add_argument("--pre-push", metavar="REMOTE", help="Read Git pre-push updates on stdin")
    args = parser.parse_args()
    if args.pre_push:
        return pre_push(args.pre_push)
    if args.commit_range:
        failed = False
        for sha in git("rev-list", args.commit_range).decode().split():
            findings = scan(entries(sha))
            for item in findings:
                print(f"BLOCK {sha[:12]} {item['path']}: {', '.join(item['reasons'])}", file=sys.stderr)
            failed |= bool(findings)
        if not failed:
            print("Repository guard: commit range passed")
        return int(failed)
    findings = scan(entries(args.ref))
    for item in findings:
        print(f"BLOCK {item['path']}: {', '.join(item['reasons'])}", file=sys.stderr)
    if not findings:
        print("Repository guard: full Git tree/index passed")
    return int(bool(findings))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (subprocess.CalledProcessError, ValueError) as exc:
        print(f"Repository guard could not complete: {type(exc).__name__}", file=sys.stderr)
        raise SystemExit(2)
