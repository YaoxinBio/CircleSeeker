#!/usr/bin/env python3
"""Check release metadata and equality to the frozen computational source."""
import argparse
import ast
import hashlib
import json
import re
import subprocess
from pathlib import Path

SOURCE_COMMIT = "44c79b8eafdb1f60f949d69c9eeab12a69ce570b"
REPOSITORY = "https://github.com/YaoxinBio/CircleSeeker"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]

    def git(*items: str) -> str:
        return subprocess.check_output(["git", "-C", str(root), *items], text=True).strip()

    if args.output and git("status", "--porcelain", "--untracked-files=normal"):
        raise SystemExit("Commit release inputs before writing source provenance")

    tree = ast.parse((root / "src/circleseeker/__version__.py").read_text())
    versions = [ast.literal_eval(node.value) for node in tree.body
                if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == "__version__"
                        for target in node.targets)]
    if versions != ["1.5.2"]:
        raise SystemExit("Unexpected release version")
    version = versions[0]
    for filename, pattern in [
        ("conda-recipe/meta.yaml", r'{% set version = "' + re.escape(version) + r'" %}'),
        ("CITATION.cff", r"^version: " + re.escape(version) + r"$"),
        ("README.md", r"version-" + re.escape(version) + r"-blue.svg"),
    ]:
        if not re.search(pattern, (root / filename).read_text(), re.MULTILINE):
            raise SystemExit(f"Version metadata mismatch: {filename}")
    for filename in ["pyproject.toml", "README.md", "CITATION.cff", "conda-recipe/meta.yaml"]:
        contents = (root / filename).read_text()
        if REPOSITORY not in contents or "github.com/leoxqy/CircleSeeker" in contents:
            raise SystemExit(f"Active repository URL mismatch: {filename}")
    current_files = set(git("ls-files", "src/circleseeker").splitlines())
    reference_files = set(git("ls-tree", "-r", "--name-only", SOURCE_COMMIT, "src/circleseeker").splitlines())
    if current_files != reference_files:
        raise SystemExit("Computational source file set changed")
    hashes = {}
    for name in sorted(current_files):
        if name == "src/circleseeker/__version__.py":
            continue
        reference = subprocess.check_output(["git", "-C", str(root), "show", f"{SOURCE_COMMIT}:{name}"])
        current = (root / name).read_bytes()
        if current != reference:
            raise SystemExit(f"Computational source changed: {name}")
        hashes[name] = hashlib.sha256(current).hexdigest()
    record = {"version": version, "repository": REPOSITORY,
              "release_commit": git("rev-parse", "HEAD"),
              "computational_source_commit": SOURCE_COMMIT,
              "excluded_metadata_file": "src/circleseeker/__version__.py",
              "algorithm_source_identical": True, "source_sha256": hashes}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2) + "\n")
    print(f"Release {version}: metadata consistent; {len(hashes)} source files match {SOURCE_COMMIT}.")


if __name__ == "__main__":
    main()
