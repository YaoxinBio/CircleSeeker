#!/usr/bin/env python3
"""Check release metadata and the diff against the analysis baseline.

1.5.2 was prepared as a packaging-only release: its computational source was
byte-identical to the commit the completed analyses ran on, and this script
asserted exactly that. That is no longer true. Two defects were found that made
deliverables vary between runs of identical code on identical input, and fixing
them - together with the scale work that made GlioSarc_P01_Tumor finish - changed
eleven source files.

Asserting identity would now be asserting something false, so the check states
the truth instead: the files that differ from the baseline must be exactly the
ones declared here. An unintended change still fails the gate; the intended ones
are named, and the written record says the source is not identical.
"""
import argparse
import ast
import hashlib
import json
import re
import subprocess
from pathlib import Path

ANALYSIS_BASELINE = "a662ad74e1cf984af0417017a36572e3a725c934"
REPOSITORY = "https://github.com/YaoxinBio/CircleSeeker"

# Files this release deliberately changes since the analysis baseline.
# See CHANGELOG.md and docs/releases/v1.5.2.md for what changed in each.
EXPECTED_CHANGES = {
    "src/circleseeker/core/pipeline.py",
    "src/circleseeker/modules/cecc_build.py",
    "src/circleseeker/modules/ecc_output_formatter.py",
    "src/circleseeker/modules/ecc_packager.py",
    "src/circleseeker/modules/ecc_summary.py",
    "src/circleseeker/modules/ecc_unify.py",
    "src/circleseeker/modules/splitreads_core.py",
    "src/circleseeker/modules/tandem_to_ring.py",
    "src/circleseeker/modules/um_classify.py",
    "src/circleseeker/modules/umc_process.py",
    "src/circleseeker/utils/read_support.py",
}


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
    reference_files = set(
        git("ls-tree", "-r", "--name-only", ANALYSIS_BASELINE, "src/circleseeker").splitlines()
    )
    if current_files != reference_files:
        added = sorted(current_files - reference_files)
        removed = sorted(reference_files - current_files)
        raise SystemExit(f"Computational source file set changed: +{added} -{removed}")
    hashes = {}
    changed = set()
    for name in sorted(current_files):
        current = (root / name).read_bytes()
        hashes[name] = hashlib.sha256(current).hexdigest()
        if name == "src/circleseeker/__version__.py":
            continue
        reference = subprocess.check_output(
            ["git", "-C", str(root), "show", f"{ANALYSIS_BASELINE}:{name}"]
        )
        if current != reference:
            changed.add(name)
    if changed != EXPECTED_CHANGES:
        unexpected = sorted(changed - EXPECTED_CHANGES)
        reverted = sorted(EXPECTED_CHANGES - changed)
        raise SystemExit(
            "Source differs from the analysis baseline in files this release does "
            f"not declare: unexpected={unexpected} no-longer-changed={reverted}"
        )
    record = {"version": version, "repository": REPOSITORY,
              "release_commit": git("rev-parse", "HEAD"),
              "analysis_baseline_commit": ANALYSIS_BASELINE,
              "algorithm_source_identical": False,
              "files_changed_since_analysis_baseline": sorted(changed),
              "source_sha256": hashes}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2) + "\n")
    print(f"Release {version}: metadata consistent; {len(hashes)} source files, "
          f"{len(changed)} changed since {ANALYSIS_BASELINE[:7]} as declared.")


if __name__ == "__main__":
    main()
