#!/usr/bin/env python3
"""Check release metadata and the diff against the analysis baseline.

This release was prepared as packaging-only: its computational source was
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

# Files this release deliberately changes since the analysis baseline, each
# pinned to the content that was reviewed. Naming the files alone would let any
# further edit inside one of them through - changing a default such as
# breakpoint_depth would still report "11 changed as declared" - so the hash is
# what the gate actually compares. Regenerate deliberately with
# --update-manifest after reviewing a change, and commit the result.
EXPECTED_CHANGES: dict[str, str] = {
    "src/circleseeker/core/pipeline.py": "b32044ee8cceceeeed9eed665f4aead172e0857b2f5368996a32dacf85045699",
    "src/circleseeker/modules/cecc_build.py": "7b673a3fe555b6edbc7ad141921feabf7eb3f6e9b0608e6974ea5b620f2b15f6",
    "src/circleseeker/modules/ecc_output_formatter.py": "c07cd8fd61935dc2c244e71a6b404bbe2fa01856d93f0b15f4ad33120c02a208",
    "src/circleseeker/modules/ecc_packager.py": "f3642ff047aabd8015403cf590fc075df7ff42d0402f3b1529ed3a0a6b0e5169",
    "src/circleseeker/modules/ecc_summary.py": "56c594bd9c3e1d0a245a35feab587f9afc30eaf3b1817b4ebf29acd105c8cd83",
    "src/circleseeker/modules/ecc_unify.py": "57b5858d73adb8ce61e16f0f2b1fa053a4fef23b05c04a30f83390105ea3e530",
    "src/circleseeker/modules/splitreads_core.py": "8a0cdbfdd5ba7aad9da7bf772ad89ce38d23a942a3e459aabad940ba936896fa",
    "src/circleseeker/modules/tandem_to_ring.py": "221bedeb1ead909c87c738763dbfe79359d0d115c5e03067a093e6795f2b6391",
    "src/circleseeker/modules/um_classify.py": "5d5576f2dd1da96e4c58284f3d07147ef0c14e2c5493c1483fe24b89b3d3743c",
    "src/circleseeker/modules/umc_process.py": "72d448cc0a249e2b4921a9f03e8df86ae3d113a46b625d874a22e0414d2743c4",
    "src/circleseeker/utils/read_support.py": "25296f2276b12688af582c7762ad2a4a685f4228ed10b9091d61202d84418669",
}


def _git(root, *items, text: bool = False):
    """Run git, avoiding macOS' /usr/bin/git shim where it cannot work.

    That shim resolves the real binary through xcrun, and xcrun cannot load its
    library across architectures - so an x86_64 interpreter under Rosetta fails
    to launch an arm64 git through it. Try real binaries first.
    """
    for exe in ("/opt/homebrew/bin/git",
                "/Library/Developer/CommandLineTools/usr/bin/git",
                "git"):
        try:
            return subprocess.check_output([exe, "-C", str(root), *items],
                                           stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            continue
    raise SystemExit("git " + " ".join(items) + " failed with every candidate binary")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--update-manifest",
        action="store_true",
        help="Rewrite the pinned hashes in this file to match the working tree. "
             "Use after reviewing a change to a declared file, then commit.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]

    def git(*items: str) -> str:
        return _git(root, *items).decode().strip()

    if args.output and git("status", "--porcelain", "--untracked-files=normal"):
        raise SystemExit("Commit release inputs before writing source provenance")

    tree = ast.parse((root / "src/circleseeker/__version__.py").read_text())
    versions = [ast.literal_eval(node.value) for node in tree.body
                if isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == "__version__"
                        for target in node.targets)]
    if versions != ["1.6.1"]:
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
        reference = _git(root, "show", f"{ANALYSIS_BASELINE}:{name}")
        if current != reference:
            changed.add(name)
    if args.update_manifest:
        me = Path(__file__)
        text = me.read_text()
        block = "\n".join(f'    "{name}": "{hashes[name]}",' for name in sorted(changed))
        text = re.sub(
            r"(EXPECTED_CHANGES: dict\[str, str\] = \{\n).*?(\n\})",
            lambda m: m.group(1) + block + m.group(2),
            text,
            count=1,
            flags=re.S,
        )
        me.write_text(text)
        print(f"Manifest updated: {len(changed)} files pinned. Review and commit.")
        return
    if changed != set(EXPECTED_CHANGES):
        unexpected = sorted(changed - set(EXPECTED_CHANGES))
        reverted = sorted(set(EXPECTED_CHANGES) - changed)
        raise SystemExit(
            "Source differs from the analysis baseline in files this release does "
            f"not declare: unexpected={unexpected} no-longer-changed={reverted}"
        )
    drifted = sorted(
        name for name, digest in EXPECTED_CHANGES.items() if hashes[name] != digest
    )
    if drifted:
        raise SystemExit(
            "Declared file changed since it was pinned - review it, then rerun with "
            f"--update-manifest: {drifted}"
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
