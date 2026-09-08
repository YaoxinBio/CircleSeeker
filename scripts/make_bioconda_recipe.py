#!/usr/bin/env python3
"""Render the local recipe with a release sdist URL and its actual checksum."""
import argparse
import hashlib
import re
import tarfile
from email.parser import BytesParser
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdist", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    version_text = (root / "src/circleseeker/__version__.py").read_text()
    match = re.search(r'^__version__ = "([^"]+)"', version_text, re.MULTILINE)
    if match is None:
        parser.error("Cannot read the source version")
    version = match.group(1)
    if args.sdist.name != f"circleseeker-{version}.tar.gz":
        parser.error("The sdist filename must match the source version")
    with tarfile.open(args.sdist, "r:gz") as archive:
        member = archive.extractfile(f"circleseeker-{version}/PKG-INFO")
        if member is None:
            parser.error("The sdist lacks PKG-INFO")
        metadata = BytesParser().parsebytes(member.read())
        if metadata["Name"] != "circleseeker" or metadata["Version"] != version:
            parser.error("The sdist metadata does not match this checkout")
    checksum = hashlib.sha256(args.sdist.read_bytes()).hexdigest()
    recipe = (root / "conda-recipe/meta.yaml").read_text()
    source = (
        "source:\n"
        f"  url: https://github.com/YaoxinBio/CircleSeeker/releases/download/v{version}/{args.sdist.name}\n"
        f"  sha256: {checksum}\n"
    )
    recipe, replacements = re.subn(r"source:\n.*?(?=\nbuild:)", source, recipe, flags=re.DOTALL)
    if replacements != 1:
        parser.error("Expected exactly one source block")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(recipe)
    print(f"Wrote {args.output}; version={version}, sha256={checksum}")


if __name__ == "__main__":
    main()
