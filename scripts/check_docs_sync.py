#!/usr/bin/env python3
"""
CircleSeeker Bilingual Documentation Sync Checker

This script checks that bilingual documentation pairs are in sync.
It verifies:
1. Each Chinese document has a corresponding English version
2. Modification times are reasonably close
3. With --interfaces, compares documented commands, defaults and output columns with the installation
"""

import argparse
import re
import shlex
import sys
from pathlib import Path
from datetime import datetime, timedelta

# Define bilingual document pairs (Chinese -> English)
BILINGUAL_PAIRS = {
    "CLI_Reference.md": "CLI_Reference_en.md",
    "Pipeline_Modules.md": "Pipeline_Modules_en.md",
    "Configuration_Reference.md": "Configuration_Reference_en.md",
    "Output_Format_Reference.md": "Output_Format_Reference_en.md",
    "Validation_Methodology.md": "Validation_Methodology_en.md",
    "UMC_Classification_Model.md": "UMC_Classification_Model_en.md",
    "SplitReads_Core_Algorithm.md": "SplitReads_Core_Algorithm_en.md",
    "Runtime_Tmp_Dir.md": "Runtime_Tmp_Dir_en.md",
}

# Maximum allowed time difference between paired documents (in days)
MAX_TIME_DIFF_DAYS = 7


def get_docs_dir() -> Path:
    """Get the docs directory path."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    return project_root / "docs"


def check_file_exists(docs_dir: Path, filename: str) -> bool:
    """Check if a file exists in the docs directory."""
    return (docs_dir / filename).exists()


def get_mtime(docs_dir: Path, filename: str) -> datetime:
    """Get the modification time of a file."""
    filepath = docs_dir / filename
    return datetime.fromtimestamp(filepath.stat().st_mtime)


def check_bilingual_sync() -> tuple[list[str], list[str], list[str]]:
    """
    Check bilingual documentation sync status.

    Returns:
        Tuple of (missing_files, out_of_sync, synced) document lists
    """
    docs_dir = get_docs_dir()
    missing_files = []
    out_of_sync = []
    synced = []

    for zh_doc, en_doc in BILINGUAL_PAIRS.items():
        missing = [name for name in (zh_doc, en_doc) if not check_file_exists(docs_dir, name)]
        if missing:
            missing_files.extend(missing)
            continue

        # Check modification time difference
        zh_mtime = get_mtime(docs_dir, zh_doc)
        en_mtime = get_mtime(docs_dir, en_doc)
        time_diff = abs(zh_mtime - en_mtime)

        if time_diff > timedelta(days=MAX_TIME_DIFF_DAYS):
            out_of_sync.append(
                f"{zh_doc} <-> {en_doc} "
                f"(diff: {time_diff.days} days)"
            )
        else:
            synced.append(f"{zh_doc} <-> {en_doc}")

    return missing_files, out_of_sync, synced


def check_cli_examples(doc: Path, cli) -> list[str]:
    """Check shell examples against Click without executing a pipeline."""
    import click

    errors = []
    blocks = re.findall(r"^```(?:bash|sh)\n(.*?)^```$", doc.read_text(), re.M | re.S)
    for line in "\n".join(blocks).replace("\\\n", " ").splitlines():
        line = line.strip()
        if not re.match(r"circleseeker\b", line, re.I):
            continue
        program, *args = shlex.split(line, comments=True)
        command = cli.commands.get(args[0]) if args and not args[0].startswith("-") else cli
        if program != "circleseeker" or command is None:
            errors.append(f"{doc.name}: unknown command in {line}")
            continue
        flags = {
            flag for param in command.get_params(click.Context(command))
            if isinstance(param, click.Option) for flag in param.opts
        }
        for arg in args:
            if arg.startswith("--") and arg.partition("=")[0] not in flags:
                errors.append(f"{doc.name}: unknown option {arg}")
    return errors


def check_interfaces(docs_dir: Path) -> list[str]:
    """Compare the two references with the installed CLI, defaults and CSV headers."""
    import pandas as pd
    from click.testing import CliRunner
    from circleseeker.cli.main import cli
    from circleseeker.config import CeccBuildConfig
    from circleseeker.modules.ecc_output_formatter import (
        generate_reads_table, generate_regions_table, generate_summary_table,
    )

    errors = []
    defaults = CeccBuildConfig()
    public_options = {
        flag for opt in cli.params if not opt.hidden for flag in opt.opts if flag.startswith("--")
    }
    sample = pd.DataFrame([{
        "eccDNA_id": "UeccDNA0001", "eccDNA_type": "UeccDNA",
        "chr": "chr1", "start0": 2, "end0": 6, "read_name": "r1",
    }])
    regions = generate_regions_table(uecc_df=sample)
    tables = (generate_summary_table(sample, regions), regions, generate_reads_table(uecc_df=sample))
    headers = {",".join(table.columns) for table in tables}
    for suffix in ("", "_en"):
        output_doc = docs_dir / f"Output_Format_Reference{suffix}.md"
        documented_headers = set(re.findall(r"^```csv\n([^\n]+)\n```$", output_doc.read_text(), re.M))
        for header in sorted(headers - documented_headers):
            errors.append(f"{output_doc.name}: missing current CSV header {header}")

        cli_doc = docs_dir / f"CLI_Reference{suffix}.md"
        cli_text = cli_doc.read_text()
        documented_flags = set(re.findall(r"--[a-z][a-z-]*", cli_text))
        for flag in sorted(public_options - documented_flags):
            errors.append(f"{cli_doc.name}: missing public option {flag}")
        for name in cli.commands:
            if f"`{name}`" not in cli_text:
                errors.append(f"{cli_doc.name}: missing subcommand {name}")

        config_doc = docs_dir / f"Configuration_Reference{suffix}.md"
        for doc in (cli_doc, config_doc):
            errors.extend(check_cli_examples(doc, cli))

        model_doc = docs_dir / f"UMC_Classification_Model{suffix}.md"
        for doc in (config_doc, model_doc):
            rows = [line.strip().strip("|").split("|") for line in doc.read_text().splitlines()
                    if line.startswith("| `")]
            for key in ("min_match_degree", "position_tolerance", "edge_tolerance"):
                expected = str(getattr(defaults, key))
                values = {
                    row[2].strip().strip("`").removesuffix(" bp") for row in rows
                    if len(row) >= 3 and row[0].strip() == f"`{key}`"
                }
                if values != {expected}:
                    errors.append(f"{doc.name}: expected {key} default {expected}")

        template_cmd = "circleseeker init-config --output-file config.yaml"
        if template_cmd not in config_doc.read_text().splitlines():
            errors.append(f"{config_doc.name}: missing template command {template_cmd}")

    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli, ["init-config", "--output-file", "config.yaml"])
        if result.exit_code or not Path("config.yaml").is_file():
            errors.append(f"init-config failed to generate config.yaml: {result.output}")
    return errors


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interfaces", action="store_true",
                        help="Also check the public interfaces. Requires an installed CircleSeeker.")
    args = parser.parse_args()
    print("=" * 60)
    print("CircleSeeker Bilingual Documentation Sync Check")
    print("=" * 60)
    print()

    docs_dir = get_docs_dir()
    if not docs_dir.exists():
        print(f"ERROR: docs directory not found: {docs_dir}")
        return 1

    missing_files, out_of_sync, synced = check_bilingual_sync()

    if synced:
        print("Synced documents:")
        for doc in synced:
            print(f"  [OK] {doc}")
        print()

    # Report out-of-sync documents (warning)
    if out_of_sync:
        print("WARNING: Out-of-sync documents (modification time differs):")
        for doc in out_of_sync:
            print(f"  [WARN] {doc}")
        print()

    if missing_files:
        print("ERROR: Missing reference files:")
        for doc in missing_files:
            print(f"  [MISSING] {doc}")
        print()

    interface_errors = check_interfaces(docs_dir) if args.interfaces and not missing_files else []
    for error in interface_errors:
        print(f"  [ERROR] {error}")
    if args.interfaces and not missing_files and not interface_errors:
        print("  [OK] CLI examples, CeccBuild defaults and packaged CSV headers")

    # Summary
    print("-" * 60)
    print(f"Summary: {len(synced)} synced, {len(out_of_sync)} warnings, "
          f"{len(missing_files) + len(interface_errors)} errors")

    if missing_files or interface_errors:
        print("\nFAILED: Documentation checks failed.")
        return 1

    if out_of_sync:
        print("\nPASSED with warnings: Some documents may need sync update.")
        return 0

    print("\nPASSED: All bilingual documentation is in sync!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
