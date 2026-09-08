"""Preserve candidate-level RCA support before grouping reads or genomic segments."""

from __future__ import annotations

import json
import math
from typing import Any

import pandas as pd


def _names(row: pd.Series) -> list[str]:
    for col in ("reads", "read_name", "eReads"):
        value = row.get(col)
        if isinstance(value, str) and value.strip() not in ("", "nan", "NA"):
            return list(dict.fromkeys(x.strip() for x in value.split(";") if x.strip()))
    return []


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (ValueError, TypeError):
        return None
    return number if math.isfinite(number) and number >= 0 else None


def collect_support(df: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Union support by candidate, counting its repeated segment rows only once.

    Different repeat candidates from one read contribute separately, matching
    the candidate-support sum used by the pipeline. Repeated representations
    of one candidate contribute once. Missing historical per-read values stay
    unknown; a multi-read aggregate is never assigned to each of its reads.
    """
    support: dict[str, dict[str, Any]] = {}

    def add(key: str, read: str, copies: Any) -> None:
        entry = {"read": read, "copy_number": _number(copies)}
        if key in support and support[key] != entry:
            raise ValueError(f"Conflicting RCA support for candidate {key}")
        support[key] = entry

    for _, row in df.iterrows():
        encoded = row.get("candidate_support")
        if isinstance(encoded, str) and encoded.strip():
            for key, item in json.loads(encoded).items():
                add(key, str(item["read"]), item["copy_number"])
            continue
        reads = _names(row)
        if not reads:
            continue
        encoded = row.get("per_read_copy_number")
        values = encoded.split(";") if isinstance(encoded, str) else []
        if values and len(values) != len(reads):
            raise ValueError("RCA values do not align with supporting read names")
        query = row.get("query_id")
        has_query = isinstance(query, str) and query.strip() not in ("", "NA", "nan")
        candidate = str(query if has_query else row.get("eccDNA_id", "")).removesuffix("|circular")
        measured = next(
            (
                _number(row.get(col))
                for col in ("copy_number", "repeat_number", "eRepeatNum")
                if _number(row.get(col)) is not None
            ),
            None,
        )
        for i, read in enumerate(reads):
            copies = values[i] if values else (measured if len(reads) == 1 else None)
            key = candidate if has_query and len(reads) == 1 else f"{candidate}::{read}"
            add(key, read, copies)
    return support


def support_fields(df: pd.DataFrame, reads: str | None = None) -> dict[str, Any]:
    """Return serialized provenance and aligned per-read/aggregate copy counts."""
    support = collect_support(df)
    by_read: dict[str, list[float | None]] = {}
    for item in support.values():
        by_read.setdefault(item["read"], []).append(item["copy_number"])
    ordered = [x.strip() for x in reads.split(";") if x.strip()] if reads else sorted(by_read)
    values = {
        read: (
            sum(v for v in numbers if v is not None)
            if all(v is not None for v in numbers)
            else None
        )
        for read, numbers in by_read.items()
    }

    def encode(value: float | None) -> str:
        return "" if value is None else format(value, ".15g")

    total = (
        sum(v for v in values.values() if v is not None)
        if values and all(v is not None for v in values.values())
        else float("nan")
    )
    return {
        "candidate_support": json.dumps(support, sort_keys=True, separators=(",", ":")),
        "per_read_copy_number": ";".join(encode(values.get(read)) for read in ordered),
        "copy_number": total,
    }
