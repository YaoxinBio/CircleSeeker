"""Directed circular segment identity and oriented breakpoint coordinates."""

from __future__ import annotations

from typing import Any

import pandas as pd
import edlib

ORDER_COLUMNS = ("segment_in_circle", "segment_order", "seg_index", "region_idx")


def ordered_segments(df: pd.DataFrame) -> pd.DataFrame:
    """Preserve the supplied traversal; genomic coordinates are not a path order."""
    for col in ORDER_COLUMNS:
        if col in df.columns:
            order = pd.to_numeric(df[col], errors="coerce")
            if order.isna().any() or order.duplicated().any():
                raise ValueError(f"Missing or duplicate circular segment order: {col}")
            return df.loc[order.sort_values(kind="stable").index]
    return df


def cycle(df: pd.DataFrame) -> tuple[tuple[str, int, int, str], ...]:
    ordered = ordered_segments(df)
    return tuple(
        (str(r.chr), int(r.start0), int(r.end0), str(getattr(r, "strand", ".")))
        for r in ordered.itertuples()
    )


def reverse_cycle(segments: tuple) -> tuple:
    return tuple(
        (ch, start, end, {"+": "-", "-": "+"}.get(strand, "."))
        for ch, start, end, strand in reversed(segments)
    )


def canonical_cycle(segments: tuple) -> tuple:
    if not segments:
        return ()
    return min(
        order[i:] + order[:i]
        for order in (segments, reverse_cycle(segments))
        for i in range(len(order))
    )


def same_cycle(a: tuple, b: tuple, tolerance: int = 0) -> bool:
    """Require matching directed order up to ring rotation / whole-ring reversal."""
    if len(a) != len(b) or not a or any(s[3] not in ("+", "-") for s in a + b):
        return False
    for order in (b, reverse_cycle(b)):
        for offset in range(len(order)):
            rotated = order[offset:] + order[:offset]
            if all(
                x[0] == y[0]
                and x[3] == y[3]
                and abs(x[1] - y[1]) <= tolerance
                and abs(x[2] - y[2]) <= tolerance
                for x, y in zip(a, rotated)
            ):
                return True
    return False


def same_circular_sequence(a: str, b: str) -> bool:
    """Require at most 1% whole-circle edits, allowing rotation and reversal.

    A strict equality gate splits noisy consensus copies when a base error
    changes their canonical origin. Align one unit against a doubled target;
    charge unaligned/extra target span too, so a matching subfragment cannot
    establish whole-molecule identity. This gate is used only after directed
    coordinate compatibility has independently been established.
    """
    if not isinstance(a, str) or not isinstance(b, str) or not a or not b:
        return False
    a, b = a.upper(), b.upper()
    if a == "NAN" or b == "NAN":
        return False
    if set(a + b) - set("ACGTN") or not (set(a) & set("ACGT")):
        return False
    a, b = sorted((a, b), key=lambda seq: (len(seq), seq))
    reverse_b = b.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    if len(a) == len(b) and (a in b + b or a in reverse_b + reverse_b):
        return True
    if "N" in a or "N" in b:
        return False
    budget = max(len(a), len(b)) // 100
    if abs(len(a) - len(b)) > budget:
        return False
    for target in (b, reverse_b):
        result = edlib.align(a, target + target, mode="HW", task="locations", k=budget)
        distance = result["editDistance"]
        if distance < 0:
            continue
        for start, end in result["locations"]:
            if start < len(target) and distance + abs(end - start + 1 - len(target)) <= budget:
                return True
    return False


def junction_end(start: Any, end: Any, strand: str, *, outgoing: bool) -> tuple[int, int]:
    """One-base half-open interval at the traversal's entry or exit boundary."""
    start, end = int(start), int(end)
    if start < 0 or end <= start or strand not in ("+", "-"):
        raise ValueError("Invalid segment coordinates or strand for circular junction")
    right = (strand == "+") == outgoing
    return (end - 1, end) if right else (start, start + 1)


def validate_sequences(
    df: pd.DataFrame, *, id_col: str = "eccDNA_id", sequence_col: str = "eSeq"
) -> None:
    """Reject missing/conflicting sequence and length before publishing outputs."""
    if sequence_col not in df.columns:
        raise ValueError("Cecc output requires its representative sequence")
    for eid, group in df.groupby(id_col, sort=False):
        lengths = pd.to_numeric(group["length"], errors="coerce")
        sequences = group[sequence_col].dropna().unique()
        if lengths.isna().any() or lengths.nunique() != 1 or len(sequences) != 1:
            raise ValueError(f"Inconsistent Cecc sequence or length for {eid}")
        length = float(lengths.iloc[0])
        seq = sequences[0]
        if (
            length <= 0
            or not length.is_integer()
            or not isinstance(seq, str)
            or len(seq) != int(length)
        ):
            raise ValueError(f"Cecc FASTA length disagrees with metadata for {eid}")
