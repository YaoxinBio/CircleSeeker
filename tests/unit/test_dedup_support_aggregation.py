"""RCA support semantics that must survive bulk cluster writeback."""

import json

import pandas as pd
import pytest

from circleseeker.modules.ecc_dedup import EccDedup


def test_support_preserves_segment_order_and_read_order():
    source = pd.DataFrame([
        dict(cluster_id="b", eccDNA_id="B1", query_id="x|rep0|100|3", reads="x", copy_number=3),
        dict(cluster_id="b", eccDNA_id="B1", query_id="x|rep0|100|3", reads="x", copy_number=3),
        dict(cluster_id="a", eccDNA_id="A1", query_id="z|rep0|100|2", reads="z", copy_number=2),
        dict(cluster_id="b", eccDNA_id="B2", query_id="y|rep0|100|5", reads="y", copy_number=5),
        dict(cluster_id="a", eccDNA_id="A2", query_id="z|rep1|100|4", reads="z", copy_number=4),
    ])
    result = pd.DataFrame({
        "cluster_id": ["b", "a", "b"],
        "reads": ["y;x", "z", "y;x"],
        "copy_number": pd.array([99, 99, 99], dtype="Int64"),
        "repeat_number": [13, 14, 15],
        "strand": ["-", "+", "+"],
        "segment_in_circle": [2, 1, 1],
        "candidate_support": pd.array(["old", "old", "old"], dtype="string"),
        "per_read_copy_number": pd.array(["old", "old", "old"], dtype="string"),
    }, index=[4, 1, 7])
    fixed = EccDedup()._attach_cluster_support(source, result, "Cecc")
    assert fixed.index.tolist() == [4, 1, 7]
    assert fixed.strand.tolist() == ["-", "+", "+"]
    assert fixed.segment_in_circle.tolist() == [2, 1, 1]
    assert fixed.copy_number.tolist() == [8, 6, 8]
    assert str(fixed.copy_number.dtype) == "Int64"
    assert fixed.repeat_number.tolist() == [13, 14, 15]
    assert fixed.per_read_copy_number.tolist() == ["5;3", "6", "5;3"]
    assert len(json.loads(fixed.loc[4, "candidate_support"])) == 2
    assert set(json.loads(fixed.loc[1, "candidate_support"])) == {"z|rep0|100|2", "z|rep1|100|4"}


def test_unknown_support_retains_representative_values():
    source = pd.DataFrame([
        dict(cluster_id="unknown", eccDNA_id="U1", reads="a;b", copy_number=8),
        dict(cluster_id="known", eccDNA_id="U2", query_id="c|rep0|100|5", reads="c", copy_number=5),
    ])
    result = pd.DataFrame({
        "cluster_id": ["known", "unknown"],
        "reads": ["c", "b;a"],
        "copy_number": pd.array([100, 8], dtype="Int64"),
        "repeat_number": [100, 11],
    })
    fixed = EccDedup()._attach_cluster_support(source, result, "Uecc", update_repeat_number=True)
    assert fixed.copy_number.tolist() == [5, 8]
    assert fixed.repeat_number.tolist() == [5, 11]
    assert fixed.per_read_copy_number.tolist() == ["5", ";"]
    assert all(v["copy_number"] is None for v in json.loads(fixed.loc[1, "candidate_support"]).values())


def test_support_without_representative_read_column_uses_sorted_names():
    source = pd.DataFrame([
        dict(cluster_id="b", eccDNA_id="B1", query_id="z|rep0|100|3", reads="z", copy_number=3),
        dict(cluster_id="b", eccDNA_id="B2", query_id="a|rep0|100|5", reads="a", copy_number=5),
    ])
    result = pd.DataFrame([dict(cluster_id="b", copy_number=3)])
    fixed = EccDedup()._attach_cluster_support(source, result, "Uecc", update_repeat_number=True)
    assert fixed.per_read_copy_number.tolist() == ["5;3"]
    assert fixed.copy_number.tolist() == [8]
    assert fixed.repeat_number.tolist() == [8]


def test_conflicting_candidate_support_is_still_rejected():
    source = pd.DataFrame([
        dict(cluster_id="b", eccDNA_id="B1", query_id="x|rep0|100|3", reads="x", copy_number=3),
        dict(cluster_id="b", eccDNA_id="B1", query_id="x|rep0|100|3", reads="x", copy_number=5),
    ])
    result = pd.DataFrame([dict(cluster_id="b", reads="x", copy_number=3)])
    with pytest.raises(ValueError, match="Conflicting RCA support"):
        EccDedup()._attach_cluster_support(source, result, "Uecc")


def test_empty_support_does_not_invent_a_copy_total():
    source = pd.DataFrame([dict(cluster_id="empty", eccDNA_id="U1", reads="")])
    result = pd.DataFrame([dict(cluster_id="empty", reads="")])
    fixed = EccDedup()._attach_cluster_support(source, result, "Uecc", update_repeat_number=True)
    assert fixed.candidate_support.tolist() == ["{}"]
    assert fixed.per_read_copy_number.tolist() == [""]
    assert "copy_number" not in fixed and "repeat_number" not in fixed
