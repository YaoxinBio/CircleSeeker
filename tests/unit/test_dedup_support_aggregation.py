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


class TestSupportCollectionAvoidsIterrows:
    """collect_support must not build a Series per row.

    support_fields is called once per cluster from ecc_dedup and once per
    location signature from umc_process - about 1.1 million times each on
    GlioSarc_P01_Tumor.  iterrows constructs an object Series for every row,
    which fable measured at 50-100 us before any of the per-row logic runs.
    Every access in this function is `row.get(...)`, which a dict answers
    identically.
    """

    @staticmethod
    def _mixed_frame():
        return pd.DataFrame(
            [
                dict(reads="a;b", per_read_copy_number="2;3", query_id="q1|circular",
                     copy_number=5, eccDNA_id="U1"),
                dict(reads="c", per_read_copy_number=None, query_id=None,
                     copy_number=float("nan"), eccDNA_id="U2", repeat_number=7),
                dict(reads="d;d", per_read_copy_number=None, query_id="q3",
                     copy_number=1.5, eccDNA_id="U3"),
            ]
        )

    def test_collect_support_does_not_call_iterrows(self, monkeypatch):
        from circleseeker.utils.read_support import collect_support

        def boom(self):
            raise AssertionError("collect_support still walks rows with iterrows")

        monkeypatch.setattr(pd.DataFrame, "iterrows", boom)
        support = collect_support(self._mixed_frame())
        assert support  # and it still produced something

    def test_values_survive_mixed_dtypes(self):
        from circleseeker.utils.read_support import collect_support

        support = collect_support(self._mixed_frame())

        assert support["q1::a"] == {"read": "a", "copy_number": 2.0}
        assert support["q1::b"] == {"read": "b", "copy_number": 3.0}
        # single read, no query_id: the key is built from eccDNA_id and read,
        # and the copy number falls back to repeat_number
        assert support["U2::c"] == {"read": "c", "copy_number": 7.0}
        # duplicate read names collapse to one entry
        assert support["q3"] == {"read": "d", "copy_number": 1.5}

    def test_encoded_candidate_support_still_short_circuits(self):
        from circleseeker.utils.read_support import collect_support

        df = pd.DataFrame(
            [dict(candidate_support='{"k1": {"read": "r1", "copy_number": 4}}', reads="ignored")]
        )
        assert collect_support(df) == {"k1": {"read": "r1", "copy_number": 4.0}}
