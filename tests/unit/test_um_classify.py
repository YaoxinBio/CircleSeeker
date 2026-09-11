from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.um_classify import UMeccClassifier

def test_um_classify_subset(tmp_path):
    blast_tsv = tmp_path / "step3_subset.tsv"

    q_uecc = "readU|rep1|100|1"
    q_mecc = "readM|rep1|100|2"
    q_unclassified = "readX|rep1|100|1"

    rows = [
        # Single high-quality alignment -> Uecc
        [q_uecc, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        # Two non-overlapping, full-length alignments -> Mecc
        [q_mecc, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        [q_mecc, "chr1", 98.5, 100, 0, 0, 1, 100, 400, 499, 0.0, 190.0, "plus"],
        # High-quality but not full-length coverage -> unclassified
        [q_unclassified, "chr2", 99.0, 90, 0, 0, 1, 90, 100, 189, 0.0, 180.0, "plus"],
        [q_unclassified, "chr2", 98.0, 90, 0, 0, 1, 90, 300, 389, 0.0, 170.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier()
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert set(uecc_df["query_id"].unique()) == {q_uecc}
    assert set(mecc_df["query_id"].unique()) == {q_mecc}
    assert set(unclassified_df["query_id"].unique()) == {q_unclassified}

    assert set(uecc_df["eccdna_type"].unique()) == {"Uecc"}
    assert set(mecc_df["eccdna_type"].unique()) == {"Mecc"}

    assert len(unclassified_df) == 2
    assert {"unclass_reason", "quality_category"}.issubset(unclassified_df.columns)
    assert set(unclassified_df["quality_category"].unique()) == {"High_quality"}


def test_um_classify_invalid_column_count(tmp_path):
    classifier = UMeccClassifier()
    alignment_df = pd.DataFrame([["query1", "subject1"]])

    with pytest.raises(ValueError):
        classifier.classify_alignment_results(alignment_df, tmp_path / "invalid")


def test_um_classify_overlap_components_mecc(tmp_path):
    classifier = UMeccClassifier()
    query_id = "readA|rep1|100|2"

    alignment_df = pd.DataFrame(
        [
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 100, 199, 0, 0, "plus"],
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 150, 249, 0, 0, "plus"],
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 400, 499, 0, 0, "plus"],
            [query_id, "chr1", 99.0, 100, 1, 0, 1, 100, 450, 549, 0, 0, "plus"],
        ]
    )

    uecc_df, mecc_df, unclassified_df = classifier.classify_alignment_results(
        alignment_df, tmp_path / "um_classify"
    )

    assert uecc_df.empty
    assert mecc_df["query_id"].nunique() == 1
    assert len(mecc_df) == 2
    assert unclassified_df.empty


def test_um_classify_uecc_requires_low_second_locus(tmp_path):
    """If a strong 2nd locus exists, the query should not be labeled as Uecc."""
    blast_tsv = tmp_path / "step3_second_locus.tsv"

    query_id = "readU2|rep1|100|1"

    rows = [
        # Best locus: full ring coverage
        [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        # Second locus: substantial (90%) ring coverage, but below theta_m=0.95
        [query_id, "chr1", 99.0, 90, 0, 0, 1, 90, 500, 589, 0.0, 190.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert mecc_df.empty
    assert set(unclassified_df["query_id"].unique()) == {query_id}


def test_um_classify_mecc_can_use_lower_theta_m(tmp_path):
    """Mecc can be recovered by lowering theta_m while keeping U strict."""
    blast_tsv = tmp_path / "step3_theta_m.tsv"

    query_id = "readM2|rep1|100|1"

    rows = [
        # Locus A: full ring coverage
        [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus"],
        # Locus B: 90% ring coverage (should count toward Mecc when theta_m<=0.9)
        [query_id, "chr1", 99.0, 90, 0, 0, 1, 90, 500, 589, 0.0, 190.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.85, theta_u2_max=0.05)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert set(mecc_df["query_id"].unique()) == {query_id}
    assert mecc_df["query_id"].nunique() == 1
    assert len(mecc_df) == 2
    assert unclassified_df.empty


def test_um_classify_uecc_vetoes_significant_secondary_mapping(tmp_path):
    """A strong single-locus coverage can still be vetoed by significant secondary evidence.

    Note: After the adaptive threshold improvements, small secondary mappings (<5% of primary)
    are allowed. This test uses a 10% secondary mapping which exceeds the relative threshold.
    """
    blast_tsv = tmp_path / "step3_secondary_evidence.tsv"

    query_id = "readC|rep1|10000|1"

    rows = [
        # Primary locus: ~90% ring coverage
        [query_id, "chr1", 99.0, 9000, 0, 0, 1, 9000, 1000, 9999, 0.0, 200.0, "plus"],
        # Secondary locus: 10% ring coverage on another chromosome (exceeds relative threshold)
        [query_id, "chr5", 99.0, 1000, 0, 0, 9001, 10000, 5000, 5999, 0.0, 50.0, "plus"],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    # 10% secondary / 90% primary = 11.1% > 5% threshold => vetoed
    assert uecc_df.empty
    assert mecc_df.empty
    assert set(unclassified_df["query_id"].unique()) == {query_id}


def test_um_classify_uecc_can_be_gated_by_mapq(tmp_path):
    """When enabled, mapq_u_min should veto Uecc calls in low-uniqueness mappings."""
    blast_tsv = tmp_path / "step3_mapq_gate.tsv"

    query_id = "readU_mapq|rep1|100|1"

    rows = [
        # Full ring coverage, but mapq is 0 (ambiguous / repetitive mapping)
        [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus", 0],
    ]

    pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)

    classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05, mapq_u_min=20)
    uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)

    assert uecc_df.empty
    assert mecc_df.empty
    assert set(unclassified_df["query_id"].unique()) == {query_id}


def test_um_classify_confidence_score_is_monotonic_in_mapq(tmp_path):
    """For the same read, confidence_score should increase as MAPQ increases."""
    query_id = "readU_conf|rep1|100|1"

    def run_with_mapq(mapq: int) -> float:
        blast_tsv = tmp_path / f"step3_confidence_mapq_{mapq}.tsv"
        rows = [
            [query_id, "chr1", 99.0, 100, 0, 0, 1, 100, 100, 199, 0.0, 200.0, "plus", mapq],
        ]
        pd.DataFrame(rows).to_csv(blast_tsv, sep="\t", header=False, index=False)
        classifier = UMeccClassifier(theta_u=0.95, theta_m=0.95, theta_u2_max=0.05, mapq_u_min=0)
        uecc_df, mecc_df, unclassified_df = classifier.run(blast_tsv)
        assert set(uecc_df["query_id"].unique()) == {query_id}
        assert mecc_df.empty
        assert unclassified_df.empty
        assert "confidence_score" in uecc_df.columns
        return float(uecc_df["confidence_score"].iloc[0])

    score_0 = run_with_mapq(0)
    score_20 = run_with_mapq(20)
    score_60 = run_with_mapq(60)
    assert score_0 <= score_20 <= score_60


# ---------------------------------------------------------------------------
# New unit-test classes for individual static/class/instance methods
# ---------------------------------------------------------------------------

import numpy as np
from circleseeker.utils.column_standards import ColumnStandard


# ---- helpers ---------------------------------------------------------------

def _make_preprocessed_df(rows, columns=None):
    """Build a DataFrame that looks like the output of _preprocess_alignment_df."""
    if columns is None:
        columns = UMeccClassifier.ALIGNMENT_COLUMNS
    return pd.DataFrame(rows, columns=columns)


def _make_alignment_row(
    query_id="read1|rep1|100|1",
    subject_id="chr1",
    identity=99.0,
    alignment_length=100,
    mismatches=0,
    gap_opens=0,
    q_start=1,
    q_end=100,
    s_start=100,
    s_end=199,
    evalue=0.0,
    bit_score=200.0,
    sstrand="plus",
):
    return [
        query_id, subject_id, identity, alignment_length,
        mismatches, gap_opens, q_start, q_end,
        s_start, s_end, evalue, bit_score, sstrand,
    ]


# ---- 1. _as_fraction -------------------------------------------------------

class TestAsFraction:
    def test_fraction_unchanged(self):
        assert UMeccClassifier._as_fraction(0.5) == 0.5

    def test_percent_converted(self):
        assert UMeccClassifier._as_fraction(95) == 0.95

    def test_none_returns_zero(self):
        assert UMeccClassifier._as_fraction(None) == 0.0

    def test_non_numeric_returns_zero(self):
        assert UMeccClassifier._as_fraction("abc") == 0.0

    def test_zero(self):
        assert UMeccClassifier._as_fraction(0) == 0.0

    def test_one(self):
        assert UMeccClassifier._as_fraction(1.0) == 1.0

    def test_boundary_just_above_one(self):
        assert UMeccClassifier._as_fraction(1.01) == pytest.approx(0.0101)


# ---- 2. _clamp01 -----------------------------------------------------------

class TestClamp01:
    def test_in_range(self):
        assert UMeccClassifier._clamp01(0.5) == 0.5

    def test_negative(self):
        assert UMeccClassifier._clamp01(-0.5) == 0.0

    def test_above_one(self):
        assert UMeccClassifier._clamp01(1.5) == 1.0

    def test_boundary_zero(self):
        assert UMeccClassifier._clamp01(0.0) == 0.0

    def test_boundary_one(self):
        assert UMeccClassifier._clamp01(1.0) == 1.0

    def test_non_numeric(self):
        assert UMeccClassifier._clamp01("xyz") == 0.0


# ---- 3. _norm_mapq ---------------------------------------------------------

class TestNormMapq:
    def test_max_mapq(self):
        assert UMeccClassifier._norm_mapq(60) == pytest.approx(1.0)

    def test_zero_mapq(self):
        assert UMeccClassifier._norm_mapq(0) == pytest.approx(0.0)

    def test_half_mapq(self):
        assert UMeccClassifier._norm_mapq(30) == pytest.approx(0.5)

    def test_non_numeric(self):
        assert UMeccClassifier._norm_mapq("bad") == 0.0


# ---- 4. _norm_identity -----------------------------------------------------

class TestNormIdentity:
    def test_100_returns_one(self):
        assert UMeccClassifier._norm_identity(100) == pytest.approx(1.0)

    def test_90_returns_zero(self):
        assert UMeccClassifier._norm_identity(90) == pytest.approx(0.0)

    def test_95_returns_half(self):
        assert UMeccClassifier._norm_identity(95) == pytest.approx(0.5)

    def test_below_90_clamped(self):
        assert UMeccClassifier._norm_identity(85) == 0.0

    def test_non_numeric(self):
        assert UMeccClassifier._norm_identity("bad") == 0.0


# ---- 5. _geom_mean ---------------------------------------------------------

class TestGeomMean:
    def test_single_value(self):
        assert UMeccClassifier._geom_mean([0.8]) == pytest.approx(0.8)

    def test_two_values(self):
        expected = (0.64 * 0.81) ** 0.5
        assert UMeccClassifier._geom_mean([0.64, 0.81]) == pytest.approx(expected)

    def test_with_zero(self):
        assert UMeccClassifier._geom_mean([0.5, 0.0]) == 0.0

    def test_empty(self):
        assert UMeccClassifier._geom_mean([]) == 0.0


# ---- 6. _infer_query_coords_style ------------------------------------------

class TestInferQueryCoordsStyle:
    def test_0based_detected(self):
        # alignment_length == q_end - q_start  -->  0-based
        df = pd.DataFrame({
            "q_start": [0, 10],
            "q_end": [100, 60],
            "alignment_length": [100, 50],
        })
        assert UMeccClassifier._infer_query_coords_style(df) == "0based"

    def test_blast_detected(self):
        # alignment_length == q_end - q_start + 1  -->  blast (1-based)
        df = pd.DataFrame({
            "q_start": [1, 11],
            "q_end": [100, 60],
            "alignment_length": [100, 50],
        })
        assert UMeccClassifier._infer_query_coords_style(df) == "blast"

    def test_empty_df_returns_blast(self):
        df = pd.DataFrame(columns=["q_start", "q_end", "alignment_length"])
        assert UMeccClassifier._infer_query_coords_style(df) == "blast"


# ---- 7. _query_interval_0based_half_open -----------------------------------

class TestQueryInterval0BasedHalfOpen:
    def test_0based_no_adjustment(self):
        assert UMeccClassifier._query_interval_0based_half_open(0, 100, "0based") == (0, 100)

    def test_blast_subtracts_one_from_start(self):
        # blast style only subtracts 1 from q_start; q_end is unchanged
        assert UMeccClassifier._query_interval_0based_half_open(1, 100, "blast") == (0, 100)

    def test_reversed_coordinates(self):
        # When q_end < q_start after conversion, should be swapped
        result = UMeccClassifier._query_interval_0based_half_open(100, 1, "0based")
        assert result == (1, 100)


# ---- 8. _project_query_interval_to_ring ------------------------------------

class TestProjectQueryIntervalToRing:
    def test_non_wrapping(self):
        # q_start=1, q_end=50, cons_len=100, blast -> 0-based [0,49) on ring of 100
        segs = UMeccClassifier._project_query_interval_to_ring(0, 50, 100, "0based")
        assert segs == [(0, 50)]

    def test_wrapping(self):
        # Interval [80, 120) on ring of 100 -> wraps: [80,100) + [0,20)
        segs = UMeccClassifier._project_query_interval_to_ring(80, 120, 100, "0based")
        assert sorted(segs) == [(0, 20), (80, 100)]

    def test_full_coverage(self):
        segs = UMeccClassifier._project_query_interval_to_ring(0, 100, 100, "0based")
        assert segs == [(0, 100)]

    def test_zero_cons_len(self):
        assert UMeccClassifier._project_query_interval_to_ring(0, 50, 0, "0based") == []

    def test_empty_interval(self):
        assert UMeccClassifier._project_query_interval_to_ring(50, 50, 100, "0based") == []


# ---- 9. _union_len ---------------------------------------------------------

class TestUnionLen:
    def test_non_overlapping(self):
        assert UMeccClassifier._union_len([(0, 10), (20, 30)]) == 20

    def test_overlapping(self):
        assert UMeccClassifier._union_len([(0, 15), (10, 30)]) == 30

    def test_adjacent(self):
        assert UMeccClassifier._union_len([(0, 10), (10, 20)]) == 20

    def test_empty(self):
        assert UMeccClassifier._union_len([]) == 0

    def test_single(self):
        assert UMeccClassifier._union_len([(5, 15)]) == 10

    def test_nested(self):
        assert UMeccClassifier._union_len([(0, 100), (10, 50)]) == 100


# ---- 10. _coverage_fraction_for_alignments ---------------------------------

class TestCoverageFractionForAlignments:
    def test_basic(self):
        df = pd.DataFrame({"q_start": [0, 50], "q_end": [50, 100]})
        result = UMeccClassifier._coverage_fraction_for_alignments(df, 100, "0based")
        assert result == pytest.approx(1.0)

    def test_empty_df(self):
        df = pd.DataFrame(columns=["q_start", "q_end"])
        assert UMeccClassifier._coverage_fraction_for_alignments(df, 100, "0based") == 0.0

    def test_zero_cons_len(self):
        df = pd.DataFrame({"q_start": [0], "q_end": [50]})
        assert UMeccClassifier._coverage_fraction_for_alignments(df, 0, "0based") == 0.0


# ---- 11. _ensure_alignment_columns -----------------------------------------

class TestEnsureAlignmentColumns:
    def test_correct_13_columns(self):
        clf = UMeccClassifier()
        data = [_make_alignment_row()]
        df = pd.DataFrame(data)
        result = clf._ensure_alignment_columns(df, "test")
        assert list(result.columns) == clf.ALIGNMENT_COLUMNS

    def test_correct_14_columns(self):
        clf = UMeccClassifier()
        row = _make_alignment_row() + [60]
        df = pd.DataFrame([row])
        result = clf._ensure_alignment_columns(df, "test")
        assert list(result.columns) == clf.ALIGNMENT_COLUMNS_WITH_MAPQ

    def test_wrong_count_raises(self):
        clf = UMeccClassifier()
        df = pd.DataFrame([["a", "b", "c"]])
        with pytest.raises(ValueError):
            clf._ensure_alignment_columns(df, "test")


# ---- 12. _preprocess_alignment_df ------------------------------------------

class TestPreprocessAlignmentDf:
    def test_basic_parse(self):
        clf = UMeccClassifier()
        row = _make_alignment_row()
        df = pd.DataFrame([row])
        result = clf._preprocess_alignment_df(df, "test")
        assert ColumnStandard.CHR in result.columns
        assert ColumnStandard.START0 in result.columns
        assert ColumnStandard.END0 in result.columns
        assert ColumnStandard.READS in result.columns
        assert ColumnStandard.LENGTH in result.columns
        assert len(result) == 1

    def test_strand_conversion_minus(self):
        clf = UMeccClassifier()
        row = _make_alignment_row(sstrand="minus")
        df = pd.DataFrame([row])
        result = clf._preprocess_alignment_df(df, "test")
        assert result["strand"].iloc[0] == "-"

    def test_query_id_parsing(self):
        clf = UMeccClassifier()
        row = _make_alignment_row(query_id="myread|consensus|200|3")
        df = pd.DataFrame([row])
        result = clf._preprocess_alignment_df(df, "test")
        assert result[ColumnStandard.READS].iloc[0] == "myread"
        assert result[ColumnStandard.LENGTH].iloc[0] == 200
        assert result[ColumnStandard.COPY_NUMBER].iloc[0] == 3.0


# ---- 13. _cluster_loci -----------------------------------------------------

class TestClusterLoci:
    def _make_preprocessed_group(self, clf, rows_raw):
        df = pd.DataFrame(rows_raw)
        df = clf._preprocess_alignment_df(df, "test")
        return df

    def test_single_locus(self):
        clf = UMeccClassifier()
        row = _make_alignment_row()
        df = self._make_preprocessed_group(clf, [row])
        loci = clf._cluster_loci(df)
        assert len(loci) == 1

    def test_two_distant_loci(self):
        clf = UMeccClassifier()
        rows = [
            _make_alignment_row(s_start=100, s_end=199),
            _make_alignment_row(s_start=50000, s_end=50099),
        ]
        df = self._make_preprocessed_group(clf, rows)
        loci = clf._cluster_loci(df)
        assert len(loci) == 2

    def test_rows_with_a_missing_chromosome_never_cluster_together(self):
        """`groupby([chr, strand])` drops NaN keys - that is its default - so a
        malformed row with an empty subject_id was never clustered with
        anything. Bucketing by hand would instead collect every such row into
        one bucket keyed (nan, strand) and union them into a single locus.
        """
        clf = UMeccClassifier()
        rows = [
            _make_alignment_row(s_start=100, s_end=199),
            _make_alignment_row(s_start=110, s_end=209),
            _make_alignment_row(s_start=120, s_end=219),
        ]
        df = self._make_preprocessed_group(clf, rows)
        chr_col = df.columns[list(df.columns).index(ColumnStandard.CHR)]
        df.loc[df.index[:2], chr_col] = np.nan

        loci = clf._cluster_loci(df)

        # two NaN rows that overlap each other, plus one real row
        assert len(loci) == 3
        assert all(len(members) == 1 for members in loci.values())

    def test_overlapping_same_cluster(self):
        clf = UMeccClassifier()
        rows = [
            _make_alignment_row(s_start=100, s_end=199),
            _make_alignment_row(s_start=110, s_end=209),
        ]
        df = self._make_preprocessed_group(clf, rows)
        loci = clf._cluster_loci(df)
        assert len(loci) == 1


# ---- 14. _is_full_length_repeat --------------------------------------------

class TestIsFullLengthRepeat:
    def _make_group(self, rlengths, cons_len=100):
        rows = []
        for rl in rlengths:
            rows.append({
                ColumnStandard.LENGTH: cons_len,
                "Rlength": rl,
            })
        return pd.DataFrame(rows)

    def test_two_full_length(self):
        clf = UMeccClassifier(min_full_length_coverage=95.0)
        group = self._make_group([98, 97], cons_len=100)
        assert clf._is_full_length_repeat(group) is True

    def test_one_short(self):
        clf = UMeccClassifier(min_full_length_coverage=95.0)
        group = self._make_group([98, 50], cons_len=100)
        assert clf._is_full_length_repeat(group) is False

    def test_empty(self):
        clf = UMeccClassifier()
        group = pd.DataFrame(columns=[ColumnStandard.LENGTH, "Rlength"])
        assert clf._is_full_length_repeat(group) is False


# ---- 15. _passes_identity_gap ----------------------------------------------

class TestPassesIdentityGap:
    def test_within_gap(self):
        clf = UMeccClassifier(max_identity_gap_for_mecc=5.0)
        group = pd.DataFrame({"identity": [99.0, 97.0, 95.0]})
        assert clf._passes_identity_gap(group) is True

    def test_exceeds_gap(self):
        clf = UMeccClassifier(max_identity_gap_for_mecc=1.0)
        group = pd.DataFrame({"identity": [99.0, 90.0]})
        assert clf._passes_identity_gap(group) is False

    def test_none_max_identity_gap(self):
        clf = UMeccClassifier(max_identity_gap_for_mecc=None)
        group = pd.DataFrame({"identity": [99.0, 50.0]})
        assert clf._passes_identity_gap(group) is True

    def test_empty_group(self):
        clf = UMeccClassifier(max_identity_gap_for_mecc=5.0)
        group = pd.DataFrame(columns=["identity"])
        assert clf._passes_identity_gap(group) is False


# ---- 16. format_outputs ----------------------------------------------------

class TestFormatOutputs:
    def test_adds_match_degree(self):
        clf = UMeccClassifier()
        uecc_df = pd.DataFrame({
            "query_id": ["q1"],
            ColumnStandard.GAP_PERCENTAGE: [5.0],
        })
        mecc_df = pd.DataFrame({
            "query_id": ["q2"],
            ColumnStandard.GAP_PERCENTAGE: [3.0],
        })
        u_out, m_out = clf.format_outputs(uecc_df, mecc_df)
        assert ColumnStandard.MATCH_DEGREE in u_out.columns
        assert u_out[ColumnStandard.MATCH_DEGREE].iloc[0] == pytest.approx(95.0)
        assert ColumnStandard.MATCH_DEGREE in m_out.columns
        assert m_out[ColumnStandard.MATCH_DEGREE].iloc[0] == pytest.approx(97.0)

    def test_handles_empty_uecc(self):
        clf = UMeccClassifier()
        uecc_df = pd.DataFrame()
        mecc_df = pd.DataFrame({
            "query_id": ["q2"],
            ColumnStandard.GAP_PERCENTAGE: [3.0],
        })
        u_out, m_out = clf.format_outputs(uecc_df, mecc_df)
        assert u_out.empty
        assert ColumnStandard.MATCH_DEGREE in m_out.columns

    def test_handles_empty_mecc(self):
        clf = UMeccClassifier()
        uecc_df = pd.DataFrame({
            "query_id": ["q1"],
            ColumnStandard.GAP_PERCENTAGE: [5.0],
        })
        mecc_df = pd.DataFrame()
        u_out, m_out = clf.format_outputs(uecc_df, mecc_df)
        assert m_out.empty
        assert ColumnStandard.MATCH_DEGREE in u_out.columns

    def test_handles_both_empty(self):
        clf = UMeccClassifier()
        u_out, m_out = clf.format_outputs(pd.DataFrame(), pd.DataFrame())
        assert u_out.empty
        assert m_out.empty


# ---- 17. _u_has_significant_secondary_mapping ------------------------------

class TestUHasSignificantSecondaryMapping:
    def _build_df(self, rows):
        """Build a small DataFrame with the columns the method needs."""
        return pd.DataFrame(rows, columns=[
            ColumnStandard.CHR, ColumnStandard.START0, ColumnStandard.END0,
            "q_start", "q_end",
        ])

    def test_no_secondary_returns_false(self):
        clf = UMeccClassifier()
        # All alignments on same chr, within gap_bp of best locus
        df = self._build_df([
            ["chr1", 100, 200, 0, 100],
        ])
        result = clf._u_has_significant_secondary_mapping(
            df, best_chr="chr1", best_start0=100, best_end0=200,
            cons_len=100, q_style="0based",
        )
        assert result is False

    def test_small_secondary_returns_false(self):
        clf = UMeccClassifier()
        # Primary on chr1 and a tiny secondary on chr2.
        # The secondary covers 2/1000 = 0.2% of the ring, well below the
        # adaptive relative ratio (secondary_ratio = 0.002/1.0 = 0.2% < 5%).
        df = self._build_df([
            ["chr1", 100, 1100, 0, 1000],
            ["chr2", 500, 502, 998, 1000],
        ])
        result = clf._u_has_significant_secondary_mapping(
            df, best_chr="chr1", best_start0=100, best_end0=1100,
            cons_len=1000, q_style="0based", u_cov=1.0,
        )
        assert result is False

    def test_large_secondary_returns_true(self):
        clf = UMeccClassifier()
        # Primary on chr1 (100bp); secondary on chr2 covers 50% of ring
        df = self._build_df([
            ["chr1", 100, 200, 0, 100],
            ["chr2", 500, 600, 0, 50],
        ])
        result = clf._u_has_significant_secondary_mapping(
            df, best_chr="chr1", best_start0=100, best_end0=200,
            cons_len=100, q_style="0based", u_cov=1.0,
        )
        assert result is True

    def test_empty_returns_false(self):
        clf = UMeccClassifier()
        df = pd.DataFrame(columns=[
            ColumnStandard.CHR, ColumnStandard.START0, ColumnStandard.END0,
            "q_start", "q_end",
        ])
        result = clf._u_has_significant_secondary_mapping(
            df, best_chr="chr1", best_start0=0, best_end0=100,
            cons_len=100, q_style="0based",
        )
        assert result is False


# ---- 18. extract_unclassified -----------------------------------------------

class TestExtractUnclassified:
    def _make_df(self):
        clf = UMeccClassifier()
        rows = [
            _make_alignment_row(query_id="readA|rep1|100|1"),
            _make_alignment_row(query_id="readB|rep1|100|1"),
            _make_alignment_row(query_id="readC|rep1|100|1"),
        ]
        df = pd.DataFrame(rows)
        df = clf._preprocess_alignment_df(df, "test")
        return clf, df

    def test_basic_extraction(self):
        clf, df = self._make_df()
        classified = {"readA|rep1|100|1"}
        result = clf.extract_unclassified(df, classified)
        result_queries = set(result["query_id"].unique())
        assert "readA|rep1|100|1" not in result_queries
        assert "readB|rep1|100|1" in result_queries
        assert "readC|rep1|100|1" in result_queries

    def test_all_classified_returns_empty(self):
        clf, df = self._make_df()
        all_qids = set(df["query_id"].unique())
        result = clf.extract_unclassified(df, all_qids)
        assert result.empty

    def test_quality_category_column_added(self):
        clf, df = self._make_df()
        result = clf.extract_unclassified(df, set())
        assert "quality_category" in result.columns
        assert "unclass_reason" in result.columns


class TestUnclassifiedBreakdownIsVectorised:
    """The HQ/LQ split is a DEBUG line; it must not cost a groupby loop.

    Both report sites iterated `unclassified_df.groupby("query_id")` and ran
    `(group["quality_category"] == "High_quality").any()` per group, purely to
    fill in one debug message.  The group count is unclassified queries, which
    grows with the sample - fable measured 32 us per group, i.e. 32-64 s at
    1-2 million groups.
    """

    @staticmethod
    def _frame():
        return pd.DataFrame(
            {
                "query_id": ["q1", "q1", "q2", "q3", "q3", "q4"],
                "quality_category": [
                    "Low_quality",
                    "High_quality",   # q1 has one HQ -> counts as HQ
                    "Low_quality",    # q2 is LQ only
                    "Low_quality",
                    "Low_quality",    # q3 is LQ only
                    "High_quality",   # q4 is HQ
                ],
            }
        )

    @staticmethod
    def _reference(frame):
        hq = lq = 0
        for _, group in frame.groupby("query_id"):
            if (group["quality_category"] == "High_quality").any():
                hq += 1
            else:
                lq += 1
        return hq, lq

    def test_counts_match_the_per_group_loop(self):
        from circleseeker.modules.um_classify import UMeccClassifier

        frame = self._frame()
        assert UMeccClassifier._unclassified_quality_breakdown(frame) == self._reference(frame)

    def test_empty_frame_reports_zero(self):
        from circleseeker.modules.um_classify import UMeccClassifier

        empty = pd.DataFrame({"query_id": [], "quality_category": []})
        assert UMeccClassifier._unclassified_quality_breakdown(empty) == (0, 0)

    def test_missing_column_is_treated_as_low_quality(self):
        from circleseeker.modules.um_classify import UMeccClassifier

        frame = pd.DataFrame({"query_id": ["a", "b", "b"]})
        assert UMeccClassifier._unclassified_quality_breakdown(frame) == (0, 2)

    def test_breakdown_does_not_iterate_groups(self, monkeypatch):
        from circleseeker.modules.um_classify import UMeccClassifier

        frame = self._frame()
        original = pd.core.groupby.generic.DataFrameGroupBy.__iter__

        def boom(self):
            raise AssertionError("breakdown still iterates groups")

        monkeypatch.setattr(pd.core.groupby.generic.DataFrameGroupBy, "__iter__", boom)
        try:
            assert UMeccClassifier._unclassified_quality_breakdown(frame) == (2, 2)
        finally:
            monkeypatch.setattr(pd.core.groupby.generic.DataFrameGroupBy, "__iter__", original)


class TestPerQueryAlignmentLookup:
    """Looking up one query's alignments must not build an Index per query.

    `df.groupby("query_id", sort=False).groups` materialises a pandas Index
    object for every group.  At ~2 million distinct query_ids that is millions
    of objects and hundreds of MB resident; fable measured 11.2 s vs 2.7 s for
    `.indices` on a 4M-row / 2M-group frame.
    """

    @staticmethod
    def _frame():
        return pd.DataFrame(
            {
                "query_id": ["q1", "q2", "q1", "q3", "q2"],
                "mapq": [60, 10, 55, 30, 20],
            },
            index=[10, 11, 12, 13, 14],
        )

    def test_lookup_matches_a_boolean_selection(self):
        from circleseeker.modules.um_classify import UMeccClassifier

        frame = self._frame()
        positions = UMeccClassifier._positions_by_query(frame)

        for query in ("q1", "q2", "q3"):
            expected = frame[frame["query_id"] == query]
            got = UMeccClassifier._alignments_for_query(frame, positions, query, None)
            pd.testing.assert_frame_equal(got, expected)

    def test_unknown_query_returns_the_fallback(self):
        from circleseeker.modules.um_classify import UMeccClassifier

        frame = self._frame()
        positions = UMeccClassifier._positions_by_query(frame)
        sentinel = frame.head(1)

        got = UMeccClassifier._alignments_for_query(frame, positions, "absent", sentinel)
        assert got is sentinel

    def test_lookup_does_not_use_the_groups_mapping(self, monkeypatch):
        from circleseeker.modules.um_classify import UMeccClassifier

        def boom(self):
            raise AssertionError("built .groups instead of .indices")

        monkeypatch.setattr(
            pd.core.groupby.generic.DataFrameGroupBy, "groups", property(boom)
        )
        positions = UMeccClassifier._positions_by_query(self._frame())
        assert set(positions) == {"q1", "q2", "q3"}


class TestPerLocusWorkAvoidsFrameSlices:
    """Per-locus work must not slice a whole DataFrame.

    classify_uecc_mecc did `group.loc[idxs]` once per locus for coverage, and
    again inside locus_evidence. Each slice materialises all 11+ columns while
    the callers read only q_start/q_end, or mapq/identity. cProfile on a real
    2% subset attributed 3,961,886 calls to pandas' _take_nd_ndarray and as
    many to maybe_promote, with frame.__getitem__ at 22% cumulative.
    """

    @staticmethod
    def _group():
        return pd.DataFrame(
            {
                "q_start": [0, 100, 250, 400],
                "q_end": [120, 260, 410, 500],
                "mapq": [60, 55, 30, 60],
                "identity": [99.5, 98.0, 95.0, 99.9],
            },
            index=[11, 12, 13, 14],
        )

    def test_coverage_from_arrays_matches_the_frame_version(self):
        clf = UMeccClassifier()
        group = self._group()
        positions = [0, 2, 3]
        idxs = [group.index[p] for p in positions]

        from_frame = clf._coverage_fraction_for_alignments(group.loc[idxs], 500, "0-based")
        from_arrays = clf._coverage_fraction_from_arrays(
            group["q_start"].to_numpy()[positions],
            group["q_end"].to_numpy()[positions],
            500,
            "0-based",
        )
        assert from_arrays == from_frame

    def test_coverage_handles_empty_and_nan(self):
        import numpy as np

        clf = UMeccClassifier()
        assert clf._coverage_fraction_from_arrays(np.array([]), np.array([]), 500, "0-based") == 0.0
        starts = np.array([0.0, float("nan"), 200.0])
        ends = np.array([100.0, 150.0, float("nan")])
        frame = pd.DataFrame({"q_start": starts, "q_end": ends})
        assert clf._coverage_fraction_from_arrays(starts, ends, 500, "0-based") == (
            clf._coverage_fraction_for_alignments(frame, 500, "0-based")
        )

    def test_locus_evidence_from_arrays_matches_the_frame_version(self):
        clf = UMeccClassifier()
        group = self._group()
        positions = [1, 2]
        idxs = [group.index[p] for p in positions]
        sliced = group.loc[idxs]

        expected = (
            int(sliced["mapq"].max()), int(sliced["mapq"].min()),
            float(sliced["identity"].max()), float(sliced["identity"].min()),
        )
        got = clf._locus_evidence_from_arrays(
            group["mapq"].to_numpy()[positions], group["identity"].to_numpy()[positions]
        )
        assert got == expected

    def test_locus_evidence_on_empty_positions(self):
        import numpy as np

        clf = UMeccClassifier()
        assert clf._locus_evidence_from_arrays(np.array([]), np.array([])) == (0, 0, 0.0, 0.0)


class TestLociClusteringFromArrays:
    """Clustering must not need a DataFrame, nor an inner groupby.

    _cluster_loci took the group frame, ran `groupby([chr, strand])` inside it
    - one more pandas groupby per query, millions of them - and keyed union-find
    on index labels. The same work over plain arrays is the last piece needed to
    drop DataFrame slicing from the per-query loop; measured separately, just
    iterating groupby and touching three columns costs 7.09s on a 220k-row / 78k
    group subset, against 0.29s for sorted boundary slices (24x).

    Two orderings are load-bearing and must survive: the inner grouping is
    first-appearance order (union order decides which member becomes the root),
    and loci are renumbered by `str()` of the root's index LABEL, which in
    coverage ties decides which locus wins.
    """

    @staticmethod
    def _classifier():
        return UMeccClassifier(pos_tol_bp=50, theta_locus=0.95)

    @staticmethod
    def _group(index):
        return pd.DataFrame(
            {
                ColumnStandard.CHR: ["chr1", "chr1", "chr2", "chr1", "chr2"],
                ColumnStandard.STRAND: ["+", "+", "-", "-", "-"],
                ColumnStandard.START0: [1000, 1020, 5000, 8000, 5010],
                ColumnStandard.END0: [2000, 2010, 6000, 9000, 6020],
            },
            index=index,
        )

    def _compare(self, index):
        clf = self._classifier()
        group = self._group(index)
        expected = clf._cluster_loci(group)

        labels = list(group.index)
        got = clf._cluster_loci_from_arrays(
            labels,
            group[ColumnStandard.CHR].to_numpy(),
            group[ColumnStandard.STRAND].to_numpy(),
            group[ColumnStandard.START0].to_numpy(),
            group[ColumnStandard.END0].to_numpy(),
        )
        # array version returns positions; map back to labels to compare
        as_labels = {lid: [labels[p] for p in pos] for lid, pos in got.items()}
        assert as_labels == expected, (as_labels, expected)

    def test_matches_the_frame_version_with_plain_indices(self):
        self._compare([0, 1, 2, 3, 4])

    def test_matches_with_labels_whose_string_order_differs(self):
        # str() ordering puts "10" before "9"; positional ordering would not
        self._compare([9, 10, 100, 11, 2])

    def test_matches_with_non_contiguous_labels(self):
        self._compare([57, 3, 900, 12, 41])

    def test_empty_group(self):
        clf = self._classifier()
        import numpy as np

        assert clf._cluster_loci_from_arrays([], np.array([]), np.array([]),
                                             np.array([]), np.array([])) == {}

    def test_single_alignment(self):
        clf = self._classifier()
        group = self._group([0, 1, 2, 3, 4]).iloc[:1]
        labels = list(group.index)
        got = clf._cluster_loci_from_arrays(
            labels,
            group[ColumnStandard.CHR].to_numpy(),
            group[ColumnStandard.STRAND].to_numpy(),
            group[ColumnStandard.START0].to_numpy(),
            group[ColumnStandard.END0].to_numpy(),
        )
        assert {lid: [labels[p] for p in pos] for lid, pos in got.items()} == clf._cluster_loci(group)


def _span_over_positions(frame, positions):
    """The reductions classify_uecc_mecc runs inline over a locus' positions.

    The production code does not call a helper for these - it keeps them inline
    because the mapq and the coordinate reductions need different fallbacks when
    a column is missing. This mirrors the arithmetic so the property that made
    the frame slice unnecessary stays under test.
    """
    return (
        int(np.nanmax(frame["mapq"].to_numpy()[positions])),
        str(frame[ColumnStandard.CHR].to_numpy()[positions][0]),
        int(np.nanmin(frame[ColumnStandard.START0].to_numpy()[positions])),
        int(np.nanmax(frame[ColumnStandard.END0].to_numpy()[positions])),
    )


class TestLocusStatsFromArrays:
    """The Uecc branch reads six values off a locus; none needs a frame slice.

    `locus_df = group.loc[idxs]` was taken per classified query - and 95.6% of
    queries classify - only to read mapq.max(), chr.iloc[0], start0.min(),
    end0.max() and the alignment_length argmax. Those are array reductions.
    """

    @staticmethod
    def _frame():
        return pd.DataFrame(
            {
                "mapq": [30, 60, 45],
                ColumnStandard.CHR: ["chr7", "chr7", "chr7"],
                ColumnStandard.START0: [500, 300, 900],
                ColumnStandard.END0: [1500, 1200, 1800],
                "alignment_length": [400, 900, 900],
            },
            index=[5, 6, 7],
        )

    def test_reductions_match_the_frame_version(self):
        clf = UMeccClassifier()
        frame = self._frame()
        positions = [0, 1, 2]

        got = _span_over_positions(frame, positions)
        assert got == (
            int(frame["mapq"].max()),
            str(frame[ColumnStandard.CHR].iloc[0]),
            int(frame[ColumnStandard.START0].min()),
            int(frame[ColumnStandard.END0].max()),
        )

    def test_representative_is_the_first_maximum(self):
        clf = UMeccClassifier()
        frame = self._frame()
        positions = [0, 1, 2]

        # idxmax returns the FIRST maximum: alignment_length 900 appears twice
        expected_label = frame["alignment_length"].idxmax()
        rel = clf._argmax_position(frame["alignment_length"].to_numpy(), positions)
        assert frame.index[rel] == expected_label

    def test_reductions_skip_nan_like_pandas(self):
        import numpy as np

        clf = UMeccClassifier()
        frame = pd.DataFrame(
            {
                "mapq": [float("nan"), 60.0, 45.0],
                ColumnStandard.CHR: ["chr7", "chr7", "chr7"],
                ColumnStandard.START0: [float("nan"), 300.0, 900.0],
                ColumnStandard.END0: [1500.0, float("nan"), 1800.0],
                "alignment_length": [float("nan"), 900.0, 400.0],
            }
        )
        positions = [0, 1, 2]
        got = _span_over_positions(frame, positions)
        assert got == (
            int(frame["mapq"].max()), str(frame[ColumnStandard.CHR].iloc[0]),
            int(frame[ColumnStandard.START0].min()), int(frame[ColumnStandard.END0].max()),
        )
        rel = clf._argmax_position(frame["alignment_length"].to_numpy(), positions)
        assert frame.index[rel] == frame["alignment_length"].idxmax()

    def test_subset_of_positions(self):
        clf = UMeccClassifier()
        frame = self._frame()
        positions = [1, 2]
        sliced = frame.iloc[positions]

        got = _span_over_positions(frame, positions)
        assert got == (
            int(sliced["mapq"].max()), str(sliced[ColumnStandard.CHR].iloc[0]),
            int(sliced[ColumnStandard.START0].min()), int(sliced[ColumnStandard.END0].max()),
        )
        rel = clf._argmax_position(frame["alignment_length"].to_numpy(), positions)
        assert frame.index[rel] == sliced["alignment_length"].idxmax()


class TestMissingGroupKeysMatchGroupby:
    """`x != x` only recognises float NaN.

    An object column carries None, for which `x != x` is False - those rows
    would bucket together and union into one locus - and pd.NA, for which it is
    neither True nor False and raises when used as a condition. groupby's
    default dropna=True drops all three.
    """

    @staticmethod
    def _buckets(missing):
        frame = pd.DataFrame(
            {"chr": [missing, missing, "chr1"], "strand": ["+", "+", "+"], "i": [0, 1, 2]}
        )
        expected = {
            str(key): list(sub["i"])
            for key, sub in frame.groupby(["chr", "strand"], sort=False)
        }
        chrom_vals = frame["chr"].to_numpy()
        strand_vals = frame["strand"].to_numpy()
        got: dict = {}
        for i in range(len(frame)):
            chrom, strand = chrom_vals[i], strand_vals[i]
            if pd.isna(chrom) or pd.isna(strand):
                continue
            got.setdefault((chrom, strand), []).append(i)
        return expected, {str(k): v for k, v in got.items()}

    @pytest.mark.parametrize("missing", [np.nan, None, pd.NA], ids=["nan", "None", "pd.NA"])
    def test_every_missing_value_is_dropped_like_groupby(self, missing):
        expected, got = self._buckets(missing)
        assert got == expected

    @pytest.mark.parametrize("missing", [None, pd.NA], ids=["None", "pd.NA"])
    def test_the_identity_check_would_not_have(self, missing):
        """Guard the guard: `x != x` fails on exactly these two."""
        try:
            skipped = bool(missing != missing)
        except TypeError:
            return  # pd.NA raises, which is the failure this guards against
        assert skipped is False
