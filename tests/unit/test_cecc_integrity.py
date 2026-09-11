"""Regression cases for errors found by the 2026-09-08 historical-data audit."""

import json
from pathlib import Path

import pandas as pd
import pytest

from circleseeker.modules.cecc_build import CeccBuild, _metadata_for_candidate
from circleseeker.modules.ecc_dedup import CDHitClusters, EccDedup
from circleseeker.modules.ecc_output_formatter import (
    generate_cecc_bedpe,
    generate_reads_table,
    generate_fasta_files,
    format_output,
)
from circleseeker.modules.umc_process import (
    CeccProcessor,
    SequenceLibrary,
    UMCProcessConfig,
    extract_ring_sequence,
)
from circleseeker.utils.circular_structure import cycle, same_cycle
from circleseeker.utils.read_support import support_fields


def segments(eid="C1", read="readA", copies=3, order=(1, 2, 3, 4), strands=("+", "-", "+", "-")):
    starts = {1: 1000, 2: 9000, 3: 100, 4: 3000}
    return pd.DataFrame(
        [
            dict(
                eccDNA_id=eid,
                query_id=f"{read}|rep0|400|{copies}|circular",
                reads=read,
                chr=f"chr{locus}",
                start0=starts[locus],
                end0=starts[locus] + 100,
                strand=strand,
                length=400,
                copy_number=copies,
                segment_in_circle=i,
                num_segments=4,
                match_degree=100,
                confidence_score=0.9,
                eSeq="ACGT" * 100,
            )
            for i, (locus, strand) in enumerate(zip(order, strands), 1)
        ]
    )


def clusters(*groups):
    c = CDHitClusters()
    c.members = {str(i): list(ids) for i, ids in enumerate(groups)}
    c.rep_of = {str(i): ids[0] for i, ids in enumerate(groups)}
    return c


def test_candidate_does_not_borrow_another_repeat_metadata():
    q = "read|rep0|363|8|circular"
    meta = {"read": {"reads": "read", "length": 2760, "copy_number": 3}}
    assert _metadata_for_candidate(q, meta) == {"reads": "read", "length": 363, "copy_number": 8}
    with pytest.raises(ValueError, match="disagree"):
        _metadata_for_candidate(q, {q.removesuffix("|circular"): meta["read"]})


def test_missing_candidate_cannot_use_read_alias(tmp_path):
    fasta = tmp_path / "rings.fa"
    fasta.write_text(">read|rep0|4|2|circular\nACGTACGT\n>read|rep1|3|5|circular\nACAACA\n")
    lib = SequenceLibrary()
    lib.load_fasta(fasta)
    assert lib.find_sequence("read|rep1|3|5") == "ACAACA"
    assert lib.find_sequence("read") is None
    assert lib.find_sequence("read|rep9|3|5") is None


@pytest.mark.parametrize("sequence,length", [("ACGTACGT", 5), ("ACGTACGA", 4)])
def test_invalid_doubled_consensus_is_rejected(sequence, length):
    with pytest.raises(ValueError):
        extract_ring_sequence(sequence, 1, length)


@pytest.mark.parametrize("kind", ["Cecc", "Mecc"])
def test_dedup_preserves_each_segment_direction(kind):
    source = segments()
    result = EccDedup().process_mecc_cecc(source, clusters(["C1"]), kind)
    assert result.sort_values("segment_in_circle").strand.tolist() == ["+", "-", "+", "-"]


def test_export_preserves_cycle_and_each_read_support(tmp_path):
    source = pd.concat([segments(), segments("C2", "readB", 5)], ignore_index=True)
    d = EccDedup()
    result = d.process_mecc_cecc(source, clusters(["C1", "C2"]), "Cecc")
    result = d.merge_cecc_by_tolerance(result)
    result = d.dedupe_cecc_segments(result)
    result = d.renumber_eccdna_ids(result, "Cecc")
    d.write_cecc_outputs(result, tmp_path, "test")
    core = pd.read_csv(tmp_path / "test_CeccSegments.core.csv")
    assert core.sort_values("seg_index").chr.tolist() == ["chr1", "chr2", "chr3", "chr4"]
    reads = generate_reads_table(cecc_df=core)
    assert dict(zip(reads.read_name, reads.copy_number)) == {"readA": 3, "readB": 5}
    assert set(reads.eccDNA_copy_number) == {8}
    assert len(pd.read_csv(tmp_path / "test_CeccJunctions.bedpe", sep="\t", header=None)) == 4
    d._generate_unified_confirmed_table({"Cecc": result}, tmp_path, "test")
    packaged = tmp_path / "packaged"
    format_output(
        unified_csv=tmp_path / "test_eccDNA_Confirmed.csv",
        uecc_core_csv=None,
        mecc_sites_csv=None,
        cecc_segments_csv=tmp_path / "test_CeccSegments.core.csv",
        uecc_fasta=None,
        mecc_fasta=None,
        cecc_fasta=tmp_path / "test_CeccDNA_C.fasta",
        output_dir=packaged,
    )
    assert ">CeccDNA0001|" in (packaged / "CeccDNA/cecc.fasta").read_text()


@pytest.mark.parametrize("alignment_difference", ["position", "strand", "order", "segment_count"])
def test_cdhit_cluster_remains_one_call_with_representative_structure(
    tmp_path, alignment_difference
):
    representative = segments()
    member = segments("C2", "readB", 5)
    # Similar sequences can have different alignment representations in repeats.
    # CD-HIT has already assigned both candidates to the same sequence cluster.
    member["eSeq"] = "T" + member.iloc[0].eSeq[1:]
    if alignment_difference == "position":
        member["start0"] += 1000
        member["end0"] += 1000
    elif alignment_difference == "strand":
        member.loc[0, "strand"] = "-"
    elif alignment_difference == "order":
        member["segment_in_circle"] = [1, 3, 2, 4]
    else:
        member = member.iloc[:3].copy()
        member["num_segments"] = 3

    # Put the nonrepresentative first to detect accidental row/metadata mixing.
    source = pd.concat([member, representative], ignore_index=True)
    dedup = EccDedup()
    result = dedup.process_mecc_cecc(source, clusters(["C1", "C2"]), "Cecc")
    assert result.eccDNA_id.nunique() == 1
    assert set(result.query_id) == set(representative.query_id)
    assert cycle(result) == cycle(representative)
    assert set(result.eSeq) == set(representative.eSeq)
    assert set(result.length) == {400}
    assert set(result.num_segments) == {4}
    assert set(result.num_merged) == {2}
    support = json.loads(result.iloc[0].candidate_support)
    assert set(support) == {"readA|rep0|400|3", "readB|rep0|400|5"}
    assert support["readA|rep0|400|3"]["copy_number"] == 3
    assert support["readB|rep0|400|5"]["copy_number"] == 5
    dedup.write_cecc_outputs(result, tmp_path, "cluster")
    core = pd.read_csv(tmp_path / "cluster_CeccSegments.core.csv")
    reads = generate_reads_table(cecc_df=core)
    assert dict(zip(reads.read_name, reads.copy_number)) == {"readA": 3, "readB": 5}
    assert set(reads.eccDNA_copy_number) == {8}


def test_old_multiread_aggregate_does_not_become_per_read_value():
    core = pd.DataFrame([dict(eccDNA_id="C1", read_name="a;b", copy_number=8)])
    result = generate_reads_table(cecc_df=core)
    assert result.copy_number.isna().all()
    assert set(result.eccDNA_copy_number) == {8}


def test_same_candidate_is_not_counted_again_for_each_segment_or_merge():
    first = segments()
    fields = support_fields(first)
    assert fields["copy_number"] == 3
    duplicate = first.assign(**fields)
    combined = support_fields(pd.concat([duplicate, first], ignore_index=True))
    assert combined["copy_number"] == 3
    distinct = first.copy()
    distinct["query_id"] = "readA|rep1|400|5|circular"
    distinct["copy_number"] = 5
    assert support_fields(pd.concat([first, distinct]))["per_read_copy_number"] == "8"


def test_coordinate_merge_keeps_different_directed_topologies():
    a = segments(strands=("+",) * 4)
    b = segments("C2", "readB", order=(1, 3, 2, 4), strands=("+",) * 4)
    # A coordinate-only merge across separate sequence clusters still requires
    # compatible directed structures; it must not mix their segment paths.
    result = EccDedup().merge_cecc_by_tolerance(pd.concat([a, b], ignore_index=True))
    assert result.eccDNA_id.nunique() == 2


def test_coordinate_merge_preserves_the_baseline_structural_unit():
    a = segments()
    b = segments("C2", "readB")
    b["eSeq"] = "AGCT" * 100
    result = EccDedup().merge_cecc_by_tolerance(pd.concat([a, b], ignore_index=True))
    assert result.eccDNA_id.nunique() == 1
    assert set(result.copy_number) == {6}
    assert all(len(seq) == 400 for seq in result.eSeq)


def test_coordinate_merge_combines_noisy_rotated_copies_but_preserves_support():
    import random

    rng = random.Random(491)
    seq = "".join(rng.choices("ACGT", k=400))
    noisy = seq[:101] + ("A" if seq[101] != "A" else "C") + seq[102:]
    a = segments()
    b = segments("C2", "readB", copies=5)
    a["eSeq"] = seq
    b["eSeq"] = noisy[213:] + noisy[:213]
    result = EccDedup().merge_cecc_by_tolerance(pd.concat([a, b], ignore_index=True))
    assert result.eccDNA_id.nunique() == 1
    assert set(result.copy_number) == {8}


def test_reverse_complement_cycle_is_equivalent_but_single_strand_flip_is_not():
    a = cycle(segments())
    b = tuple(
        (ch, start, end, "-" if strand == "+" else "+") for ch, start, end, strand in reversed(a)
    )
    assert same_cycle(a, b)
    changed = (a[0][:3] + ("-",),) + a[1:]
    assert not same_cycle(a, changed)


def test_early_structure_grouping_preserves_representative_and_all_support():
    a = segments()
    b = segments("C2", "readB", 5)
    c = segments("C3", "readC", 7)
    c["eSeq"] = "T" + a.iloc[0].eSeq[1:]
    p = CeccProcessor(SequenceLibrary(), UMCProcessConfig())
    result = p.cluster_by_signature(pd.concat([a, b, c], ignore_index=True))
    assert result.query_id.nunique() == 1
    assert set(result.query_id) == set(a.query_id)
    assert cycle(result) == cycle(a)
    assert set(result.eSeq) == set(a.eSeq)
    assert set(result.length) == {400}
    combined = result.iloc[0]
    assert combined.copy_number == 15
    assert combined.per_read_copy_number == "3;5;7"
    assert set(json.loads(combined.candidate_support)) == {
        "readA|rep0|400|3", "readB|rep0|400|5", "readC|rep0|400|7"
    }


def test_repeated_locus_in_different_circle_positions_is_not_deleted():
    frame = pd.DataFrame(
        [
            dict(eccDNA_id="C1", chr="chr1", start0=10, end0=20, strand="+", segment_in_circle=1),
            dict(eccDNA_id="C1", chr="chr1", start0=10, end0=20, strand="-", segment_in_circle=2),
        ]
    )
    assert len(EccDedup().dedupe_cecc_segments(frame)) == 2


def test_negative_strand_endpoints_and_closing_edge(tmp_path):
    regions = pd.DataFrame(
        [
            dict(
                eccDNA_id="C1",
                region_idx=1,
                chr="chrA",
                start=100,
                end=200,
                strand="-",
                role="head",
            ),
            dict(
                eccDNA_id="C1",
                region_idx=2,
                chr="chrB",
                start=300,
                end=400,
                strand="+",
                role="tail",
            ),
        ]
    )
    out = tmp_path / "junctions.bedpe"
    generate_cecc_bedpe(regions, out)
    records = pd.read_csv(out, sep="\t", header=None).values.tolist()
    assert records[0][:6] == ["chrA", 100, 101, "chrB", 300, 301]
    assert records[1][:6] == ["chrB", 399, 400, "chrA", 199, 200]


def test_identity_comes_from_aligned_bases(tmp_path):
    maf = tmp_path / "test.maf"
    maf.write_text(
        "a score=100\ns chr1 0 8 + 100 ACGTACGT\ns perfect 0 8 + 8 ACGTACGT\n\na score=40\ns chr1 0 8 + 100 ACGTACGT\ns half 0 8 + 8 ACGTTGCA\n\n"
    )
    parsed = CeccBuild()._parse_last_split_maf(maf)
    assert parsed["perfect"][0].identity == 100
    assert parsed["half"][0].identity == 50


def test_invalid_fasta_is_rejected_before_any_cecc_output(tmp_path):
    frame = segments()
    frame["eSeq"] = "ACGT"
    with pytest.raises(ValueError, match="length"):
        EccDedup().write_cecc_outputs(frame, tmp_path, "bad")
    assert list(tmp_path.iterdir()) == []


def test_final_packager_rejects_sequence_length_conflict(tmp_path):
    summary = pd.DataFrame(
        [dict(eccDNA_id="C1", type="Cecc", state="Confirmed", length=400, location="chr1:0-400(+)")]
    )
    with pytest.raises(ValueError, match="length"):
        generate_fasta_files({"C1": "ACGT"}, tmp_path, summary)
    assert list(tmp_path.iterdir()) == []


def test_format_output_validates_before_writing_anything(tmp_path):
    """A bad confirmed Cecc must stop the export with the output still empty.

    generate_fasta_files trusts `validated=True`, so format_output's own call is
    the only thing standing between a wrong sequence and the written tables.
    Calling the validation helper directly would not test that: moving the call
    below the `to_csv` lines leaves such a test passing, and leaves the tables
    on disk. This drives the real entry point and looks at the directory.
    """
    source = pd.concat([segments(), segments("C2", "readB", 5)], ignore_index=True)
    d = EccDedup()
    result = d.process_mecc_cecc(source, clusters(["C1", "C2"]), "Cecc")
    result = d.merge_cecc_by_tolerance(result)
    result = d.dedupe_cecc_segments(result)
    result = d.renumber_eccdna_ids(result, "Cecc")
    d.write_cecc_outputs(result, tmp_path, "test")
    d._generate_unified_confirmed_table({"Cecc": result}, tmp_path, "test")

    # truncate every ring so its length disagrees with the summary
    fasta = tmp_path / "test_CeccDNA_C.fasta"
    fasta.write_text(
        "".join(
            line if line.startswith(">") else line[:4]
            for line in fasta.read_text().splitlines(keepends=True)
        )
    )

    packaged = tmp_path / "packaged"
    with pytest.raises(ValueError):
        format_output(
            unified_csv=tmp_path / "test_eccDNA_Confirmed.csv",
            cecc_segments_csv=tmp_path / "test_CeccSegments.core.csv",
            cecc_fasta=fasta,
            output_dir=packaged,
        )

    written = [p for p in packaged.rglob("*") if p.is_file()] if packaged.exists() else []
    assert written == [], f"files written despite the failure: {written}"
