"""Regression tests for inferred direction and contradictory cycle evidence."""

from pathlib import Path
from itertools import permutations

import networkx as nx
import pandas as pd
import pysam
import pytest
from Bio import SeqIO

from circleseeker.modules.cecc_build import CeccBuild, LastAlignment
from circleseeker.modules.iecc_curator import curate_ecc_tables, generate_fasta_sequences
from circleseeker.modules.splitreads_adapter import write_overview_tsv
from circleseeker.modules.splitreads_core import SplitReadsCore, check_breakpoint_direction
from circleseeker.utils.circular_structure import same_cycle


def test_inferred_adapter_preserves_mixed_segment_directions(tmp_path: Path):
    raw = tmp_path / "eccDNA_final.txt"
    raw.write_text(
        "id\tmerge_region\tmerge_len\tnum_region\tnumreads\tcoverage\tctc\n"
        "simple\tchr1:1-120_-\t120\t1\t5\t5\tTrue\n"
        "mixed\tchr2:1-60_-;chr1:121-180_+\t120\t2\t5\t5\tTrue\n"
    )
    overview = write_overview_tsv(raw, tmp_path / "overview.tsv")
    simple, chimeric = curate_ecc_tables(overview)
    assert simple.iloc[0].strand == "-"
    assert (simple.iloc[0].start0, simple.iloc[0].end0) == (0, 120)
    assert chimeric.sort_values("seg_index").strand.tolist() == ["-", "+"]
    assert chimeric.sort_values("seg_index").chr.tolist() == ["chr2", "chr1"]


def test_inferred_fasta_respects_each_segment_direction(tmp_path: Path):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\n" + "AACG" * 50 + "\n>chr2\n" + "TTGC" * 50 + "\n")
    pysam.faidx(str(reference))
    simple = pd.DataFrame([
        dict(eccDNA_id="S", chr="chr1", start0=0, end0=120, strand="-", length=120)
    ])
    # Deliberately reverse dataframe order; traversal must follow seg_index.
    chimeric = pd.DataFrame([
        dict(eccDNA_id="C", chr="chr1", start0=120, end0=180, strand="+", seg_index=2),
        dict(eccDNA_id="C", chr="chr2", start0=0, end0=60, strand="-", seg_index=1),
    ])
    simple_path, chimeric_path = generate_fasta_sequences(
        simple, chimeric, reference, tmp_path / "result"
    )
    assert str(next(SeqIO.parse(simple_path, "fasta")).seq) == "CGTT" * 30
    assert str(next(SeqIO.parse(chimeric_path, "fasta")).seq) == "GCAA" * 15 + "AACG" * 15


def test_single_base_inferred_segment_keeps_its_explicit_direction(tmp_path: Path):
    raw = tmp_path / "eccDNA_final.txt"
    raw.write_text(
        "id\tmerge_region\tmerge_len\tnum_region\tnumreads\tcoverage\tctc\n"
        "mixed\tchr1:1-1_-;chr2:1-119_+\t120\t2\t5\t5\tTrue\n"
    )
    overview = write_overview_tsv(raw, tmp_path / "overview.tsv")
    _, chimeric = curate_ecc_tables(overview)
    first = chimeric.sort_values("seg_index").iloc[0]
    assert (first.start0, first.end0, first.strand) == (0, 1, "-")


def alignment(chrom, start, end, qstart, qend, length=400):
    return LastAlignment(
        chrom=chrom, ref_start=start, ref_end=end, query_start=qstart,
        query_end=qend, query_len=2 * length, score=100, identity=99.9, strand="+",
    )


def test_full_consensus_alignment_cannot_be_augmented_into_a_chimeric_cycle():
    # A partial A, short B, full-period A, and unrelated tail can create a
    # spurious A-B-A locus loop. A already covers every base of this consensus.
    alns = [alignment("chr1", 1080, 1380, 0, 300),
            alignment("chr2", 2000, 2020, 300, 320),
            alignment("chr1", 1000, 1400, 320, 720),
            alignment("chr3", 3000, 3080, 720, 800)]
    builder = CeccBuild()
    loci, _ = builder._assign_locus_id(alns)
    assert builder._find_cycle_in_doubled_sequence(alns, loci, 400) is None


def test_partial_start_and_complete_repeat_preserve_a_supported_cycle():
    alns = [alignment("chr1", 1100, 1200, 0, 100),
            alignment("chr2", 2000, 2200, 100, 300),
            alignment("chr1", 1000, 1200, 300, 500),
            alignment("chr2", 2000, 2200, 500, 700),
            alignment("chr1", 1000, 1100, 700, 800)]
    builder = CeccBuild()
    loci, _ = builder._assign_locus_id(alns)
    result = builder._find_cycle_in_doubled_sequence(alns, loci, 400)
    assert result is not None
    path, _, coordinates = result
    assert len(set(path)) == 2
    assert set(coordinates.values()) == {(300, 500), (100, 300)}


def test_opposite_reads_supply_the_same_oriented_breakpoint():
    forward = pd.DataFrame([
        dict(mergeid="A", strand=1, q_start=0, q_end=100, ovl_5end=1, ovl_3end=1),
        dict(mergeid="B", strand=1, q_start=100, q_end=200, ovl_5end=1, ovl_3end=1),
    ])
    reverse = forward.iloc[::-1].copy().reset_index(drop=True)
    reverse["strand"] = -1
    reverse["q_start"] = [0, 100]
    reverse["q_end"] = [100, 200]
    assert check_breakpoint_direction(reverse) == check_breakpoint_direction(forward)


@pytest.mark.parametrize("insertion_order", list(permutations(("Chr1", "Chr4", "Chr5"))))
def test_inferred_open_path_follows_supported_edges(insertion_order):
    graph = nx.MultiDiGraph()
    nodes = {chrom: f"{chrom}_100_200" for chrom in insertion_order}
    graph.add_nodes_from(nodes.values())
    graph.add_edge(nodes["Chr5"], nodes["Chr4"])
    graph.add_edge(nodes["Chr4"], nodes["Chr1"])
    pairs = {(nodes["Chr5"], nodes["Chr4"]): {"+_+"},
             (nodes["Chr4"], nodes["Chr1"]): {"+_+"}}
    majority = {nodes["Chr1"]: "+", nodes["Chr4"]: "+", nodes["Chr5"]: "-"}
    result = SplitReadsCore._resolve_component_regions(1, set(nodes.values()), graph, pairs, majority)
    actual = []
    for segment in result[1].split(","):
        chrom, start, end, strand = segment.rsplit("_", 3)
        actual.append((chrom, int(start), int(end), strand))
    expected = tuple((chrom, 100, 200, "+") for chrom in ["Chr1", "Chr5", "Chr4"])
    assert same_cycle(tuple(actual), expected)
    assert result[3] is True  # The observed path is resolved.
    assert result[5] is False  # The missing closing edge remains inferred.


def test_cycle_must_have_a_directionally_compatible_closing_edge():
    graph = nx.MultiDiGraph([("A", "B"), ("B", "C"), ("C", "A")])
    pairs = {("A", "B"): {"+_+"}, ("B", "C"): {"+_+"}, ("C", "A"): {"-_+"}}
    result = SplitReadsCore._resolve_component_regions(1, {"A", "B", "C"}, graph, pairs, {})
    assert result[1] == ""
    assert result[3] is False
    assert result[5] is False


def test_supported_cycle_includes_its_closing_direction():
    graph = nx.MultiDiGraph([("A", "B"), ("B", "C"), ("C", "A")])
    pairs = {("A", "B"): {"+_-"}, ("B", "C"): {"-_+"}, ("C", "A"): {"+_+"}}
    result = SplitReadsCore._resolve_component_regions(1, {"A", "B", "C"}, graph, pairs, {})
    assert result[1] == "A_+,B_-,C_+"
    assert result[3] is True
    assert result[5] is True


@pytest.mark.parametrize("pairs,closed", [
    ({("A", "B"): {"+_+"}}, False),
    ({("A", "B"): {"+_+"}, ("B", "A"): {"+_+"}}, True),
    # Opposite observations of one mixed-strand junction are still one edge.
    ({("A", "B"): {"+_-"}, ("B", "A"): {"+_-"}}, False),
    ({("A", "B"): {"+_-"}, ("B", "A"): {"-_+"}}, True),
])
def test_two_segment_closure_requires_both_distinct_junctions(pairs, closed):
    graph = nx.MultiDiGraph()
    graph.add_edges_from(pairs)
    result = SplitReadsCore._resolve_component_regions(1, {"A", "B"}, graph, pairs, {})
    assert result[3] is True
    assert result[5] is closed


@pytest.mark.parametrize("nodes", [("A", "B"), ("A", "B", "C")])
def test_edge_removed_by_support_filter_cannot_close_an_inferred_path(nodes):
    graph = nx.MultiDiGraph()
    graph.add_edges_from(zip(nodes, nodes[1:]))
    pairs = {(left, right): {"+_+"} for left, right in graph.edges()}
    # Evidence dictionaries retain observations below the graph's depth cutoff.
    pairs[nodes[-1], nodes[0]] = {"+_+"}
    result = SplitReadsCore._resolve_component_regions(1, set(nodes), graph, pairs, {})
    assert result[3] is True
    assert result[5] is False
