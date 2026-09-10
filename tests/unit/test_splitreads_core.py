from pathlib import Path
import sys

import pandas as pd
import pytest
import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from circleseeker.modules.splitreads_core import (
    SplitReadsConfig,
    Region,
    PdRegion,
    is_overlapping,
    any_overlapping_range,
    combine_region_strand,
    make_group_order,
    split_list_step_window,
    get_sum_pattern,
    get_idx_longest_pattern,
    len_loci,
    format_merge_region,
    majority_strand,
    reverse_strand,
    check_breakpoint_direction,
    chk_circular_subgraph,
    SplitReadsCore,
    _process_read_hits,
    _sorted_hits_frame,
    PdRegion,
    combine_region_strand,
)

import logging


# ============================================================================
# TestSplitReadsConfig
# ============================================================================


class TestSplitReadsConfig:
    def test_from_dict_basic(self):
        cfg = SplitReadsConfig.from_dict({"preset": "map-hifi", "mapq": 30})
        assert cfg.preset == "map-hifi"
        assert cfg.mapq == 30

    def test_from_dict_defaults(self):
        cfg = SplitReadsConfig.from_dict({})
        assert cfg.preset == "map-hifi"
        assert cfg.mapq == 20
        assert cfg.exclude_chrs == ""
        assert cfg.allow_gap == 10
        assert cfg.allow_overlap == 10
        assert cfg.min_region_size == 200
        assert cfg.overlap_check_size == 50
        assert cfg.breakpoint_depth == 5
        assert cfg.average_depth == 5.0
        assert cfg.threads == 0
        assert cfg.skip_variant is True

    def test_from_dict_extra_keys_ignored(self):
        cfg = SplitReadsConfig.from_dict({"foo": "bar", "preset": "map-hifi"})
        assert cfg.preset == "map-hifi"

    def test_from_dict_partial(self):
        cfg = SplitReadsConfig.from_dict({"mapq": 50})
        assert cfg.mapq == 50
        assert cfg.preset == "map-hifi"
        assert cfg.threads == 0
        assert cfg.allow_gap == 10

    def test_default_threads_zero(self):
        cfg = SplitReadsConfig()
        assert cfg.threads == 0


# ============================================================================
# TestRegion
# ============================================================================


class TestRegion:
    def test_region_attributes(self):
        r = Region("chr1", 100, 200)
        assert r.chrom == "chr1"
        assert r.start == 100
        assert r.end == 200

    def test_get_str_region(self):
        r = Region("chr1", 100, 200)
        assert r.get_str_region() == "chr1_100_200"

    def test_region_different_chrom(self):
        r = Region("chrX", 0, 500)
        assert r.chrom == "chrX"
        assert r.start == 0
        assert r.end == 500


# ============================================================================
# TestPdRegion
# ============================================================================


class TestPdRegion:
    @pytest.fixture()
    def sample_series(self):
        return pd.Series({
            "readid": "read1", "q_len": 1000, "q_start": 100, "q_end": 300,
            "ref": "chr1", "r_start": 1000, "r_end": 1500,
            "matchLen": 190, "blockLen": 200, "mapq": 60, "strand": 1,
        })

    def test_calc_q_len(self, sample_series):
        pr = PdRegion(sample_series)
        assert pr.calc_q_len() == 200

    def test_cal_r_len(self, sample_series):
        pr = PdRegion(sample_series)
        assert pr.cal_r_len() == 500

    def test_get_list_attr(self, sample_series):
        pr = PdRegion(sample_series)
        attrs = pr.get_list_attr()
        assert len(attrs) == 11
        assert attrs == ["read1", 1000, 100, 300, "chr1", 1000, 1500, 190, 200, 60, 1]

    def test_pdregion_stores_all_attributes(self, sample_series):
        pr = PdRegion(sample_series)
        assert pr.readid == "read1"
        assert pr.q_len == 1000
        assert pr.q_start == 100
        assert pr.q_end == 300
        assert pr.ref == "chr1"
        assert pr.r_start == 1000
        assert pr.r_end == 1500
        assert pr.matchLen == 190
        assert pr.blockLen == 200
        assert pr.mapq == 60
        assert pr.strand == 1


# ============================================================================
# TestIsOverlapping
# ============================================================================


class TestIsOverlapping:
    def test_exact_match(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 0, "chr1", 100, 200) is True

    def test_within_offset(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 10, "chr1", 105, 195) is True

    def test_outside_offset(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 5, "chr1", 200, 300) is False

    def test_different_chrom(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 100, "chr2", 100, 200) is False

    def test_zero_offset_exact(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 0, "chr1", 100, 200) is True

    def test_boundary_at_offset(self):
        r = Region("chr1", 100, 200)
        assert is_overlapping(r, 10, "chr1", 110, 210) is True


# ============================================================================
# TestAnyOverlappingRange
# ============================================================================


class TestAnyOverlappingRange:
    def test_complete_overlap(self):
        assert any_overlapping_range(100, 200, 100, 200) is True

    def test_partial_overlap(self):
        assert any_overlapping_range(100, 200, 150, 250) is True

    def test_adjacent_touch(self):
        assert any_overlapping_range(100, 200, 200, 300) is True

    def test_no_overlap(self):
        assert any_overlapping_range(100, 200, 300, 400) is False

    def test_contained(self):
        assert any_overlapping_range(100, 400, 200, 300) is True


# ============================================================================
# TestCombineRegionStrand
# ============================================================================


class TestCombineRegionStrand:
    def test_positive_strand(self):
        s = pd.Series({"ref": "chr1", "r_start": 100, "r_end": 200, "strand": "+"})
        assert combine_region_strand(s) == "chr1_100_200_+"

    def test_negative_strand(self):
        s = pd.Series({"ref": "chr1", "r_start": 100, "r_end": 200, "strand": "-"})
        assert combine_region_strand(s) == "chr1_100_200_-"


# ============================================================================
# TestMakeGroupOrder
# ============================================================================


class TestMakeGroupOrder:
    def test_single_region(self):
        result = make_group_order(["chr1_100_200_+"], 50)
        assert result == ["chr1_100_200_+"]

    def test_overlapping_regions(self):
        result = make_group_order(["chr1_100_200_+", "chr1_110_210_-"], 50)
        # Second region overlaps first within offset=50, so it maps to first region's coords + its own strand
        assert result[0] == "chr1_100_200_+"
        assert result[1] == "chr1_100_200_-"

    def test_non_overlapping(self):
        result = make_group_order(["chr1_100_200_+", "chr1_500_600_+"], 50)
        assert len(result) == 2
        assert result[0] == "chr1_100_200_+"
        assert result[1] == "chr1_500_600_+"

    def test_empty_list(self):
        result = make_group_order([], 50)
        assert result == []

    def test_cross_chromosome(self):
        result = make_group_order(["chr1_100_200_+", "chr2_100_200_+"], 50)
        assert len(result) == 2
        assert result[0] == "chr1_100_200_+"
        assert result[1] == "chr2_100_200_+"


# ============================================================================
# TestSplitListStepWindow
# ============================================================================


class TestSplitListStepWindow:
    def test_basic_split(self):
        result = split_list_step_window([1, 2, 3, 4, 5, 6], step=1, window=3)
        assert result == [[1, 2, 3], [3, 4, 5], [5, 6]]

    def test_no_overlap(self):
        result = split_list_step_window([1, 2, 3, 4], step=0, window=2)
        assert result == [[1, 2], [3, 4]]

    def test_single_element(self):
        result = split_list_step_window([1], step=0, window=1)
        assert result == [[1]]

    def test_empty_list(self):
        result = split_list_step_window([], step=0, window=2)
        assert result == []


# ============================================================================
# TestGetSumPattern
# ============================================================================


class TestGetSumPattern:
    def test_repeated_single_pattern(self):
        result = get_sum_pattern(["A", "A", "A"])
        assert result == [True, "A", 3]

    def test_no_repeat(self):
        result = get_sum_pattern(["A", "B", "C"])
        assert result == [False, "", 0]

    def test_ctc_pattern(self):
        result = get_sum_pattern(["A", "B", "A", "B"])
        assert result[0] is True

    def test_single_element(self):
        result = get_sum_pattern(["A"])
        assert result == [False, "", 0]

    def test_repeated_multi_pattern(self):
        result = get_sum_pattern(["A", "B", "A", "B", "A", "B"])
        assert result[0] is True


# ============================================================================
# TestGetIdxLongestPattern
# ============================================================================


class TestGetIdxLongestPattern:
    def test_simple_match(self):
        result = get_idx_longest_pattern(["A", "B", "A", "B"], ["A", "B"])
        assert len(result) > 0
        # Should match all four indices as two consecutive AB patterns merge
        assert result == [0, 1, 2, 3]

    def test_no_match(self):
        result = get_idx_longest_pattern(["A", "B", "C"], ["X"])
        assert result == []

    def test_intermittent_match(self):
        result = get_idx_longest_pattern(["A", "B", "C", "A", "B"], ["A", "B"])
        # Two matches at [0,1] and [3,4] but not contiguous, return longest
        assert len(result) == 2

    def test_consecutive_merge(self):
        result = get_idx_longest_pattern(["A", "B", "A", "B", "A", "B"], ["A", "B"])
        # Three consecutive AB matches -> all merged into one range
        assert result == [0, 1, 2, 3, 4, 5]


# ============================================================================
# TestLenLoci
# ============================================================================


class TestLenLoci:
    def test_single_region(self):
        assert len_loci("chr1_100_200_+") == 100

    def test_multiple_regions(self):
        assert len_loci("chr1_100_200_+,chr2_300_500_-") == 300

    def test_zero_length(self):
        assert len_loci("chr1_100_100_+") == 0


# ============================================================================
# TestFormatMergeRegion
# ============================================================================


class TestFormatMergeRegion:
    def test_single_region(self):
        assert format_merge_region("chr1_100_200_+") == "chr1:101-200_+"

    def test_multiple_regions(self):
        result = format_merge_region("chr1_100_200_+,chr2_300_500_-")
        assert result == "chr1:101-200_+,chr2:301-500_-"


# ============================================================================
# TestMajorityStrand
# ============================================================================


class TestMajorityStrand:
    def test_majority_positive(self):
        df = pd.DataFrame({"strand": [1, 1, -1]})
        assert majority_strand(df) == "+"

    def test_majority_negative(self):
        df = pd.DataFrame({"strand": [-1, -1, 1]})
        assert majority_strand(df) == "-"

    def test_all_same(self):
        df = pd.DataFrame({"strand": [1, 1]})
        assert majority_strand(df) == "+"


# ============================================================================
# TestReverseStrand
# ============================================================================


class TestReverseStrand:
    def test_plus_to_minus(self):
        assert reverse_strand("+") == "-"

    def test_minus_to_plus(self):
        assert reverse_strand("-") == "+"

    def test_mixed(self):
        assert reverse_strand("+-") == "-+"


# ============================================================================
# TestCheckBreakpointDirection
# ============================================================================


class TestCheckBreakpointDirection:
    def test_plus_plus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [1, 1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [0, 1],
            "ovl_3end": [1, 0],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "+_+"

    def test_minus_minus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [-1, -1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [1, 0],
            "ovl_3end": [0, 1],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        # The nodes are reversed, so both directions must be complemented.
        assert result[0] == ("regionB", "regionA", "+_+", True)

    def test_minus_plus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [-1, 1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [1, 1],
            "ovl_3end": [0, 0],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "-_+"

    def test_plus_minus_pattern(self):
        df = pd.DataFrame({
            "q_start": [100, 150],
            "q_end": [200, 250],
            "strand": [1, -1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [0, 0],
            "ovl_3end": [1, 1],
        })
        result = check_breakpoint_direction(df)
        assert len(result) == 1
        assert result[0][2] == "+_-"

    def test_no_overlap_breaks(self):
        df = pd.DataFrame({
            "q_start": [100, 500],
            "q_end": [200, 600],
            "strand": [1, 1],
            "mergeid": ["regionA", "regionB"],
            "ovl_5end": [0, 1],
            "ovl_3end": [1, 0],
        })
        result = check_breakpoint_direction(df)
        assert result == []


# ============================================================================
# TestChkCircularSubgraph
# ============================================================================


class TestChkCircularSubgraph:
    def test_self_loop(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "A")
        subgraph = G.to_undirected().subgraph(["A"])
        dict_pair_strand = {}
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        # Single node with self-loop
        assert result[3] is True  # sl (contains self-loop)
        assert result[4] is True  # cyclic (cycle_basis on single node with self-loop)

    def test_two_node_cycle(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "B")
        G.add_edge("B", "A")
        subgraph = G.to_undirected().subgraph(["A", "B"])
        dict_pair_strand = {
            ("A", "B"): {"+_+"},
            ("B", "A"): {"+_+"},
        }
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        assert result[4] is True  # cyclic

    def test_non_cycle(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "B")
        subgraph = G.to_undirected().subgraph(["A", "B"])
        # Only one direction, strand pairs don't form a cycle match
        dict_pair_strand = {("A", "B"): {"+_-"}}
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        assert result[4] is False  # not cyclic

    def test_three_node_cycle(self):
        G = nx.MultiDiGraph()
        G.add_edge("A", "B")
        G.add_edge("B", "C")
        G.add_edge("C", "A")
        subgraph = G.to_undirected().subgraph(["A", "B", "C"])
        dict_pair_strand = {}
        result = chk_circular_subgraph(G, subgraph, dict_pair_strand)
        assert result[4] is True  # cyclic (triangle forms cycle)


# ============================================================================
# TestWriteEmptyOutput
# ============================================================================


class TestWriteEmptyOutput:
    def _make_core(self):
        core = SplitReadsCore.__new__(SplitReadsCore)
        core.logger = logging.getLogger("test")
        return core

    def test_creates_file(self, tmp_path):
        core = self._make_core()
        result_path = core._write_empty_output(tmp_path)
        assert result_path.exists()
        assert result_path.name == "eccDNA_final.txt"

    def test_correct_columns(self, tmp_path):
        core = self._make_core()
        result_path = core._write_empty_output(tmp_path)
        df = pd.read_csv(result_path, sep="\t")
        expected_columns = [
            "id", "merge_region", "merge_len", "num_region",
            "ctc", "numreads", "totalbase", "coverage",
        ]
        assert list(df.columns) == expected_columns


class TestComponentResolutionDoesNotCopyTheWholeGraph:
    """The breakpoint graph must not be converted once per component.

    `_resolve_component_regions` used to run `G.to_undirected().subgraph(...)`,
    which materialized the complete breakpoint graph for every connected
    component.  On GlioSarc_P01_Tumor (19,125 components) that turned a
    seconds-long loop into 110 minutes.  Restricting first and converting the
    induced graph as a view is graph-theoretically identical, because inducing
    on a node set and dropping edge direction commute.
    """

    def test_component_regions_reuse_the_prebuilt_conversion(self):
        """Given the hoisted graph, a component must not convert G again.

        The loop builds one undirected graph and induces every component from
        it.  Converting per component was O(components x full_graph): 19,125
        components on GlioSarc_P01_Tumor spent 110 minutes here.
        """
        graph = nx.MultiDiGraph()
        graph.add_edge("target", "target")
        for index in range(100):
            node = f"unrelated_{index}"
            graph.add_edge(node, node)

        undirected = graph.to_undirected()

        def fail_full_graph_conversion(*_args, **_kwargs):
            pytest.fail("component resolution converted the complete graph again")

        graph.to_undirected = fail_full_graph_conversion

        result = SplitReadsCore._resolve_component_regions(
            1,
            {"target"},
            graph,
            {("target", "target"): {"+_+"}},
            {"target": "+"},
            undirected_graph=undirected,
        )

        assert result[1] == "target_+"
        assert result[4] is True
        assert result[5] is True

    def test_selection_equals_the_per_component_conversion(self):
        """The hoisted selection must return what the old per-component one did."""
        graph = nx.MultiDiGraph()
        for left, right in [
            ("a", "b"), ("b", "c"), ("c", "a"), ("b", "d"), ("d", "a"),
        ]:
            graph.add_edge(left, right, weight=3)
        pair_strand = {}
        for left, right in graph.edges():
            pair_strand[(left, right)] = {"+_+"}
            pair_strand[(right, left)] = {"+_+"}
        majority = {n: "+" for n in graph.nodes()}
        undirected = graph.to_undirected()

        for component in nx.connected_components(undirected):
            legacy = SplitReadsCore._resolve_component_regions(
                3, component, graph, pair_strand, majority
            )
            hoisted = SplitReadsCore._resolve_component_regions(
                3, component, graph, pair_strand, majority, undirected_graph=undirected
            )
            assert tuple(legacy) == tuple(hoisted)

    def test_component_selection_is_reproducible_and_complete(self):
        """Node order is fixed by construction; edges still match the view.

        The order deliberately differs from `subgraph()`'s: that one is the
        set-iteration order and moves with PYTHONHASHSEED. What must hold is
        that ours is the same on every run and describes the same graph.
        """
        graph = nx.MultiDiGraph()
        chord_edges = [
            ("n1", "n2"), ("n2", "n3"), ("n3", "n4"),
            ("n4", "n5"), ("n5", "n6"), ("n6", "n1"),
            ("n2", "n5"),
        ]
        for left, right in chord_edges:
            graph.add_edge(left, right, weight=3)
        component = set(graph.nodes())
        undirected = graph.to_undirected()

        first = SplitReadsCore._undirected_component(undirected, component)
        second = SplitReadsCore._undirected_component(undirected, set(component))
        view = undirected.subgraph(component)

        assert list(first.nodes()) == list(second.nodes()) == sorted(component)
        for node in sorted(component):
            assert list(first.adj[node]) == list(second.adj[node])
        assert set(first.nodes()) == set(view.nodes())
        assert sorted(map(sorted, first.edges())) == sorted(map(sorted, view.edges()))
        assert nx.cycle_basis(nx.Graph(first)) == nx.cycle_basis(nx.Graph(second))


    def test_region_order_matches_whole_graph_conversion(self):
        """Node order feeds the region string, so both selections must agree."""
        graph = nx.MultiDiGraph()
        edges = [
            ("chr1_100_200", "chr1_300_400"),
            ("chr1_300_400", "chr1_500_600"),
            ("chr1_500_600", "chr1_100_200"),
            ("chr2_10_20", "chr2_30_40"),
            ("chr2_30_40", "chr2_10_20"),
        ]
        for left, right in edges:
            graph.add_edge(left, right, weight=3)

        pair_strand = {}
        for left, right in edges:
            pair_strand[(left, right)] = {"+_+"}
            pair_strand[(right, left)] = {"+_+"}
        majority = {node: "+" for node in graph.nodes()}

        for component in nx.connected_components(graph.to_undirected()):
            optimized = SplitReadsCore._resolve_component_regions(
                7, component, graph, pair_strand, majority
            )
            # The node order deliberately no longer matches the view's: that one
            # is set-iteration order and moves with PYTHONHASHSEED.
            induced = SplitReadsCore._undirected_component(graph.to_undirected(), component)
            assert list(induced.nodes()) == sorted(component)
            assert set(induced.nodes()) == set(graph.to_undirected().subgraph(component).nodes())
            assert optimized[2] == len(component)


class TestEccdnaStatsIndexesTheReadTable:
    """Per-node lookups must not rescan the whole read/region table.

    The loop ran `read_merged_ins_df[read_merged_ins_df["mergeid"] == node]`
    once per component node.  On GlioSarc_P01_Tumor that is 24,629 scans of a
    321,724-row object column, the same quadratic shape already fixed in
    ecc_dedup and in the v2 packager.
    """

    @staticmethod
    def _read_table(tally_cls, mergeids, rows_per_id=4):
        records = []
        for position, mergeid in enumerate(mergeids):
            for repeat in range(rows_per_id):
                records.append(
                    {
                        "mergeid": tally_cls(mergeid),
                        "readid": f"read{position}_{repeat}",
                        "q_start": repeat * 10,
                        "q_end": repeat * 10 + 7,
                    }
                )
        return pd.DataFrame(records)

    def test_stats_do_not_rescan_the_table_for_every_node(self):
        class Tally(str):
            calls = 0

            def __eq__(self, other):
                Tally.calls += 1
                return str.__eq__(self, other)

            def __hash__(self):
                return str.__hash__(self)

        mergeids = [f"chr1_{start}_{start + 50}" for start in range(0, 400, 50)]
        read_table = self._read_table(Tally, mergeids)
        summary = pd.DataFrame(
            [
                {
                    "id": f"ec{index}",
                    "regions": f"{mergeid}_+",
                    "num_nodes": 1,
                    "is_cyclic": True,
                }
                for index, mergeid in enumerate(mergeids, start=1)
            ]
        )

        Tally.calls = 0
        SplitReadsCore._collect_eccdna_stats(summary, read_table)

        # Building the index hashes the column a couple of times (isin, then
        # groupby), so the count is a small multiple of the row count and its
        # exact value moves with PYTHONHASHSEED. What must not happen is the
        # per-node scan, which is rows x nodes.
        per_node_scan = len(read_table) * len(mergeids)
        assert Tally.calls < per_node_scan / 3, (
            f"{Tally.calls} comparisons for {len(read_table)} rows and "
            f"{len(mergeids)} nodes; a per-node scan would be {per_node_scan}"
        )

    def test_stats_values_match_the_per_node_scan(self):
        mergeids = ["chr1_100_200", "chr1_300_400", "chr5_10_60"]
        read_table = self._read_table(str, mergeids, rows_per_id=3)
        summary = pd.DataFrame(
            [
                {
                    "id": "ec1",
                    "regions": "chr1_100_200_+,chr1_300_400_-",
                    "num_nodes": 2,
                    "is_cyclic": True,
                },
                {
                    "id": "ec2",
                    "regions": "chr5_10_60_+",
                    "num_nodes": 1,
                    "is_cyclic": False,
                },
            ]
        )

        results = SplitReadsCore._collect_eccdna_stats(summary, read_table)

        # ec1 spans two 100 bp regions and draws 3 reads from each; every read
        # contributes q_end - q_start = 7 bases.
        assert results[0][:5] == ("ec1", "chr1_100_200_+,chr1_300_400_-", 200, 2, "True")
        assert results[0][5] == 6
        assert results[0][6] == 42
        assert results[0][7] == "0.21"
        assert results[1][:5] == ("ec2", "chr5_10_60_+", 50, 1, "False")
        assert results[1][5] == 3
        assert results[1][6] == 21
        assert results[1][7] == "0.42"


class TestSingletonComponentCyclicity:
    """A one-node component is cyclic exactly when it carries a self-loop.

    The single-node branch of `chk_circular_subgraph` rebuilt a directed and
    then an undirected graph just to run `cycle_basis`, and the undirected
    conversion deep-copies graph-level attributes.  Most GlioSarc components
    are singletons, so this ran ~19k times per sample.  These cases pin the
    behaviour the cheap self-loop test has to reproduce.
    """

    @staticmethod
    def _single_node_graph(with_self_loop: bool):
        graph = nx.MultiGraph()
        graph.add_node("chr1_100_200")
        if with_self_loop:
            graph.add_edge("chr1_100_200", "chr1_100_200")
        return graph

    def test_self_loop_singleton_is_cyclic(self):
        graph = self._single_node_graph(with_self_loop=True)
        result = chk_circular_subgraph(
            nx.MultiDiGraph(graph),
            graph,
            {("chr1_100_200", "chr1_100_200"): {"+_+"}},
        )
        assert result[1] == 1
        assert result[3] is True
        assert result[4] is True

    def test_bare_singleton_is_not_cyclic(self):
        graph = self._single_node_graph(with_self_loop=False)
        result = chk_circular_subgraph(nx.MultiDiGraph(graph), graph, {})
        assert result[1] == 1
        assert result[3] is False
        assert result[4] is False

    def test_empty_component_is_not_cyclic(self):
        graph = nx.MultiGraph()
        result = chk_circular_subgraph(nx.MultiDiGraph(), graph, {})
        assert result[1] == 0
        assert result[3] is False
        assert result[4] is False


class TestCycleDetectionKeepsTheDirectedRoundTrip:
    """cycle_basis must keep seeing the adjacency the stored outputs came from.

    Replacing `nx.DiGraph(g).to_undirected()` with `nx.Graph(g)` looked
    equivalent - both collapse parallel edges, and on small graphs the
    adjacency order matches. It is not: an end-to-end rerun of
    ara_UMC_5200_rep1_10X_HiFi against its stored 44c79b8 output changed the
    traversal of a 9-segment inferred CeccDNA, which reordered
    eccDNA_regions.csv and renumbered two CeccDNA entries. The round trip
    builds two directed edges per undirected edge and merges them back, and
    the merged order is not in general the original one.

    cycle_basis walks that adjacency and its first cycle becomes the region
    string, so this conversion stays as it was.
    """

    def test_traversal_uses_the_round_trip(self):
        graph = nx.MultiGraph()
        for left, right in [("c", "a"), ("a", "b"), ("b", "c"), ("b", "d"), ("d", "a")]:
            graph.add_edge(left, right)

        result = chk_circular_subgraph(nx.MultiDiGraph(graph), graph, {})

        # print_nodes comes from cycle_basis over the round-trip graph
        expected = nx.cycle_basis(nx.DiGraph(graph).to_undirected())
        assert expected, "test graph must contain a cycle"
        if len(expected[0]) == len(list(graph.nodes())):
            assert result[0] == ",".join(expected[0])


class TestEccdnaStatsIndexesTheReadTable:
    """Per-node lookups must not rescan the whole read/region table.

    The loop ran `read_merged_ins_df[read_merged_ins_df["mergeid"] == node]`
    once per component node.  On GlioSarc_P01_Tumor that is 24,629 scans of a
    321,724-row object column, the same quadratic shape already fixed in
    ecc_dedup and in the v2 packager.
    """

    @staticmethod
    def _read_table(tally_cls, mergeids, rows_per_id=4):
        records = []
        for position, mergeid in enumerate(mergeids):
            for repeat in range(rows_per_id):
                records.append(
                    {
                        "mergeid": tally_cls(mergeid),
                        "readid": f"read{position}_{repeat}",
                        "q_start": repeat * 10,
                        "q_end": repeat * 10 + 7,
                    }
                )
        return pd.DataFrame(records)

    def test_stats_do_not_rescan_the_table_for_every_node(self):
        class Tally(str):
            calls = 0

            def __eq__(self, other):
                Tally.calls += 1
                return str.__eq__(self, other)

            def __hash__(self):
                return str.__hash__(self)

        mergeids = [f"chr1_{start}_{start + 50}" for start in range(0, 400, 50)]
        read_table = self._read_table(Tally, mergeids)
        summary = pd.DataFrame(
            [
                {
                    "id": f"ec{index}",
                    "regions": f"{mergeid}_+",
                    "num_nodes": 1,
                    "is_cyclic": True,
                }
                for index, mergeid in enumerate(mergeids, start=1)
            ]
        )

        Tally.calls = 0
        SplitReadsCore._collect_eccdna_stats(summary, read_table)

        # Building the index hashes the column a couple of times (isin, then
        # groupby), so the count is a small multiple of the row count and its
        # exact value moves with PYTHONHASHSEED. What must not happen is the
        # per-node scan, which is rows x nodes.
        per_node_scan = len(read_table) * len(mergeids)
        assert Tally.calls < per_node_scan / 3, (
            f"{Tally.calls} comparisons for {len(read_table)} rows and "
            f"{len(mergeids)} nodes; a per-node scan would be {per_node_scan}"
        )

    def test_stats_values_match_the_per_node_scan(self):
        mergeids = ["chr1_100_200", "chr1_300_400", "chr5_10_60"]
        read_table = self._read_table(str, mergeids, rows_per_id=3)
        summary = pd.DataFrame(
            [
                {
                    "id": "ec1",
                    "regions": "chr1_100_200_+,chr1_300_400_-",
                    "num_nodes": 2,
                    "is_cyclic": True,
                },
                {
                    "id": "ec2",
                    "regions": "chr5_10_60_+",
                    "num_nodes": 1,
                    "is_cyclic": False,
                },
            ]
        )

        results = SplitReadsCore._collect_eccdna_stats(summary, read_table)

        # ec1 spans two 100 bp regions and draws 3 reads from each; every read
        # contributes q_end - q_start = 7 bases.
        assert results[0][:5] == ("ec1", "chr1_100_200_+,chr1_300_400_-", 200, 2, "True")
        assert results[0][5] == 6
        assert results[0][6] == 42
        assert results[0][7] == "0.21"
        assert results[1][:5] == ("ec2", "chr5_10_60_+", 50, 1, "False")
        assert results[1][5] == 3
        assert results[1][6] == 21
        assert results[1][7] == "0.42"


class TestSingletonComponentCyclicity:
    """A one-node component is cyclic exactly when it carries a self-loop.

    The single-node branch of `chk_circular_subgraph` rebuilt a directed and
    then an undirected graph just to run `cycle_basis`, and the undirected
    conversion deep-copies graph-level attributes.  Most GlioSarc components
    are singletons, so this ran ~19k times per sample.  These cases pin the
    behaviour the cheap self-loop test has to reproduce.
    """

    @staticmethod
    def _single_node_graph(with_self_loop: bool):
        graph = nx.MultiGraph()
        graph.add_node("chr1_100_200")
        if with_self_loop:
            graph.add_edge("chr1_100_200", "chr1_100_200")
        return graph

    def test_self_loop_singleton_is_cyclic(self):
        graph = self._single_node_graph(with_self_loop=True)
        result = chk_circular_subgraph(
            nx.MultiDiGraph(graph),
            graph,
            {("chr1_100_200", "chr1_100_200"): {"+_+"}},
        )
        assert result[1] == 1
        assert result[3] is True
        assert result[4] is True

    def test_bare_singleton_is_not_cyclic(self):
        graph = self._single_node_graph(with_self_loop=False)
        result = chk_circular_subgraph(nx.MultiDiGraph(graph), graph, {})
        assert result[1] == 1
        assert result[3] is False
        assert result[4] is False

    def test_empty_component_is_not_cyclic(self):
        graph = nx.MultiGraph()
        result = chk_circular_subgraph(nx.MultiDiGraph(), graph, {})
        assert result[1] == 0
        assert result[3] is False
        assert result[4] is False


class TestEccdnaStatsIndexesTheReadTable:
    """Per-node lookups must not rescan the whole read/region table.

    The loop ran `read_merged_ins_df[read_merged_ins_df["mergeid"] == node]`
    once per component node.  On GlioSarc_P01_Tumor that is 24,629 scans of a
    321,724-row object column, the same quadratic shape already fixed in
    ecc_dedup and in the v2 packager.
    """

    @staticmethod
    def _read_table(tally_cls, mergeids, rows_per_id=4):
        records = []
        for position, mergeid in enumerate(mergeids):
            for repeat in range(rows_per_id):
                records.append(
                    {
                        "mergeid": tally_cls(mergeid),
                        "readid": f"read{position}_{repeat}",
                        "q_start": repeat * 10,
                        "q_end": repeat * 10 + 7,
                    }
                )
        return pd.DataFrame(records)

    def test_stats_do_not_rescan_the_table_for_every_node(self):
        class Tally(str):
            calls = 0

            def __eq__(self, other):
                Tally.calls += 1
                return str.__eq__(self, other)

            def __hash__(self):
                return str.__hash__(self)

        mergeids = [f"chr1_{start}_{start + 50}" for start in range(0, 400, 50)]
        read_table = self._read_table(Tally, mergeids)
        summary = pd.DataFrame(
            [
                {
                    "id": f"ec{index}",
                    "regions": f"{mergeid}_+",
                    "num_nodes": 1,
                    "is_cyclic": True,
                }
                for index, mergeid in enumerate(mergeids, start=1)
            ]
        )

        Tally.calls = 0
        SplitReadsCore._collect_eccdna_stats(summary, read_table)

        # Building the index hashes the column a couple of times (isin, then
        # groupby), so the count is a small multiple of the row count and its
        # exact value moves with PYTHONHASHSEED. What must not happen is the
        # per-node scan, which is rows x nodes.
        per_node_scan = len(read_table) * len(mergeids)
        assert Tally.calls < per_node_scan / 3, (
            f"{Tally.calls} comparisons for {len(read_table)} rows and "
            f"{len(mergeids)} nodes; a per-node scan would be {per_node_scan}"
        )

    def test_stats_values_match_the_per_node_scan(self):
        mergeids = ["chr1_100_200", "chr1_300_400", "chr5_10_60"]
        read_table = self._read_table(str, mergeids, rows_per_id=3)
        summary = pd.DataFrame(
            [
                {
                    "id": "ec1",
                    "regions": "chr1_100_200_+,chr1_300_400_-",
                    "num_nodes": 2,
                    "is_cyclic": True,
                },
                {
                    "id": "ec2",
                    "regions": "chr5_10_60_+",
                    "num_nodes": 1,
                    "is_cyclic": False,
                },
            ]
        )

        results = SplitReadsCore._collect_eccdna_stats(summary, read_table)

        # ec1 spans two 100 bp regions and draws 3 reads from each; every read
        # contributes q_end - q_start = 7 bases.
        assert results[0][:5] == ("ec1", "chr1_100_200_+,chr1_300_400_-", 200, 2, "True")
        assert results[0][5] == 6
        assert results[0][6] == 42
        assert results[0][7] == "0.21"
        assert results[1][:5] == ("ec2", "chr5_10_60_+", 50, 1, "False")
        assert results[1][5] == 3
        assert results[1][6] == 21
        assert results[1][7] == "0.42"


class TestSingletonComponentCyclicity:
    """A one-node component is cyclic exactly when it carries a self-loop.

    The single-node branch of `chk_circular_subgraph` rebuilt a directed and
    then an undirected graph just to run `cycle_basis`, and the undirected
    conversion deep-copies graph-level attributes.  Most GlioSarc components
    are singletons, so this ran ~19k times per sample.  These cases pin the
    behaviour the cheap self-loop test has to reproduce.
    """

    @staticmethod
    def _single_node_graph(with_self_loop: bool):
        graph = nx.MultiGraph()
        graph.add_node("chr1_100_200")
        if with_self_loop:
            graph.add_edge("chr1_100_200", "chr1_100_200")
        return graph

    def test_self_loop_singleton_is_cyclic(self):
        graph = self._single_node_graph(with_self_loop=True)
        result = chk_circular_subgraph(
            nx.MultiDiGraph(graph),
            graph,
            {("chr1_100_200", "chr1_100_200"): {"+_+"}},
        )
        assert result[1] == 1
        assert result[3] is True
        assert result[4] is True

    def test_bare_singleton_is_not_cyclic(self):
        graph = self._single_node_graph(with_self_loop=False)
        result = chk_circular_subgraph(nx.MultiDiGraph(graph), graph, {})
        assert result[1] == 1
        assert result[3] is False
        assert result[4] is False

    def test_empty_component_is_not_cyclic(self):
        graph = nx.MultiGraph()
        result = chk_circular_subgraph(nx.MultiDiGraph(), graph, {})
        assert result[1] == 0
        assert result[3] is False
        assert result[4] is False


class TestStatsIndexKeepsTheOldMatchingRules:
    """Indexing must not widen what the per-node scan used to match.

    The old loop only ever touched rows whose mergeid appeared in a region
    string, and compared them as written. Building the index over the whole
    table converts every row - including ones the old code never looked at -
    and stringifying the group key makes a numeric mergeid column match a
    textual node name that previously never matched.
    """

    @staticmethod
    def _summary():
        return pd.DataFrame(
            [{"id": "ec1", "regions": "chr1_100_200_+", "num_nodes": 1, "is_cyclic": True}]
        )

    def test_unreferenced_rows_are_not_converted(self):
        table = pd.DataFrame(
            {
                "mergeid": ["chr1_100_200", "chr9_0_0"],
                "readid": ["r1", "r2"],
                "q_start": ["0", "not-a-number"],   # only the second is unusable
                "q_end": ["7", "also-bad"],
            }
        )

        results = SplitReadsCore._collect_eccdna_stats(self._summary(), table)

        assert results[0][6] == 7      # totalbase from the referenced row only
        assert results[0][5] == 1

    def test_numeric_mergeid_column_still_does_not_match_a_text_node(self):
        table = pd.DataFrame(
            {
                "mergeid": [100, 200],           # numeric, as bedtools may infer
                "readid": ["r1", "r2"],
                "q_start": [0, 0],
                "q_end": [7, 7],
            }
        )

        results = SplitReadsCore._collect_eccdna_stats(self._summary(), table)

        # "chr1_100_200" is not any of these keys, so nothing is contributed -
        # exactly what `df["mergeid"] == "chr1_100_200"` produced.
        assert results[0][6] == 0
        assert results[0][5] == 0


class TestComponentNodeOrderIsDeterministic:
    """A component's node order must not depend on PYTHONHASHSEED.

    networkx's subgraph view iterates nodes through FilterAtlas.__iter__:

        node_ok_shorter = 2 * len(self.NODE_OK.nodes) < len(self._atlas)
        if node_ok_shorter:
            return (n for n in self.NODE_OK.nodes if n in self._atlas)

    `NODE_OK.nodes` is the set that show_nodes stored, so whenever a component
    holds fewer than half of the graph's nodes - which is every component on a
    real sample - the order is the set's iteration order, i.e. hash-dependent.

    That order reaches `nodes = list(subgraph.nodes())` and decides where a
    multi-segment cycle starts and how it is walked. Two runs of identical code
    on identical input produced different structures for an inferred 9-segment
    CeccDNA because of this; with PYTHONHASHSEED fixed they agreed.
    """

    @staticmethod
    def _graph_and_component():
        graph = nx.MultiGraph()
        for index in range(200):
            graph.add_edge(f"chr1_{index}_{index + 100}", f"chr2_{index}_{index + 100}")
        component = {f"chr1_{i}_{i + 100}" for i in range(3)} | {
            f"chr2_{i}_{i + 100}" for i in range(3)
        }
        return graph, component

    def test_component_nodes_follow_the_parent_graph_not_the_hash(self):
        """The order is the parent's own node order, filtered to the component.

        Not sorted: the parent's order is what every result until now was
        computed from, and it is a property of the graph rather than of the
        process. Only the hash-seeded branch is bypassed.
        """
        graph, component = self._graph_and_component()
        subgraph = SplitReadsCore._undirected_component(graph, component)

        expected = [node for node in graph if node in component]
        assert list(subgraph.nodes()) == expected
        # and that is not the sorted order, so this test would catch a
        # reintroduction of sorting as well as of set iteration
        assert expected != sorted(component)

    def test_matches_the_view_where_the_view_is_deterministic(self):
        """A component holding at least half the graph takes networkx's other
        branch, which filters the parent's order and is already deterministic.
        Everything an eventual result is read from must match it exactly.
        """
        graph = nx.MultiGraph()
        for left, right in [("n0", "n9"), ("n2", "n1"), ("n1", "n3"),
                            ("n0", "n2"), ("n3", "n9"), ("n2", "n3")]:
            graph.add_edge(left, right, weight=1)
        component = set(graph.nodes())
        assert 2 * len(component) >= len(graph)

        view = graph.subgraph(component)
        ours = SplitReadsCore._undirected_component(graph, component)

        assert list(ours.nodes()) == list(view.nodes())
        assert list(ours.edges(keys=True)) == list(view.edges(keys=True))
        assert dict(ours.degree()) == dict(view.degree())
        assert nx.cycle_basis(nx.DiGraph(ours).to_undirected()) == nx.cycle_basis(
            nx.DiGraph(view).to_undirected()
        )

    def test_edges_and_membership_are_unchanged(self):
        graph, component = self._graph_and_component()
        view = graph.subgraph(component)
        ours = SplitReadsCore._undirected_component(graph, component)

        assert set(ours.nodes()) == set(view.nodes())
        assert sorted(map(sorted, ours.edges())) == sorted(map(sorted, view.edges()))
        assert ours.number_of_edges() == view.number_of_edges()
        for node in component:
            assert ours.degree[node] == view.degree[node]

    def test_parallel_edges_and_self_loops_survive(self):
        graph = nx.MultiGraph()
        graph.add_edge("a", "b")
        graph.add_edge("a", "b")          # parallel
        graph.add_edge("a", "a")          # self loop
        graph.add_edge("c", "d")
        ours = SplitReadsCore._undirected_component(graph, {"a", "b"})

        assert ours.number_of_edges() == 3
        assert list(nx.selfloop_edges(ours)) == [("a", "a")]
        assert list(ours.nodes()) == ["a", "b"]


class TestReadHitsFrameIsBuiltOnce:
    """One frame per read, not one per helper.

    _check_read_pattern and _get_merge_all each did
    `pd.DataFrame(list_hit, columns=header)` + `sort_values("q_start")` +
    `reset_index()` over the same list_hit. fable measured 0.26-0.84 ms per
    read against 1.5-3 million reads with alignments, i.e. 10-25 minutes, half
    of it this duplicated construction.

    The sort stays a pandas sort: `sort_values` uses numpy quicksort, which is
    unstable above 16 elements, so replacing it with Python's stable `sorted`
    would change the order of reads whose q_start ties.
    """

    HEADER = ["readid", "q_len", "q_start", "q_end", "ref", "r_start", "r_end",
              "matchLen", "blockLen", "mapq", "strand"]

    @classmethod
    def _hits(cls, n=4):
        rows = []
        for index in range(n):
            rows.append([
                "read1", 5000, 1000 - index * 100, 1400 - index * 100,
                "chr1", 20000 + index * 500, 20400 + index * 500,
                400, 400, 60, "+",
            ])
        return rows

    def test_shared_frame_matches_a_freshly_built_one(self):
        hits = self._hits()
        shared = _sorted_hits_frame(hits)

        expected = pd.DataFrame(hits, columns=self.HEADER)
        expected = expected.sort_values(by="q_start")
        expected.reset_index(drop=True, inplace=True)

        pd.testing.assert_frame_equal(shared, expected)

    def test_results_are_unchanged_when_the_frame_is_passed_in(self):
        hits = self._hits()
        kwargs = dict(allow_gap=10, allow_overlap=10, ref_merge_distance=1000)

        without = _process_read_hits("read1", hits, **kwargs)
        with_frame = _process_read_hits("read1", hits, frame=_sorted_hits_frame(hits), **kwargs)

        assert without == with_frame

    def test_single_hit_and_empty_input(self):
        kwargs = dict(allow_gap=10, allow_overlap=10, ref_merge_distance=1000)
        one = self._hits(1)

        assert _process_read_hits("read1", one, **kwargs) == _process_read_hits(
            "read1", one, frame=_sorted_hits_frame(one), **kwargs
        )
        assert _process_read_hits("read1", [], **kwargs) == []

    def test_ties_keep_the_pandas_ordering(self):
        # 20 hits sharing q_start: above 16 elements numpy quicksort is unstable,
        # so the shared frame must come from the same pandas call, not sorted()
        hits = [["read1", 5000, 700, 900, f"chr{i}", i * 100, i * 100 + 200,
                 200, 200, 60, "+"] for i in range(20)]
        shared = _sorted_hits_frame(hits)
        expected = pd.DataFrame(hits, columns=self.HEADER).sort_values(by="q_start")
        expected.reset_index(drop=True, inplace=True)

        assert list(shared["ref"]) == list(expected["ref"])




class TestStrandChoiceIsDeterministic:
    """`dict_pair_strand` values are sets of strings like "+_+" and "-_+".

    A region pair carrying both tandem and inverted evidence holds more than
    one, and `next(iter(...))` / `list(...)` over a set of strings returns a
    different element per process because str hashing is seeded. The chosen
    strand is written into merge_region, so the choice must be defined.
    """

    @staticmethod
    def _one_node_regions(strands):
        graph = nx.MultiDiGraph()
        graph.add_edge("chrA_10_200", "chrA_10_200", weight=5)
        undirected = graph.to_undirected()
        result = SplitReadsCore._resolve_component_regions(
            1,
            {"chrA_10_200"},
            graph,
            {("chrA_10_200", "chrA_10_200"): set(strands)},
            {"chrA_10_200": "+"},
            undirected_graph=undirected,
        )
        return result[1]

    def test_a_pair_with_both_orientations_picks_the_same_one_every_time(self):
        # both insertion orders of the same set must give the same answer
        assert self._one_node_regions(["+_+", "-_+"]) == self._one_node_regions(
            ["-_+", "+_+"]
        )

    def test_the_choice_is_the_smallest_strand_string(self):
        assert self._one_node_regions(["+_+", "-_+"]) == "chrA_10_200_+"
        assert self._one_node_regions(["-_-", "-_+"]) == "chrA_10_200_-"

    def test_a_single_orientation_is_unaffected(self):
        assert self._one_node_regions(["-_-"]) == "chrA_10_200_-"
