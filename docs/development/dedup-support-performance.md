# Cluster support aggregation performance

The candidate-support repair following v1.5.1 added a per-cluster scan of the
entire representative table. With K clusters and R representative rows, that
writeback required O(K * R) comparisons. Uecc keeps one row per cluster, so this
became approximately quadratic on million-cluster samples.

`EccDedup._attach_cluster_support` now indexes representative read names once,
calls the existing `support_fields` once per cluster, and maps the resulting
columns back in bulk. Lookup and writeback scale with K + R; candidate parsing,
grouping and sorting retain their own costs. The same helper supports one-row
Uecc representatives and multi-row Mecc/Cecc representatives.

The change preserves candidate identity, representative selection, row order,
segment structure, supporting-read order, serialized candidate support and RCA
copy accounting. Unknown support totals do not replace the representative's
existing copy or repeat value. Conflicting candidate evidence still raises an
error. CD-HIT clustering, thresholds and detection/inference rules are unchanged.

The helper logs completed cluster counts every 10,000 clusters and reports its
elapsed time. It is still a serial aggregation routine; the global thread count
does not automatically parallelize it.

Validation on 2026-09-10:

- All 1,296 selected local tests passed; two LAST-dependent tests were skipped,
  and one test was deselected. Five new support-accounting tests cover repeated
  segments, multiple candidates from one read, read order, missing totals and
  conflicting evidence.
- Bounded Uecc, Mecc and Cecc subsets from the actual GS intermediate tables
  produced exactly equal DataFrames, including dtypes/order, and identical CSV
  serialization against frozen source `44c79b8`.
- Synthetic 1,000- and 4,000-cluster inputs exactly matched the support loop
  extracted from `44c79b8`. Timings at 4,000 clusters were 7.129 seconds before
  and 0.884 seconds after the change.
- The optimized support loop took 26.399 seconds for 100,000 clusters and
  281.471 seconds for 1,000,000 clusters on fat2. These timings cover support
  aggregation and writeback, not a complete GS pipeline run. The million-cluster
  legacy loop was not rerun.

GS resumes from its completed CD-HIT checkpoint using the optimized source.
Its upstream eight steps remain attributed to `44c79b8`; downstream execution
must record the new commit and source hashes. Previous results retain their
original provenance. Final resumed output requires the usual integrity checks.
