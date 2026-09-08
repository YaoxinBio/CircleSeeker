# Cecc integrity fixes after the v1.5.1 audit

This development branch starts at tag `v1.5.1`, commit
`87ae7f09c11a4795532a7524db1958c0c11a5b2d` (2026-08-16). It is a software
repair branch. The selected publication version is 1.5.2; final publication and
replacement of archived analyses require separate acceptance.

Its internal version is `1.5.1+ceccfix.20260908`, so new logs identify the patch
instead of reporting an unmodified historical `1.5.1`. The repair branch's
Conda recipe builds the local checkout with this same internal version; the
historical archive remains accessible at the unchanged `v1.5.1` tag. This branch
is evaluated from its own Git commit and Python source before any release work.

The 2026-09-08 audit replayed 59 supporting raw reads for 52 flagged Cecc
records. All 52 reproduced the historical wrong sequence and coordinates while
borrowing length/copy metadata from another repeat candidate of the same read.
An earlier audit had established three more examples. The tests in
`tests/unit/test_cecc_integrity.py` cover this failure and downstream defects.

## Data contracts

* Candidate identity is the full `read|repeat|length|copies` token. The optional
  `|circular` suffix is presentation only. A missing candidate must never fall
  back to another repeat under the base read name. Exact legacy query names can
  use their explicitly supplied metadata; they do not create aliases.
* A candidate sequence is two identical copies of its declared consensus.
  Extracting one circle validates that invariant and returns exactly its length.
* Segment direction belongs to each representative segment (or each Mecc locus).
  It is not a majority vote across different segments or loci.
* Segment order follows the representative query traversal. The equivalent
  opposite traversal reverses order **and flips every strand**. Sequence FASTA
  canonicalization can choose another origin/strand; it represents the same
  circular sequence, not a claim that FASTA position zero is the first locus.
* Cecc coordinate signatures preserve directed traversal and the baseline's
  structural catalogue unit. Early grouping does not require exact equality of
  consensus sequences. Every candidate's sequence and metadata are validated
  before grouping; one complete representative is retained and support is
  aggregated without replacing its sequence, length or segment directions.
* CD-HIT defines the subsequent sequence clusters using the existing thresholds.
  Alignment differences among members do not split those clusters into more
  molecule calls. Repeats can give even identical circular sequences different
  genomic placements. Keep the chosen member's complete directed structure,
  rather than requiring every supporting alignment to share its placement.
* The secondary coordinate-only merge across sequence clusters still requires
  equivalent directed cycles within the original 10 bp tolerance. Different
  directions or orders cannot merge merely because unordered coordinates match.
  This pass adds no nucleotide-variant threshold. Each retained representative's
  own sequence, length and directed path remain linked and validated.
* Overlapping or opposite-strand occurrences within a circle remain distinct.
  Only exact repeated segment records are removed; the old proximity rule could
  remove genuine occurrences.
* `candidate_support` records a union of candidate ID to read and RCA copies.
  Repeated segment rows and repeated representations of the **same candidate**
  count once. Different repeat candidates from one read contribute separately,
  as in the pipeline's candidate-support sum. Per-read copies sum those candidate
  contributions; the molecule aggregate sums reads. This is an RCA support
  quantity, not genomic copy number. Unknown historical per-read decomposition
  remains missing and must not be filled with the molecule aggregate.
* Cecc BEDPE describes the full directed circle: N segments produce N junctions,
  including last-to-first. Endpoints follow each segment's strand, using one-base
  zero-based half-open intervals at its entry/exit. A one-segment input has its
  self-closing junction. This is an explicit output contract change from N−1.
* LAST identity is measured from aligned MAF strings: exact A/C/G/T matches divided
  by alignment columns with at least one non-gap character, including gaps and
  ambiguous bases in the denominator. Unknown/non-ACGT pairs do not count as
  confirmed matches. LAST MAPQ is unavailable, so its fields and `low_mapq` are
  missing. The previous confidence formula relied on a fabricated MAPQ; its
  confidence is also missing rather than replaced by an unvalidated score.
* Confirmed Cecc output requires one valid representative sequence of the declared
  length. Missing or contradictory sequences fail before the Cecc writer creates
  files. Final packaged FASTA receives the same length check.
  The pipeline packager checks this before emitting final tables or BED files;
  its errors propagate as a failed step instead of publishing a fallback output.
  Duplicate normalized FASTA IDs are rejected rather than silently overwritten.
* Final FASTA excludes inferred candidates removed by unification; retained
  sequence IDs follow the same final catalogue as the tables and BED files.

## Validation and use

The baseline's 278 selected module tests passed before changes. Some old tests
explicitly expected bugs (constant identity, base-read fallback, dropping reverse
segments, or omitted closure). Their expectations have been corrected. Other
fixtures have been corrected to supply genuinely doubled sequences of their
stated length; invalid inputs are now tested for rejection separately.

The initial repair at `055110f` introduced a specificity regression. It split
existing CD-HIT clusters by directed alignment and 10 bp coordinate agreement,
and required exact consensus equality in early structural grouping. Neither
condition is required to keep the selected representative internally consistent.
Both additional conditions are removed in this candidate; full candidate identity,
representative structure and support accounting remain enforced.

Controlled replay of the same 19 completed benchmark scenarios gave:

| Source / controlled change | TP | FP |
|---|---:|---:|
| Historical v1.5.1 | 153,809 | 1,274 |
| 055110f repair, reproduced control | 153,815 | 1,624 |
| Remove CD-HIT cluster partition only | 153,811 | 1,283 |
| Also remove early exact-sequence grouping | 153,809 | 1,276 |

All 19 frozen controls reproduced the failed repair metrics. After both changes,
17 scenarios match all historical metrics exactly. The remaining scenarios are
ara_UMC_5200_rep3_30X_HiFi and human_U_10000_rep2_30X_HiFi, each with one additional
FP. These are net differences under a fixed evaluator, not a one-to-one attribution
of every newly numbered record. The cohort replays saved Cecc intermediates through
clustering, deduplication, unification and packaging. It reuses U/M and inference
only after checking that the supporting-read partition is unchanged; it is not
19 complete raw-read runs.

The corrected candidate passed 1,273 tests (one deselected) with 72.26% coverage,
exceeding the configured 50% threshold, and mypy passed all 59 source files.
Checks ran on `fat2`, explicitly selecting the isolated checkout. Four cluster
regression cases and the early-grouping case fail against 055110f and pass with
this repair. AST checks link the tested methods to the controlled ablations.
Regression tests also cover the actual pipeline packager and propagation of
sequence-validation failures.

All 55 previously flagged real records pass replay from their saved original
alignments through Cecc construction, U/M/C processing, CD-HIT, deduplication and
formatting. Every final sequence and directed path must belong to its selected
representative. In the HeLa_T7_HiFi / CeccDNA0055 case, one support candidate has an
alternative genomic placement; its identical canonical sequence and actual shared
CD-HIT cluster were independently verified. Requiring every support alignment to
have the representative's placement would recreate the faulty partition rule.
This scoped replay is not a full-library rebuild or a final catalogue replacement.

A complete raw-input run of ara_UMC_5200_rep1_30X_HiFi also passed: 16 stages,
16 threads, verified turbo storage under `/dev/shm`, 5,054 final records,
TP 5,016 and FP 38, with all type metrics matching historical results. Sequence
lengths, supporting-copy sums, IDs and directed closure checks passed.

Integrity alone is insufficient for release acceptance. A separate benchmark gate
checks TP, recall, precision, F1 and FP; replay rejects all 19 failed repair outputs
and accepts unchanged baseline controls. It was also applied to the raw-input run.
The original batch remains cancelled and the v1.5.2 draft remains on hold. The
residual scenarios, a separately observed loss of strand in inferred-output
conversion, and full cohort acceptance remain open; no claim of complete software
or scientific-result validation is made.

Before replacing historical outputs, validate the changed software, replay real
candidate cases, re-evaluate affected benchmarks, then rebuild the relevant
library outputs and old-to-new ID mappings. No claim that the corrected final
catalogue still contains the historical number of records is justified until
sequence/structure deduplication is rerun for those libraries. Scientific
figures, source data and manuscript fields must then be reconciled to that result.
