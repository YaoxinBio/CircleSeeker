# Cecc integrity fixes after the v1.5.1 audit

This development branch starts at tag `v1.5.1`, commit
`87ae7f09c11a4795532a7524db1958c0c11a5b2d` (2026-08-16). It is a software
repair branch, not the planned 1.0.0 release or a replacement for the archived
historical analysis.

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
* Cecc coordinate signatures preserve directed traversal. Early grouping also
  requires equal canonical sequence. CD-HIT clusters are split when their
  directed structures disagree. A second coordinate-tolerance merge requires
  equivalent directed cycles within the original 10 bp coordinate tolerance and
  compatible whole-circle sequences. Its first draft required exact equality;
  full raw-input validation exposed duplicate noisy consensus copies with
  identical directed coordinates but different canonical origins. The corrected
  gate permits at most 1% whole-circle edits, including insertions/deletions,
  across rotation and reverse complement. This retains the original CD-HIT
  99% similarity intent; the added whole-circle edit bound is an explicit new
  validation criterion, not an assertion that the two similarity formulas are
  identical. Different directions/orders or substantially different sequences
  still cannot merge. Unknown sequence/direction cannot establish equivalence.
  The added dependency is [Edlib](https://pypi.org/project/edlib/) (tested with
  1.3.9.post1), using bounded infix alignment against a doubled sequence and
  charging unaligned/extra target span to the edit budget.
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

The repair candidate passed 1,267 tests (one deselected) with 72.29% coverage,
exceeding the configured 50% threshold, and mypy passed all 59 source files.
Checks ran on `fat2`, explicitly selecting this checkout with `PYTHONPATH=src`.
The CLI reported `1.5.1+ceccfix.20260908`. Regression tests include the actual
pipeline packager and verify that sequence validation failures cannot publish
fallback results. Separately, 55 previously flagged real records passed the
candidate replay through Cecc construction, U/M/C processing, CD-HIT,
deduplication and standalone formatting. That replay is not a full raw-input
pipeline rerun or a final catalogue replacement.

Before replacing historical outputs, validate the changed software, replay real
candidate cases, re-evaluate affected benchmarks, then rebuild the relevant
library outputs and old-to-new ID mappings. No claim that the corrected final
catalogue still contains the historical number of records is justified until
sequence/structure deduplication is rerun for those libraries. Scientific
figures, source data and manuscript fields must then be reconciled to that result.
