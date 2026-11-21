# Code provenance and maintenance

CircleSeeker is maintained at https://github.com/YaoxinBio/CircleSeeker.
Its original authors are Yaoxin Zhang and Leo Xinqi Yu. Original copyright,
contributor identities, commit hashes and license terms are retained.

Development up to v1.5.1 took place in a repository that is being retired. Its
releases - v1.0.0, v1.1.1, v1.1.2, v1.5.0 and v1.5.1 - remain historical releases
and their tags are preserved here; this repository continues that version
sequence and does not redefine those versions. Releases published to Bioconda
under those version numbers were built from those same tags.

v1.5.0 is the version benchmarked in the accompanying manuscript. v1.6.1 is the
first release made from this repository; it corrects two defects that made
deliverables vary between runs of identical code on identical input, so results
produced with v1.6.1 are reproducible where results from earlier versions were
not. See CHANGELOG.md.

## Version correspondence

| Identity | Role |
|---|---|
| v1.5.1 / `87ae7f09c11a4795532a7524db1958c0c11a5b2d` | Historical analysis baseline |
| `055110f531a9e08936f8c526da5c71950a6bbcee` | Rejected preliminary repair; its dataset batch and prior release draft are obsolete |
| `1.5.1+ceccfix.20260908` / `44c79b8eafdb1f60f949d69c9eeab12a69ce570b` | Frozen repair used for 59 completed runs and the first eight GS steps |
| `1.5.1+ceccfix.20260908` / `a662ad74e1cf984af0417017a36572e3a725c934` | Performance-only support aggregation change; GS resumes at step 9 from its completed CD-HIT checkpoint |
| `e68fd7b30f70cf02af34c88a6120f133b7f01a1c` | Earlier release preparation. The source at this commit identifies itself as 1.6.0. |
| v1.6.1 / `245048e2e64245430da7f6f63e45de067b272b18` | Published release, 2026-09-12. Compared with `a662ad7`, eleven implementation files changed. Their reviewed digests are pinned in the [tagged release check](https://github.com/YaoxinBio/CircleSeeker/blob/v1.6.1/scripts/check_release_metadata.py). |

The release comparison excludes only `src/circleseeker/__version__.py` when
checking algorithm-source equality. Rebuilding with the repair commit is not
retroactively recorded as running an already-published v1.6.1.
The preliminary and current repairs shared an internal version string, so that
string alone cannot identify an analysis. Use the full Git commit and recorded
source checksums. The superseded draft at `5b6a79f` used the rejected `055110f` implementation. Use the [published v1.6.1 release](https://github.com/YaoxinBio/CircleSeeker/releases/tag/v1.6.1) for the final source and comparison scope.

All 194 commits reachable from the captured historical branches, tags and four
pull-request heads were preserved locally. The three repair commits bring that
history to 197 commits before the initial release preparation. Two subsequent
repair commits (`2df5dc0` and `44c79b8`) are merged with that preparation history,
preserving their original identities. Historical PR heads and the
ONT development branch are retained as archival branches, not merged into the
release implementation. GitHub issues, review conversations and release asset
metadata are platform records and are not automatically recreated by Git copying.

The `a662ad7` performance commit is also merged with its original identity.
It changes only `ecc_dedup.py` in the computational implementation. The 59
completed runs are not retroactively relabelled as executions of this commit;
the GS run records separate upstream and downstream commits. Equivalence and
million-cluster timing checks are described in
[the performance report](docs/development/dedup-support-performance.md).

The distinct CircleSeeker-dev v2.1.3 development line is not the baseline for this
release or the dataset rebuild.

## Reproducibility

Record the exact software commit, reference data, configuration, input checksums
and dependency versions for each analysis. The installation environment in this
repository does not claim to recreate historical environments. See the
[release notes](docs/releases/v1.6.1.md) for the published comparison and the separately dated candidate checks.
