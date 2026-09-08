# Code provenance and maintenance

CircleSeeker is maintained at https://github.com/YaoxinBio/CircleSeeker.
Its original authors are Yaoxin Zhang and Leo Xinqi Yu. Original copyright,
contributor identities, commit hashes and license terms are retained.

The historical repository is https://github.com/leoxqy/CircleSeeker. Its existing
releases remain historical releases, including v1.0.0 and v1.5.1. The current
repository continues that version sequence; it does not redefine those versions.

## Version correspondence

| Identity | Role |
|---|---|
| v1.5.1 / `87ae7f09c11a4795532a7524db1958c0c11a5b2d` | Historical analysis baseline |
| `1.5.1+ceccfix.20260908` / `055110f531a9e08936f8c526da5c71950a6bbcee` | Frozen repair code used for dataset rebuilding |
| v1.5.2 preparation | Same algorithm source as the repair commit; updated version, packaging and documentation |

The release comparison excludes only `src/circleseeker/__version__.py` when
checking algorithm-source equality. Rebuilding with the repair commit is not
retroactively recorded as running an already-published v1.5.2.

All 194 commits reachable from the captured historical branches, tags and four
pull-request heads were preserved locally. The three repair commits bring that
history to 197 commits before release preparation. Historical PR heads and the
ONT development branch are retained as archival branches, not merged into the
release implementation. GitHub issues, review conversations and release asset
metadata are platform records and are not automatically recreated by Git copying.

The distinct CircleSeeker-dev v2.1.3 development line is not the baseline for this
release or the dataset rebuild.

## Reproducibility

Record the exact software commit, reference data, configuration, input checksums
and dependency versions for each analysis. The installation environment in this
repository does not claim to recreate historical environments. See the
[release notes](docs/releases/v1.5.2.md) for the completed checks and outstanding
dataset acceptance. A final version DOI is added only after archival succeeds.
