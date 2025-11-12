# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- Added a failsafe mechanism to deal with misaligned taxonomy (TaxID not found in `nodes.dmp`).
- Taxonomy nodes and names files are downloaded together with the rest of the database #36.
- `tidytable` works with relative path #30.
- `covplot` has been slightly slightly changed to fix #33 and #28. It's also faster.

### Changed
- `covplot` reports more information on sequence length and covered bases, and features the Economist theme.
- `update` now performs some tests on the integrity of `viral_database.fasta`.
- `orgs_list` now includes species name and coverage information #32.
- The viral database now includes all reference viral genomes #34.
- VirMet now saves its version and the information on the viral database into the output directory #26.
- `virmet update` now has some output even when no new sequences are added #19.

## [1.1.1] - 2017-05-22
### Fixed
- Failure to create output directory when `wolfpack` was run on a single file with slashes in the path.
- A bug in `tidytable` prevented to run if the directory was given with a trailing slash.
- `order` (a method deprecated in pandas) replaced with `sort_values`.

## [1.1] - 2017-04-26
### Added
- Subcommand `tidytable` creates two tables with summary information for samples/runs.

## [2.0.0] - 2025-11-12
### Changed
- Subcommand `fetch` is optimized and uses different software. The resulting database is also slightly different.
- Subcommand `upgrade` has been modified to match the modifications in `fetch`.
- Subcommand `index` has been modified to match the modifications in `fetch` and `upgrade`.
- Subcommand `wolfpack` has been changed substantially. It uses state-of-the-art software and is much faster.
- Subcommand `covplot` is now integrated into `wolfpack` but it can also be called separately.
- Subcommand `tidytable` is now integrated into `wolfpack` and is not a subcommand itself.