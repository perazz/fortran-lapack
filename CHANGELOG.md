# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- The BLAS sources are generated from kind-templated fypp topic modules
  instead of six per-kind modules: `la_blas_aux` plus `la_blas_level1`,
  `la_blas_level2_{ban,gen,pac,sym,tri}` and `la_blas_level3_{gen,sym,tri}`.
  The `la_blas` generic interfaces are unchanged.
- Five BLAS specific names that kept the kind letter of the precision they
  were copied from are corrected: `la_qsdot`, `la_qzasum`, `la_qznrm2`,
  `la_wdrot` and `la_wdscal` become `la_qddot`, `la_qwasum`, `la_qwnrm2`,
  `la_wqrot` and `la_wqscal`. None of them is reachable through a generic
  interface.
- `scripts/fypp_deploy.py` replaces `scripts/preprocess.sh` for regenerating
  the committed sources, and has a `--check` mode.

### Fixed

- `la_dgesdd` no longer writes `path 5t` to standard output on one of its
  execution paths.
- `la_claunhr_col_getrfnp2` declares its `cabs1` statement function
  `real(sp)`, matching the precision of the data it works on.
- Removed `la_ilaqiag`, an unreachable misnamed copy of `la_iladiag`.
- Removed `la_qlag2q` and `la_wlag2w`, unreachable converters between a kind
  and itself.
- The external BLAS interfaces for `snrm2` and `dnrm2` declared no result
  type, so any build with `LA_EXTERNAL_BLAS` failed to compile.
- The test suite seeds its random-number generator, so its results are
  reproducible from one run to the next.
