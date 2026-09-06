# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
