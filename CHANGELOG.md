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
- The LAPACK auxiliary and base routines are generated from kind-templated
  fypp topic modules instead of the six per-kind monoliths: `la_lapack_aux`
  and `la_lapack_auxiliary`, `la_lapack_blas_like_{base,l1,l2,l3,scalar,mnorm}`,
  `la_lapack_solve_aux`, `la_lapack_givens_jacobi_rot` and
  `la_lapack_householder_reflectors`. The `la_lapack` generic interfaces are
  unchanged.
- Six LAPACK specific names that kept the kind letter of the precision they
  were copied from are corrected: `la_qlag2s`, `la_qlat2s`, `la_qzsum1`,
  `la_wdrscl`, `la_wlag2c` and `la_wlat2c` become `la_qlag2d`, `la_qlat2d`,
  `la_qwsum1`, `la_wqrscl`, `la_wlag2z` and `la_wlat2z`. None of them is
  reachable through a generic interface.
- The LAPACK solve family is generated from kind-templated fypp topic modules
  instead of the six per-kind monoliths: `la_lapack_solve_tri_comp`,
  `la_lapack_solve_{chol,lu}_comp`, `la_lapack_solve_lu`,
  `la_lapack_solve_ldl_comp{,2,3,4}`, `la_lapack_solve_{chol,ldl}` and
  `la_lapack_others_sm`. The `la_lapack` generic interfaces are unchanged.
- Four mixed-precision LAPACK drivers whose name kept the kind letter of the
  precision they were copied from are corrected: `la_qsgesv`, `la_qsposv`,
  `la_wcgesv` and `la_wcposv` become `la_qdgesv`, `la_qdposv`, `la_wzgesv` and
  `la_wzposv`. None of them is reachable through a generic interface.
- The LAPACK orthogonal factorizations, cosine-sine decomposition and
  bidiagonal singular-value kernels are generated from kind-templated fypp
  topic modules instead of the six per-kind monoliths:
  `la_lapack_orthogonal_factors_{qr,rz,ql}`, `la_lapack_cosine_sine`,
  `la_lapack_lsq_constrained`, `la_lapack_svd_comp2`,
  `la_lapack_svd_bidiag_qr`, `la_lapack_eigv_gen_aux` and
  `la_lapack_eigv_gen_hess`. The `la_lapack` generic interfaces are unchanged.
- The LAPACK eigenvalue, singular-value and least-squares drivers are
  generated from kind-templated fypp topic modules instead of the six per-kind
  monoliths: `la_lapack_eigv_sym{,_comp}`, `la_lapack_eigv_tridiag{,2,3}`,
  `la_lapack_eigv_comp{,2}`, `la_lapack_eigv_gen{,2,3}`,
  `la_lapack_eigv_svd_{drivers,drivers2,bidiag_dc}`, `la_lapack_svd_comp`,
  `la_lapack_lsq` and `la_lapack_lsq_aux`. The `la_lapack` generic interfaces
  are unchanged.
- Removed the six per-kind LAPACK modules `la_lapack_{s,d,q,c,z,w}`. Every
  routine they held now lives in one of the 47 kind-templated topic modules,
  which the `la_lapack` umbrella imports directly.
- The `la_lapack` umbrella is generated from `include/la_lapack_interfaces.fypp`,
  a data table of its 498 generic interfaces and their 1524 external-library
  stubs, the way `la_blas` already was. `src/la_lapack.f90` and
  `src/la_lapack_aux.f90` are renamed `src/la_lapack.F90` and
  `src/la_lapack_aux.F90`: both carry cpp directives, so they now follow the
  extension rule the rest of the generated tree follows. Build files that list
  either source by name need the new spelling.
- The `generated-sources` continuous-integration job runs
  `scripts/fypp_deploy.py --check` and fails on any drift between a template
  and its committed output.

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
