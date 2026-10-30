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
- The BLAS and LAPACK backend is `pure`, apart from the 26 LAPACK functions
  that take an `intent(out)` argument, which a pure function cannot have, and
  `?lacon`, which keeps its state in `save` variables and has no callers
  (LAPACK replaced it with `?lacn2`). The 26 functions keep their interface;
  every caller inside the library uses their pure subroutine twins instead.
  With `LA_EXTERNAL_LAPACK`, the external interfaces of the routines that
  became pure are declared `pure` as well.
- The high-level procedures that waited for a pure backend are `pure`: the
  subroutines `eig`, `eigh`, `svd`, `pseudoinvert`, `invert`, `schur`,
  `schur_space`, `solve_lstsq`, `solve_weighted_lstsq`,
  `solve_constrained_lstsq` and `constrained_lstsq_space`, the forms of
  `eigvals` and `eigvalsh` without `err`, and the operators `.inv.` and
  `.pinv.`. `mnorm` calls the pure `lange_sub`. A function that takes an `err`
  argument stays impure, since a pure function cannot.
- The generalized `eigvals(a, b)` takes `b` as `intent(in)`. It was declared
  `intent(inout)` but never modified, so a constant or an expression can now
  be passed.
  
### Added

- `pseudoinvert(a, pinva [, rtol] [, overwrite_a] [, err])`: with
  `overwrite_a=.true.` the SVD works on `a` itself instead of a copy, as in
  `svd`, `eig` and `eigh`. `a` is left unchanged by default, as it always was.
  `.pinv.` uses it on its own copy of the operand, so it copies once. A call
  that passed `err` by position now needs `err=`.
- A `pure` subroutine twin for each LAPACK function with an `intent(out)`
  argument: the matrix norms `?langb`, `?lange`, `?lanhb`, `?lanhe`, `?lanhf`,
  `?lanhp`, `?lanhs`, `?lansb`, `?lansf`, `?lansp`, `?lansy`, `?lantb`,
  `?lantp` and `?lantr`, and the extra-precise refinement helpers
  `?la_gbrcond`, `?la_gercond`, `?la_porcond`, `?la_syrcond`, their `_c`
  variants, `?la_hercond_c`, `?la_herpvgrw`, `?la_porpvgrw` and
  `?la_syrpvgrw`. `la_?lange_sub(norm, m, n, a, lda, work, res)` takes the
  arguments of `la_?lange` and returns its value in `res`. `la_lapack` exposes
  them as the generics `lange_sub` and so on, bound to the bundled code in
  every build, `LA_EXTERNAL_LAPACK` included.
- Linear-algebra test suite ported from the Fortran standard library
  (test-drive dev-dependency).
- `solve_lu(a, b, x [, pivot] [, overwrite_a] [, err])`, a `pure` subroutine form
  of `solve` that writes into a caller-provided solution array and can reuse a
  caller-provided pivot array, so a repeated solve of the same size allocates
  nothing.
- `operator(.det.)`: `.det. a` is the `pure` determinant of a square matrix, with
  no `overwrite_a` flag and no error argument.
- `invert(a, inva [, pivot] [, err])`, which writes the inverse into a second
  matrix instead of overwriting `a`, and an optional `pivot` argument on both
  forms of `invert`.
- `la_lapack` exposes seven more kind-agnostic generics, whose specifics were
  already there for all six kinds: `geqp3`, `ggev3`, `gges3`, `gghd3`, `geevx`,
  `geesx` and `gesvx`. As for every other generic, the quadruple-precision
  procedures stay outside the `LA_EXTERNAL_LAPACK` selection.
- `solve_chol(a, b, x [, lower] [, overwrite_a] [, err])` factorizes a symmetric
  or Hermitian positive definite matrix and solves in one call, and
  `solve_lower_chol(l, b, x [, err])` / `solve_upper_chol(u, b, x [, err])` solve
  from a Cholesky factor computed earlier. All three are `pure` and take one or
  many right-hand sides.
- `scripts/templatize.py --add-lapack-generic <name>...` builds the interface
  table entry of a generic from the declarations of its specifics and inserts it
  in alphabetical order, so `include/la_lapack_interfaces.fypp` never needs to be
  written by hand.
- Seven matrix property predicates in the new `la_matrix_property_checks`
  module: `is_square`, `is_diagonal`, `is_symmetric`, `is_skew_symmetric`,
  `is_hermitian`, `is_triangular(a, uplo)` and `is_hessenberg(a, uplo)`.
- `trace(a)`, `outer_product(u, v)`, `cross_product(a, b)`,
  `kronecker_product(a, b)` and `hermitian(a)`, all `pure`, for every real and
  complex kind.
- `diag` gains three overloads: `diag(v, k)` places a vector on the `k`-th
  diagonal of a larger matrix, and `diag(a)` / `diag(a, k)` extract a diagonal
  of a matrix as a vector.
- `lstsq_space(a, b, lrwork, liwork [, lcwork])`, the workspace query of the
  least-squares solver, which used to be a private routine.
- `solve_lstsq(a, b, x [, real_storage] [, int_storage] [, cmpl_storage]
  [, cond] [, singvals] [, overwrite_a] [, rank] [, err])`, the subroutine form
  of `lstsq`. Given the working arrays it allocates nothing, and it can return
  the singular values.
- `weighted_lstsq(w, a, b, ...)` and `solve_weighted_lstsq(w, a, b, x, ...)`,
  which minimize `||D(b - Ax)||` with `D = diag(sqrt(w))`.
- `constrained_lstsq(a, b, c, d, ...)` and
  `solve_constrained_lstsq(a, b, c, d, x, ...)`, which minimize `||b - Ax||`
  subject to `Cx = d` through LAPACK's `gglse`, with the matching
  `constrained_lstsq_space(a, c, lwork [, err])` workspace query.
- The optional precisions are fpm features. `quad` carries the 128-bit kinds
  `qp`/`w` and is part of the `default` profile, so a plain `fpm build` and
  every consumer that asks for nothing keep the kinds they had. `--profile lean`
  drops them, and with them a third of the library's source lines. `xdp` carries
  the 80-bit extended
  kinds `xdp`/`x`/`y`; it is never a default, because 80-bit reals exist on x86
  and x86_64 only. `--profile allkinds` turns both on. See the new
  "Precision kinds and fpm features" section of the README.
- `external-blas` and `external-lapack` are features too, and `--profile
  external` replaces the `--flag "-DLA_EXTERNAL_..."` pair the continuous
  integration used to spell out.
- 80-bit extended-precision instances of every BLAS and LAPACK routine and of
  the whole high-level API, initials `x` for real and `y` for complex, behind
  `LA_WITH_XDP`. `xdp` sits above `dp` as a branch of its own, so `qp` keeps
  `dp` below it and `la_qdgesv`, `la_qlag2d`, `la_dlag2q`, `la_wzgesv`,
  `la_zlag2w` and `la_wlag2z` keep their names and their meaning.
- `la_xddot`, the extended-precision accumulation dot product, is public from
  `la_blas_level1` but is not part of the `sdot` generic: its dummy arguments
  are `real(dp)`, exactly those of `la_qddot`, so a generic interface could not
  tell the two apart.
- `la_constants` exports `xdp` alongside `sp`, `dp` and `qp`, and the logical
  parameters `la_with_qp` and `la_with_xdp`. An optional kind that a build left
  out is `-1`, so code that imports the kind number still compiles.
- `scripts/guard_kinds.py` writes the cpp fences into the templates and has a
  `--check` mode the `generated-sources` job runs.
- `la_state` messages format array arguments of any rank from 2 to 14: the elements in
  array element order, then the shape, as in `[1 2 3 4] (2x2)`. Only scalars and vectors
  were formatted before.
- The matrix exponential, ported from the Fortran standard library, in the new
  `la_matrix_functions` module: `expm(a [, order])` is `pure`, `expm(a [, order], err)`
  returns the error state instead, and `matrix_exp(a [, order] [, work] [, iwork] [, err])` /
  `matrix_exp(a, e [, order] [, work] [, iwork] [, err])` are `pure` subroutines that
  overwrite `a` or write into `e`. Given `work` and `iwork`, sized with
  `expm_space(a, lwork, liwork)`, they allocate no work arrays. Every real and complex
  kind; `order` is the Padé order, 10 by default.
- Generalized least squares, ported from the Fortran standard library: minimize
  `(b - Ax)^H W^-1 (b - Ax)` for observations with correlated errors, whose covariance
  matrix `W` is symmetric (Hermitian) positive definite, through LAPACK's `ggglm`.
  `generalized_lstsq(w, a, b [, prefactored_w])` is `pure`, and
  `generalized_lstsq(w, a, b [, prefactored_w], err)` returns the error state instead.
  `solve_generalized_lstsq(w, a, b, x [, storage] [, prefactored_w] [, overwrite_a]
  [, overwrite_w] [, err])` is the `pure` subroutine form. Given `storage`, sized with
  `generalized_lstsq_space(a, lwork [, err])`, and both overwrite flags, it allocates
  nothing. With `prefactored_w=.true.`, `w` is any square root `B` of the covariance,
  `W = B B^H`, instead of the covariance itself. Every real and complex kind.

### Changed

- The `fypp` templates of the high-level modules regenerate their committed
  sources byte for byte again.
- `scripts/fypp_deploy.py` deploys every template under `fypp/src` and
  `fypp/test`. The exclusion list that held the ten high-level modules back is
  gone, so `--check` now covers the whole generated tree.
- Every generated source under `src/` is renamed from `.f90` to `.F90`: they
  all carry cpp directives now, so they follow the extension rule the tree
  already used for `la_blas.F90` and `la_lapack.F90`. Build files that list the
  sources by name need the new spelling; fpm consumers list nothing.
- `real(qp)` and `complex(qp)` code compiles only when `LA_WITH_QP` is defined,
  which the `quad` feature does. A build that defines neither macro carries
  `sp` and `dp` alone.

### Fixed

- `eig` and `eigvals` no longer leak memory on a real matrix: the real and
  imaginary parts of the eigenvalues were allocated as pointers and never
  freed, 2n values per call.
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
- `qr` now writes the whole orthogonal factor: `?orgqr`/`?ungqr` is asked for
  `size(q,2)` columns instead of `size(a,2)`, so columns `n+1:m` of a full `Q` are
  no longer left undefined. `qr_space` sizes its query for the same column count.
- `mnorm(a, order, dim)` handles any pair of dimensions of a rank-3 or higher
  array. The permuted copy is allocated before it is written and is built with the
  inverse permutation, so `dim(1)/=1` no longer writes through a null pointer and
  `dim=[1,k]` with `k>2` no longer collapses the wrong axes.
- **Behaviour change:** `lstsq` returns a solution of size `size(a,2)`, as its
  interface documents, instead of one the size of the right-hand side. The values
  are unchanged; code that read the leading `size(a,2)` entries keeps working, code
  that sized a receiving array from `b` has to be updated.
- `lstsq` frees its internal copy of the coefficient matrix, and no longer
  deallocates the caller's matrix when `overwrite_a` is set.
- `la_eye_s` is part of the `eye` generic interface, so `eye(m, mold=0.0_sp)`
  resolves.
- `mnorm` accepts `fro` as a spelling of the Frobenius (Euclidean) order.
- `test_la_cholesky` compiles, and `la_tests` runs it.
