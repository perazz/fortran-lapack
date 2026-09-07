# fortran-lapack
This package provides precision-agnostic, high-level linear algebra APIs for `real` and `complex` arguments in Modern Fortran. The APIs are similar to NumPy/SciPy operations, and leverage a Modern Fortran implementation of the [Reference-LAPACK](http://github.com/reference-LAPACK) library.

A full and standardized implementation of the present library has been integrated into the [Fortran Standard Library](http://stdlib.fortran-lang.org/), and as such, most users should seek to access the functionality from `stdlib`. The present library is kept in place for those who seek a compact implementation of it.

# Browse API

All procedures work with all types (`real`, `complex`) and kinds (32, 64, 128-bit floats).

## [chol](@ref la_cholesky::chol) - Cholesky factorization of a matrix (function).

### Syntax

`c = chol(a [, lower] [, other_zeroed])`

### Description

This function computes the Cholesky factorization of a real symmetric or complex Hermitian matrix \f$ A \f$:

\f[
A = L L^T = U^T U
\f]

where \f$ L \f$ is a lower triangular matrix and \f$ U \f$ is an upper triangular matrix. 
The function returns the factorized matrix as a new allocation, without modifying the input matrix.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [n,n] \f$, representing the symmetric/Hermitian input matrix.
- `lower` (optional): A logical flag indicating whether the lower (\f$ L \f$) or upper (\f$ U \f$) triangular factor should be computed. Defaults to `lower = .true.`.
- `other_zeroed` (optional): A logical flag determining whether the unused half of the returned matrix should be explicitly zeroed. Defaults to `other_zeroed = .true.`.

### Return value

- `c`: A `real` or `complex` matrix of size \f$ [n,n] \f$, containing the Cholesky factors. The returned matrix is triangular (upper or lower, as selected).

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the input matrix ihas invalid size.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if numerical instability prevents factorization.
- If error handling is not provided, exceptions will trigger an `error stop`.

### Notes

- The function is based on LAPACK's [POTRF](@ref la_lapack::potrf) routines.
- This function allocates a new matrix to store the factorization. For an in-place version, use [cholesky](@ref la_cholesky::cholesky).

## [cholesky](@ref la_cholesky::cholesky) - Cholesky factorization of a matrix (subroutine).

### Syntax

`call cholesky(a [, c] [, lower] [, other_zeroed])`

### Description

This subroutine computes the Cholesky factorization of a real symmetric or complex Hermitian matrix \f$ A \f$:

\f[
A = L L^T = U^T U
\f]

where \f$ L \f$ is a lower triangular matrix and \f$ U \f$ is an upper triangular matrix. The factorization is performed in-place, modifying the input matrix `a`, 
or on a pre-allocated matrix `c` with the same type and kind as `a`.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [n,n] \f$, representing the symmetric/Hermitian input matrix. if `c` is not provided, on return it contains the Cholesky factorization.
- `c` (optional): A matrix of size \f$ [n,n] \f$, of the same type and kind as `a`, containing the Cholesky factorization. If provided, `a` is unchanged.
- `lower` (optional): A logical flag indicating whether the lower (\f$ L \f$) or upper (\f$ U \f$) triangular factor should be computed. Defaults to `lower = .true.`.
- `other_zeroed` (optional): A logical flag determining whether the unused half of the matrix should be explicitly zeroed. Defaults to `other_zeroed = .true.`.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the input matrix is not positive definite.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if numerical instability prevents factorization.
- If error handling is not provided, exceptions will trigger an `error stop`.

### Notes

- The returned Cholesky factorization matrix is triangular (upper or lower, as selected).
- The subroutine is based on LAPACK's [POTRF](@ref la_lapack::potrf) routines.
- This subroutine modifies the input matrix in-place. For a version that returns a newly allocated matrix, use [chol](@ref la_cholesky::chol).

## [eig](@ref la_eig::eig) - Eigendecomposition of a square matrix.

### Syntax

`call eig(a [, b] [, lambda] [, right] [, left] [, overwrite_a] [, overwrite_b] [, err])`

### Description

This interface provides methods for computing the eigenvalues and eigenvectors of a real or complex matrix. 
It supports both standard and generalized eigenvalue problems, allowing for the decomposition of a matrix `A` alone or a pair of matrices `(A, B)` in the generalized case.

Given a square matrix \f$ A \f$, this routine computes its eigenvalues \f$ \lambda \f$ and, optionally, its right or left eigenvectors:

\f[
A v = \lambda v
\f]

where \f$ v \f$ represents an eigenvector corresponding to eigenvalue \f$ \lambda \f$.

In the generalized eigenvalue problem case, the routine solves:

\f[
A v = \lambda B v
\f]

The computation supports both `real` and `complex` matrices. If requested, eigenvectors are returned as additional output arguments. The function provides options to allow in-place modification of `A` and `B` for performance optimization.

**Note:** The solution is based on LAPACK's [GEEV](@ref la_lapack::geev) and [GGEV](@ref la_lapack::ggev) routines.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$[n,n]\f$, representing the input matrix to be decomposed.
- `b` (optional): A `real` or `complex` matrix of size \f$[n,n]\f$, same type and kind as `a`, representing the second matrix in the generalized eigenvalue problem.
- `lambda`: A `complex` or `real` array of length \f$ n \f$, containing the computed eigenvalues.
- `right` (optional): A `complex` matrix of size \f$[n,n]\f$ containing the right eigenvectors as columns.
- `left` (optional): A `complex` matrix of size \f$[n,n]\f$ containing the left eigenvectors as columns.
- `overwrite_a` (optional): A logical flag indicating whether `A` can be overwritten for performance optimization.
- `overwrite_b` (optional): A logical flag indicating whether `B` can be overwritten (only in the generalized case).
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable to handle errors. If not provided, execution will stop on errors.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrices have invalid/incompatible sizes.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the eigendecomposition fails.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- The computed eigenvectors are normalized.
- If computing real eigenvalues, an error is returned if eigenvalues have nonzero imaginary parts.
- This routine is based on LAPACK's [GEEV](@ref la_lapack::geev) and [GGEV](@ref la_lapack::ggev) routines.
- Overwriting `A` or `B` can improve performance but destroys the original matrix data.

## [eigh](@ref la_eig::eigh) - Eigendecomposition of a real symmetric or complex Hermitian matrix.

### Syntax

`call eigh(a, lambda [, vectors] [, upper_a] [, overwrite_a] [, err])`

### Description

This interface provides methods for computing the eigenvalues and optionally the eigenvectors of a real symmetric or complex Hermitian matrix.

Given a real symmetric or complex Hermitian matrix \f$ A \f$, this routine computes its eigenvalues \f$ \lambda \f$ and, optionally, its right eigenvectors:

\f[
A v = \lambda v
\f]

where \f$ v \f$ represents an eigenvector corresponding to eigenvalue \f$ \lambda \f$.

The computation supports both real and complex matrices, and the eigenvectors, if requested, are returned as orthonormal vectors.

**Note:** The solution is based on LAPACK's [SYEVD](@ref la_lapack::syevd) and [HEEVD](@ref la_lapack::heevd) routines.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$[n,n]\f$, representing the input matrix to be decomposed. The matrix is overwritten with the eigenvalues on output.
- `lambda`: A `real` array of length \f$ n \f$, with the same kind as `a`, containing the computed eigenvalues.
- `vectors` (optional): A matrix of size \f$[n,n]\f$, with the same type and kind as `a`, containing the right eigenvectors stored as columns.
- `overwrite_a` (optional): A logical flag indicating whether the matrix `A` can be overwritten for performance optimization.
- `upper_a` (optional): A logical flag indicating whether the upper half of matrix `A` should be used for computation (default is the lower half).
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable to handle errors. If not provided, execution will stop on errors.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrices have invalid/incompatible sizes.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the eigendecomposition fails.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- The computed eigenvectors are orthonormal.
- Eigenvalues are real for symmetric or Hermitian matrices.
- This routine is based on LAPACK's [SYEVD](@ref la_lapack::syevd) and [HEEVD](@ref la_lapack::heevd) routines.
- Overwriting the matrix `A` can improve performance but destroys the original matrix data.

## [eigvals](@ref la_eig::eigvals) - Eigenvalues of a square matrix (interface).

### Syntax

`lambda = eigvals(a [, b] [, err])`

### Description

This interface provides methods for computing the eigenvalues of a real or complex square matrix. 
It supports both standard and generalized eigenvalue problems.

- In the standard eigenvalue problem, the function computes the eigenvalues \f$\lambda\f$ of the matrix \f$A\f$ such that:

\f[
A v = \lambda v
\f]

where \f$v\f$ is the eigenvector corresponding to eigenvalue \f$\lambda\f$.

- In the generalized eigenvalue problem, the function solves:

\f[
A v = \lambda B v
\f]

where \f$A\f$ and \f$B\f$ are the input matrices and \f$\lambda\f$ is the eigenvalue.

The function returns an array of eigenvalues computed for the input matrix \f$A\f$, and optionally the matrix \f$B\f$ for the generalized case.

**Note:** The solution is based on LAPACK's [GEEV](@ref la_lapack::geev) and [GGEV](@ref la_lapack::ggev) routines.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$[n,n]\f$, representing the input matrix to be decomposed.
- `b` (optional, generalized case): A matrix of size \f$[n,n]\f$ and same type and kind as `a`, representing the second matrix in the generalized eigenvalue problem.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable to handle errors. If not provided, execution will stop on errors.

### Return value

- `lambda`: A `complex` array of eigenvalues, computed from the input matrix \f$A\f$ (and \f$B\f$ if in the generalized case).

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrices have invalid/incompatible sizes.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the eigendecomposition fails.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- The eigenvalues are returned in a complex array, even for real matrices.
- For the generalized eigenvalue problem, matrix \f$B\f$ must be provided and will be modified in-place.
- This routine is based on LAPACK's [GEEV](@ref la_lapack::geev) and [GGEV](@ref la_lapack::ggev) routines.

## [eigvalsh](@ref la_eig::eigvalsh) - Eigenvalues of a real symmetric or complex Hermitian matrix.

### Syntax

`lambda = eigvalsh(a [, upper_a] [, err])`

### Description

This interface provides methods for computing the eigenvalues of a real symmetric or complex Hermitian matrix. 
The function computes the eigenvalues of the matrix \f$A\f$, and returns them in an array. 
The user can specify whether to use the upper or lower half of the matrix for computation.

- The function solves the eigenvalue problem:

\f[
A v = \lambda v
\f]

where \f$v\f$ is the eigenvector corresponding to eigenvalue \f$\lambda\f$.

- The computation supports both real and complex matrices. Regardless, due to symmetry the eigenvalues are returned as an array of `real` values.

The user can specify whether to use the upper or lower half of the matrix \f$A\f$ for the computation (default: lower half).

**Note:** The solution is based on LAPACK's [SYEV](@ref la_lapack::syev) and [HEEV](@ref la_lapack::heev) routines.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$[n,n]\f$, representing the real symmetric or complex Hermitian matrix to be decomposed.
- `upper_a` (optional): A logical flag indicating whether to use the upper half (`.true.`) or the lower half (`.false.`) of \f$A\f$ for the computation. The default is lower.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable to handle errors. If not provided, execution will stop on errors.

### Return value

- `lambda`: A `real` array containing the computed eigenvalues of the matrix \f$A\f$.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrix has invalid size.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the eigendecomposition fails.
- If `err` is not provided, execution will stop on errors.

### Notes

- The eigenvalues are returned in a `real` array, with the same kind as the input matrix `a`.
- This routine is based on LAPACK's [SYEV](@ref la_lapack::syev) and [HEEV](@ref la_lapack::heev) routines.

## [solve](@ref la_solve::solve) - Solve a linear matrix equation or a linear system of equations.

### Syntax

`x = solve(a, b [, overwrite_a] [, err])`  

### Description

Solve linear systems - one (`b(:)`) or many (`b(:,:)`).  

### Arguments

- `a`: A `real` or `complex` coefficient matrix. If `overwrite_a=.true.`, it is destroyed by the call.
- `b`: A rank-1 (one system) or rank-2 (many systems) array of the same kind as `a`, containing the right-hand-side vector(s).
- `overwrite_a` (optional, default = `.false.`): If `.true.`, input matrix `a` will be used as temporary storage and overwritten, to avoid internal data allocation.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable. 

### Return value

For a full-rank matrix, returns an array value that represents the solution to the linear system of equations.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the matrix is singular to working precision.
- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrix and rhs vectors have invalid/incompatible sizes.
- If `err` is not present, exceptions trigger an `error stop`.

## [solve_lu](@ref la_solve::solve_lu) - Solve a linear system into a pre-allocated array.

### Syntax

`call solve_lu(a, b, x [, pivot] [, overwrite_a] [, err])`

### Description

Solve linear systems - one (`b(:)`) or many (`b(:,:)`) - writing the result into the caller's array `x` instead of returning a new one. Storage for the pivot indices may be supplied as well: when `x` and `pivot` are both provided and `overwrite_a=.true.`, the call performs no internal allocation, which makes it suited to a loop over many systems of the same size. The routine is `pure`.

### Arguments

- `a`: A `real` or `complex` coefficient matrix of size \f$ [n,n] \f$. If `overwrite_a=.true.`, it is destroyed by the call.
- `b`: A rank-1 (one system) or rank-2 (many systems) array of the same kind as `a`, containing the right-hand-side vector(s).
- `x`: An array of the same shape and kind as `b`. On output it holds the solution.
- `pivot` (optional): An `integer(ilp)` array of size `n` that receives the diagonal pivot indices of the LU factorization.
- `overwrite_a` (optional, default = `.false.`): If `.true.`, input matrix `a` will be used as temporary storage and overwritten, to avoid internal data allocation.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the matrix is singular to working precision.
- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if `a`, `b`, `x` or `pivot` have invalid/incompatible sizes.
- If `err` is not present, exceptions trigger an `error stop`.

### Notes

- This subroutine is based on LAPACK's LU decomposition solvers [GESV](@ref la_lapack::gesv).
- [solve](@ref la_solve::solve) is the function form; it allocates and returns the solution instead of writing into `x`.

## [solve_chol](@ref la_solve::solve_chol) - Solve a Hermitian positive definite system.

### Syntax

`call solve_chol(a, b, x [, lower] [, overwrite_a] [, err])`

### Description

Factorize a `real` symmetric or `complex` Hermitian positive definite matrix and solve \f$ A x = b \f$ in one call, for one (`b(:)`) or many (`b(:,:)`) right-hand sides. Only the triangle `lower` selects is read. The result is written into the caller's array `x`. The routine is `pure`.

### Arguments

- `a`: A `real` symmetric or `complex` Hermitian positive definite matrix of size \f$ [n,n] \f$. If `overwrite_a=.true.`, it is overwritten with its Cholesky factor.
- `b`: A rank-1 (one system) or rank-2 (many systems) array of the same kind as `a`, containing the right-hand-side vector(s).
- `x`: An array of the same shape and kind as `b`. On output it holds the solution.
- `lower` (optional, default = `.true.`): If `.true.`, the lower triangle of `a` is read and the factorization is \f$ A = L L^H \f$; otherwise the upper triangle is read and the factorization is \f$ A = U^H U \f$.
- `overwrite_a` (optional, default = `.false.`): If `.true.`, input matrix `a` will be used as temporary storage and overwritten, to avoid internal data allocation.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if `a` is not positive definite.
- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if `a`, `b` or `x` have invalid/incompatible sizes.
- If `err` is not present, exceptions trigger an `error stop`.

### Notes

- This subroutine is based on LAPACK's [POSV](@ref la_lapack::posv) drivers.
- To reuse a factorization across several right-hand sides, call [cholesky](@ref la_cholesky::cholesky) once and then [solve_lower_chol](@ref la_solve::solve_lower_chol) or [solve_upper_chol](@ref la_solve::solve_upper_chol).

## [solve_lower_chol](@ref la_solve::solve_lower_chol), [solve_upper_chol](@ref la_solve::solve_upper_chol) - Solve from a Cholesky factor.

### Syntax

`call solve_lower_chol(l, b, x [, err])`

`call solve_upper_chol(u, b, x [, err])`

### Description

Solve \f$ A x = b \f$ for one or many right-hand sides from a Cholesky factor computed earlier, without factorizing again. Each call costs two triangular solves. Both routines are `pure`.

### Arguments

- `l` / `u`: The lower or upper Cholesky factor of size \f$ [n,n] \f$, as returned by [cholesky](@ref la_cholesky::cholesky) with `lower=.true.` or `lower=.false.`.
- `b`: A rank-1 (one system) or rank-2 (many systems) array of the same kind as the factor, containing the right-hand-side vector(s).
- `x`: An array of the same shape and kind as `b`. On output it holds the solution.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the factor, `b` or `x` have invalid/incompatible sizes.
- If `err` is not present, exceptions trigger an `error stop`.

### Notes

- Both routines are based on LAPACK's [POTRS](@ref la_lapack::potrs) routines.
- The factor is taken as given: a matrix that is not a Cholesky factor produces a wrong answer, not an error.

## [lstsq](@ref la_least_squares::lstsq) - Compute a least squares solution to a system of linear equations.

### Syntax

`x = lstsq(a, b [, cond] [, overwrite_a] [, rank] [, err])`

### Description

Solves the least-squares problem for the system \f$ A \cdot x = b \f$, where \f$ A \f$ is a square matrix of size \f$ n \times n \f$ and \f$ b \f$ is either a vector of size \f$ n \f$ or a matrix of size \f$ n \times nrhs \f$. The function minimizes the 2-norm \f$ \|b - A \cdot x\| \f$ by solving for \f$ x \f$. 

The result \f$ x \f$ is returned as an allocatable array, and it is either a vector (for a single right-hand side) or a matrix (for multiple right-hand sides).

### Arguments

- `a`: A `real` matrix of size \f$ n \times n \f$ representing the coefficient matrix. If `overwrite_a = .true.`, the contents of `a` may be modified during the computation.
- `b`: A `real` vector or matrix representing the right-hand side. The size should be \f$ n \f$ (for a single right-hand side) or \f$ n \times nrhs \f$ (for multiple right-hand sides).
- `cond` (optional): A cutoff for rank evaluation. Singular values \f$ s(i) \f$ such that \f$ s(i) \leq \text{cond} \cdot \max(s) \f$ are considered zero. 
- `overwrite_a` (optional, default = `.false.`): If `.true.`, both `a` and `b` may be overwritten and destroyed during computation. 
- `rank` (optional): An integer variable that returns the rank of the matrix \f$ A \f$.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If `err` is not provided, the function will stop execution on error.

### Return value

Returns the solution array \f$ x \f$ with size \f$ n \f$ (for a single right-hand side) or \f$ n \times nrhs \f$ (for multiple right-hand sides).

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the matrix \f$ A \f$ is singular to working precision.
- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrix `a` and the right-hand side `b` have incompatible sizes.
- If `err` is not provided, the function stops execution on error.

### Notes

- This function relies on LAPACK's least-squares solvers, such as [GELSS](@ref la_lapack::gelss).
- If `overwrite_a` is enabled, the original contents of `a` and `b` may be lost.

## [det](@ref la_determinant::det) - Determinant of a scalar or rectangular matrix.

### Syntax

`d = det(a [, overwrite_a] [, err])`

### Description

This function computes the determinant of a square matrix \f$ A \f$. The matrix must be a real matrix of size \f$ [m, n] \f$, and the determinant is computed using an efficient factorization method (e.g., LU decomposition).

### Arguments

- `a`: A real matrix of size \f$ [m, n] \f$, representing the rectangular matrix for which the determinant is calculated. If `overwrite_a`, it is an `inout` argument and may be modified during computation.
- `overwrite_a` (optional, default = `.false.`): A logical flag that determines whether the input matrix `a` can be overwritten. If `.true.`, the matrix `a` may be destroyed and modified in place to save memory.
- `err` (optional): A state return flag of  [type(la_state)](@ref la_state_type::la_state). If an error occurs and `err` is not provided, the function will stop execution.

### Return value

The function returns a `real` scalar value representing the determinant of the input matrix \f$ A \f$, with the same kind as \f$ A \f$.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the matrix `a` is not square.
- If `err` is not provided, the function will stop execution on errors.

### Notes

- The determinant of the matrix is computed using the LAPACK [getrf](@ref la_lapack::getrf) backend.
- If `overwrite_a` is enabled, the input matrix `a` will be destroyed during the computation process.



## [operator(.det.)](@ref la_determinant::operator(.det.)) - Determinant of a square matrix.

### Syntax

```fortran
d = .det. A
```

### Description

This operator computes the determinant of a square real or complex matrix \f$ A \f$ from its LU factorization, in the same way as [det](@ref la_determinant::det). It is `pure`, so it can be used inside `pure` procedures and `do concurrent` blocks; it takes no `overwrite_a` flag and never modifies its operand, which is copied internally.

### Arguments

- `A`: A `real` or `complex` square matrix of size \f$ [n,n] \f$.

### Return value

A scalar of the same type and kind as `A`, holding its determinant.

### Errors

- Unlike [det](@ref la_determinant::det), this operator **does not provide explicit error handling**: it has no `err` argument, so a non-square or singular matrix triggers an `error stop`.

### Notes

- The determinant is computed through the LAPACK [getrf](@ref la_lapack::getrf) backend.
- If error handling is required, use [det](@ref la_determinant::det) with its `err` argument instead.

## [inv](@ref la_inverse::inv) - Inverse of a square matrix.

### Syntax

`inv_a = inv(a [, err])`

### Description

This function computes the inverse \f$ A^{-1} \f$ of a real or complex square matrix \f$ A \f$, provided that \f$ A \f$ is non-singular. 
The inverse of a matrix is defined as:

\f[
A A^{-1} = A^{-1} A = I
\f]

where \f$ I \f$ is the identity matrix of the same size as \f$ A \f$. The inverse exists only if \f$ A \f$ is square and has full rank (i.e., all its singular values are nonzero).

The computation is performed using LU decomposition.

### Arguments

- `a`: A `real` or `complex` square matrix of size \f$ [n,n] \f$, representing the matrix to be inverted.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

- `inv_a`: A `real` or `complex` square matrix of size \f$ [n,n] \f$, representing the inverse of `a`.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if `a` is singular or has invalid size.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This function computes the inverse using LAPACK's LU decomposition routine [GETRF](@ref la_lapack::getrf) followed by [GETRI](@ref la_lapack::getri).
- The inverse should be used with caution in numerical computations. For solving linear systems, using [solve](@ref la_solve::solve) is usually more stable and efficient than explicitly computing the inverse.

## [invert](@ref la_inverse::invert) - Matrix inversion (subroutine).

### Syntax

`call invert(a [, pivot] [, err])`

`call invert(a, inva [, pivot] [, err])`

### Description

This subroutine computes the inverse \\( A^{-1} \\) of a real or complex square matrix \\( A \\). The first form works **in-place**, modifying `a` directly; the second writes the inverse into a second matrix `inva` of the same shape and leaves `a` untouched. Both use the LU decomposition method via LAPACK's [GETRF](@ref la_lapack::getrf) and [GETRI](@ref la_lapack::getri) routines.

Given a square matrix \\( A \\), the LU decomposition factorizes it as:

\f[
A = P L U
\f]

where:
- \\( P \\) is a permutation matrix,
- \\( L \\) is a lower triangular matrix with unit diagonal,
- \\( U \\) is an upper triangular matrix.

The inverse is then obtained by solving \\( A X = I \\) using the LU factors.

### Arguments

- `a`: A `real` or `complex` square matrix of size \\( [n,n] \\). In the in-place form it is replaced with its inverse \\( A^{-1} \\) on output; in the split form it is read only.
- `inva` (split form only): A matrix of the same shape and kind as `a`, which receives the inverse \\( A^{-1} \\).
- `pivot` (optional): An `integer(ilp)` array of size at least `n` that receives the diagonal pivot indices of the LU factorization. Supplying it avoids the internal allocation of the pivot array.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the matrix is singular.
- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if `a` has invalid size, if `inva` does not match the shape of `a`, or if `pivot` is shorter than `n`.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- The in-place form modifies `a`. If the original matrix needs to be preserved, use the split form or [inv](@ref la_inverse::inv) instead.
- The determinant of `a` can be computed before inversion using [det](@ref la_determinant::det) to check for singularity.
- The computational complexity is \\( O(n^3) \\), making it expensive for large matrices.
- It is recommended to use matrix factorizations (e.g., LU or QR) for solving linear systems instead of computing the inverse explicitly, as it is numerically more stable and efficient.


## [operator(.inv.)](@ref la_inverse::operator(.inv.)) - Compute the inverse of a square matrix.

### Syntax

```fortran
invA = .inv. A
```

### Description

This operator computes the inverse \f$ A^{-1} \f$ of a square, non-singular real or complex matrix \f$ A \f$ using an LU decomposition. The inversion satisfies:

\f[
A A^{-1} = I
\f]

where \f$ I \f$ is the identity matrix of appropriate size.

This operator is functionally equivalent to [inv](@ref la_inverse::inv) but provides a more convenient syntax. It supports operator chaining, allowing multiple inversions within expressions:

### Arguments

- `A`: A `real` or `complex` square matrix of size \f$ [n,n] \f$, representing the input matrix to be inverted.

### Return value

- `invA`: A `real` or `complex` square matrix of size \f$ [n,n] \f$, and same kind as `A` representing its inverse.  
- If `A` is singular or the inversion fails, an **empty matrix** (size \f$ [0,0] \f$) is returned instead of raising an error.

### Errors

- Unlike [inv](@ref la_inverse::inv), this operator **does not provide explicit error handling**.
- If `A` is singular or an error occurs during inversion, the function **returns an empty matrix** (size \f$ [0,0] \f$) instead of raising an exception.
- The caller should check the size of the returned matrix to determine if inversion was successful.

### Notes

- This operator internally calls LAPACK's LU decomposition routine [GETRF](@ref la_lapack::getrf) followed by [GETRI](@ref la_lapack::getri).
- The chaining property allows for concise expressions but requires caution: if any intermediate inversion fails, subsequent operations may propagate errors due to empty matrix results.
- If strict error handling is required, use [inv](@ref la_inverse::inv) instead.

## [pinv](@ref la_pseudoinverse::pinv) - Moore-Penrose pseudo-inverse of a matrix (function).

### Syntax

`pinva = pinv(a [, rtol] [, err])`

### Description

This function computes the Moore-Penrose pseudo-inverse \f$ A^+ \f$ of a real or complex matrix \f$ A \f$ using Singular Value Decomposition (SVD). The pseudo-inverse provides a generalization of the inverse for non-square and singular matrices, making it useful for solving least-squares problems and underdetermined systems.

The computation is based on the singular value decomposition (SVD):

\f[
A = U \Sigma V^T
\f]

where \f$ U \f$ and \f$ V \f$ are orthogonal matrices, and \f$ \Sigma \f$ is a diagonal matrix containing the singular values. The pseudo-inverse is computed as:

\f[
A^+ = V \Sigma^+ U^T
\f]

where \f$ \Sigma^+ \f$ is obtained by inverting the nonzero singular values.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,n] \f$, representing the input matrix to be inverted.
- `rtol` (optional): A real scalar specifying the relative tolerance for singular value truncation. Singular values smaller than `rtol * max(singular_values(A))` are set to zero. If not provided, a default machine-precision-based tolerance is used.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

- `pinva`: A `real` or `complex` matrix of size \f$ [n,m] \f$, representing the pseudo-inverse of `a`.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the SVD decomposition fails or the input matrix has invalid dimensions.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if numerical instability prevents inversion.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This function computes the pseudo-inverse using LAPACK's SVD decomposition routine [`*GESVD`](@ref la_lapack::gesvd).
- The choice of `rtol` affects numerical stability and rank estimation: setting it too high may result in an inaccurate inverse, while setting it too low may amplify numerical noise.
- This function returns a newly allocated matrix. For an in-place version, use [`pseudoinvert`](@ref la_pseudoinverse::pseudoinvert).

## [operator(.pinv.)](@ref la_pseudoinverse::operator(.pinv.)) - Compute the Moore-Penrose pseudo-inverse of a matrix.

### Syntax

`pinva = .pinv. a`

### Description

This operator computes the Moore-Penrose pseudo-inverse \f$ A^+ \f$ of a real or complex matrix \f$ A \f$ using Singular Value Decomposition (SVD). The pseudo-inverse is useful for solving least-squares problems and handling singular or underdetermined systems.

Given the singular value decomposition (SVD):

\f[
A = U \Sigma V^T
\f]

the pseudo-inverse is computed as:

\f[
A^+ = V \Sigma^+ U^T
\f]

where \f$ \Sigma^+ \f$ is the inverse of the nonzero singular values in \f$ \Sigma \f$.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,n] \f$, representing the input matrix to be inverted.

### Return value

- `pinva`: A `real` or `complex` matrix of size \f$ [n,m] \f$, representing the pseudo-inverse of `a`.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the SVD decomposition fails or the input matrix has invalid dimensions.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if numerical instability prevents inversion.
- If an error occurs, execution will stop.

### Notes

- This operator internally calls [pinv](@ref la_pseudoinverse::pinv) and behaves identically.
- The pseudo-inverse is computed using LAPACK's SVD decomposition routine [GESVD](@ref la_lapack::gesvd).
- This operator is a convenient shorthand for calling the functional interface `pinv(a)`.

## [svd](@ref la_svd::svd) - Singular Value Decomposition (SVD) of a matrix

### Syntax

`call svd(a, s [, u] [, vt] [, overwrite_a] [, full_matrices] [, err])`

### Description

This subroutine computes the Singular Value Decomposition (SVD) of a matrix \f$ A \f$:

\f[
A = U \cdot S \cdot V^T
\f]

where:
- \f$ A \f$ is the input matrix of size \f$ [m,n] \f$.
- \f$ U \f$ is an orthogonal matrix of size \f$ [m,m] \f$ (or \f$ [m,k] \f$ for the reduced problem), containing the left singular vectors of \f$ A \f$.
- \f$ S \f$ is a diagonal matrix containing the singular values of size \f$ [k,k] \f$ with \f$ k = \min(m,n) \f$.
- \f$ V^T \f$ is an orthogonal matrix of size \f$ [n,n] \f$ (or \f$ [k,n] \f$ for the reduced problem), containing the right singular vectors of \f$ A^T \f$.

The singular values are returned in the array \f$ S \f$, and optionally, the matrices \f$ U \f$ and \f$ V^T \f$ are computed and returned.

### Arguments

- `a`: A `real` matrix of size \f$ [m,n] \f$ representing the input matrix \f$ A \f$. If `overwrite_a = .true.`, this matrix may be modified during computation. This is an `inout` argument.
- `s`: A `real` array of size \f$ k = \min(m,n) \f$, containing the singular values of \f$ A \f$. This is an output argument.
- `u`: An optional `real` matrix of the same type and kind as `a`, representing the left singular vectors of \f$ A \f$. This has shape \f$ [m,m] \f$ for the full problem or \f$ [m,k] \f$ for the reduced problem. This is an output argument.
- `vt`: An optional `real` matrix of the same type and kind as `a`, representing the right singular vectors of \f$ A^T \f$. This has shape \f$ [n,n] \f$ for the full problem or \f$ [k,n] \f$ for the reduced problem. This is an output argument.
- `overwrite_a`: (Optional, default = `.false.`) A logical flag indicating whether the input matrix `a` may be overwritten during computation. If `.true.`, `a` is overwritten to avoid additional memory allocation.
- `full_matrices`: (Optional, default = `.true.`) A logical flag that determines whether to compute full-sized matrices \f$ U \f$ and \f$ V^T \f$ (shape \f$ [m,m] \f$ and \f$ [n,n] \f$). If `.false.`, computes reduced matrices of shape \f$ [m,k] \f$ and \f$ [k,n] \f$.
- `err`: (Optional) A [type(la_state)](@ref la_state_type::la_state) variable to capture the error state. If not provided, the function will stop execution on error.

### Return value

The SVD of matrix \f$ A \f$ is returned in the corresponding output arguments, with the singular values in `s`, and optionally, the matrices \f$ U \f$ and \f$ V^T \f$.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the sizes of the matrices are incompatible with the full/reduced problem.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if there is insufficient storage space.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine computes the Singular Value Decomposition using LAPACK's [GESDD](@ref la_lapack::gesdd) algorithm.
- If `overwrite_a` is enabled, the input matrix `a` may be overwritten during computation.

## [svdvals](@ref la_svd::svdvals) - Singular Values Computation (function).

### Syntax

`s = svdvals(a [, err])`

### Description

This function computes the singular values of a `real` or `complex` matrix \f$ A \f$ and returns them in a vector \f$ s \f$, where \f$ s \f$ is an array of size \f$ k = \min(m, n) \f$. 
Singular values are non-negative values that provide important insights into the properties of the matrix, such as its rank and conditioning.

This function does not compute the full Singular Value Decomposition ([SVD](@ref la_svd::svd)); instead, it directly calculates and returns only the singular values of matrix \f$ A \f$.

The Singular Value Decomposition of a matrix \f$ A \f$ is expressed as:

\f[
A = U \cdot S \cdot V^T
\f]

where:
- \f$ A \f$ is the input matrix of size \f$ [m,n] \f$,
- \f$ U \f$ and \f$ V \f$ are orthogonal matrices,
- \f$ S \f$ is a diagonal matrix with the singular values of \f$ A \f$.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,n] \f$, representing the input matrix whose singular values are to be computed.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

- `s`: A `real` array containing the singular values of the matrix \f$ A \f$, with the same type and kind as the input matrix. The size of the array is \f$ k = \min(m, n) \f$.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the input matrix has invalid dimensions or if the SVD computation fails.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This function only computes the singular values and does not compute the full SVD (i.e., matrices \f$ U \f$ and \f$ V \f$ are not computed).
- The singular values are returned as a vector, sorted in decreasing order.
- The function uses LAPACK's [GESDD](@ref la_lapack::gesdd) routine for singular value computation.

## [diag](@ref la_eye::diag) - Diagonal matrix.

### Syntax

`d = diag(n, source [, err])` for scalar input
`d = diag(source(:) [, err])` for array input

### Description

This function generates a square diagonal matrix where the diagonal elements are populated either by a scalar value or an array of values. The size of the matrix is determined by the input parameter \f$n\f$ or the size of the input array. 
If a scalar is provided, the diagonal elements are all set to the same value. If an array is provided, its length determines the size of the matrix, and its elements are placed along the diagonal.

### Arguments

- `n`: The size of the square matrix (only used if a scalar is provided for the diagonal).
- `source`: 
  - If a scalar, this value is used to populate all the diagonal elements of the matrix.
  - If an array, the elements of the array are used to populate the diagonal of the matrix. The size of the array determines the matrix size.
- `err` (optional): A state return flag of [type(la_state)](@ref la_state_type::la_state). If an error occurs and `err` is not provided, the function will stop execution.

### Return value

The function returns a matrix of size \f$n \times n\f$, where the diagonal elements are either all equal to the scalar `source` or populated by the values from the input array.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the dimensions of the matrix are invalid or if the array size does not match the expected matrix size.
- If `err` is not provided, the function will stop execution on errors.

### Notes

- The diagonal elements are set to the specified scalar or the array values in the order they appear in the input.
- If the `err` parameter is provided, the error state of the function will be returned.


## [eye](@ref la_eye::eye) - Identity matrix.

### Syntax

`eye = eye(m [, n] [, mold] [, err])`

### Description

This function constructs an identity matrix of size \f$m \times n\f$, where the diagonal elements are set to 1 and all off-diagonal elements are set to 0. If only the number of rows \f$m\f$ is provided, a square matrix of size \f$m \times m\f$ is returned. The matrix is populated with a real data type, by default `real(real64)`, or a type specified by the user.

### Arguments

- `m`: The number of rows of the identity matrix.
- `n` (optional): The number of columns of the identity matrix. If omitted, the matrix is square (\f$m \times m\f$).
- `mold` (optional): The data type to define the return type. Defaults to `real(real64)`. 
- `err` (optional): A state return flag of [type(la_state)](@ref la_state_type::la_state). If an error occurs and `err` is not provided, the function will stop execution.

### Return value

The function returns a matrix of size \f$m \times n\f$ (or \f$m \times m\f$ if \f$n\f$ is omitted) with diagonal elements set to 1 and all other elements set to 0.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the dimensions of the matrix are invalid (e.g., negative values).
- If `err` is not provided, the function will stop execution on errors.

### Notes

- The identity matrix is constructed with the specified data type, which defaults to `real(real64)` if no type is specified.
- The `mold` scalar is used to provide a function return type. 
- If the `err` parameter is provided, the error state of the function will be returned.

## [qr](@ref la_qr::qr) - QR factorization of a matrix.

### Syntax

`call qr(a, q, r [, overwrite_a] [, storage] [, err])`

### Description

This subroutine computes the QR factorization of a `real` or `complex` matrix \f$ A = Q \cdot R \f$, where \f$ Q \f$ is orthonormal and \f$ R \f$ is upper-triangular. The matrix \f$ A \f$ has size \f$ [m,n] \f$ with \f$ m \ge n \f$. The result is returned in the output matrices \f$ Q \f$ and \f$ R \f$, which have the same type and kind as \f$ A \f$. 

Given \f$ k = \min(m, n) \f$, the matrix \f$ A \f$ can be written as:

\f[
A = \left( \begin{array}{cc} Q_1 & Q_2 \end{array} \right) \cdot \left( \begin{array}{cc} R_1 & 0 \end{array} \right)
\f]

Because the lower rows of \f$ R \f$ are zeros, a reduced problem \f$ A = Q_1 R_1 \f$ can be solved. The size of the input matrices determines which problem is solved:
- For full matrices (`shape(Q) == [m,m]`, `shape(R) == [m,n]`), the full problem is solved.
- For reduced matrices (`shape(Q) == [m,k]`, `shape(R) == [k,n]`), the reduced problem is solved.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,n] \f$, representing the coefficient matrix. If `overwrite_a = .false.`, this is an input argument. If `overwrite_a = .true.`, it is an `inout` argument and is overwritten upon return.
- `q`: A rank-2 array of the same type and kind as `a`, representing the orthonormal matrix \f$ Q \f$. This is an output argument with shape \f$ [m,m] \f$ (for the full problem) or \f$ [m,k] \f$ (for the reduced problem).
- `r`: A rank-2 array of the same type and kind as `a`, representing the upper-triangular matrix \f$ R \f$. This is an output argument with shape \f$ [m,n] \f$ (for the full problem) or \f$ [k,n] \f$ (for the reduced problem).
- `storage` (optional): A rank-1 array of the same type and kind as `a`, providing working storage for the solver. Its minimum size can be determined by a call to [qr_space](@ref la_qr::qr_space). This is an output argument.
- `overwrite_a` (optional, default = `.false.`): A logical flag that determines whether the input matrix `a` can be overwritten. If `.true.`, the matrix `a` is used as temporary storage and overwritten to avoid internal memory allocation. This is an input argument.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

The QR factorization matrices \f$ Q \f$ and \f$ R \f$ are returned in the corresponding arguments.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the sizes of the matrices are incompatible with the full/reduced problem.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if there is insufficient storage space.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine computes the QR factorization using LAPACK's QR decomposition algorithm [`*GEQRF`](@ref la_lapack::geqrf).
- If `overwrite_a` is enabled, the input matrix `a` will be modified during computation.


## [qr_space](@ref la_qr::qr_space) - Workspace size for QR operations.

### Syntax

`call qr_space(a, lwork [, err])`

### Description

This subroutine computes the minimum workspace size required for performing QR factorization. The size of the workspace array needed for both QR factorization and solving the reduced problem is determined based on the input matrix \f$ A \f$.

The input matrix \f$ A \f$ has size \f$ [m,n] \f$, and the output value \f$ lwork \f$ represents the minimum size of the workspace array that should be allocated for QR operations.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,n] \f$, representing the input matrix used to determine the required workspace size.
- `lwork`: An integer variable that will return the minimum workspace size required for QR factorization.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

The workspace size \f$ lwork \f$ that should be allocated before calling the QR factorization routine is returned.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if there is an issue determining the required workspace size.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine is useful for preallocating memory for QR factorization in large systems.
- It is important to ensure that the workspace size is correctly allocated before proceeding with QR factorization to avoid memory issues.

## [schur](@ref la_schur::schur) - Schur decomposition of a matrix.

### Syntax

`call schur(a, t, z [, eigvals] [, overwrite_a] [, storage] [, err])`

### Description

This subroutine computes the Schur decomposition of a `real` or `complex` matrix \f$ A = Z T Z^H \f$, where \f$ Z \f$ is an orthonormal/unitary matrix, and \f$ T \f$ is an upper-triangular or quasi-upper-triangular matrix. The matrix \f$ A \f$ has size \f$ [m,m] \f$. 

The decomposition produces:
- \f$ T \f$, which is upper-triangular for `complex` matrices and quasi-upper-triangular for `real` matrices (with possible \f$ 2 \times 2 \f$ blocks on the diagonal).
- \f$ Z \f$, the transformation matrix, which is optional.
- Optionally, the eigenvalues corresponding to the diagonal elements of \f$ T \f$.

If a pre-allocated workspace is provided, no internal memory allocations take place.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,m] \f$. If `overwrite_a = .false.`, this is an input argument. If `overwrite_a = .true.`, it is an `inout` argument and is overwritten upon return.
- `t`: A rank-2 array of the same type and kind as `a`, representing the Schur form of `a`. This is an output argument with shape \f$ [m,m] \f$.
- `z` (optional): A rank-2 array of the same type and kind as `a`, representing the unitary/orthonormal transformation matrix \f$ Z \f$. This is an output argument with shape \f$ [m,m] \f$.
- `eigvals` (optional): A complex array of size \f$ [m] \f$, representing the eigenvalues that appear on the diagonal of \f$ T \f$. This is an output argument.
- `storage` (optional): A rank-1 array of the same type and kind as `a`, providing working storage for the solver. Its minimum size can be determined by a call to [schur_space](@ref la_schur::schur_space). This is an input argument.
- `overwrite_a` (optional, default = `.false.`): A logical flag that determines whether the input matrix `a` can be overwritten. If `.true.`, the matrix `a` is used as temporary storage and overwritten to avoid internal memory allocation. This is an input argument.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

The Schur decomposition matrices \f$ T \f$ and optionally \f$ Z \f$ are returned in the corresponding arguments.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the sizes of the matrices are incompatible.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the algorithm did not converge.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine computes the Schur decomposition using LAPACK's Schur decomposition routines ([GEES](@ref la_lapack::gees)).
- Sorting options for eigenvalues can be requested, utilizing LAPACK's eigenvalue sorting mechanism.
- If `overwrite_a` is enabled, the input matrix `a` will be modified during computation.


## [schur_space](@ref la_schur::schur_space) - Workspace size for Schur decomposition.

### Syntax

`call schur_space(a, lwork [, err])`

### Description

This subroutine computes the minimum workspace size required for performing Schur decomposition. The size of the workspace array needed is determined based on the input matrix \f$ A \f$.

The input matrix \f$ A \f$ has size \f$ [m,m] \f$, and the output value \f$ lwork \f$ represents the minimum size of the workspace array that should be allocated for Schur decomposition operations.

### Arguments

- `a`: A `real` or `complex` matrix of size \f$ [m,m] \f$, representing the input matrix used to determine the required workspace size.
- `lwork`: An integer variable that will return the minimum workspace size required for Schur decomposition.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

The workspace size \f$ lwork \f$ that should be allocated before calling the Schur decomposition routine is returned.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if there is an issue determining the required workspace size.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine is useful for preallocating memory for Schur decomposition in large systems.
- It is important to ensure that the workspace size is correctly allocated before proceeding with Schur decomposition to avoid memory issues.



# BLAS, LAPACK
Modern Fortran modules with full explicit typing features are available as modules `la_blas` and `la_lapack`. 
The reference Fortran-77 library, forked from Release 3.10.1, was automatically processed and modernized.
The following refactorings are applied: 
- All datatypes and accuracy constants standardized into a module (`stdlib`-compatible names)
- Both libraries available for 32, 64 and 128-bit floats
- Free format, lower-case style
- `implicit none(type, external)` everywhere
- all `pure` procedures where possible
- `intent` added to all procedure arguments
- Removed `DO 10 .... 10 CONTINUE`, replaced with `do..end do` loops or labelled `loop_10: do ... cycle loop_10 ... end do loop_10` in case control statements are present
- BLAS split into ten kind-templated topic modules, LAPACK into 47, each holding every precision of the routines of one topic
- All procedures prefixed (with `stdlib_`, currently).
- F77-style `parameter`s removed, and numeric constants moved to the top of each module.
- Ambiguity in single vs. double precision constants (`0.0`, `0.d0`, `(1.0,0.0)`) removed
- preprocessor-based OpenMP directives retained.

Grouping every precision of a topic in one module hopefully allows for cross-procedural inlining which is otherwise impossible without link-time optimization.

# Building
An automated build is currently available via the [Fortran Package Manager](https://fpm.fortran-lang.org).
To add fortran-lapack to your project, simply add it as a dependency: 

```
[dependencies]
fortran-lapack = { git="https://github.com/perazz/fortran-lapack.git" }
```

`fortran-lapack` is compatible with the LAPACK API. If high-performance external BLAS/LAPACK libraries are available, it is sufficient to define macros

```
[dependencies]
fortran-lapack = { git="https://github.com/perazz/fortran-lapack.git", preprocess.cpp.macros=["LA_EXTERNAL_BLAS", "LA_EXTERNAL_LAPACK"] }
```

# Extension to external BLAS/LAPACK libraries

Generic interfaces to most BLAS/LAPACK functions are exposed to modules `la_blas` and `la_lapack`. These interfaces drop the initial letter to wrap a precision-agnostic version. For example, `axpy` is a precision-agnostic interface to `saxpy`, `daxpy`, `caxpy`, `zaxpy`, `qaxpy`, `waxpy`. 
The naming convention is: 

Type     | 32-bit | 64-bit | 128-bit
---      | ---    | ---    | --- 
real     | `s`    | `d`    | `q`
complex  | `c`    | `z`    | `w`

All public interfaces in `la_blas` and `la_lapack` allow seamless linking against external libraries via a simple pre-processor flag. 
When an external library is available, just define macros `LA_EXTERNAL_BLAS` and `LA_EXTERNAL_LAPACK`. The kind-agnostic interface
will just point to the external function. All such interfaces follow this template:  

```fortran  
interface axpy
#ifdef LA_EXTERNAL_BLAS
    ! Use external library
    pure subroutine saxpy(n, a, x, incx, y, incy)
      import :: ik, sp
      integer, parameter :: wp = sp
      integer(ik), intent(in) :: n
      real(wp), intent(in) :: a
      real(wp), intent(in) :: x(*)
      integer(ik), intent(in) :: incx
      real(wp), intent(inout) :: y(*)
      integer(ik), intent(in) :: incy
    end subroutine saxpy
#else
    ! Use internal implementation
    module procedure la_saxpy
#endif
end interface
```

# Regenerating sources

The Fortran under `src/` and `test/` is generated from the kind-templated [fypp](https://fypp.readthedocs.io) sources
under `fypp/src/` and `fypp/test/`, with the shared kind algebra in `include/`. The generated files are committed, so
building or installing the package never needs fypp. After editing a template, regenerate with

```bash
python3 scripts/fypp_deploy.py           # rewrite src/ and test/ from the templates
python3 scripts/fypp_deploy.py --check   # verify the committed tree matches the templates
```

`--check` is what continuous integration runs; it prints the templates it does not own yet and the reason for each.
`la_blas` and `la_lapack` are umbrella modules that re-export 57 topic modules, ten for BLAS and 47 for LAPACK; the
generic interfaces they publish are data tables under `include/`, regenerated with `python3 scripts/templatize.py
--blas-interfaces` and `python3 scripts/templatize.py --lapack-interfaces`.
Two further scripts support that layout:
`scripts/templatize.py` converts per-kind Fortran into one template per topic, driven by `scripts/la_modules.tsv`,
and `scripts/check_generated.py` compares the regenerated tree against a git reference routine by routine.

# Licensing

LAPACK is a freely-available software package. It is available from [netlib](https://www.netlib.org/lapack/) via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software packages (and has been). Credit for the library should be given to the [LAPACK authors](https://www.netlib.org/lapack/contributor-list.html).
The license used for the software is the [modified BSD license](https://www.netlib.org/lapack/LICENSE.txt).
According to the original license, we changed the name of the routines and commented the changes made to the original.

# Acknowledgments
Part of this work was supported by the [Sovereign Tech Fund](https://www.sovereigntechfund.de).
