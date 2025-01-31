# fortran-lapack
This package provides precision-agnostic, high-level linear algebra APIs for `real` and `complex` arguments in Modern Fortran. The APIs are similar to NumPy/SciPy operations, and leverage a Modern Fortran implementation of the [Reference-LAPACK](http://github.com/reference-LAPACK) library.

A full and standardized implementation of the present library has been integrated into the [Fortran Standard Library](http://stdlib.fortran-lang.org/), and as such, most users should seek to access the functionality from `stdlib`. The present library is kept in place for those who seek a compact implementation of it.

# Browse API

All procedures work with all types (`real`, `complex`) and kinds (32, 64, 128-bit floats).

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

## `inv(A)`
**Type**: Function  
**Description**: Inverse of a scalar or square matrix.  
**Optional arguments**:  
- `err`: Return state handler.

## `invert(A)`
**Type**: Subroutine  
**Description**: In-place inverse of a scalar or square matrix.  
**Optional arguments**:  
- `err`: Return state handler.  

**Usage**: `call invert(A, err=err)` where `A` is replaced with $A^{-1}$.

## `.inv.A`
**Type**: Operator  
**Description**: Inverse of a scalar or square matrix.  

**Effect**: `A` is replaced with $A^{-1}$.

## [pseudoinvert](@ref la_pseudoinverse::pseudoinvert) - Moore-Penrose pseudo-inverse of a matrix.

### Syntax

`call pseudoinvert(a, pinva [, rtol] [, err])`

### Description

This subroutine computes the Moore-Penrose pseudo-inverse \f$ A^+ \f$ of a real or complex matrix \f$ A \f$ using Singular Value Decomposition (SVD). The pseudo-inverse is a generalization of the matrix inverse that can be computed for non-square and singular matrices. It is particularly useful for solving least-squares problems and underdetermined systems.

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
- `pinva`: A `real` or `complex` matrix of size \f$ [n,m] \f$, representing the pseudo-inverse of `a`. This is an output argument.
- `rtol` (optional): A real scalar specifying the relative tolerance for singular value truncation. Singular values smaller than `rtol * max(singular_values(A))` are set to zero. If not provided, a default machine-precision-based tolerance is used.
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Return value

The pseudo-inverse of the input matrix `a` is returned in `pinva`.

### Errors

- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if the dimensions of `pinva` do not match the expected output size.
- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the SVD decomposition fails.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine computes the pseudo-inverse using LAPACK’s SVD decomposition routine [`*GESVD`](@ref la_lapack::gesvd).
- The choice of `rtol` affects numerical stability and rank estimation: setting it too high may result in an inaccurate inverse, while setting it too low may amplify numerical noise.
- This version requires `pinva` to be pre-allocated. To obtain the required size before allocation, use [`pseudoinvert_space`](@ref la_pseudoinvert::pseudoinvert_space).


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

- This function computes the pseudo-inverse using LAPACK’s SVD decomposition routine [`*GESVD`](@ref la_lapack::gesvd).
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
- The pseudo-inverse is computed using LAPACK’s SVD decomposition routine [GESVD](@ref la_lapack::gesvd).
- This operator is a convenient shorthand for calling the functional interface `pinv(a)`.


## `svd(A)`
**Type**: Subroutine  
**Description**: Singular value decomposition of $A = U S V^t$.  
**Optional arguments**:  
- `s`: Singular values.  
- `u`: Left singular vectors.  
- `vt`: Right singular vectors.  
- `full_matrices`: Defaults to `.false.`.  
- `err`: State handler.  

**Usage**: `call svd(A, s, u, vt, full_matrices=.false., err=state)`.

## `svdvals(A)`
**Type**: Function  
**Description**: Singular values $S$ from $A = U S V^t$.  
**Usage**: `s = svdvals(A)` where `s` is a real array with the same precision as `A`.

## `eye(m)`
**Type**: Function  
**Description**: Identity matrix of size `m`.  
**Optional arguments**:  
- `n`: Optional column size.  
- `mold`: Optional datatype (default: real64).  
- `err`: Error handler.

## `eigvals(A)`
**Type**: Function  
**Description**: Eigenvalues of matrix $A$.  
**Optional arguments**:  
- `err`: State handler.

## `eig(A, lambda)`
**Type**: Subroutine  
**Description**: Eigenproblem of matrix $A`.  
**Optional arguments**:  
- `left`: Output left eigenvector matrix.  
- `right`: Output right eigenvector matrix.  
- `overwrite_a`: Option to let A be destroyed.  
- `err`: Return state handler.

## `eigvalsh(A)`
**Type**: Function  
**Description**: Eigenvalues of symmetric or Hermitian matrix $A$.  
**Optional arguments**:  
- `upper_a`: Choose to use upper or lower triangle.  
- `err`: State handler.

## `eigh(A, lambda)`
**Type**: Subroutine  
**Description**: Eigenproblem of symmetric or Hermitian matrix $A`.  
**Optional arguments**:  
- `vector`: Output eigenvectors.  
- `upper_a`: Choose to use upper or lower triangle.  
- `overwrite_a`: Option to let A be destroyed.  
- `err`: Return state handler.

## `diag(n, source)`
**Type**: Function  
**Description**: Diagonal matrix from scalar input value.  
**Optional arguments**:  
- `err`: Error handler.

## `diag(source)`
**Type**: Function  
**Description**: Diagonal matrix from array input values.  
**Optional arguments**:  
- `err`: Error handler.

## [qr](@ref la_qr::qr) - Compute the QR factorization of a matrix.

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
- BLAS modularized into a single-file module
- LAPACK modularized into a single-file module
- All procedures prefixed (with `stdlib_`, currently).
- F77-style `parameter`s removed, and numeric constants moved to the top of each module.
- Ambiguity in single vs. double precision constants (`0.0`, `0.d0`, `(1.0,0.0)`) removed
- preprocessor-based OpenMP directives retained.

The single-source module structure hopefully allows for cross-procedural inlining which is otherwise impossible without link-time optimization.

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

# Licensing

LAPACK is a freely-available software package. It is available from [netlib](https://www.netlib.org/lapack/) via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software packages (and has been). Credit for the library should be given to the [LAPACK authors](https://www.netlib.org/lapack/contributor-list.html).
The license used for the software is the [modified BSD license](https://www.netlib.org/lapack/LICENSE.txt).
According to the original license, we changed the name of the routines and commented the changes made to the original.

# Acknowledgments
Part of this work was supported by the [Sovereign Tech Fund](https://www.sovereigntechfund.de).
