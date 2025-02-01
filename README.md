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

## [invert](@ref la_inverse::invert) - In-place matrix inversion

### Syntax

`call invert(a [, err])`

### Description

This subroutine computes the inverse \\( A^{-1} \\) of a real or complex square matrix \\( A \\) **in-place**, modifying `a` directly. It uses the LU decomposition method via LAPACK's [GETRF](@ref la_lapack::getrf) and [GETRI](@ref la_lapack::getri) routines.

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

- `a`: A `real` or `complex` square matrix of size \\( [n,n] \\). On output, it is replaced with its inverse \\( A^{-1} \\).
- `err` (optional): A [type(la_state)](@ref la_state_type::la_state) variable that returns the error state. If not provided, the function will stop execution on error.

### Errors

- Raises [LINALG_ERROR](@ref la_state_type::linalg_error) if the matrix is singular.
- Raises [LINALG_VALUE_ERROR](@ref la_state_type::linalg_value_error) if `a` has invalid size.
- If `err` is not provided, exceptions will trigger an `error stop`.

### Notes

- This subroutine modifies `a` in-place. If the original matrix needs to be preserved, use [inv](@ref la_inverse::inv) instead.
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

- This subroutine computes the Schur decomposition using LAPACK's Schur decomposition routines ([GEES](@ref la_lapack:gees)).
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
