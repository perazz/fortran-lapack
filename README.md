# *** WORK IN PROGRESS ***

# Current API
Procedure   | Type | Description | Optional arguments
---        | ---         | --- | ---
`solve(A,b)` | function | Solve linear systems - one (`b(:)`) or many (`b(:,:)`) | `solve(A,b,overwrite_a,err)`: option to let A be destroyed, return state handler `err`
`lstsq(A,b)` | function | Solve non-square systems in a least squares sense - one (`b(:)`) or many (`b(:,:)`) | `lstsq(A,b,cond,overwrite_a,rank,err)`: `cond` is optional SVD cutoff; `rank` to return matrix rank, `err` to return state handler
`det(A)` | function | Determinant of a scalar or square matrix | `det(A,overwrite_a,err=err)`: option to let A be destroyed, return state handler `err`
`inv(A)` | function | Inverse of a scalar or square matrix | `inv(A,err=err)`: A is not destroyed; return state handler `err`
`pinv(A)` | function | Moore-Penrose Pseudo-Inverse of a matrix | `pinv(A,rtol,err=err)`: A is not destroyed; optional singular value threshold `rtol`; return state handler `err`
`invert(A)` | subroutine | In-place inverse of a scalar or square matrix | `call invert(A,err=err)`: A is replaced with $A^{-1}$, return state handler `err`
`.inv.A` | operator | Inverse of a scalar or square matrix | A is replaced with $A^{-1}$
`.pinv.A` | operator | Moore-Penrose Pseudo-Inverse | A is replaced with $A^{-1}$
`svd(A)` | subroutine | Singular value decomposition of $A = U S V^t$ | `call svd(A,s,u,vt,full_matrices=.false.,err=state)`, all optional arguments but `A,s`
`svdvals(A)` | function | Singular values $S$ from $A = U S V^t$ | `s = svdvals(A)`, real array with same precision as `A`
`eye(m)` | function | Identity matrix of size `m` | `eye(m,n,mold,err)`: Optional column size `n`, datatype `dtype` (default: real64), error handler
`eigvals(A)` | function | Eigenvalues of matrix $A$ | `eigvals(A,err)`: Optional state handler `err`
`eig(A,lambda)` | subroutine | Eigenproblem of matrix $A$ | `eig(A,lambda,left,right,overwrite_a,err)`: optional output eigenvector matrices (left and/or right)
`eigvalsh(A)` | function | Eigenvalues of symmetric or hermitian matrix $A$ | `eigvalsh(A,upper_a,err)`: Choose to use upper or lower triangle; optional state handler `err`
`eigh(A,lambda)` | subroutine | Eigenproblem of symmetric or hermitianmatrix $A$ | `eigh(A,lambda,vector,upper_a,overwrite_a,err)`: optional output eigenvectors 
`diag(n,source)` | function | Diagonal matrix from scalar input value | `diag(n,source,err)`: Optional error handler
`diag(source)` | function | Diagonal matrix from array input values | `diag(source,err)`: Optional error handler
`qr(A,Q,R)` | subroutine | QR factorization | `qr(A,Q,R,storage=work,err=err)`: Optional pre-allocated working storage, error handler
`qr_space(A,lwork)` | subroutine | QR Working space size | `qr_space(A,lwork,err)`: Optional error handler

All procedures work with all types (`real`, `complex`) and kinds (32, 64, 128-bit floats).

# fortran-lapack
This package contains a Modern Fortran implementation of the [Reference-LAPACK](http://github.com/reference-LAPACK) library.
The reference Fortran-77 library is automatically downloaded from its master repository, and processed to create Modern Fortran modules with full explicit typing features. 
Release 3.10.1 is currently targeted. Function interfaces are unchanged from the original implementation, and allow future extension to handle its usage through external implementations.
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
# Extension to external BLAS/LAPACK libraries

This task is in progress. The names of all procedures have been prefixed not to pollute the original BLAS/LAPACK namespace, so that handling of external libraries can be accomplished via a preprocessor flag. For example:

```fortran  
#ifdef EXTERNAL_BLAS
interface 
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
end interface
#else
interface saxpy
    module procedure stdlib_saxpy
end interface
#endif
```

# Licensing

LAPACK is a freely-available software package. It is available from [netlib](https://www.netlib.org/lapack/) via anonymous ftp and the World Wide Web. Thus, it can be included in commercial software packages (and has been). Credit for the library should be given to the [LAPACK authors](https://www.netlib.org/lapack/contributor-list.html).
The license used for the software is the [modified BSD license](https://www.netlib.org/lapack/LICENSE.txt).
According to the original license, we changed the name of the routines and commented the changes made to the original.

# Acknowledgments
The development of this package is supported by the [Sovereign Tech Fund](https://www.sovereigntechfund.de).
