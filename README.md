# *** WORK IN PROGRESS ***

# fortran-lapack
This package contains a Modern Fortran implementation of the [Reference-LAPACK](http://github.com/reference-LAPACK) library.
The reference Fortran-77 library is automatically downloaded from its master repository, and processed to create Modern Fortran modules with full explicit typing features. 
Release 3.10.1 is currently targeted. Function interfaces are unchanged from the original implementation, and allow future extension to handle its usage through external implementations.
The following refactorings are applied: 
- All datatypes and accuracy constants standardized into a module (`stdlib`-compatible names)
- Free format, lower-case style
- `implicit none(type, external)` everywhere
- BLAS modularized into a single-file module
- LAPACK modularized into a single-file module
- All procedures prefixed (with `stdlib_`, currently).
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
