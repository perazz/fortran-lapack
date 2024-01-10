# fortran-lapack
This package contains a Modern Fortran implementation of the [Reference-LAPACK](http://github.com/reference-LAPACK) library.
The reference Fortran-77 library is automatically downloaded from its master repository, and processed to create Modern Fortran modules with full explicit typing features. 
Function interfaces are unchanged from the original implementation, and allow future extension to handle its usage through external implementations.

# Building
An automated build is currently available via the [Fortran Package Manager](https://fpm.fortran-lang.org).
To add fortran-lapack to your project, simply add it as a dependency: 

```
[dependencies]
fortran-lapack = { git="https://github.com/perazz/fortran-lapack.git" }
```

# Acknowledgments
The development of this package is supported by the [Sovereign Tech Fund](https://www.sovereigntechfund.de).
