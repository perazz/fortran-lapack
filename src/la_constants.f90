module la_constants
     use iso_fortran_env,only:real32,real64,real128,int32,int64
     use,intrinsic :: ieee_arithmetic,only:ieee_is_nan
#if defined(_OPENMP)
     use omp_lib
#endif
     implicit none(type,external)
     public

     integer,parameter :: sp = real32
     integer,parameter :: dp = real64
     integer,parameter :: qp = real128
     integer,parameter :: lk = kind(.true.)
     ! Integer size support for ILP64 builds should be done here
     integer,parameter :: ilp = int32
     private :: int32,int64

end module la_constants
