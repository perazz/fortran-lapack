module stdlib_linalg_constants
     use iso_fortran_env,only:int32,int64
     use,intrinsic :: ieee_arithmetic,only:ieee_is_nan
#if defined(_OPENMP)
     use omp_lib
#endif
     implicit none(type,external)
     public

     integer,parameter :: sp = selected_real_kind(6)
     integer,parameter :: dp = selected_real_kind(15)
     integer,parameter :: qp = selected_real_kind(33)
     integer,parameter :: lk = kind(.true.)
     ! Integer size support for ILP64 builds should be done here
     integer,parameter :: ilp = int32
     private :: int32,int64

end module stdlib_linalg_constants
