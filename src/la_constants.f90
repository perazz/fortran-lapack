!> Supported kind parameters
module la_constants
     use iso_fortran_env,only:real32,real64,real128,int32,int64
     use,intrinsic :: ieee_arithmetic,only:ieee_is_nan
#if defined(_OPENMP)
     use omp_lib
#endif
     implicit none(type,external)
     public

     !> Single-precision floats
     integer,parameter :: sp = real32
     
     !> Double-precision floats
     integer,parameter :: dp = real64
     
     !> Quadruple-precision floats
     integer,parameter :: qp = real128
     
     !> Internal logical kind
     integer,parameter :: lk = kind(.true.)
     
     !> 32-bit integer size type 
     integer,parameter :: ilp = int32
     
     private :: int32,int64

end module la_constants
