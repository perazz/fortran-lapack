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
     
     !> 80-bit extended-precision floats; -1 when the build did not ask for them
#ifdef LA_WITH_XDP
     integer,parameter :: xdp = selected_real_kind(18)
#else
     integer,parameter :: xdp = -1
#endif

     !> Quadruple-precision floats; -1 when the build did not ask for them
#ifdef LA_WITH_QP
     integer,parameter :: qp = real128
#else
     integer,parameter :: qp = -1
#endif

     !> Whether the extended- and quadruple-precision procedures are part of this build
#ifdef LA_WITH_XDP
     logical,parameter :: la_with_xdp = .true.
#else
     logical,parameter :: la_with_xdp = .false.
#endif
#ifdef LA_WITH_QP
     logical,parameter :: la_with_qp = .true.
#else
     logical,parameter :: la_with_qp = .false.
#endif

     !> Internal logical kind
     integer,parameter :: lk = kind(.true.)
     
     !> 32-bit integer size type
     integer,parameter :: ilp = int32
     
#ifdef LA_WITH_XDP
     integer,parameter,private :: xdp_is_a_kind_of_its_own = &
        1/merge(1,0,selected_real_kind(18) /= selected_real_kind(33)) ! LA_WITH_XDP needs a target with 80-bit reals
#endif
     private :: int32,int64

end module la_constants
