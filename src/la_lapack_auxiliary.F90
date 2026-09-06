!> LAPACK auxiliary: machine parameters, safe division, band scaling
module la_lapack_auxiliary
     use la_constants
     use la_blas_aux
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slabad
     public :: la_sladiv2
     public :: la_slamch
     public :: la_slamc3
     public :: la_slaqsb
     public :: la_scsum1
     public :: la_sladiv1
     public :: la_dlabad
     public :: la_dladiv2
     public :: la_dlamch
     public :: la_dlamc3
     public :: la_dlaqsb
     public :: la_dzsum1
     public :: la_dladiv1
#ifdef LA_WITH_XDP
     public :: la_xlabad
     public :: la_xladiv2
     public :: la_xlamch
     public :: la_xlamc3
     public :: la_xlaqsb
     public :: la_xysum1
     public :: la_xladiv1
#endif
#ifdef LA_WITH_QP
     public :: la_qlabad
     public :: la_qladiv2
     public :: la_qlamch
     public :: la_qlamc3
     public :: la_qlaqsb
     public :: la_qwsum1
     public :: la_qladiv1
#endif
     public :: la_claqsb
     public :: la_crot
     public :: la_zlaqsb
     public :: la_zrot
#ifdef LA_WITH_XDP
     public :: la_ylaqsb
     public :: la_yrot
#endif
#ifdef LA_WITH_QP
     public :: la_wlaqsb
     public :: la_wrot
#endif

     contains

     !> SLABAD: takes as input the values computed by SLAMCH for underflow and
     !> overflow, and returns the square root of each of these values if the
     !> log of LARGE is sufficiently large.  This subroutine is intended to
     !> identify machines with a large exponent range, such as the Crays, and
     !> redefine the underflow and overflow limits to be the square roots of
     !> the values computed by SLAMCH.  This subroutine is needed because
     !> SLAMCH does not compensate for poor arithmetic in the upper half of
     !> the exponent range, as is found on a Cray.

     pure subroutine la_slabad(small,large)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(inout) :: large,small
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: log10,sqrt
           ! Executable Statements
           ! if it looks like we're on a cray, take the square root of
           ! small and large to avoid overflow and underflow problems.
           if (log10(large) > 2000.) then
              small = sqrt(small)
              large = sqrt(large)
           end if
           return
     end subroutine la_slabad
     !> DLABAD: takes as input the values computed by DLAMCH for underflow and
     !> overflow, and returns the square root of each of these values if the
     !> log of LARGE is sufficiently large.  This subroutine is intended to
     !> identify machines with a large exponent range, such as the Crays, and
     !> redefine the underflow and overflow limits to be the square roots of
     !> the values computed by DLAMCH.  This subroutine is needed because
     !> DLAMCH does not compensate for poor arithmetic in the upper half of
     !> the exponent range, as is found on a Cray.

     pure subroutine la_dlabad(small,large)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(inout) :: large,small
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: log10,sqrt
           ! Executable Statements
           ! if it looks like we're on a cray, take the square root of
           ! small and large to avoid overflow and underflow problems.
           if (log10(large) > 2000._dp) then
              small = sqrt(small)
              large = sqrt(large)
           end if
           return
     end subroutine la_dlabad
#ifdef LA_WITH_XDP
     !> XLABAD: takes as input the values computed by XLAMCH for underflow and
     !> overflow, and returns the square root of each of these values if the
     !> log of LARGE is sufficiently large.  This subroutine is intended to
     !> identify machines with a large exponent range, such as the Crays, and
     !> redefine the underflow and overflow limits to be the square roots of
     !> the values computed by XLAMCH.  This subroutine is needed because
     !> XLAMCH does not compensate for poor arithmetic in the upper half of
     !> the exponent range, as is found on a Cray.

     pure subroutine la_xlabad(small,large)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(inout) :: large,small
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: log10,sqrt
           ! Executable Statements
           ! if it looks like we're on a cray, take the square root of
           ! small and large to avoid overflow and underflow problems.
           if (log10(large) > 2000._xdp) then
              small = sqrt(small)
              large = sqrt(large)
           end if
           return
     end subroutine la_xlabad
#endif
#ifdef LA_WITH_QP
     !> QLABAD: takes as input the values computed by QLAMCH for underflow and
     !> overflow, and returns the square root of each of these values if the
     !> log of LARGE is sufficiently large.  This subroutine is intended to
     !> identify machines with a large exponent range, such as the Crays, and
     !> redefine the underflow and overflow limits to be the square roots of
     !> the values computed by QLAMCH.  This subroutine is needed because
     !> QLAMCH does not compensate for poor arithmetic in the upper half of
     !> the exponent range, as is found on a Cray.

     pure subroutine la_qlabad(small,large)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(inout) :: large,small
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: log10,sqrt
           ! Executable Statements
           ! if it looks like we're on a cray, take the square root of
           ! small and large to avoid overflow and underflow problems.
           if (log10(large) > 2000._qp) then
              small = sqrt(small)
              large = sqrt(large)
           end if
           return
     end subroutine la_qlabad
#endif

     pure real(sp) function la_sladiv2(a,b,c,d,r,t)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: a,b,c,d,r,t
        ! =====================================================================

           ! Local Scalars
           real(sp) :: br
           ! Executable Statements
           if (r /= zero) then
              br = b*r
              if (br /= zero) then
                 la_sladiv2 = (a + br)*t
              else
                 la_sladiv2 = a*t + (b*t)*r
              end if
           else
              la_sladiv2 = (a + d*(b/c))*t
           end if
           return
     end function la_sladiv2
     pure real(dp) function la_dladiv2(a,b,c,d,r,t)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: a,b,c,d,r,t
        ! =====================================================================

           ! Local Scalars
           real(dp) :: br
           ! Executable Statements
           if (r /= zero) then
              br = b*r
              if (br /= zero) then
                 la_dladiv2 = (a + br)*t
              else
                 la_dladiv2 = a*t + (b*t)*r
              end if
           else
              la_dladiv2 = (a + d*(b/c))*t
           end if
           return
     end function la_dladiv2
#ifdef LA_WITH_XDP
     pure real(xdp) function la_xladiv2(a,b,c,d,r,t)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: a,b,c,d,r,t
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: br
           ! Executable Statements
           if (r /= zero) then
              br = b*r
              if (br /= zero) then
                 la_xladiv2 = (a + br)*t
              else
                 la_xladiv2 = a*t + (b*t)*r
              end if
           else
              la_xladiv2 = (a + d*(b/c))*t
           end if
           return
     end function la_xladiv2
#endif
#ifdef LA_WITH_QP
     pure real(qp) function la_qladiv2(a,b,c,d,r,t)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: a,b,c,d,r,t
        ! =====================================================================

           ! Local Scalars
           real(qp) :: br
           ! Executable Statements
           if (r /= zero) then
              br = b*r
              if (br /= zero) then
                 la_qladiv2 = (a + br)*t
              else
                 la_qladiv2 = a*t + (b*t)*r
              end if
           else
              la_qladiv2 = (a + d*(b/c))*t
           end if
           return
     end function la_qladiv2
#endif

     !> SLAMCH: determines single precision machine parameters.

     pure real(sp) function la_slamch(cmach)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: cmach
       ! =====================================================================

           ! Local Scalars
           real(sp) :: rnd,eps,sfmin,small,rmach
           ! Intrinsic Functions
           intrinsic :: digits,epsilon,huge,maxexponent,minexponent,radix,tiny
           ! Executable Statements
           ! assume rounding, not chopping. always.
           rnd = one
           if (one == rnd) then
              eps = epsilon(zero)*0.5
           else
              eps = epsilon(zero)
           end if
           if (la_lsame(cmach,'E')) then
              rmach = eps
           else if (la_lsame(cmach,'S')) then
              sfmin = tiny(zero)
              small = one/huge(zero)
              if (small >= sfmin) then
                 ! use small plus a bit, to avoid the possibility of rounding
                 ! causing overflow when computing  1/sfmin.
                 sfmin = small*(one + eps)
              end if
              rmach = sfmin
           else if (la_lsame(cmach,'B')) then
              rmach = radix(zero)
           else if (la_lsame(cmach,'P')) then
              rmach = eps*radix(zero)
           else if (la_lsame(cmach,'N')) then
              rmach = digits(zero)
           else if (la_lsame(cmach,'R')) then
              rmach = rnd
           else if (la_lsame(cmach,'M')) then
              rmach = minexponent(zero)
           else if (la_lsame(cmach,'U')) then
              rmach = tiny(zero)
           else if (la_lsame(cmach,'L')) then
              rmach = maxexponent(zero)
           else if (la_lsame(cmach,'O')) then
              rmach = huge(zero)
           else
              rmach = zero
           end if
           la_slamch = rmach
           return
     end function la_slamch
     !> DLAMCH: determines double precision machine parameters.

     pure real(dp) function la_dlamch(cmach)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: cmach
       ! =====================================================================

           ! Local Scalars
           real(dp) :: rnd,eps,sfmin,small,rmach
           ! Intrinsic Functions
           intrinsic :: digits,epsilon,huge,maxexponent,minexponent,radix,tiny
           ! Executable Statements
           ! assume rounding, not chopping. always.
           rnd = one
           if (one == rnd) then
              eps = epsilon(zero)*0.5
           else
              eps = epsilon(zero)
           end if
           if (la_lsame(cmach,'E')) then
              rmach = eps
           else if (la_lsame(cmach,'S')) then
              sfmin = tiny(zero)
              small = one/huge(zero)
              if (small >= sfmin) then
                 ! use small plus a bit, to avoid the possibility of rounding
                 ! causing overflow when computing  1/sfmin.
                 sfmin = small*(one + eps)
              end if
              rmach = sfmin
           else if (la_lsame(cmach,'B')) then
              rmach = radix(zero)
           else if (la_lsame(cmach,'P')) then
              rmach = eps*radix(zero)
           else if (la_lsame(cmach,'N')) then
              rmach = digits(zero)
           else if (la_lsame(cmach,'R')) then
              rmach = rnd
           else if (la_lsame(cmach,'M')) then
              rmach = minexponent(zero)
           else if (la_lsame(cmach,'U')) then
              rmach = tiny(zero)
           else if (la_lsame(cmach,'L')) then
              rmach = maxexponent(zero)
           else if (la_lsame(cmach,'O')) then
              rmach = huge(zero)
           else
              rmach = zero
           end if
           la_dlamch = rmach
           return
     end function la_dlamch
#ifdef LA_WITH_XDP
     !> XLAMCH: determines extended precision machine parameters.

     pure real(xdp) function la_xlamch(cmach)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: cmach
       ! =====================================================================

           ! Local Scalars
           real(xdp) :: rnd,eps,sfmin,small,rmach
           ! Intrinsic Functions
           intrinsic :: digits,epsilon,huge,maxexponent,minexponent,radix,tiny
           ! Executable Statements
           ! assume rounding, not chopping. always.
           rnd = one
           if (one == rnd) then
              eps = epsilon(zero)*0.5
           else
              eps = epsilon(zero)
           end if
           if (la_lsame(cmach,'E')) then
              rmach = eps
           else if (la_lsame(cmach,'S')) then
              sfmin = tiny(zero)
              small = one/huge(zero)
              if (small >= sfmin) then
                 ! use small plus a bit, to avoid the possibility of rounding
                 ! causing overflow when computing  1/sfmin.
                 sfmin = small*(one + eps)
              end if
              rmach = sfmin
           else if (la_lsame(cmach,'B')) then
              rmach = radix(zero)
           else if (la_lsame(cmach,'P')) then
              rmach = eps*radix(zero)
           else if (la_lsame(cmach,'N')) then
              rmach = digits(zero)
           else if (la_lsame(cmach,'R')) then
              rmach = rnd
           else if (la_lsame(cmach,'M')) then
              rmach = minexponent(zero)
           else if (la_lsame(cmach,'U')) then
              rmach = tiny(zero)
           else if (la_lsame(cmach,'L')) then
              rmach = maxexponent(zero)
           else if (la_lsame(cmach,'O')) then
              rmach = huge(zero)
           else
              rmach = zero
           end if
           la_xlamch = rmach
           return
     end function la_xlamch
#endif
#ifdef LA_WITH_QP
     !> QLAMCH: determines quad precision machine parameters.

     pure real(qp) function la_qlamch(cmach)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: cmach
       ! =====================================================================

           ! Local Scalars
           real(qp) :: rnd,eps,sfmin,small,rmach
           ! Intrinsic Functions
           intrinsic :: digits,epsilon,huge,maxexponent,minexponent,radix,tiny
           ! Executable Statements
           ! assume rounding, not chopping. always.
           rnd = one
           if (one == rnd) then
              eps = epsilon(zero)*0.5
           else
              eps = epsilon(zero)
           end if
           if (la_lsame(cmach,'E')) then
              rmach = eps
           else if (la_lsame(cmach,'S')) then
              sfmin = tiny(zero)
              small = one/huge(zero)
              if (small >= sfmin) then
                 ! use small plus a bit, to avoid the possibility of rounding
                 ! causing overflow when computing  1/sfmin.
                 sfmin = small*(one + eps)
              end if
              rmach = sfmin
           else if (la_lsame(cmach,'B')) then
              rmach = radix(zero)
           else if (la_lsame(cmach,'P')) then
              rmach = eps*radix(zero)
           else if (la_lsame(cmach,'N')) then
              rmach = digits(zero)
           else if (la_lsame(cmach,'R')) then
              rmach = rnd
           else if (la_lsame(cmach,'M')) then
              rmach = minexponent(zero)
           else if (la_lsame(cmach,'U')) then
              rmach = tiny(zero)
           else if (la_lsame(cmach,'L')) then
              rmach = maxexponent(zero)
           else if (la_lsame(cmach,'O')) then
              rmach = huge(zero)
           else
              rmach = zero
           end if
           la_qlamch = rmach
           return
     end function la_qlamch
#endif

     pure real(sp) function la_slamc3(a,b)
        ! -- lapack auxiliary routine --
           ! univ. of tennessee, univ. of california berkeley and nag ltd..
           ! Scalar Arguments
           real(sp),intent(in) :: a,b
       ! =====================================================================
           ! Executable Statements
           la_slamc3 = a + b
           return
     end function la_slamc3
     pure real(dp) function la_dlamc3(a,b)
        ! -- lapack auxiliary routine --
           ! univ. of tennessee, univ. of california berkeley and nag ltd..
           ! Scalar Arguments
           real(dp),intent(in) :: a,b
       ! =====================================================================
           ! Executable Statements
           la_dlamc3 = a + b
           return
     end function la_dlamc3
#ifdef LA_WITH_XDP
     pure real(xdp) function la_xlamc3(a,b)
        ! -- lapack auxiliary routine --
           ! univ. of tennessee, univ. of california berkeley and nag ltd..
           ! Scalar Arguments
           real(xdp),intent(in) :: a,b
       ! =====================================================================
           ! Executable Statements
           la_xlamc3 = a + b
           return
     end function la_xlamc3
#endif
#ifdef LA_WITH_QP
     pure real(qp) function la_qlamc3(a,b)
        ! -- lapack auxiliary routine --
           ! univ. of tennessee, univ. of california berkeley and nag ltd..
           ! Scalar Arguments
           real(qp),intent(in) :: a,b
       ! =====================================================================
           ! Executable Statements
           la_qlamc3 = a + b
           return
     end function la_qlamc3
#endif

     !> SLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_slaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(sp),intent(in) :: amax,scond
           ! Array Arguments
           real(sp),intent(inout) :: ab(ldab,*)
           real(sp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: thresh = 0.1e+0_sp

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_slaqsb
     !> DLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_dlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(dp),intent(in) :: amax,scond
           ! Array Arguments
           real(dp),intent(inout) :: ab(ldab,*)
           real(dp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: thresh = 0.1e+0_dp

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_dlaqsb
#ifdef LA_WITH_XDP
     !> XLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_xlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(xdp),intent(in) :: amax,scond
           ! Array Arguments
           real(xdp),intent(inout) :: ab(ldab,*)
           real(xdp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: thresh = 0.1e+0_xdp

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_xlamch('SAFE MINIMUM')/la_xlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_xlaqsb
#endif
#ifdef LA_WITH_QP
     !> QLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_qlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(qp),intent(in) :: amax,scond
           ! Array Arguments
           real(qp),intent(inout) :: ab(ldab,*)
           real(qp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: thresh = 0.1e+0_qp

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_qlaqsb
#endif

     !> SCSUM1: takes the sum of the absolute values of a complex
     !> vector and returns a single precision result.
     !> Based on SCASUM from the Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure real(sp) function la_scsum1(n,cx,incx)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           real(sp) :: stemp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_scsum1 = zero
           stemp = zero
           if (n <= 0) return
           if (incx == 1) go to 20
           ! code for increment not equal to 1
           nincx = n*incx
           do i = 1,nincx,incx
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_scsum1 = stemp
           return
           ! code for increment equal to 1
           20 continue
           do i = 1,n
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_scsum1 = stemp
           return
     end function la_scsum1
     !> DZSUM1: takes the sum of the absolute values of a complex
     !> vector and returns a double precision result.
     !> Based on DZASUM from the Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure real(dp) function la_dzsum1(n,cx,incx)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           real(dp) :: stemp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_dzsum1 = zero
           stemp = zero
           if (n <= 0) return
           if (incx == 1) go to 20
           ! code for increment not equal to 1
           nincx = n*incx
           do i = 1,nincx,incx
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_dzsum1 = stemp
           return
           ! code for increment equal to 1
           20 continue
           do i = 1,n
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_dzsum1 = stemp
           return
     end function la_dzsum1
#ifdef LA_WITH_XDP
     !> XYSUM1: takes the sum of the absolute values of a complex
     !> vector and returns a extended precision result.
     !> Based on XYASUM from the Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure real(xdp) function la_xysum1(n,cx,incx)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(xdp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           real(xdp) :: stemp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_xysum1 = zero
           stemp = zero
           if (n <= 0) return
           if (incx == 1) go to 20
           ! code for increment not equal to 1
           nincx = n*incx
           do i = 1,nincx,incx
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_xysum1 = stemp
           return
           ! code for increment equal to 1
           20 continue
           do i = 1,n
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_xysum1 = stemp
           return
     end function la_xysum1
#endif
#ifdef LA_WITH_QP
     !> QWSUM1: takes the sum of the absolute values of a complex
     !> vector and returns a quad precision result.
     !> Based on QWASUM from the Level 1 BLAS.
     !> The change is to use the 'genuine' absolute value.

     pure real(qp) function la_qwsum1(n,cx,incx)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           real(qp) :: stemp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           la_qwsum1 = zero
           stemp = zero
           if (n <= 0) return
           if (incx == 1) go to 20
           ! code for increment not equal to 1
           nincx = n*incx
           do i = 1,nincx,incx
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_qwsum1 = stemp
           return
           ! code for increment equal to 1
           20 continue
           do i = 1,n
              ! next line modified.
              stemp = stemp + abs(cx(i))
           end do
           la_qwsum1 = stemp
           return
     end function la_qwsum1
#endif

     pure subroutine la_sladiv1(a,b,c,d,p,q)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(inout) :: a
           real(sp),intent(in) :: b,c,d
           real(sp),intent(out) :: p,q
        ! =====================================================================

           ! Local Scalars
           real(sp) :: r,t
           ! Executable Statements
           r = d/c
           t = one/(c + d*r)
           p = la_sladiv2(a,b,c,d,r,t)
           a = -a
           q = la_sladiv2(b,a,c,d,r,t)
           return
     end subroutine la_sladiv1
     pure subroutine la_dladiv1(a,b,c,d,p,q)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(inout) :: a
           real(dp),intent(in) :: b,c,d
           real(dp),intent(out) :: p,q
        ! =====================================================================

           ! Local Scalars
           real(dp) :: r,t
           ! Executable Statements
           r = d/c
           t = one/(c + d*r)
           p = la_dladiv2(a,b,c,d,r,t)
           a = -a
           q = la_dladiv2(b,a,c,d,r,t)
           return
     end subroutine la_dladiv1
#ifdef LA_WITH_XDP
     pure subroutine la_xladiv1(a,b,c,d,p,q)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(inout) :: a
           real(xdp),intent(in) :: b,c,d
           real(xdp),intent(out) :: p,q
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: r,t
           ! Executable Statements
           r = d/c
           t = one/(c + d*r)
           p = la_xladiv2(a,b,c,d,r,t)
           a = -a
           q = la_xladiv2(b,a,c,d,r,t)
           return
     end subroutine la_xladiv1
#endif
#ifdef LA_WITH_QP
     pure subroutine la_qladiv1(a,b,c,d,p,q)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(inout) :: a
           real(qp),intent(in) :: b,c,d
           real(qp),intent(out) :: p,q
        ! =====================================================================

           ! Local Scalars
           real(qp) :: r,t
           ! Executable Statements
           r = d/c
           t = one/(c + d*r)
           p = la_qladiv2(a,b,c,d,r,t)
           a = -a
           q = la_qladiv2(b,a,c,d,r,t)
           return
     end subroutine la_qladiv1
#endif

     !> CLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_claqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(sp),intent(in) :: amax,scond
           ! Array Arguments
           real(sp),intent(in) :: s(*)
           complex(sp),intent(inout) :: ab(ldab,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: thresh = 0.1e+0_sp

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_claqsb
     !> ZLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_zlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(dp),intent(in) :: amax,scond
           ! Array Arguments
           real(dp),intent(in) :: s(*)
           complex(dp),intent(inout) :: ab(ldab,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: thresh = 0.1e+0_dp

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_zlaqsb
#ifdef LA_WITH_XDP
     !> YLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_ylaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(xdp),intent(in) :: amax,scond
           ! Array Arguments
           real(xdp),intent(in) :: s(*)
           complex(xdp),intent(inout) :: ab(ldab,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: thresh = 0.1e+0_xdp

           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_xlamch('SAFE MINIMUM')/la_xlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_ylaqsb
#endif
#ifdef LA_WITH_QP
     !> WLAQSB: equilibrates a symmetric band matrix A using the scaling
     !> factors in the vector S.

     pure subroutine la_wlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: kd,ldab,n
           real(qp),intent(in) :: amax,scond
           ! Array Arguments
           real(qp),intent(in) :: s(*)
           complex(qp),intent(inout) :: ab(ldab,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: thresh = 0.1e+0_qp

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: cj,large,small
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored in band format.
                 do j = 1,n
                    cj = s(j)
                    do i = max(1,j - kd),j
                       ab(kd + 1 + i - j,j) = cj*s(i)*ab(kd + 1 + i - j,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,min(n,j + kd)
                       ab(1 + i - j,j) = cj*s(i)*ab(1 + i - j,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_wlaqsb
#endif

     !> CROT:   applies a plane rotation, where the cos (C) is real and the
     !> sin (S) is complex, and the vectors CX and CY are complex.

     pure subroutine la_crot(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(sp),intent(in) :: c
           complex(sp),intent(in) :: s
           ! Array Arguments
           complex(sp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(sp) :: stemp
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
           ix = 1
           iy = 1
           if (incx < 0) ix = (-n + 1)*incx + 1
           if (incy < 0) iy = (-n + 1)*incy + 1
           do i = 1,n
              stemp = c*cx(ix) + s*cy(iy)
              cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
              cx(ix) = stemp
              ix = ix + incx
              iy = iy + incy
           end do
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              stemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - conjg(s)*cx(i)
              cx(i) = stemp
           end do
           return
     end subroutine la_crot
     !> ZROT:   applies a plane rotation, where the cos (C) is real and the
     !> sin (S) is complex, and the vectors CX and CY are complex.

     pure subroutine la_zrot(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(dp),intent(in) :: c
           complex(dp),intent(in) :: s
           ! Array Arguments
           complex(dp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(dp) :: stemp
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
           ix = 1
           iy = 1
           if (incx < 0) ix = (-n + 1)*incx + 1
           if (incy < 0) iy = (-n + 1)*incy + 1
           do i = 1,n
              stemp = c*cx(ix) + s*cy(iy)
              cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
              cx(ix) = stemp
              ix = ix + incx
              iy = iy + incy
           end do
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              stemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - conjg(s)*cx(i)
              cx(i) = stemp
           end do
           return
     end subroutine la_zrot
#ifdef LA_WITH_XDP
     !> YROT:   applies a plane rotation, where the cos (C) is real and the
     !> sin (S) is complex, and the vectors CX and CY are complex.

     pure subroutine la_yrot(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(xdp),intent(in) :: c
           complex(xdp),intent(in) :: s
           ! Array Arguments
           complex(xdp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(xdp) :: stemp
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
           ix = 1
           iy = 1
           if (incx < 0) ix = (-n + 1)*incx + 1
           if (incy < 0) iy = (-n + 1)*incy + 1
           do i = 1,n
              stemp = c*cx(ix) + s*cy(iy)
              cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
              cx(ix) = stemp
              ix = ix + incx
              iy = iy + incy
           end do
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              stemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - conjg(s)*cx(i)
              cx(i) = stemp
           end do
           return
     end subroutine la_yrot
#endif
#ifdef LA_WITH_QP
     !> WROT:   applies a plane rotation, where the cos (C) is real and the
     !> sin (S) is complex, and the vectors CX and CY are complex.

     pure subroutine la_wrot(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(qp),intent(in) :: c
           complex(qp),intent(in) :: s
           ! Array Arguments
           complex(qp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(qp) :: stemp
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
           ix = 1
           iy = 1
           if (incx < 0) ix = (-n + 1)*incx + 1
           if (incy < 0) iy = (-n + 1)*incy + 1
           do i = 1,n
              stemp = c*cx(ix) + s*cy(iy)
              cy(iy) = c*cy(iy) - conjg(s)*cx(ix)
              cx(ix) = stemp
              ix = ix + incx
              iy = iy + incy
           end do
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              stemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - conjg(s)*cx(i)
              cx(i) = stemp
           end do
           return
     end subroutine la_wrot
#endif

end module la_lapack_auxiliary
