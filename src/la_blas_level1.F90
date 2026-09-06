!> BLAS level 1: vector operations
module la_blas_level1
     use la_constants
     use la_blas_aux
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sasum
     public :: la_saxpy
     public :: la_scopy
     public :: la_sdot
     public :: la_snrm2
     public :: la_srot
     public :: la_srotg
     public :: la_srotm
     public :: la_srotmg
     public :: la_sscal
     public :: la_sswap
     public :: la_scasum
     public :: la_scnrm2
     public :: la_dasum
     public :: la_daxpy
     public :: la_dcopy
     public :: la_ddot
     public :: la_dnrm2
     public :: la_drot
     public :: la_drotg
     public :: la_drotm
     public :: la_drotmg
     public :: la_dscal
     public :: la_dsdot
     public :: la_dswap
     public :: la_dzasum
     public :: la_dznrm2
#ifdef LA_WITH_XDP
     public :: la_xasum
     public :: la_xaxpy
     public :: la_xcopy
     public :: la_xdot
     public :: la_xnrm2
     public :: la_xrot
     public :: la_xrotg
     public :: la_xrotm
     public :: la_xrotmg
     public :: la_xscal
     public :: la_xddot
     public :: la_xswap
     public :: la_xyasum
     public :: la_xynrm2
#endif
#ifdef LA_WITH_QP
     public :: la_qasum
     public :: la_qaxpy
     public :: la_qcopy
     public :: la_qdot
     public :: la_qnrm2
     public :: la_qrot
     public :: la_qrotg
     public :: la_qrotm
     public :: la_qrotmg
     public :: la_qscal
     public :: la_qddot
     public :: la_qswap
     public :: la_qwasum
     public :: la_qwnrm2
#endif
     public :: la_caxpy
     public :: la_ccopy
     public :: la_cdotc
     public :: la_cdotu
     public :: la_csrot
     public :: la_csscal
     public :: la_crotg
     public :: la_cscal
     public :: la_cswap
     public :: la_zaxpy
     public :: la_zcopy
     public :: la_zdotc
     public :: la_zdotu
     public :: la_zdrot
     public :: la_zdscal
     public :: la_zrotg
     public :: la_zscal
     public :: la_zswap
#ifdef LA_WITH_XDP
     public :: la_yaxpy
     public :: la_ycopy
     public :: la_ydotc
     public :: la_ydotu
     public :: la_yxrot
     public :: la_yxscal
     public :: la_yrotg
     public :: la_yscal
     public :: la_yswap
#endif
#ifdef LA_WITH_QP
     public :: la_waxpy
     public :: la_wcopy
     public :: la_wdotc
     public :: la_wdotu
     public :: la_wqrot
     public :: la_wqscal
     public :: la_wrotg
     public :: la_wscal
     public :: la_wswap
#endif
     public :: la_sdsdot

     contains

     !> SASUM: takes the sum of the absolute values.
     !> uses unrolled loops for increment equal to one.

     pure real(sp) function la_sasum(n,sx,incx)
        use la_constants_sp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: stemp
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: abs,mod
           la_sasum = zero
           stemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,6)
              if (m /= 0) then
                 do i = 1,m
                    stemp = stemp + abs(sx(i))
                 end do
                 if (n < 6) then
                    la_sasum = stemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,6
                 stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2)) + abs(sx(i + 3)) + abs(sx(i + &
                           4)) + abs(sx(i + 5))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + abs(sx(i))
              end do
           end if
           la_sasum = stemp
           return
     end function la_sasum
     !> DASUM: takes the sum of the absolute values.

     pure real(dp) function la_dasum(n,dx,incx)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(dp),intent(in) :: dx(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dtemp
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: abs,mod
           la_dasum = zero
           dtemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,6)
              if (m /= 0) then
                 do i = 1,m
                    dtemp = dtemp + abs(dx(i))
                 end do
                 if (n < 6) then
                    la_dasum = dtemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,6
                 dtemp = dtemp + abs(dx(i)) + abs(dx(i + 1)) + abs(dx(i + 2)) + abs(dx(i + 3)) + abs(dx(i + &
                           4)) + abs(dx(i + 5))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 dtemp = dtemp + abs(dx(i))
              end do
           end if
           la_dasum = dtemp
           return
     end function la_dasum
#ifdef LA_WITH_XDP
     !> XASUM: takes the sum of the absolute values.

     pure real(xdp) function la_xasum(n,xx,incx)
        use la_constants_xdp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(xdp),intent(in) :: xx(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: xtemp
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: abs,mod
           la_xasum = zero
           xtemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,6)
              if (m /= 0) then
                 do i = 1,m
                    xtemp = xtemp + abs(xx(i))
                 end do
                 if (n < 6) then
                    la_xasum = xtemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,6
                 xtemp = xtemp + abs(xx(i)) + abs(xx(i + 1)) + abs(xx(i + 2)) + abs(xx(i + 3)) + abs(xx(i + &
                           4)) + abs(xx(i + 5))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 xtemp = xtemp + abs(xx(i))
              end do
           end if
           la_xasum = xtemp
           return
     end function la_xasum
#endif
#ifdef LA_WITH_QP
     !> QASUM: takes the sum of the absolute values.

     pure real(qp) function la_qasum(n,qx,incx)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(qp),intent(in) :: qx(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: qtemp
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: abs,mod
           la_qasum = zero
           qtemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,6)
              if (m /= 0) then
                 do i = 1,m
                    qtemp = qtemp + abs(qx(i))
                 end do
                 if (n < 6) then
                    la_qasum = qtemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,6
                 qtemp = qtemp + abs(qx(i)) + abs(qx(i + 1)) + abs(qx(i + 2)) + abs(qx(i + 3)) + abs(qx(i + &
                           4)) + abs(qx(i + 5))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 qtemp = qtemp + abs(qx(i))
              end do
           end if
           la_qasum = qtemp
           return
     end function la_qasum
#endif

     !> SAXPY: constant times a vector plus a vector.
     !> uses unrolled loops for increments equal to one.

     pure subroutine la_saxpy(n,sa,sx,incx,sy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: sa
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*)
           real(sp),intent(inout) :: sy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (sa == 0.0_sp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,4)
              if (m /= 0) then
                 do i = 1,m
                    sy(i) = sy(i) + sa*sx(i)
                 end do
              end if
              if (n < 4) return
              mp1 = m + 1
              do i = mp1,n,4
                 sy(i) = sy(i) + sa*sx(i)
                 sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
                 sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
                 sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
               sy(iy) = sy(iy) + sa*sx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
     end subroutine la_saxpy
     !> DAXPY: constant times a vector plus a vector.
     !> uses unrolled loops for increments equal to one.

     pure subroutine la_daxpy(n,da,dx,incx,dy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: da
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: dx(*)
           real(dp),intent(inout) :: dy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (da == 0.0_dp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,4)
              if (m /= 0) then
                 do i = 1,m
                    dy(i) = dy(i) + da*dx(i)
                 end do
              end if
              if (n < 4) return
              mp1 = m + 1
              do i = mp1,n,4
                 dy(i) = dy(i) + da*dx(i)
                 dy(i + 1) = dy(i + 1) + da*dx(i + 1)
                 dy(i + 2) = dy(i + 2) + da*dx(i + 2)
                 dy(i + 3) = dy(i + 3) + da*dx(i + 3)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
               dy(iy) = dy(iy) + da*dx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
     end subroutine la_daxpy
#ifdef LA_WITH_XDP
     !> XAXPY: constant times a vector plus a vector.
     !> uses unrolled loops for increments equal to one.

     pure subroutine la_xaxpy(n,xa,xx,incx,xy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: xa
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(xdp),intent(in) :: xx(*)
           real(xdp),intent(inout) :: xy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (xa == 0.0_xdp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,4)
              if (m /= 0) then
                 do i = 1,m
                    xy(i) = xy(i) + xa*xx(i)
                 end do
              end if
              if (n < 4) return
              mp1 = m + 1
              do i = mp1,n,4
                 xy(i) = xy(i) + xa*xx(i)
                 xy(i + 1) = xy(i + 1) + xa*xx(i + 1)
                 xy(i + 2) = xy(i + 2) + xa*xx(i + 2)
                 xy(i + 3) = xy(i + 3) + xa*xx(i + 3)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
               xy(iy) = xy(iy) + xa*xx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
     end subroutine la_xaxpy
#endif
#ifdef LA_WITH_QP
     !> QAXPY: constant times a vector plus a vector.
     !> uses unrolled loops for increments equal to one.

     pure subroutine la_qaxpy(n,qa,qx,incx,qy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: qa
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(qp),intent(in) :: qx(*)
           real(qp),intent(inout) :: qy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (qa == 0.0_qp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,4)
              if (m /= 0) then
                 do i = 1,m
                    qy(i) = qy(i) + qa*qx(i)
                 end do
              end if
              if (n < 4) return
              mp1 = m + 1
              do i = mp1,n,4
                 qy(i) = qy(i) + qa*qx(i)
                 qy(i + 1) = qy(i + 1) + qa*qx(i + 1)
                 qy(i + 2) = qy(i + 2) + qa*qx(i + 2)
                 qy(i + 3) = qy(i + 3) + qa*qx(i + 3)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
               qy(iy) = qy(iy) + qa*qx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
     end subroutine la_qaxpy
#endif

     !> SCOPY: copies a vector, x, to a vector, y.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_scopy(n,sx,incx,sy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*)
           real(sp),intent(out) :: sy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,7)
              if (m /= 0) then
                 do i = 1,m
                    sy(i) = sx(i)
                 end do
                 if (n < 7) return
              end if
              mp1 = m + 1
              do i = mp1,n,7
                 sy(i) = sx(i)
                 sy(i + 1) = sx(i + 1)
                 sy(i + 2) = sx(i + 2)
                 sy(i + 3) = sx(i + 3)
                 sy(i + 4) = sx(i + 4)
                 sy(i + 5) = sx(i + 5)
                 sy(i + 6) = sx(i + 6)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 sy(iy) = sx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_scopy
     !> DCOPY: copies a vector, x, to a vector, y.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_dcopy(n,dx,incx,dy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: dx(*)
           real(dp),intent(out) :: dy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,7)
              if (m /= 0) then
                 do i = 1,m
                    dy(i) = dx(i)
                 end do
                 if (n < 7) return
              end if
              mp1 = m + 1
              do i = mp1,n,7
                 dy(i) = dx(i)
                 dy(i + 1) = dx(i + 1)
                 dy(i + 2) = dx(i + 2)
                 dy(i + 3) = dx(i + 3)
                 dy(i + 4) = dx(i + 4)
                 dy(i + 5) = dx(i + 5)
                 dy(i + 6) = dx(i + 6)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 dy(iy) = dx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_dcopy
#ifdef LA_WITH_XDP
     !> XCOPY: copies a vector, x, to a vector, y.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_xcopy(n,xx,incx,xy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(xdp),intent(in) :: xx(*)
           real(xdp),intent(out) :: xy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,7)
              if (m /= 0) then
                 do i = 1,m
                    xy(i) = xx(i)
                 end do
                 if (n < 7) return
              end if
              mp1 = m + 1
              do i = mp1,n,7
                 xy(i) = xx(i)
                 xy(i + 1) = xx(i + 1)
                 xy(i + 2) = xx(i + 2)
                 xy(i + 3) = xx(i + 3)
                 xy(i + 4) = xx(i + 4)
                 xy(i + 5) = xx(i + 5)
                 xy(i + 6) = xx(i + 6)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 xy(iy) = xx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_xcopy
#endif
#ifdef LA_WITH_QP
     !> QCOPY: copies a vector, x, to a vector, y.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_qcopy(n,qx,incx,qy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(qp),intent(in) :: qx(*)
           real(qp),intent(out) :: qy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,7)
              if (m /= 0) then
                 do i = 1,m
                    qy(i) = qx(i)
                 end do
                 if (n < 7) return
              end if
              mp1 = m + 1
              do i = mp1,n,7
                 qy(i) = qx(i)
                 qy(i + 1) = qx(i + 1)
                 qy(i + 2) = qx(i + 2)
                 qy(i + 3) = qx(i + 3)
                 qy(i + 4) = qx(i + 4)
                 qy(i + 5) = qx(i + 5)
                 qy(i + 6) = qx(i + 6)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 qy(iy) = qx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_qcopy
#endif

     !> SDOT: forms the dot product of two vectors.
     !> uses unrolled loops for increments equal to one.

     pure real(sp) function la_sdot(n,sx,incx,sy,incy)
        use la_constants_sp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*),sy(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: stemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           stemp = zero
           la_sdot = zero
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    stemp = stemp + sx(i)*sy(i)
                 end do
                 if (n < 5) then
                    la_sdot = stemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,5
               stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + &
                         sx(i + 4)*sy(i + 4)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 stemp = stemp + sx(ix)*sy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_sdot = stemp
           return
     end function la_sdot
     !> DDOT: forms the dot product of two vectors.
     !> uses unrolled loops for increments equal to one.

     pure real(dp) function la_ddot(n,dx,incx,dy,incy)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: dx(*),dy(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dtemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           la_ddot = zero
           dtemp = zero
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    dtemp = dtemp + dx(i)*dy(i)
                 end do
                 if (n < 5) then
                    la_ddot = dtemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,5
               dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + &
                         dx(i + 4)*dy(i + 4)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 dtemp = dtemp + dx(ix)*dy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_ddot = dtemp
           return
     end function la_ddot
#ifdef LA_WITH_XDP
     !> XDOT: forms the dot product of two vectors.
     !> uses unrolled loops for increments equal to one.

     pure real(xdp) function la_xdot(n,xx,incx,xy,incy)
        use la_constants_xdp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(xdp),intent(in) :: xx(*),xy(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: xtemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           la_xdot = zero
           xtemp = zero
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    xtemp = xtemp + xx(i)*xy(i)
                 end do
                 if (n < 5) then
                    la_xdot = xtemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,5
               xtemp = xtemp + xx(i)*xy(i) + xx(i + 1)*xy(i + 1) + xx(i + 2)*xy(i + 2) + xx(i + 3)*xy(i + 3) + &
                         xx(i + 4)*xy(i + 4)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 xtemp = xtemp + xx(ix)*xy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_xdot = xtemp
           return
     end function la_xdot
#endif
#ifdef LA_WITH_QP
     !> QDOT: forms the dot product of two vectors.
     !> uses unrolled loops for increments equal to one.

     pure real(qp) function la_qdot(n,qx,incx,qy,incy)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(qp),intent(in) :: qx(*),qy(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: qtemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           la_qdot = zero
           qtemp = zero
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    qtemp = qtemp + qx(i)*qy(i)
                 end do
                 if (n < 5) then
                    la_qdot = qtemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,5
               qtemp = qtemp + qx(i)*qy(i) + qx(i + 1)*qy(i + 1) + qx(i + 2)*qy(i + 2) + qx(i + 3)*qy(i + 3) + &
                         qx(i + 4)*qy(i + 4)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 qtemp = qtemp + qx(ix)*qy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_qdot = qtemp
           return
     end function la_qdot
#endif

     !> !
     !>
     !> SNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> SNRM2 := sqrt( x'*x ).

     pure function la_snrm2(n,x,incx)
        use la_constants_sp
        real(sp) :: la_snrm2
        ! -- reference blas level1 routine (version 3.9.1_sp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._sp)
        real(sp),parameter :: maxn = huge(0.0_sp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        real(sp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(sp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_snrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_snrm2 = scl*sqrt(sumsq)
        return
     end function la_snrm2
     !> !
     !>
     !> DNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> DNRM2 := sqrt( x'*x )

     pure function la_dnrm2(n,x,incx)
        use la_constants_dp
        real(dp) :: la_dnrm2
        ! -- reference blas level1 routine (version 3.9.1_dp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._dp)
        real(dp),parameter :: maxn = huge(0.0_dp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        real(dp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(dp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_dnrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_dnrm2 = scl*sqrt(sumsq)
        return
     end function la_dnrm2
#ifdef LA_WITH_XDP
     !> !
     !>
     !> XNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> XNRM2 := sqrt( x'*x )

     pure function la_xnrm2(n,x,incx)
        use la_constants_xdp
        real(xdp) :: la_xnrm2
        ! -- reference blas level1 routine (version 3.9.1_xdp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._xdp)
        real(xdp),parameter :: maxn = huge(0.0_xdp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        real(xdp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(xdp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_xnrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_xnrm2 = scl*sqrt(sumsq)
        return
     end function la_xnrm2
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> QNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> QNRM2 := sqrt( x'*x )

     pure function la_qnrm2(n,x,incx)
        use la_constants_qp
        real(qp) :: la_qnrm2
        ! -- reference blas level1 routine (version 3.9.1_qp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._qp)
        real(qp),parameter :: maxn = huge(0.0_qp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        real(qp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(qp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_qnrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_qnrm2 = scl*sqrt(sumsq)
        return
     end function la_qnrm2
#endif

     !> applies a plane rotation.

     pure subroutine la_srot(n,sx,incx,sy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: c,s
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(inout) :: sx(*),sy(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: stemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 stemp = c*sx(i) + s*sy(i)
                 sy(i) = c*sy(i) - s*sx(i)
                 sx(i) = stemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 stemp = c*sx(ix) + s*sy(iy)
                 sy(iy) = c*sy(iy) - s*sx(ix)
                 sx(ix) = stemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_srot
     !> DROT: applies a plane rotation.

     pure subroutine la_drot(n,dx,incx,dy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: c,s
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(inout) :: dx(*),dy(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dtemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 dtemp = c*dx(i) + s*dy(i)
                 dy(i) = c*dy(i) - s*dx(i)
                 dx(i) = dtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 dtemp = c*dx(ix) + s*dy(iy)
                 dy(iy) = c*dy(iy) - s*dx(ix)
                 dx(ix) = dtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_drot
#ifdef LA_WITH_XDP
     !> XROT: applies a plane rotation.

     pure subroutine la_xrot(n,xx,incx,xy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: c,s
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(xdp),intent(inout) :: xx(*),xy(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: xtemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 xtemp = c*xx(i) + s*xy(i)
                 xy(i) = c*xy(i) - s*xx(i)
                 xx(i) = xtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 xtemp = c*xx(ix) + s*xy(iy)
                 xy(iy) = c*xy(iy) - s*xx(ix)
                 xx(ix) = xtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_xrot
#endif
#ifdef LA_WITH_QP
     !> QROT: applies a plane rotation.

     pure subroutine la_qrot(n,qx,incx,qy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: c,s
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(qp),intent(inout) :: qx(*),qy(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: qtemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 qtemp = c*qx(i) + s*qy(i)
                 qy(i) = c*qy(i) - s*qx(i)
                 qx(i) = qtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 qtemp = c*qx(ix) + s*qy(iy)
                 qy(iy) = c*qy(iy) - s*qx(ix)
                 qx(ix) = qtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_qrot
#endif

     !> !
     !>
     !> The computation uses the formulas
     !> sigma = sgn(a)    if |a| >  |b|
     !> = sgn(b)    if |b| >= |a|
     !> r = sigma*sqrt( a**2 + b**2 )
     !> c = 1; s = 0      if r = 0
     !> c = a/r; s = b/r  if r != 0
     !> The subroutine also computes
     !> z = s    if |a| > |b|,
     !> = 1/c  if |b| >= |a| and c != 0
     !> = 1    if c = 0
     !> This allows c and s to be reconstructed from z as follows:
     !> If z = 1, set c = 0, s = 1.
     !> If |z| < 1, set c = sqrt(1 - z**2) and s = z.
     !> If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).

     pure subroutine la_srotg(a,b,c,s)
        use la_constants_sp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._sp)
        ! Scaling Constants
        ! Scalar Arguments
        real(sp),intent(inout) :: a,b
        real(sp),intent(out) :: c,s
        ! Local Scalars
        real(sp) :: anorm,bnorm,scl,sigma,r,z
        anorm = abs(a)
        bnorm = abs(b)
        if (bnorm == zero) then
           c = one
           s = zero
           b = zero
        else if (anorm == zero) then
           c = zero
           s = one
           a = b
           b = one
        else
           scl = min(safmax,max(safmin,anorm,bnorm))
           if (anorm > bnorm) then
              sigma = sign(one,a)
           else
              sigma = sign(one,b)
           end if
           r = sigma*(scl*sqrt((a/scl)**2 + (b/scl)**2))
           c = a/r
           s = b/r
           if (anorm > bnorm) then
              z = s
           else if (c /= zero) then
              z = one/c
           else
              z = one
           end if
           a = r
           b = z
        end if
        return
     end subroutine la_srotg
     !> !
     !>
     !> The computation uses the formulas
     !> sigma = sgn(a)    if |a| >  |b|
     !> = sgn(b)    if |b| >= |a|
     !> r = sigma*sqrt( a**2 + b**2 )
     !> c = 1; s = 0      if r = 0
     !> c = a/r; s = b/r  if r != 0
     !> The subroutine also computes
     !> z = s    if |a| > |b|,
     !> = 1/c  if |b| >= |a| and c != 0
     !> = 1    if c = 0
     !> This allows c and s to be reconstructed from z as follows:
     !> If z = 1, set c = 0, s = 1.
     !> If |z| < 1, set c = sqrt(1 - z**2) and s = z.
     !> If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).

     pure subroutine la_drotg(a,b,c,s)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._dp)
        ! Scaling Constants
        ! Scalar Arguments
        real(dp),intent(inout) :: a,b
        real(dp),intent(out) :: c,s
        ! Local Scalars
        real(dp) :: anorm,bnorm,scl,sigma,r,z
        anorm = abs(a)
        bnorm = abs(b)
        if (bnorm == zero) then
           c = one
           s = zero
           b = zero
        else if (anorm == zero) then
           c = zero
           s = one
           a = b
           b = one
        else
           scl = min(safmax,max(safmin,anorm,bnorm))
           if (anorm > bnorm) then
              sigma = sign(one,a)
           else
              sigma = sign(one,b)
           end if
           r = sigma*(scl*sqrt((a/scl)**2 + (b/scl)**2))
           c = a/r
           s = b/r
           if (anorm > bnorm) then
              z = s
           else if (c /= zero) then
              z = one/c
           else
              z = one
           end if
           a = r
           b = z
        end if
        return
     end subroutine la_drotg
#ifdef LA_WITH_XDP
     !> !
     !>
     !> The computation uses the formulas
     !> sigma = sgn(a)    if |a| >  |b|
     !> = sgn(b)    if |b| >= |a|
     !> r = sigma*sqrt( a**2 + b**2 )
     !> c = 1; s = 0      if r = 0
     !> c = a/r; s = b/r  if r != 0
     !> The subroutine also computes
     !> z = s    if |a| > |b|,
     !> = 1/c  if |b| >= |a| and c != 0
     !> = 1    if c = 0
     !> This allows c and s to be reconstructed from z as follows:
     !> If z = 1, set c = 0, s = 1.
     !> If |z| < 1, set c = sqrt(1 - z**2) and s = z.
     !> If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).

     pure subroutine la_xrotg(a,b,c,s)
        use la_constants_xdp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._xdp)
        ! Scaling Constants
        ! Scalar Arguments
        real(xdp),intent(inout) :: a,b
        real(xdp),intent(out) :: c,s
        ! Local Scalars
        real(xdp) :: anorm,bnorm,scl,sigma,r,z
        anorm = abs(a)
        bnorm = abs(b)
        if (bnorm == zero) then
           c = one
           s = zero
           b = zero
        else if (anorm == zero) then
           c = zero
           s = one
           a = b
           b = one
        else
           scl = min(safmax,max(safmin,anorm,bnorm))
           if (anorm > bnorm) then
              sigma = sign(one,a)
           else
              sigma = sign(one,b)
           end if
           r = sigma*(scl*sqrt((a/scl)**2 + (b/scl)**2))
           c = a/r
           s = b/r
           if (anorm > bnorm) then
              z = s
           else if (c /= zero) then
              z = one/c
           else
              z = one
           end if
           a = r
           b = z
        end if
        return
     end subroutine la_xrotg
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> The computation uses the formulas
     !> sigma = sgn(a)    if |a| >  |b|
     !> = sgn(b)    if |b| >= |a|
     !> r = sigma*sqrt( a**2 + b**2 )
     !> c = 1; s = 0      if r = 0
     !> c = a/r; s = b/r  if r != 0
     !> The subroutine also computes
     !> z = s    if |a| > |b|,
     !> = 1/c  if |b| >= |a| and c != 0
     !> = 1    if c = 0
     !> This allows c and s to be reconstructed from z as follows:
     !> If z = 1, set c = 0, s = 1.
     !> If |z| < 1, set c = sqrt(1 - z**2) and s = z.
     !> If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).

     pure subroutine la_qrotg(a,b,c,s)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._qp)
        ! Scaling Constants
        ! Scalar Arguments
        real(qp),intent(inout) :: a,b
        real(qp),intent(out) :: c,s
        ! Local Scalars
        real(qp) :: anorm,bnorm,scl,sigma,r,z
        anorm = abs(a)
        bnorm = abs(b)
        if (bnorm == zero) then
           c = one
           s = zero
           b = zero
        else if (anorm == zero) then
           c = zero
           s = one
           a = b
           b = one
        else
           scl = min(safmax,max(safmin,anorm,bnorm))
           if (anorm > bnorm) then
              sigma = sign(one,a)
           else
              sigma = sign(one,b)
           end if
           r = sigma*(scl*sqrt((a/scl)**2 + (b/scl)**2))
           c = a/r
           s = b/r
           if (anorm > bnorm) then
              z = s
           else if (c /= zero) then
              z = one/c
           else
              z = one
           end if
           a = r
           b = z
        end if
        return
     end subroutine la_qrotg
#endif

     !> APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
     !> (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
     !> (SX**T)
     !> SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
     !> LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
     !> WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> SFLAG=-1._sp     SFLAG=0._sp        SFLAG=1._sp     SFLAG=-2.E0
     !> (SH11  SH12)    (1._sp  SH12)    (SH11  1._sp)    (1._sp  0._sp)
     !> H=(          )    (          )    (          )    (          )
     !> (SH21  SH22),   (SH21  1._sp),   (-1._sp SH22),   (0._sp  1._sp).
     !> SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.

     pure subroutine la_srotm(n,sx,incx,sy,incy,sparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: sparam(5)
           real(sp),intent(inout) :: sx(*),sy(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: sflag,sh11,sh12,sh21,sh22,two,w,z,zero
           integer(ilp) :: i,kx,ky,nsteps
           ! Data Statements
           zero = 0.0_sp
           two = 2.0_sp
           sflag = sparam(1)
           if (n <= 0 .or. (sflag + two == zero)) return
           if (incx == incy .and. incx > 0) then
              nsteps = n*incx
              if (sflag < zero) then
                 sh11 = sparam(2)
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 sh22 = sparam(5)
                 do i = 1,nsteps,incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z*sh12
                    sy(i) = w*sh21 + z*sh22
                 end do
              else if (sflag == zero) then
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 do i = 1,nsteps,incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w + z*sh12
                    sy(i) = w*sh21 + z
                 end do
              else
                 sh11 = sparam(2)
                 sh22 = sparam(5)
                 do i = 1,nsteps,incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z
                    sy(i) = -w + sh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              if (sflag < zero) then
                 sh11 = sparam(2)
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 sh22 = sparam(5)
                 do i = 1,n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w*sh11 + z*sh12
                    sy(ky) = w*sh21 + z*sh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (sflag == zero) then
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 do i = 1,n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w + z*sh12
                    sy(ky) = w*sh21 + z
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else
                  sh11 = sparam(2)
                  sh22 = sparam(5)
                  do i = 1,n
                     w = sx(kx)
                     z = sy(ky)
                     sx(kx) = w*sh11 + z
                     sy(ky) = -w + sh22*z
                     kx = kx + incx
                     ky = ky + incy
                 end do
              end if
           end if
           return
     end subroutine la_srotm
     !> APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
     !> (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
     !> (DY**T)
     !> DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
     !> LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
     !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> DFLAG=-1._dp     DFLAG=0._dp        DFLAG=1._dp     DFLAG=-2.D0
     !> (DH11  DH12)    (1._dp  DH12)    (DH11  1._dp)    (1._dp  0._dp)
     !> H=(          )    (          )    (          )    (          )
     !> (DH21  DH22),   (DH21  1._dp),   (-1._dp DH22),   (0._dp  1._dp).
     !> SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.

     pure subroutine la_drotm(n,dx,incx,dy,incy,dparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: dparam(5)
           real(dp),intent(inout) :: dx(*),dy(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dflag,dh11,dh12,dh21,dh22,two,w,z,zero
           integer(ilp) :: i,kx,ky,nsteps
           ! Data Statements
           zero = 0.0_dp
           two = 2.0_dp
           dflag = dparam(1)
           if (n <= 0 .or. (dflag + two == zero)) return
           if (incx == incy .and. incx > 0) then
              nsteps = n*incx
              if (dflag < zero) then
                 dh11 = dparam(2)
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 dh22 = dparam(5)
                 do i = 1,nsteps,incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w*dh11 + z*dh12
                    dy(i) = w*dh21 + z*dh22
                 end do
              else if (dflag == zero) then
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 do i = 1,nsteps,incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w + z*dh12
                    dy(i) = w*dh21 + z
                 end do
              else
                 dh11 = dparam(2)
                 dh22 = dparam(5)
                 do i = 1,nsteps,incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w*dh11 + z
                    dy(i) = -w + dh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              if (dflag < zero) then
                 dh11 = dparam(2)
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 dh22 = dparam(5)
                 do i = 1,n
                    w = dx(kx)
                    z = dy(ky)
                    dx(kx) = w*dh11 + z*dh12
                    dy(ky) = w*dh21 + z*dh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (dflag == zero) then
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 do i = 1,n
                    w = dx(kx)
                    z = dy(ky)
                    dx(kx) = w + z*dh12
                    dy(ky) = w*dh21 + z
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else
                  dh11 = dparam(2)
                  dh22 = dparam(5)
                  do i = 1,n
                     w = dx(kx)
                     z = dy(ky)
                     dx(kx) = w*dh11 + z
                     dy(ky) = -w + dh22*z
                     kx = kx + incx
                     ky = ky + incy
                 end do
              end if
           end if
           return
     end subroutine la_drotm
#ifdef LA_WITH_XDP
     !> APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
     !> (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
     !> (DY**T)
     !> DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
     !> LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
     !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> DFLAG=-1._xdp     DFLAG=0._xdp        DFLAG=1._xdp     DFLAG=-2.D0
     !> (DH11  DH12)    (1._xdp  DH12)    (DH11  1._xdp)    (1._xdp  0._xdp)
     !> H=(          )    (          )    (          )    (          )
     !> (DH21  DH22),   (DH21  1._xdp),   (-1._xdp DH22),   (0._xdp  1._xdp).
     !> SEE XROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.

     pure subroutine la_xrotm(n,xx,incx,xy,incy,xparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(xdp),intent(in) :: xparam(5)
           real(xdp),intent(inout) :: xx(*),xy(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: xflag,xh11,xh12,xh21,xh22,two,w,z,zero
           integer(ilp) :: i,kx,ky,nsteps
           ! Data Statements
           zero = 0.0_xdp
           two = 2.0_xdp
           xflag = xparam(1)
           if (n <= 0 .or. (xflag + two == zero)) return
           if (incx == incy .and. incx > 0) then
              nsteps = n*incx
              if (xflag < zero) then
                 xh11 = xparam(2)
                 xh12 = xparam(4)
                 xh21 = xparam(3)
                 xh22 = xparam(5)
                 do i = 1,nsteps,incx
                    w = xx(i)
                    z = xy(i)
                    xx(i) = w*xh11 + z*xh12
                    xy(i) = w*xh21 + z*xh22
                 end do
              else if (xflag == zero) then
                 xh12 = xparam(4)
                 xh21 = xparam(3)
                 do i = 1,nsteps,incx
                    w = xx(i)
                    z = xy(i)
                    xx(i) = w + z*xh12
                    xy(i) = w*xh21 + z
                 end do
              else
                 xh11 = xparam(2)
                 xh22 = xparam(5)
                 do i = 1,nsteps,incx
                    w = xx(i)
                    z = xy(i)
                    xx(i) = w*xh11 + z
                    xy(i) = -w + xh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              if (xflag < zero) then
                 xh11 = xparam(2)
                 xh12 = xparam(4)
                 xh21 = xparam(3)
                 xh22 = xparam(5)
                 do i = 1,n
                    w = xx(kx)
                    z = xy(ky)
                    xx(kx) = w*xh11 + z*xh12
                    xy(ky) = w*xh21 + z*xh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (xflag == zero) then
                 xh12 = xparam(4)
                 xh21 = xparam(3)
                 do i = 1,n
                    w = xx(kx)
                    z = xy(ky)
                    xx(kx) = w + z*xh12
                    xy(ky) = w*xh21 + z
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else
                  xh11 = xparam(2)
                  xh22 = xparam(5)
                  do i = 1,n
                     w = xx(kx)
                     z = xy(ky)
                     xx(kx) = w*xh11 + z
                     xy(ky) = -w + xh22*z
                     kx = kx + incx
                     ky = ky + incy
                 end do
              end if
           end if
           return
     end subroutine la_xrotm
#endif
#ifdef LA_WITH_QP
     !> APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
     !> (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
     !> (DY**T)
     !> DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
     !> LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
     !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> DFLAG=-1._qp     DFLAG=0._qp        DFLAG=1._qp     DFLAG=-2.D0
     !> (DH11  DH12)    (1._qp  DH12)    (DH11  1._qp)    (1._qp  0._qp)
     !> H=(          )    (          )    (          )    (          )
     !> (DH21  DH22),   (DH21  1._qp),   (-1._qp DH22),   (0._qp  1._qp).
     !> SEE QROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.

     pure subroutine la_qrotm(n,qx,incx,qy,incy,qparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(qp),intent(in) :: qparam(5)
           real(qp),intent(inout) :: qx(*),qy(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: qflag,qh11,qh12,qh21,qh22,two,w,z,zero
           integer(ilp) :: i,kx,ky,nsteps
           ! Data Statements
           zero = 0.0_qp
           two = 2.0_qp
           qflag = qparam(1)
           if (n <= 0 .or. (qflag + two == zero)) return
           if (incx == incy .and. incx > 0) then
              nsteps = n*incx
              if (qflag < zero) then
                 qh11 = qparam(2)
                 qh12 = qparam(4)
                 qh21 = qparam(3)
                 qh22 = qparam(5)
                 do i = 1,nsteps,incx
                    w = qx(i)
                    z = qy(i)
                    qx(i) = w*qh11 + z*qh12
                    qy(i) = w*qh21 + z*qh22
                 end do
              else if (qflag == zero) then
                 qh12 = qparam(4)
                 qh21 = qparam(3)
                 do i = 1,nsteps,incx
                    w = qx(i)
                    z = qy(i)
                    qx(i) = w + z*qh12
                    qy(i) = w*qh21 + z
                 end do
              else
                 qh11 = qparam(2)
                 qh22 = qparam(5)
                 do i = 1,nsteps,incx
                    w = qx(i)
                    z = qy(i)
                    qx(i) = w*qh11 + z
                    qy(i) = -w + qh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              if (qflag < zero) then
                 qh11 = qparam(2)
                 qh12 = qparam(4)
                 qh21 = qparam(3)
                 qh22 = qparam(5)
                 do i = 1,n
                    w = qx(kx)
                    z = qy(ky)
                    qx(kx) = w*qh11 + z*qh12
                    qy(ky) = w*qh21 + z*qh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (qflag == zero) then
                 qh12 = qparam(4)
                 qh21 = qparam(3)
                 do i = 1,n
                    w = qx(kx)
                    z = qy(ky)
                    qx(kx) = w + z*qh12
                    qy(ky) = w*qh21 + z
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else
                  qh11 = qparam(2)
                  qh22 = qparam(5)
                  do i = 1,n
                     w = qx(kx)
                     z = qy(ky)
                     qx(kx) = w*qh11 + z
                     qy(ky) = -w + qh22*z
                     kx = kx + incx
                     ky = ky + incy
                 end do
              end if
           end if
           return
     end subroutine la_qrotm
#endif

     !> CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
     !> THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)    SY2)**T.
     !> WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> SFLAG=-1._sp     SFLAG=0._sp        SFLAG=1._sp     SFLAG=-2.E0
     !> (SH11  SH12)    (1._sp  SH12)    (SH11  1._sp)    (1._sp  0._sp)
     !> H=(          )    (          )    (          )    (          )
     !> (SH21  SH22),   (SH21  1._sp),   (-1._sp SH22),   (0._sp  1._sp).
     !> LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
     !> RESPECTIVELY. (VALUES OF 1._sp, -1._sp, OR 0._sp IMPLIED BY THE
     !> VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
     !> THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
     !> INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
     !> OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.

     pure subroutine la_srotmg(sd1,sd2,sx1,sy1,sparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(inout) :: sd1,sd2,sx1
           real(sp),intent(in) :: sy1
           ! Array Arguments
           real(sp),intent(out) :: sparam(5)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: gam,gamsq,one,rgamsq,sflag,sh11,sh12,sh21,sh22,sp1,sp2,sq1,sq2, &
                      stemp,su,two,zero
           ! Intrinsic Functions
           intrinsic :: abs
           ! Data Statements
           zero = 0.0_sp
           one = 1.0_sp
           two = 2.0_sp
           gam = 4096.0_sp
           gamsq = 1.67772e7_sp
           rgamsq = 5.96046e-8_sp
           if (sd1 < zero) then
              ! go zero-h-d-and-sx1..
              sflag = -one
              sh11 = zero
              sh12 = zero
              sh21 = zero
              sh22 = zero
              sd1 = zero
              sd2 = zero
              sx1 = zero
           else
              ! case-sd1-nonnegative
              sp2 = sd2*sy1
              if (sp2 == zero) then
                 sflag = -two
                 sparam(1) = sflag
                 return
              end if
              ! regular-case..
              sp1 = sd1*sx1
              sq2 = sp2*sy1
              sq1 = sp1*sx1
              if (abs(sq1) > abs(sq2)) then
                 sh21 = -sy1/sx1
                 sh12 = sp2/sp1
                 su = one - sh12*sh21
                if (su > zero) then
                  sflag = zero
                  sd1 = sd1/su
                  sd2 = sd2/su
                  sx1 = sx1*su
                else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  sflag = -one
                  sh11 = zero
                  sh12 = zero
                  sh21 = zero
                  sh22 = zero
                  sd1 = zero
                  sd2 = zero
                  sx1 = zero
                end if
              else
                 if (sq2 < zero) then
                    ! go zero-h-d-and-sx1..
                    sflag = -one
                    sh11 = zero
                    sh12 = zero
                    sh21 = zero
                    sh22 = zero
                    sd1 = zero
                    sd2 = zero
                    sx1 = zero
                 else
                    sflag = one
                    sh11 = sp1/sp2
                    sh22 = sx1/sy1
                    su = one + sh11*sh22
                    stemp = sd2/su
                    sd2 = sd1/su
                    sd1 = stemp
                    sx1 = sy1*su
                 end if
              end if
           ! procedure..scale-check
              if (sd1 /= zero) then
                 do while ((sd1 <= rgamsq) .or. (sd1 >= gamsq))
                    if (sflag == zero) then
                       sh11 = one
                       sh22 = one
                       sflag = -one
                    else
                       sh21 = -one
                       sh12 = one
                       sflag = -one
                    end if
                    if (sd1 <= rgamsq) then
                       sd1 = sd1*gam**2
                       sx1 = sx1/gam
                       sh11 = sh11/gam
                       sh12 = sh12/gam
                    else
                       sd1 = sd1/gam**2
                       sx1 = sx1*gam
                       sh11 = sh11*gam
                       sh12 = sh12*gam
                    end if
                 end do
              end if
              if (sd2 /= zero) then
                 do while ((abs(sd2) <= rgamsq) .or. (abs(sd2) >= gamsq))
                    if (sflag == zero) then
                       sh11 = one
                       sh22 = one
                       sflag = -one
                    else
                       sh21 = -one
                       sh12 = one
                       sflag = -one
                    end if
                    if (abs(sd2) <= rgamsq) then
                       sd2 = sd2*gam**2
                       sh21 = sh21/gam
                       sh22 = sh22/gam
                    else
                       sd2 = sd2/gam**2
                       sh21 = sh21*gam
                       sh22 = sh22*gam
                    end if
                 end do
              end if
           end if
           if (sflag < zero) then
              sparam(2) = sh11
              sparam(3) = sh21
              sparam(4) = sh12
              sparam(5) = sh22
           else if (sflag == zero) then
              sparam(3) = sh21
              sparam(4) = sh12
           else
              sparam(2) = sh11
              sparam(5) = sh22
           end if
           sparam(1) = sflag
           return
     end subroutine la_srotmg
     !> CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
     !> THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(DD1)*DX1,SQRT(DD2)    DY2)**T.
     !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> DFLAG=-1._dp     DFLAG=0._dp        DFLAG=1._dp     DFLAG=-2.D0
     !> (DH11  DH12)    (1._dp  DH12)    (DH11  1._dp)    (1._dp  0._dp)
     !> H=(          )    (          )    (          )    (          )
     !> (DH21  DH22),   (DH21  1._dp),   (-1._dp DH22),   (0._dp  1._dp).
     !> LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
     !> RESPECTIVELY. (VALUES OF 1._dp, -1._dp, OR 0._dp IMPLIED BY THE
     !> VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
     !> THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
     !> INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
     !> OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.

     pure subroutine la_drotmg(dd1,dd2,dx1,dy1,dparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(inout) :: dd1,dd2,dx1
           real(dp),intent(in) :: dy1
           ! Array Arguments
           real(dp),intent(out) :: dparam(5)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dflag,dh11,dh12,dh21,dh22,dp1,dp2,dq1,dq2,dtemp,du,gam,gamsq, &
                     one,rgamsq,two,zero
           ! Intrinsic Functions
           intrinsic :: abs
           ! Data Statements
           zero = 0.0_dp
           one = 1.0_dp
           two = 2.0_dp
           gam = 4096.0_dp
           gamsq = 16777216.0_dp
           rgamsq = 5.9604645e-8_dp
           if (dd1 < zero) then
              ! go zero-h-d-and-dx1..
              dflag = -one
              dh11 = zero
              dh12 = zero
              dh21 = zero
              dh22 = zero
              dd1 = zero
              dd2 = zero
              dx1 = zero
           else
              ! case-dd1-nonnegative
              dp2 = dd2*dy1
              if (dp2 == zero) then
                 dflag = -two
                 dparam(1) = dflag
                 return
              end if
              ! regular-case..
              dp1 = dd1*dx1
              dq2 = dp2*dy1
              dq1 = dp1*dx1
              if (abs(dq1) > abs(dq2)) then
                 dh21 = -dy1/dx1
                 dh12 = dp2/dp1
                 du = one - dh12*dh21
                if (du > zero) then
                  dflag = zero
                  dd1 = dd1/du
                  dd2 = dd2/du
                  dx1 = dx1*du
                else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  dflag = -one
                  dh11 = zero
                  dh12 = zero
                  dh21 = zero
                  dh22 = zero
                  dd1 = zero
                  dd2 = zero
                  dx1 = zero
                end if
              else
                 if (dq2 < zero) then
                    ! go zero-h-d-and-dx1..
                    dflag = -one
                    dh11 = zero
                    dh12 = zero
                    dh21 = zero
                    dh22 = zero
                    dd1 = zero
                    dd2 = zero
                    dx1 = zero
                 else
                    dflag = one
                    dh11 = dp1/dp2
                    dh22 = dx1/dy1
                    du = one + dh11*dh22
                    dtemp = dd2/du
                    dd2 = dd1/du
                    dd1 = dtemp
                    dx1 = dy1*du
                 end if
              end if
           ! procedure..scale-check
              if (dd1 /= zero) then
                 do while ((dd1 <= rgamsq) .or. (dd1 >= gamsq))
                    if (dflag == zero) then
                       dh11 = one
                       dh22 = one
                       dflag = -one
                    else
                       dh21 = -one
                       dh12 = one
                       dflag = -one
                    end if
                    if (dd1 <= rgamsq) then
                       dd1 = dd1*gam**2
                       dx1 = dx1/gam
                       dh11 = dh11/gam
                       dh12 = dh12/gam
                    else
                       dd1 = dd1/gam**2
                       dx1 = dx1*gam
                       dh11 = dh11*gam
                       dh12 = dh12*gam
                    end if
                 end do
              end if
              if (dd2 /= zero) then
                 do while ((abs(dd2) <= rgamsq) .or. (abs(dd2) >= gamsq))
                    if (dflag == zero) then
                       dh11 = one
                       dh22 = one
                       dflag = -one
                    else
                       dh21 = -one
                       dh12 = one
                       dflag = -one
                    end if
                    if (abs(dd2) <= rgamsq) then
                       dd2 = dd2*gam**2
                       dh21 = dh21/gam
                       dh22 = dh22/gam
                    else
                       dd2 = dd2/gam**2
                       dh21 = dh21*gam
                       dh22 = dh22*gam
                    end if
                 end do
              end if
           end if
           if (dflag < zero) then
              dparam(2) = dh11
              dparam(3) = dh21
              dparam(4) = dh12
              dparam(5) = dh22
           else if (dflag == zero) then
              dparam(3) = dh21
              dparam(4) = dh12
           else
              dparam(2) = dh11
              dparam(5) = dh22
           end if
           dparam(1) = dflag
           return
     end subroutine la_drotmg
#ifdef LA_WITH_XDP
     !> CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
     !> THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(DD1)*DX1,SQRT(DD2)    DY2)**T.
     !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> DFLAG=-1._xdp     DFLAG=0._xdp        DFLAG=1._xdp     DFLAG=-2.D0
     !> (DH11  DH12)    (1._xdp  DH12)    (DH11  1._xdp)    (1._xdp  0._xdp)
     !> H=(          )    (          )    (          )    (          )
     !> (DH21  DH22),   (DH21  1._xdp),   (-1._xdp DH22),   (0._xdp  1._xdp).
     !> LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
     !> RESPECTIVELY. (VALUES OF 1._xdp, -1._xdp, OR 0._xdp IMPLIED BY THE
     !> VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
     !> THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
     !> INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
     !> OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.

     pure subroutine la_xrotmg(xd1,xd2,xx1,xy1,xparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(inout) :: xd1,xd2,xx1
           real(xdp),intent(in) :: xy1
           ! Array Arguments
           real(xdp),intent(out) :: xparam(5)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: xflag,xh11,xh12,xh21,xh22,xp1,xp2,xq1,xq2,xtemp,xu,gam,gamsq, &
                     one,rgamsq,two,zero
           ! Intrinsic Functions
           intrinsic :: abs
           ! Data Statements
           zero = 0.0_xdp
           one = 1.0_xdp
           two = 2.0_xdp
           gam = 4096.0_xdp
           gamsq = 16777216.0_xdp
           rgamsq = 5.9604645e-8_xdp
           if (xd1 < zero) then
              ! go zero-h-d-and-xx1..
              xflag = -one
              xh11 = zero
              xh12 = zero
              xh21 = zero
              xh22 = zero
              xd1 = zero
              xd2 = zero
              xx1 = zero
           else
              ! case-xd1-nonnegative
              xp2 = xd2*xy1
              if (xp2 == zero) then
                 xflag = -two
                 xparam(1) = xflag
                 return
              end if
              ! regular-case..
              xp1 = xd1*xx1
              xq2 = xp2*xy1
              xq1 = xp1*xx1
              if (abs(xq1) > abs(xq2)) then
                 xh21 = -xy1/xx1
                 xh12 = xp2/xp1
                 xu = one - xh12*xh21
                if (xu > zero) then
                  xflag = zero
                  xd1 = xd1/xu
                  xd2 = xd2/xu
                  xx1 = xx1*xu
                else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  xflag = -one
                  xh11 = zero
                  xh12 = zero
                  xh21 = zero
                  xh22 = zero
                  xd1 = zero
                  xd2 = zero
                  xx1 = zero
                end if
              else
                 if (xq2 < zero) then
                    ! go zero-h-d-and-xx1..
                    xflag = -one
                    xh11 = zero
                    xh12 = zero
                    xh21 = zero
                    xh22 = zero
                    xd1 = zero
                    xd2 = zero
                    xx1 = zero
                 else
                    xflag = one
                    xh11 = xp1/xp2
                    xh22 = xx1/xy1
                    xu = one + xh11*xh22
                    xtemp = xd2/xu
                    xd2 = xd1/xu
                    xd1 = xtemp
                    xx1 = xy1*xu
                 end if
              end if
           ! procedure..scale-check
              if (xd1 /= zero) then
                 do while ((xd1 <= rgamsq) .or. (xd1 >= gamsq))
                    if (xflag == zero) then
                       xh11 = one
                       xh22 = one
                       xflag = -one
                    else
                       xh21 = -one
                       xh12 = one
                       xflag = -one
                    end if
                    if (xd1 <= rgamsq) then
                       xd1 = xd1*gam**2
                       xx1 = xx1/gam
                       xh11 = xh11/gam
                       xh12 = xh12/gam
                    else
                       xd1 = xd1/gam**2
                       xx1 = xx1*gam
                       xh11 = xh11*gam
                       xh12 = xh12*gam
                    end if
                 end do
              end if
              if (xd2 /= zero) then
                 do while ((abs(xd2) <= rgamsq) .or. (abs(xd2) >= gamsq))
                    if (xflag == zero) then
                       xh11 = one
                       xh22 = one
                       xflag = -one
                    else
                       xh21 = -one
                       xh12 = one
                       xflag = -one
                    end if
                    if (abs(xd2) <= rgamsq) then
                       xd2 = xd2*gam**2
                       xh21 = xh21/gam
                       xh22 = xh22/gam
                    else
                       xd2 = xd2/gam**2
                       xh21 = xh21*gam
                       xh22 = xh22*gam
                    end if
                 end do
              end if
           end if
           if (xflag < zero) then
              xparam(2) = xh11
              xparam(3) = xh21
              xparam(4) = xh12
              xparam(5) = xh22
           else if (xflag == zero) then
              xparam(3) = xh21
              xparam(4) = xh12
           else
              xparam(2) = xh11
              xparam(5) = xh22
           end if
           xparam(1) = xflag
           return
     end subroutine la_xrotmg
#endif
#ifdef LA_WITH_QP
     !> CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
     !> THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(DD1)*DX1,SQRT(DD2)    DY2)**T.
     !> WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     !> DFLAG=-1._qp     DFLAG=0._qp        DFLAG=1._qp     DFLAG=-2.D0
     !> (DH11  DH12)    (1._qp  DH12)    (DH11  1._qp)    (1._qp  0._qp)
     !> H=(          )    (          )    (          )    (          )
     !> (DH21  DH22),   (DH21  1._qp),   (-1._qp DH22),   (0._qp  1._qp).
     !> LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
     !> RESPECTIVELY. (VALUES OF 1._qp, -1._qp, OR 0._qp IMPLIED BY THE
     !> VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
     !> THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
     !> INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
     !> OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.

     pure subroutine la_qrotmg(qd1,qd2,qx1,qy1,qparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(inout) :: qd1,qd2,qx1
           real(qp),intent(in) :: qy1
           ! Array Arguments
           real(qp),intent(out) :: qparam(5)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: qflag,qh11,qh12,qh21,qh22,qp1,qp2,qq1,qq2,qtemp,qu,gam,gamsq, &
                     one,rgamsq,two,zero
           ! Intrinsic Functions
           intrinsic :: abs
           ! Data Statements
           zero = 0.0_qp
           one = 1.0_qp
           two = 2.0_qp
           gam = 4096.0_qp
           gamsq = 16777216.0_qp
           rgamsq = 5.9604645e-8_qp
           if (qd1 < zero) then
              ! go zero-h-d-and-qx1..
              qflag = -one
              qh11 = zero
              qh12 = zero
              qh21 = zero
              qh22 = zero
              qd1 = zero
              qd2 = zero
              qx1 = zero
           else
              ! case-qd1-nonnegative
              qp2 = qd2*qy1
              if (qp2 == zero) then
                 qflag = -two
                 qparam(1) = qflag
                 return
              end if
              ! regular-case..
              qp1 = qd1*qx1
              qq2 = qp2*qy1
              qq1 = qp1*qx1
              if (abs(qq1) > abs(qq2)) then
                 qh21 = -qy1/qx1
                 qh12 = qp2/qp1
                 qu = one - qh12*qh21
                if (qu > zero) then
                  qflag = zero
                  qd1 = qd1/qu
                  qd2 = qd2/qu
                  qx1 = qx1*qu
                else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  qflag = -one
                  qh11 = zero
                  qh12 = zero
                  qh21 = zero
                  qh22 = zero
                  qd1 = zero
                  qd2 = zero
                  qx1 = zero
                end if
              else
                 if (qq2 < zero) then
                    ! go zero-h-d-and-qx1..
                    qflag = -one
                    qh11 = zero
                    qh12 = zero
                    qh21 = zero
                    qh22 = zero
                    qd1 = zero
                    qd2 = zero
                    qx1 = zero
                 else
                    qflag = one
                    qh11 = qp1/qp2
                    qh22 = qx1/qy1
                    qu = one + qh11*qh22
                    qtemp = qd2/qu
                    qd2 = qd1/qu
                    qd1 = qtemp
                    qx1 = qy1*qu
                 end if
              end if
           ! procedure..scale-check
              if (qd1 /= zero) then
                 do while ((qd1 <= rgamsq) .or. (qd1 >= gamsq))
                    if (qflag == zero) then
                       qh11 = one
                       qh22 = one
                       qflag = -one
                    else
                       qh21 = -one
                       qh12 = one
                       qflag = -one
                    end if
                    if (qd1 <= rgamsq) then
                       qd1 = qd1*gam**2
                       qx1 = qx1/gam
                       qh11 = qh11/gam
                       qh12 = qh12/gam
                    else
                       qd1 = qd1/gam**2
                       qx1 = qx1*gam
                       qh11 = qh11*gam
                       qh12 = qh12*gam
                    end if
                 end do
              end if
              if (qd2 /= zero) then
                 do while ((abs(qd2) <= rgamsq) .or. (abs(qd2) >= gamsq))
                    if (qflag == zero) then
                       qh11 = one
                       qh22 = one
                       qflag = -one
                    else
                       qh21 = -one
                       qh12 = one
                       qflag = -one
                    end if
                    if (abs(qd2) <= rgamsq) then
                       qd2 = qd2*gam**2
                       qh21 = qh21/gam
                       qh22 = qh22/gam
                    else
                       qd2 = qd2/gam**2
                       qh21 = qh21*gam
                       qh22 = qh22*gam
                    end if
                 end do
              end if
           end if
           if (qflag < zero) then
              qparam(2) = qh11
              qparam(3) = qh21
              qparam(4) = qh12
              qparam(5) = qh22
           else if (qflag == zero) then
              qparam(3) = qh21
              qparam(4) = qh12
           else
              qparam(2) = qh11
              qparam(5) = qh22
           end if
           qparam(1) = qflag
           return
     end subroutine la_qrotmg
#endif

     !> SSCAL: scales a vector by a constant.
     !> uses unrolled loops for increment equal to 1.

     pure subroutine la_sscal(n,sa,sx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: sa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(sp),intent(inout) :: sx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    sx(i) = sa*sx(i)
                 end do
                 if (n < 5) return
              end if
              mp1 = m + 1
              do i = mp1,n,5
                 sx(i) = sa*sx(i)
                 sx(i + 1) = sa*sx(i + 1)
                 sx(i + 2) = sa*sx(i + 2)
                 sx(i + 3) = sa*sx(i + 3)
                 sx(i + 4) = sa*sx(i + 4)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 sx(i) = sa*sx(i)
              end do
           end if
           return
     end subroutine la_sscal
     !> DSCAL: scales a vector by a constant.
     !> uses unrolled loops for increment equal to 1.

     pure subroutine la_dscal(n,da,dx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: da
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(dp),intent(inout) :: dx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    dx(i) = da*dx(i)
                 end do
                 if (n < 5) return
              end if
              mp1 = m + 1
              do i = mp1,n,5
                 dx(i) = da*dx(i)
                 dx(i + 1) = da*dx(i + 1)
                 dx(i + 2) = da*dx(i + 2)
                 dx(i + 3) = da*dx(i + 3)
                 dx(i + 4) = da*dx(i + 4)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 dx(i) = da*dx(i)
              end do
           end if
           return
     end subroutine la_dscal
#ifdef LA_WITH_XDP
     !> XSCAL: scales a vector by a constant.
     !> uses unrolled loops for increment equal to 1.

     pure subroutine la_xscal(n,xa,xx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: xa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(xdp),intent(inout) :: xx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    xx(i) = xa*xx(i)
                 end do
                 if (n < 5) return
              end if
              mp1 = m + 1
              do i = mp1,n,5
                 xx(i) = xa*xx(i)
                 xx(i + 1) = xa*xx(i + 1)
                 xx(i + 2) = xa*xx(i + 2)
                 xx(i + 3) = xa*xx(i + 3)
                 xx(i + 4) = xa*xx(i + 4)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 xx(i) = xa*xx(i)
              end do
           end if
           return
     end subroutine la_xscal
#endif
#ifdef LA_WITH_QP
     !> QSCAL: scales a vector by a constant.
     !> uses unrolled loops for increment equal to 1.

     pure subroutine la_qscal(n,qa,qx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: qa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           real(qp),intent(inout) :: qx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,m,mp1,nincx
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n,5)
              if (m /= 0) then
                 do i = 1,m
                    qx(i) = qa*qx(i)
                 end do
                 if (n < 5) return
              end if
              mp1 = m + 1
              do i = mp1,n,5
                 qx(i) = qa*qx(i)
                 qx(i + 1) = qa*qx(i + 1)
                 qx(i + 2) = qa*qx(i + 2)
                 qx(i + 3) = qa*qx(i + 3)
                 qx(i + 4) = qa*qx(i + 4)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 qx(i) = qa*qx(i)
              end do
           end if
           return
     end subroutine la_qscal
#endif

     !> Compute the inner product of two vectors with extended
     !> precision accumulation and result.
     !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
     !> DSDOT: = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
     !> where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
     !> defined in a similar way using INCY.

     pure real(dp) function la_dsdot(n,sx,incx,sy,incy)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*),sy(*)
        ! authors:
        ! ========
        ! lawson, c. l., (jpl), hanson, r. j., (snla),
        ! kincaid, d. r., (u. of texas), krogh, f. t., (jpl)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,kx,ky,ns
           ! Intrinsic Functions
           intrinsic :: real
           la_dsdot = zero
           if (n <= 0) return
           if (incx == incy .and. incx > 0) then
           ! code for equal, positive, non-unit increments.
              ns = n*incx
              do i = 1,ns,incx
                 la_dsdot = la_dsdot + real(sx(i),KIND=dp)*real(sy(i),KIND=dp)
              end do
           else
           ! code for unequal or nonpositive increments.
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              do i = 1,n
                 la_dsdot = la_dsdot + real(sx(kx),KIND=dp)*real(sy(ky),KIND=dp)
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           return
     end function la_dsdot
#ifdef LA_WITH_XDP
     !> Compute the inner product of two vectors with extended
     !> precision accumulation and result.
     !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
     !> XDDOT: = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
     !> where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
     !> defined in a similar way using INCY.

     pure real(xdp) function la_xddot(n,sx,incx,sy,incy)
        use la_constants_xdp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: sx(*),sy(*)
        ! authors:
        ! ========
        ! lawson, c. l., (jpl), hanson, r. j., (snla),
        ! kincaid, d. r., (u. of texas), krogh, f. t., (jpl)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,kx,ky,ns
           ! Intrinsic Functions
           intrinsic :: real
           la_xddot = zero
           if (n <= 0) return
           if (incx == incy .and. incx > 0) then
           ! code for equal, positive, non-unit increments.
              ns = n*incx
              do i = 1,ns,incx
                 la_xddot = la_xddot + real(sx(i),KIND=xdp)*real(sy(i),KIND=xdp)
              end do
           else
           ! code for unequal or nonpositive increments.
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              do i = 1,n
                 la_xddot = la_xddot + real(sx(kx),KIND=xdp)*real(sy(ky),KIND=xdp)
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           return
     end function la_xddot
#endif
#ifdef LA_WITH_QP
     !> Compute the inner product of two vectors with extended
     !> precision accumulation and result.
     !> Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
     !> QDDOT: = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
     !> where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
     !> defined in a similar way using INCY.

     pure real(qp) function la_qddot(n,sx,incx,sy,incy)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: sx(*),sy(*)
        ! authors:
        ! ========
        ! lawson, c. l., (jpl), hanson, r. j., (snla),
        ! kincaid, d. r., (u. of texas), krogh, f. t., (jpl)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,kx,ky,ns
           ! Intrinsic Functions
           intrinsic :: real
           la_qddot = zero
           if (n <= 0) return
           if (incx == incy .and. incx > 0) then
           ! code for equal, positive, non-unit increments.
              ns = n*incx
              do i = 1,ns,incx
                 la_qddot = la_qddot + real(sx(i),KIND=qp)*real(sy(i),KIND=qp)
              end do
           else
           ! code for unequal or nonpositive increments.
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              do i = 1,n
                 la_qddot = la_qddot + real(sx(kx),KIND=qp)*real(sy(ky),KIND=qp)
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           return
     end function la_qddot
#endif

     !> SSWAP: interchanges two vectors.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_sswap(n,sx,incx,sy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(inout) :: sx(*),sy(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: stemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
             ! clean-up loop
              m = mod(n,3)
              if (m /= 0) then
                 do i = 1,m
                    stemp = sx(i)
                    sx(i) = sy(i)
                    sy(i) = stemp
                 end do
                 if (n < 3) return
              end if
              mp1 = m + 1
              do i = mp1,n,3
                 stemp = sx(i)
                 sx(i) = sy(i)
                 sy(i) = stemp
                 stemp = sx(i + 1)
                 sx(i + 1) = sy(i + 1)
                 sy(i + 1) = stemp
                 stemp = sx(i + 2)
                 sx(i + 2) = sy(i + 2)
                 sy(i + 2) = stemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 stemp = sx(ix)
                 sx(ix) = sy(iy)
                 sy(iy) = stemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_sswap
     !> DSWAP: interchanges two vectors.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_dswap(n,dx,incx,dy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(dp),intent(inout) :: dx(*),dy(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: dtemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
             ! clean-up loop
              m = mod(n,3)
              if (m /= 0) then
                 do i = 1,m
                    dtemp = dx(i)
                    dx(i) = dy(i)
                    dy(i) = dtemp
                 end do
                 if (n < 3) return
              end if
              mp1 = m + 1
              do i = mp1,n,3
                 dtemp = dx(i)
                 dx(i) = dy(i)
                 dy(i) = dtemp
                 dtemp = dx(i + 1)
                 dx(i + 1) = dy(i + 1)
                 dy(i + 1) = dtemp
                 dtemp = dx(i + 2)
                 dx(i + 2) = dy(i + 2)
                 dy(i + 2) = dtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 dtemp = dx(ix)
                 dx(ix) = dy(iy)
                 dy(iy) = dtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_dswap
#ifdef LA_WITH_XDP
     !> XSWAP: interchanges two vectors.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_xswap(n,xx,incx,xy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(xdp),intent(inout) :: xx(*),xy(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: xtemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
             ! clean-up loop
              m = mod(n,3)
              if (m /= 0) then
                 do i = 1,m
                    xtemp = xx(i)
                    xx(i) = xy(i)
                    xy(i) = xtemp
                 end do
                 if (n < 3) return
              end if
              mp1 = m + 1
              do i = mp1,n,3
                 xtemp = xx(i)
                 xx(i) = xy(i)
                 xy(i) = xtemp
                 xtemp = xx(i + 1)
                 xx(i + 1) = xy(i + 1)
                 xy(i + 1) = xtemp
                 xtemp = xx(i + 2)
                 xx(i + 2) = xy(i + 2)
                 xy(i + 2) = xtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 xtemp = xx(ix)
                 xx(ix) = xy(iy)
                 xy(iy) = xtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_xswap
#endif
#ifdef LA_WITH_QP
     !> QSWAP: interchanges two vectors.
     !> uses unrolled loops for increments equal to 1.

     pure subroutine la_qswap(n,qx,incx,qy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(qp),intent(inout) :: qx(*),qy(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: qtemp
           integer(ilp) :: i,ix,iy,m,mp1
           ! Intrinsic Functions
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
             ! clean-up loop
              m = mod(n,3)
              if (m /= 0) then
                 do i = 1,m
                    qtemp = qx(i)
                    qx(i) = qy(i)
                    qy(i) = qtemp
                 end do
                 if (n < 3) return
              end if
              mp1 = m + 1
              do i = mp1,n,3
                 qtemp = qx(i)
                 qx(i) = qy(i)
                 qy(i) = qtemp
                 qtemp = qx(i + 1)
                 qx(i + 1) = qy(i + 1)
                 qy(i + 1) = qtemp
                 qtemp = qx(i + 2)
                 qx(i + 2) = qy(i + 2)
                 qy(i + 2) = qtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 qtemp = qx(ix)
                 qx(ix) = qy(iy)
                 qy(iy) = qtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_qswap
#endif

     !> SCASUM: takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
     !> returns a single precision result.

     pure real(sp) function la_scasum(n,cx,incx)
        use la_constants_sp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: stemp
           integer(ilp) :: i,nincx
           ! Intrinsic Functions
           intrinsic :: abs,aimag,real
           la_scasum = zero
           stemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 stemp = stemp + abs(real(cx(i),KIND=sp)) + abs(aimag(cx(i)))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + abs(real(cx(i),KIND=sp)) + abs(aimag(cx(i)))
              end do
           end if
           la_scasum = stemp
           return
     end function la_scasum
     !> DZASUM: takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
     !> returns a double precision result.

     pure real(dp) function la_dzasum(n,zx,incx)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: stemp
           integer(ilp) :: i,nincx
           la_dzasum = zero
           stemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 stemp = stemp + la_dcabs1(zx(i))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + la_dcabs1(zx(i))
              end do
           end if
           la_dzasum = stemp
           return
     end function la_dzasum
#ifdef LA_WITH_XDP
     !> XYASUM: takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
     !> returns a extended precision result.

     pure real(xdp) function la_xyasum(n,yx,incx)
        use la_constants_xdp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(xdp),intent(in) :: yx(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: stemp
           integer(ilp) :: i,nincx
           la_xyasum = zero
           stemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 stemp = stemp + la_xcabs1(yx(i))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + la_xcabs1(yx(i))
              end do
           end if
           la_xyasum = stemp
           return
     end function la_xyasum
#endif
#ifdef LA_WITH_QP
     !> QWASUM: takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
     !> returns a quad precision result.

     pure real(qp) function la_qwasum(n,wx,incx)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(in) :: wx(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: stemp
           integer(ilp) :: i,nincx
           la_qwasum = zero
           stemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 stemp = stemp + la_qcabs1(wx(i))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + la_qcabs1(wx(i))
              end do
           end if
           la_qwasum = stemp
           return
     end function la_qwasum
#endif

     !> !
     !>
     !> SCNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> SCNRM2 := sqrt( x**H*x )

     pure function la_scnrm2(n,x,incx)
        use la_constants_sp
        real(sp) :: la_scnrm2
        ! -- reference blas level1 routine (version 3.9.1_sp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._sp)
        real(sp),parameter :: maxn = huge(0.0_sp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        complex(sp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(sp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_scnrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=sp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_scnrm2 = scl*sqrt(sumsq)
        return
     end function la_scnrm2
     !> !
     !>
     !> DZNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> DZNRM2 := sqrt( x**H*x )

     pure function la_dznrm2(n,x,incx)
        use la_constants_dp
        real(dp) :: la_dznrm2
        ! -- reference blas level1 routine (version 3.9.1_dp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._dp)
        real(dp),parameter :: maxn = huge(0.0_dp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        complex(dp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(dp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_dznrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=dp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_dznrm2 = scl*sqrt(sumsq)
        return
     end function la_dznrm2
#ifdef LA_WITH_XDP
     !> !
     !>
     !> XYNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> XYNRM2 := sqrt( x**H*x )

     pure function la_xynrm2(n,x,incx)
        use la_constants_xdp
        real(xdp) :: la_xynrm2
        ! -- reference blas level1 routine (version 3.9.1_xdp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._xdp)
        real(xdp),parameter :: maxn = huge(0.0_xdp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        complex(xdp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(xdp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_xynrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=xdp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_xynrm2 = scl*sqrt(sumsq)
        return
     end function la_xynrm2
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> QWNRM2: returns the euclidean norm of a vector via the function
     !> name, so that
     !> QWNRM2 := sqrt( x**H*x )

     pure function la_qwnrm2(n,x,incx)
        use la_constants_qp
        real(qp) :: la_qwnrm2
        ! -- reference blas level1 routine (version 3.9.1_qp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! Constants
        integer,parameter :: wp = kind(1._qp)
        real(qp),parameter :: maxn = huge(0.0_qp)
        ! .. blue's scaling constants ..
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        ! Array Arguments
        complex(qp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(qp) :: abig,amed,asml,ax,scl,sumsq,ymax,ymin
        ! quick return if possible
        la_qwnrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=qp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        la_qwnrm2 = scl*sqrt(sumsq)
        return
     end function la_qwnrm2
#endif

     !> CAXPY: constant times a vector plus a vector.

     pure subroutine la_caxpy(n,ca,cx,incx,cy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: ca
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*)
           complex(sp),intent(inout) :: cy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (la_scabs1(ca) == 0.0e+0_sp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 cy(i) = cy(i) + ca*cx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 cy(iy) = cy(iy) + ca*cx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_caxpy
     !> ZAXPY: constant times a vector plus a vector.

     pure subroutine la_zaxpy(n,za,zx,incx,zy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: za
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*)
           complex(dp),intent(inout) :: zy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (la_dcabs1(za) == 0.0_dp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 zy(i) = zy(i) + za*zx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 zy(iy) = zy(iy) + za*zx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_zaxpy
#ifdef LA_WITH_XDP
     !> YAXPY: constant times a vector plus a vector.

     pure subroutine la_yaxpy(n,ya,yx,incx,yy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: ya
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(xdp),intent(in) :: yx(*)
           complex(xdp),intent(inout) :: yy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (la_xcabs1(ya) == 0.0_xdp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 yy(i) = yy(i) + ya*yx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 yy(iy) = yy(iy) + ya*yx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_yaxpy
#endif
#ifdef LA_WITH_QP
     !> WAXPY: constant times a vector plus a vector.

     pure subroutine la_waxpy(n,wa,wx,incx,wy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: wa
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(qp),intent(in) :: wx(*)
           complex(qp),intent(inout) :: wy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (la_qcabs1(wa) == 0.0_qp) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 wy(i) = wy(i) + wa*wx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 wy(iy) = wy(iy) + wa*wx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_waxpy
#endif

     !> CCOPY: copies a vector x to a vector y.
     pure subroutine la_ccopy(n,cx,incx,cy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*)
           complex(sp),intent(out) :: cy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
               cy(i) = cx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 cy(iy) = cx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_ccopy
     !> ZCOPY: copies a vector, x, to a vector, y.

     pure subroutine la_zcopy(n,zx,incx,zy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*)
           complex(dp),intent(out) :: zy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
               zy(i) = zx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 zy(iy) = zx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_zcopy
#ifdef LA_WITH_XDP
     !> YCOPY: copies a vector, x, to a vector, y.

     pure subroutine la_ycopy(n,yx,incx,yy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(xdp),intent(in) :: yx(*)
           complex(xdp),intent(out) :: yy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
               yy(i) = yx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 yy(iy) = yx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_ycopy
#endif
#ifdef LA_WITH_QP
     !> WCOPY: copies a vector, x, to a vector, y.

     pure subroutine la_wcopy(n,wx,incx,wy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(qp),intent(in) :: wx(*)
           complex(qp),intent(out) :: wy(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
               wy(i) = wx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 wy(iy) = wx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_wcopy
#endif

     !> CDOTC: forms the dot product of two complex vectors
     !> CDOTC = X^H * Y

     pure complex(sp) function la_cdotc(n,cx,incx,cy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*),cy(*)
        ! =====================================================================
           ! Local Scalars
           complex(sp) :: ctemp
           integer(ilp) :: i,ix,iy
           ! Intrinsic Functions
           intrinsic :: conjg
           ctemp = (0.0_sp,0.0_sp)
           la_cdotc = (0.0_sp,0.0_sp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ctemp = ctemp + conjg(cx(i))*cy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = ctemp + conjg(cx(ix))*cy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_cdotc = ctemp
           return
     end function la_cdotc
     !> ZDOTC: forms the dot product of two complex vectors
     !> ZDOTC = X^H * Y

     pure complex(dp) function la_zdotc(n,zx,incx,zy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*),zy(*)
        ! =====================================================================
           ! Local Scalars
           complex(dp) :: ztemp
           integer(ilp) :: i,ix,iy
           ! Intrinsic Functions
           intrinsic :: conjg
           ztemp = (0.0_dp,0.0_dp)
           la_zdotc = (0.0_dp,0.0_dp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ztemp = ztemp + conjg(zx(i))*zy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ztemp = ztemp + conjg(zx(ix))*zy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_zdotc = ztemp
           return
     end function la_zdotc
#ifdef LA_WITH_XDP
     !> YDOTC: forms the dot product of two complex vectors
     !> YDOTC = X^H * Y

     pure complex(xdp) function la_ydotc(n,yx,incx,yy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(xdp),intent(in) :: yx(*),yy(*)
        ! =====================================================================
           ! Local Scalars
           complex(xdp) :: ytemp
           integer(ilp) :: i,ix,iy
           ! Intrinsic Functions
           intrinsic :: conjg
           ytemp = (0.0_xdp,0.0_xdp)
           la_ydotc = (0.0_xdp,0.0_xdp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ytemp = ytemp + conjg(yx(i))*yy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ytemp = ytemp + conjg(yx(ix))*yy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_ydotc = ytemp
           return
     end function la_ydotc
#endif
#ifdef LA_WITH_QP
     !> WDOTC: forms the dot product of two complex vectors
     !> WDOTC = X^H * Y

     pure complex(qp) function la_wdotc(n,wx,incx,wy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(qp),intent(in) :: wx(*),wy(*)
        ! =====================================================================
           ! Local Scalars
           complex(qp) :: wtemp
           integer(ilp) :: i,ix,iy
           ! Intrinsic Functions
           intrinsic :: conjg
           wtemp = (0.0_qp,0.0_qp)
           la_wdotc = (0.0_qp,0.0_qp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 wtemp = wtemp + conjg(wx(i))*wy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 wtemp = wtemp + conjg(wx(ix))*wy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_wdotc = wtemp
           return
     end function la_wdotc
#endif

     !> CDOTU: forms the dot product of two complex vectors
     !> CDOTU = X^T * Y

     pure complex(sp) function la_cdotu(n,cx,incx,cy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(sp),intent(in) :: cx(*),cy(*)
        ! =====================================================================
           ! Local Scalars
           complex(sp) :: ctemp
           integer(ilp) :: i,ix,iy
           ctemp = (0.0_sp,0.0_sp)
           la_cdotu = (0.0_sp,0.0_sp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ctemp = ctemp + cx(i)*cy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = ctemp + cx(ix)*cy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_cdotu = ctemp
           return
     end function la_cdotu
     !> ZDOTU: forms the dot product of two complex vectors
     !> ZDOTU = X^T * Y

     pure complex(dp) function la_zdotu(n,zx,incx,zy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(dp),intent(in) :: zx(*),zy(*)
        ! =====================================================================
           ! Local Scalars
           complex(dp) :: ztemp
           integer(ilp) :: i,ix,iy
           ztemp = (0.0_dp,0.0_dp)
           la_zdotu = (0.0_dp,0.0_dp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ztemp = ztemp + zx(i)*zy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ztemp = ztemp + zx(ix)*zy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_zdotu = ztemp
           return
     end function la_zdotu
#ifdef LA_WITH_XDP
     !> YDOTU: forms the dot product of two complex vectors
     !> YDOTU = X^T * Y

     pure complex(xdp) function la_ydotu(n,yx,incx,yy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(xdp),intent(in) :: yx(*),yy(*)
        ! =====================================================================
           ! Local Scalars
           complex(xdp) :: ytemp
           integer(ilp) :: i,ix,iy
           ytemp = (0.0_xdp,0.0_xdp)
           la_ydotu = (0.0_xdp,0.0_xdp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ytemp = ytemp + yx(i)*yy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ytemp = ytemp + yx(ix)*yy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_ydotu = ytemp
           return
     end function la_ydotu
#endif
#ifdef LA_WITH_QP
     !> WDOTU: forms the dot product of two complex vectors
     !> WDOTU = X^T * Y

     pure complex(qp) function la_wdotu(n,wx,incx,wy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(qp),intent(in) :: wx(*),wy(*)
        ! =====================================================================
           ! Local Scalars
           complex(qp) :: wtemp
           integer(ilp) :: i,ix,iy
           wtemp = (0.0_qp,0.0_qp)
           la_wdotu = (0.0_qp,0.0_qp)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 wtemp = wtemp + wx(i)*wy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 wtemp = wtemp + wx(ix)*wy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           la_wdotu = wtemp
           return
     end function la_wdotu
#endif

     !> CSROT: applies a plane rotation, where the cos and sin (c and s) are real
     !> and the vectors cx and cy are complex.
     !> jack dongarra, linpack, 3/11/78.

     pure subroutine la_csrot(n,cx,incx,cy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(sp),intent(in) :: c,s
           ! Array Arguments
           complex(sp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(sp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ctemp = c*cx(i) + s*cy(i)
                 cy(i) = c*cy(i) - s*cx(i)
                 cx(i) = ctemp
              end do
           else
              ! code for unequal increments or equal increments not equal
                ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = c*cx(ix) + s*cy(iy)
                 cy(iy) = c*cy(iy) - s*cx(ix)
                 cx(ix) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_csrot
     !> Applies a plane rotation, where the cos and sin (c and s) are real
     !> and the vectors cx and cy are complex.
     !> jack dongarra, linpack, 3/11/78.

     pure subroutine la_zdrot(n,zx,incx,zy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(dp),intent(in) :: c,s
           ! Array Arguments
           complex(dp),intent(inout) :: zx(*),zy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(dp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ctemp = c*zx(i) + s*zy(i)
                 zy(i) = c*zy(i) - s*zx(i)
                 zx(i) = ctemp
              end do
           else
              ! code for unequal increments or equal increments not equal
                ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = c*zx(ix) + s*zy(iy)
                 zy(iy) = c*zy(iy) - s*zx(ix)
                 zx(ix) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_zdrot
#ifdef LA_WITH_XDP
     !> Applies a plane rotation, where the cos and sin (c and s) are real
     !> and the vectors cx and cy are complex.
     !> jack dongarra, linpack, 3/11/78.

     pure subroutine la_yxrot(n,yx,incx,yy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(xdp),intent(in) :: c,s
           ! Array Arguments
           complex(xdp),intent(inout) :: yx(*),yy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(xdp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ctemp = c*yx(i) + s*yy(i)
                 yy(i) = c*yy(i) - s*yx(i)
                 yx(i) = ctemp
              end do
           else
              ! code for unequal increments or equal increments not equal
                ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = c*yx(ix) + s*yy(iy)
                 yy(iy) = c*yy(iy) - s*yx(ix)
                 yx(ix) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_yxrot
#endif
#ifdef LA_WITH_QP
     !> Applies a plane rotation, where the cos and sin (c and s) are real
     !> and the vectors cx and cy are complex.
     !> jack dongarra, linpack, 3/11/78.

     pure subroutine la_wqrot(n,wx,incx,wy,incy,c,s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           real(qp),intent(in) :: c,s
           ! Array Arguments
           complex(qp),intent(inout) :: wx(*),wy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(qp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1,n
                 ctemp = c*wx(i) + s*wy(i)
                 wy(i) = c*wy(i) - s*wx(i)
                 wx(i) = ctemp
              end do
           else
              ! code for unequal increments or equal increments not equal
                ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = c*wx(ix) + s*wy(iy)
                 wy(iy) = c*wy(iy) - s*wx(ix)
                 wx(ix) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_wqrot
#endif

     !> CSSCAL: scales a complex vector by a real constant.

     pure subroutine la_csscal(n,sa,cx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: sa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(inout) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           ! Intrinsic Functions
           intrinsic :: aimag,cmplx,real
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 cx(i) = cmplx(sa*real(cx(i),KIND=sp),sa*aimag(cx(i)),KIND=sp)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 cx(i) = cmplx(sa*real(cx(i),KIND=sp),sa*aimag(cx(i)),KIND=sp)
              end do
           end if
           return
     end subroutine la_csscal
     !> ZDSCAL: scales a vector by a constant.

     pure subroutine la_zdscal(n,da,zx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: da
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(inout) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           ! Intrinsic Functions
           intrinsic :: cmplx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 zx(i) = cmplx(da,0.0_dp,KIND=dp)*zx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 zx(i) = cmplx(da,0.0_dp,KIND=dp)*zx(i)
              end do
           end if
           return
     end subroutine la_zdscal
#ifdef LA_WITH_XDP
     !> YXSCAL: scales a vector by a constant.

     pure subroutine la_yxscal(n,xa,yx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: xa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(xdp),intent(inout) :: yx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           ! Intrinsic Functions
           intrinsic :: cmplx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 yx(i) = cmplx(xa,0.0_xdp,KIND=xdp)*yx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 yx(i) = cmplx(xa,0.0_xdp,KIND=xdp)*yx(i)
              end do
           end if
           return
     end subroutine la_yxscal
#endif
#ifdef LA_WITH_QP
     !> WQSCAL: scales a vector by a constant.

     pure subroutine la_wqscal(n,qa,wx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: qa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(inout) :: wx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           ! Intrinsic Functions
           intrinsic :: cmplx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 wx(i) = cmplx(qa,0.0_qp,KIND=qp)*wx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 wx(i) = cmplx(qa,0.0_qp,KIND=qp)*wx(i)
              end do
           end if
           return
     end subroutine la_wqscal
#endif

     !> !
     !>
     !> The computation uses the formulas
     !> |x| = sqrt( Re(x)**2 + Im(x)**2 )
     !> sgn(x) = x / |x|  if x /= 0
     !> = 1        if x  = 0
     !> c = |a| / sqrt(|a|**2 + |b|**2)
     !> s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
     !> When a and b are real and r /= 0, the formulas simplify to
     !> r = sgn(a)*sqrt(|a|**2 + |b|**2)
     !> c = a / r
     !> s = b / r
     !> the same as in SROTG when |a| > |b|.  When |b| >= |a|, the
     !> sign of c and s will be different from those computed by SROTG
     !> if the signs of a and b are not the same.

     pure subroutine la_crotg(a,b,c,s)
        use la_constants_sp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._sp)
        ! Scaling Constants
        ! Scalar Arguments
        real(sp),intent(out) :: c
        complex(sp),intent(inout) :: a
        complex(sp),intent(in) :: b
        complex(sp),intent(out) :: s
        ! Local Scalars
        real(sp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(sp) :: f,fs,g,gs,r,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(sp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=sp)**2 + aimag(t)**2
        ! Executable Statements
        f = a
        g = b
        if (g == czero) then
           c = one
           s = czero
           r = f
        else if (f == czero) then
           c = zero
           g1 = max(abs(real(g,KIND=sp)),abs(aimag(g)))
           if (g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              g2 = abssq(g)
              d = sqrt(g2)
              s = conjg(g)/d
              r = d
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              d = sqrt(g2)
              s = conjg(gs)/d
              r = d*u
           end if
        else
           f1 = max(abs(real(f,KIND=sp)),abs(aimag(f)))
           g1 = max(abs(real(g,KIND=sp)),abs(aimag(g)))
     if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              f2 = abssq(f)
              g2 = abssq(g)
              h2 = f2 + g2
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = f2*p
              s = conjg(g)*(f*p)
              r = f*(h2*p)
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,f1,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              if (f1*uu < rtmin) then
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
                 v = min(safmax,max(safmin,f1))
                 vv = one/v
                 w = v*uu
                 fs = f*vv
                 f2 = abssq(fs)
                 h2 = f2*w**2 + g2
              else
                 ! otherwise use the same scaling for f and g.
                 w = one
                 fs = f*uu
                 f2 = abssq(fs)
                 h2 = f2 + g2
              end if
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = (f2*p)*w
              s = conjg(gs)*(fs*p)
              r = (fs*(h2*p))*u
           end if
        end if
        a = r
        return
     end subroutine la_crotg
     !> !
     !>
     !> The computation uses the formulas
     !> |x| = sqrt( Re(x)**2 + Im(x)**2 )
     !> sgn(x) = x / |x|  if x /= 0
     !> = 1        if x  = 0
     !> c = |a| / sqrt(|a|**2 + |b|**2)
     !> s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
     !> When a and b are real and r /= 0, the formulas simplify to
     !> r = sgn(a)*sqrt(|a|**2 + |b|**2)
     !> c = a / r
     !> s = b / r
     !> the same as in DROTG when |a| > |b|.  When |b| >= |a|, the
     !> sign of c and s will be different from those computed by DROTG
     !> if the signs of a and b are not the same.

     pure subroutine la_zrotg(a,b,c,s)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._dp)
        ! Scaling Constants
        ! Scalar Arguments
        real(dp),intent(out) :: c
        complex(dp),intent(inout) :: a
        complex(dp),intent(in) :: b
        complex(dp),intent(out) :: s
        ! Local Scalars
        real(dp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(dp) :: f,fs,g,gs,r,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(dp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=dp)**2 + aimag(t)**2
        ! Executable Statements
        f = a
        g = b
        if (g == czero) then
           c = one
           s = czero
           r = f
        else if (f == czero) then
           c = zero
           g1 = max(abs(real(g,KIND=dp)),abs(aimag(g)))
           if (g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              g2 = abssq(g)
              d = sqrt(g2)
              s = conjg(g)/d
              r = d
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              d = sqrt(g2)
              s = conjg(gs)/d
              r = d*u
           end if
        else
           f1 = max(abs(real(f,KIND=dp)),abs(aimag(f)))
           g1 = max(abs(real(g,KIND=dp)),abs(aimag(g)))
     if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              f2 = abssq(f)
              g2 = abssq(g)
              h2 = f2 + g2
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = f2*p
              s = conjg(g)*(f*p)
              r = f*(h2*p)
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,f1,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              if (f1*uu < rtmin) then
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
                 v = min(safmax,max(safmin,f1))
                 vv = one/v
                 w = v*uu
                 fs = f*vv
                 f2 = abssq(fs)
                 h2 = f2*w**2 + g2
              else
                 ! otherwise use the same scaling for f and g.
                 w = one
                 fs = f*uu
                 f2 = abssq(fs)
                 h2 = f2 + g2
              end if
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = (f2*p)*w
              s = conjg(gs)*(fs*p)
              r = (fs*(h2*p))*u
           end if
        end if
        a = r
        return
     end subroutine la_zrotg
#ifdef LA_WITH_XDP
     !> !
     !>
     !> The computation uses the formulas
     !> |x| = sqrt( Re(x)**2 + Im(x)**2 )
     !> sgn(x) = x / |x|  if x /= 0
     !> = 1        if x  = 0
     !> c = |a| / sqrt(|a|**2 + |b|**2)
     !> s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
     !> When a and b are real and r /= 0, the formulas simplify to
     !> r = sgn(a)*sqrt(|a|**2 + |b|**2)
     !> c = a / r
     !> s = b / r
     !> the same as in XROTG when |a| > |b|.  When |b| >= |a|, the
     !> sign of c and s will be different from those computed by XROTG
     !> if the signs of a and b are not the same.

     pure subroutine la_yrotg(a,b,c,s)
        use la_constants_xdp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._xdp)
        ! Scaling Constants
        ! Scalar Arguments
        real(xdp),intent(out) :: c
        complex(xdp),intent(inout) :: a
        complex(xdp),intent(in) :: b
        complex(xdp),intent(out) :: s
        ! Local Scalars
        real(xdp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(xdp) :: f,fs,g,gs,r,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(xdp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=xdp)**2 + aimag(t)**2
        ! Executable Statements
        f = a
        g = b
        if (g == czero) then
           c = one
           s = czero
           r = f
        else if (f == czero) then
           c = zero
           g1 = max(abs(real(g,KIND=xdp)),abs(aimag(g)))
           if (g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              g2 = abssq(g)
              d = sqrt(g2)
              s = conjg(g)/d
              r = d
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              d = sqrt(g2)
              s = conjg(gs)/d
              r = d*u
           end if
        else
           f1 = max(abs(real(f,KIND=xdp)),abs(aimag(f)))
           g1 = max(abs(real(g,KIND=xdp)),abs(aimag(g)))
     if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              f2 = abssq(f)
              g2 = abssq(g)
              h2 = f2 + g2
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = f2*p
              s = conjg(g)*(f*p)
              r = f*(h2*p)
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,f1,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              if (f1*uu < rtmin) then
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
                 v = min(safmax,max(safmin,f1))
                 vv = one/v
                 w = v*uu
                 fs = f*vv
                 f2 = abssq(fs)
                 h2 = f2*w**2 + g2
              else
                 ! otherwise use the same scaling for f and g.
                 w = one
                 fs = f*uu
                 f2 = abssq(fs)
                 h2 = f2 + g2
              end if
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = (f2*p)*w
              s = conjg(gs)*(fs*p)
              r = (fs*(h2*p))*u
           end if
        end if
        a = r
        return
     end subroutine la_yrotg
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> The computation uses the formulas
     !> |x| = sqrt( Re(x)**2 + Im(x)**2 )
     !> sgn(x) = x / |x|  if x /= 0
     !> = 1        if x  = 0
     !> c = |a| / sqrt(|a|**2 + |b|**2)
     !> s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
     !> When a and b are real and r /= 0, the formulas simplify to
     !> r = sgn(a)*sqrt(|a|**2 + |b|**2)
     !> c = a / r
     !> s = b / r
     !> the same as in QROTG when |a| > |b|.  When |b| >= |a|, the
     !> sign of c and s will be different from those computed by QROTG
     !> if the signs of a and b are not the same.

     pure subroutine la_wrotg(a,b,c,s)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Constants
        integer,parameter :: wp = kind(1._qp)
        ! Scaling Constants
        ! Scalar Arguments
        real(qp),intent(out) :: c
        complex(qp),intent(inout) :: a
        complex(qp),intent(in) :: b
        complex(qp),intent(out) :: s
        ! Local Scalars
        real(qp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(qp) :: f,fs,g,gs,r,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(qp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=qp)**2 + aimag(t)**2
        ! Executable Statements
        f = a
        g = b
        if (g == czero) then
           c = one
           s = czero
           r = f
        else if (f == czero) then
           c = zero
           g1 = max(abs(real(g,KIND=qp)),abs(aimag(g)))
           if (g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              g2 = abssq(g)
              d = sqrt(g2)
              s = conjg(g)/d
              r = d
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              d = sqrt(g2)
              s = conjg(gs)/d
              r = d*u
           end if
        else
           f1 = max(abs(real(f,KIND=qp)),abs(aimag(f)))
           g1 = max(abs(real(g,KIND=qp)),abs(aimag(g)))
     if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              f2 = abssq(f)
              g2 = abssq(g)
              h2 = f2 + g2
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = f2*p
              s = conjg(g)*(f*p)
              r = f*(h2*p)
           else
              ! use scaled algorithm
              u = min(safmax,max(safmin,f1,g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              if (f1*uu < rtmin) then
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
                 v = min(safmax,max(safmin,f1))
                 vv = one/v
                 w = v*uu
                 fs = f*vv
                 f2 = abssq(fs)
                 h2 = f2*w**2 + g2
              else
                 ! otherwise use the same scaling for f and g.
                 w = one
                 fs = f*uu
                 f2 = abssq(fs)
                 h2 = f2 + g2
              end if
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = (f2*p)*w
              s = conjg(gs)*(fs*p)
              r = (fs*(h2*p))*u
           end if
        end if
        a = r
        return
     end subroutine la_wrotg
#endif

     !> CSCAL: scales a vector by a constant.

     pure subroutine la_cscal(n,ca,cx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: ca
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(inout) :: cx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 cx(i) = ca*cx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 cx(i) = ca*cx(i)
              end do
           end if
           return
     end subroutine la_cscal
     !> ZSCAL: scales a vector by a constant.

     pure subroutine la_zscal(n,za,zx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: za
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(inout) :: zx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 zx(i) = za*zx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 zx(i) = za*zx(i)
              end do
           end if
           return
     end subroutine la_zscal
#ifdef LA_WITH_XDP
     !> YSCAL: scales a vector by a constant.

     pure subroutine la_yscal(n,ya,yx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: ya
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(xdp),intent(inout) :: yx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 yx(i) = ya*yx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 yx(i) = ya*yx(i)
              end do
           end if
           return
     end subroutine la_yscal
#endif
#ifdef LA_WITH_QP
     !> WSCAL: scales a vector by a constant.

     pure subroutine la_wscal(n,wa,wx,incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: wa
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(inout) :: wx(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,nincx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1,n
                 wx(i) = wa*wx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1,nincx,incx
                 wx(i) = wa*wx(i)
              end do
           end if
           return
     end subroutine la_wscal
#endif

     !> CSWAP: interchanges two vectors.

     pure subroutine la_cswap(n,cx,incx,cy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(sp),intent(inout) :: cx(*),cy(*)
        ! =====================================================================
           ! Local Scalars
           complex(sp) :: ctemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 ctemp = cx(i)
                 cx(i) = cy(i)
                 cy(i) = ctemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ctemp = cx(ix)
                 cx(ix) = cy(iy)
                 cy(iy) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_cswap
     !> ZSWAP: interchanges two vectors.

     pure subroutine la_zswap(n,zx,incx,zy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(dp),intent(inout) :: zx(*),zy(*)
        ! =====================================================================
           ! Local Scalars
           complex(dp) :: ztemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 ztemp = zx(i)
                 zx(i) = zy(i)
                 zy(i) = ztemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ztemp = zx(ix)
                 zx(ix) = zy(iy)
                 zy(iy) = ztemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_zswap
#ifdef LA_WITH_XDP
     !> YSWAP: interchanges two vectors.

     pure subroutine la_yswap(n,yx,incx,yy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(xdp),intent(inout) :: yx(*),yy(*)
        ! =====================================================================
           ! Local Scalars
           complex(xdp) :: ytemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 ytemp = yx(i)
                 yx(i) = yy(i)
                 yy(i) = ytemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 ytemp = yx(ix)
                 yx(ix) = yy(iy)
                 yy(iy) = ytemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_yswap
#endif
#ifdef LA_WITH_QP
     !> WSWAP: interchanges two vectors.

     pure subroutine la_wswap(n,wx,incx,wy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           complex(qp),intent(inout) :: wx(*),wy(*)
        ! =====================================================================
           ! Local Scalars
           complex(qp) :: wtemp
           integer(ilp) :: i,ix,iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1,n
                 wtemp = wx(i)
                 wx(i) = wy(i)
                 wy(i) = wtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1,n
                 wtemp = wx(ix)
                 wx(ix) = wy(iy)
                 wy(iy) = wtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     end subroutine la_wswap
#endif

     !> Compute the inner product of two vectors with extended
     !> precision accumulation.
     !> Returns S.P. result with dot product accumulated in D.P.
     !> SDSDOT: = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
     !> where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
     !> defined in a similar way using INCY.

     pure real(sp) function la_sdsdot(n,sb,sx,incx,sy,incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: sb
           integer(ilp),intent(in) :: incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: sx(*),sy(*)
           ! Local Scalars
           real(dp) :: dsdot
           integer(ilp) :: i,kx,ky,ns
           ! Intrinsic Functions
           intrinsic :: real
           dsdot = sb
           if (n <= 0) then
              la_sdsdot = dsdot
              return
           end if
           if (incx == incy .and. incx > 0) then
           ! code for equal and positive increments.
              ns = n*incx
              do i = 1,ns,incx
                 dsdot = dsdot + real(sx(i),KIND=sp)*real(sy(i),KIND=sp)
              end do
           else
           ! code for unequal or nonpositive increments.
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              do i = 1,n
                 dsdot = dsdot + real(sx(kx),KIND=sp)*real(sy(ky),KIND=sp)
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           la_sdsdot = dsdot
           return
     end function la_sdsdot

end module la_blas_level1
