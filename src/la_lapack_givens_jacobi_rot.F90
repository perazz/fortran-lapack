!> Givens and Jacobi plane rotations
module la_lapack_givens_jacobi_rot
     use la_constants
     use la_blas_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_scalar
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slar2v
     public :: la_slargv
     public :: la_slartg
     public :: la_slartgp
     public :: la_slartv
     public :: la_slasr
     public :: la_dlar2v
     public :: la_dlargv
     public :: la_dlartg
     public :: la_dlartgp
     public :: la_dlartv
     public :: la_dlasr
#ifdef LA_WITH_XDP
     public :: la_xlar2v
     public :: la_xlargv
     public :: la_xlartg
     public :: la_xlartgp
     public :: la_xlartv
     public :: la_xlasr
#endif
#ifdef LA_WITH_QP
     public :: la_qlar2v
     public :: la_qlargv
     public :: la_qlartg
     public :: la_qlartgp
     public :: la_qlartv
     public :: la_qlasr
#endif
     public :: la_clacrt
     public :: la_clar2v
     public :: la_clartg
     public :: la_clartv
     public :: la_clasr
     public :: la_clargv
     public :: la_zlacrt
     public :: la_zlar2v
     public :: la_zlartg
     public :: la_zlartv
     public :: la_zlasr
     public :: la_zlargv
#ifdef LA_WITH_XDP
     public :: la_ylacrt
     public :: la_ylar2v
     public :: la_ylartg
     public :: la_ylartv
     public :: la_ylasr
     public :: la_ylargv
#endif
#ifdef LA_WITH_QP
     public :: la_wlacrt
     public :: la_wlar2v
     public :: la_wlartg
     public :: la_wlartv
     public :: la_wlasr
     public :: la_wlargv
#endif

     contains

     !> SLAR2V: applies a vector of real plane rotations from both sides to
     !> a sequence of 2-by-2 real symmetric matrices, defined by the elements
     !> of the vectors x, y and z. For i = 1,2,...,n
     !> ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
     !> ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )

     pure subroutine la_slar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(sp),intent(in) :: c(*),s(*)
           real(sp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(sp) :: ci,si,t1,t2,t3,t4,t5,t6,xi,yi,zi
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(ix)
              zi = z(ix)
              ci = c(ic)
              si = s(ic)
              t1 = si*zi
              t2 = ci*zi
              t3 = t2 - si*xi
              t4 = t2 + si*yi
              t5 = ci*xi + t1
              t6 = ci*yi - t1
              x(ix) = ci*t5 + si*t4
              y(ix) = ci*t6 - si*t3
              z(ix) = ci*t4 - si*t5
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_slar2v
     !> DLAR2V: applies a vector of real plane rotations from both sides to
     !> a sequence of 2-by-2 real symmetric matrices, defined by the elements
     !> of the vectors x, y and z. For i = 1,2,...,n
     !> ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
     !> ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )

     pure subroutine la_dlar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(dp),intent(in) :: c(*),s(*)
           real(dp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(dp) :: ci,si,t1,t2,t3,t4,t5,t6,xi,yi,zi
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(ix)
              zi = z(ix)
              ci = c(ic)
              si = s(ic)
              t1 = si*zi
              t2 = ci*zi
              t3 = t2 - si*xi
              t4 = t2 + si*yi
              t5 = ci*xi + t1
              t6 = ci*yi - t1
              x(ix) = ci*t5 + si*t4
              y(ix) = ci*t6 - si*t3
              z(ix) = ci*t4 - si*t5
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_dlar2v
#ifdef LA_WITH_XDP
     !> XLAR2V: applies a vector of real plane rotations from both sides to
     !> a sequence of 2-by-2 real symmetric matrices, defined by the elements
     !> of the vectors x, y and z. For i = 1,2,...,n
     !> ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
     !> ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )

     pure subroutine la_xlar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(xdp),intent(in) :: c(*),s(*)
           real(xdp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(xdp) :: ci,si,t1,t2,t3,t4,t5,t6,xi,yi,zi
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(ix)
              zi = z(ix)
              ci = c(ic)
              si = s(ic)
              t1 = si*zi
              t2 = ci*zi
              t3 = t2 - si*xi
              t4 = t2 + si*yi
              t5 = ci*xi + t1
              t6 = ci*yi - t1
              x(ix) = ci*t5 + si*t4
              y(ix) = ci*t6 - si*t3
              z(ix) = ci*t4 - si*t5
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_xlar2v
#endif
#ifdef LA_WITH_QP
     !> QLAR2V: applies a vector of real plane rotations from both sides to
     !> a sequence of 2-by-2 real symmetric matrices, defined by the elements
     !> of the vectors x, y and z. For i = 1,2,...,n
     !> ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
     !> ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )

     pure subroutine la_qlar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(qp),intent(in) :: c(*),s(*)
           real(qp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(qp) :: ci,si,t1,t2,t3,t4,t5,t6,xi,yi,zi
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(ix)
              zi = z(ix)
              ci = c(ic)
              si = s(ic)
              t1 = si*zi
              t2 = ci*zi
              t3 = t2 - si*xi
              t4 = t2 + si*yi
              t5 = ci*xi + t1
              t6 = ci*yi - t1
              x(ix) = ci*t5 + si*t4
              y(ix) = ci*t6 - si*t3
              z(ix) = ci*t4 - si*t5
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_qlar2v
#endif

     !> SLARGV: generates a vector of real plane rotations, determined by
     !> elements of the real vectors x and y. For i = 1,2,...,n
     !> (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
     !> ( -s(i)  c(i) ) ( y(i) ) = (   0  )

     pure subroutine la_slargv(n,x,incx,y,incy,c,incc)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(sp),intent(out) :: c(*)
           real(sp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(sp) :: f,g,t,tt
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           loop_10: do i = 1,n
              f = x(ix)
              g = y(iy)
              if (g == zero) then
                 c(ic) = one
              else if (f == zero) then
                 c(ic) = zero
                 y(iy) = one
                 x(ix) = g
              else if (abs(f) > abs(g)) then
                 t = g/f
                 tt = sqrt(one + t*t)
                 c(ic) = one/tt
                 y(iy) = t*c(ic)
                 x(ix) = f*tt
              else
                 t = f/g
                 tt = sqrt(one + t*t)
                 y(iy) = one/tt
                 c(ic) = t*y(iy)
                 x(ix) = g*tt
              end if
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_10
           return
     end subroutine la_slargv
     !> DLARGV: generates a vector of real plane rotations, determined by
     !> elements of the real vectors x and y. For i = 1,2,...,n
     !> (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
     !> ( -s(i)  c(i) ) ( y(i) ) = (   0  )

     pure subroutine la_dlargv(n,x,incx,y,incy,c,incc)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(dp),intent(out) :: c(*)
           real(dp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(dp) :: f,g,t,tt
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           loop_10: do i = 1,n
              f = x(ix)
              g = y(iy)
              if (g == zero) then
                 c(ic) = one
              else if (f == zero) then
                 c(ic) = zero
                 y(iy) = one
                 x(ix) = g
              else if (abs(f) > abs(g)) then
                 t = g/f
                 tt = sqrt(one + t*t)
                 c(ic) = one/tt
                 y(iy) = t*c(ic)
                 x(ix) = f*tt
              else
                 t = f/g
                 tt = sqrt(one + t*t)
                 y(iy) = one/tt
                 c(ic) = t*y(iy)
                 x(ix) = g*tt
              end if
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_10
           return
     end subroutine la_dlargv
#ifdef LA_WITH_XDP
     !> XLARGV: generates a vector of real plane rotations, determined by
     !> elements of the real vectors x and y. For i = 1,2,...,n
     !> (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
     !> ( -s(i)  c(i) ) ( y(i) ) = (   0  )

     pure subroutine la_xlargv(n,x,incx,y,incy,c,incc)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(xdp),intent(out) :: c(*)
           real(xdp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(xdp) :: f,g,t,tt
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           loop_10: do i = 1,n
              f = x(ix)
              g = y(iy)
              if (g == zero) then
                 c(ic) = one
              else if (f == zero) then
                 c(ic) = zero
                 y(iy) = one
                 x(ix) = g
              else if (abs(f) > abs(g)) then
                 t = g/f
                 tt = sqrt(one + t*t)
                 c(ic) = one/tt
                 y(iy) = t*c(ic)
                 x(ix) = f*tt
              else
                 t = f/g
                 tt = sqrt(one + t*t)
                 y(iy) = one/tt
                 c(ic) = t*y(iy)
                 x(ix) = g*tt
              end if
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_10
           return
     end subroutine la_xlargv
#endif
#ifdef LA_WITH_QP
     !> QLARGV: generates a vector of real plane rotations, determined by
     !> elements of the real vectors x and y. For i = 1,2,...,n
     !> (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
     !> ( -s(i)  c(i) ) ( y(i) ) = (   0  )

     pure subroutine la_qlargv(n,x,incx,y,incy,c,incc)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(qp),intent(out) :: c(*)
           real(qp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(qp) :: f,g,t,tt
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           loop_10: do i = 1,n
              f = x(ix)
              g = y(iy)
              if (g == zero) then
                 c(ic) = one
              else if (f == zero) then
                 c(ic) = zero
                 y(iy) = one
                 x(ix) = g
              else if (abs(f) > abs(g)) then
                 t = g/f
                 tt = sqrt(one + t*t)
                 c(ic) = one/tt
                 y(iy) = t*c(ic)
                 x(ix) = f*tt
              else
                 t = f/g
                 tt = sqrt(one + t*t)
                 y(iy) = one/tt
                 c(ic) = t*y(iy)
                 x(ix) = g*tt
              end if
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_10
           return
     end subroutine la_qlargv
#endif

     !> !
     !>
     !> SLARTG: generates a plane rotation so that
     !> [  C  S  ]  .  [ F ]  =  [ R ]
     !> [ -S  C  ]     [ G ]     [ 0 ]
     !> where C**2 + S**2 = 1.
     !> The mathematical formulas used for C and S are
     !> R = sign(F) * sqrt(F**2 + G**2)
     !> C = F / R
     !> S = G / R
     !> Hence C >= 0. The algorithm used to compute these quantities
     !> incorporates scaling to avoid overflow or underflow in computing the
     !> square root of the sum of squares.
     !> This version is discontinuous in R at F = 0 but it returns the same
     !> C and S as SLARTG for complex inputs (F,0) and (G,0).
     !> This is a more accurate version of the BLAS1 routine SROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0 and (G .ne. 0), then C=0 and S=sign(1,G) without doing any
     !> floating point operations (saves work in SBDSQR when
     !> there are zeros on the diagonal).
     !> If F exceeds G in magnitude, C will be positive.
     !> Below, wp=>sp stands for single precision from LA_CONSTANTS module.

     pure subroutine la_slartg(f,g,c,s,r)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(sp),intent(out) :: c,r,s
        real(sp),intent(in) :: f,g
        ! Local Scalars
        real(sp) :: d,f1,fs,g1,gs,p,u,uu
        ! Intrinsic Functions
        intrinsic :: abs,sign,sqrt
        ! Executable Statements
        f1 = abs(f)
        g1 = abs(g)
        if (g == zero) then
           c = one
           s = zero
           r = f
        else if (f == zero) then
           c = zero
           s = sign(one,g)
           r = g1
     else if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) &
               then
           d = sqrt(f*f + g*g)
           p = one/d
           c = f1*p
           s = g*sign(p,f)
           r = sign(d,f)
        else
           u = min(safmax,max(safmin,f1,g1))
           uu = one/u
           fs = f*uu
           gs = g*uu
           d = sqrt(fs*fs + gs*gs)
           p = one/d
           c = abs(fs)*p
           s = gs*sign(p,f)
           r = sign(d,f)*u
        end if
        return
     end subroutine la_slartg
     !> !
     !>
     !> DLARTG: generates a plane rotation so that
     !> [  C  S  ]  .  [ F ]  =  [ R ]
     !> [ -S  C  ]     [ G ]     [ 0 ]
     !> where C**2 + S**2 = 1.
     !> The mathematical formulas used for C and S are
     !> R = sign(F) * sqrt(F**2 + G**2)
     !> C = F / R
     !> S = G / R
     !> Hence C >= 0. The algorithm used to compute these quantities
     !> incorporates scaling to avoid overflow or underflow in computing the
     !> square root of the sum of squares.
     !> This version is discontinuous in R at F = 0 but it returns the same
     !> C and S as ZLARTG for complex inputs (F,0) and (G,0).
     !> This is a more accurate version of the BLAS1 routine DROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0 and (G .ne. 0), then C=0 and S=sign(1,G) without doing any
     !> floating point operations (saves work in DBDSQR when
     !> there are zeros on the diagonal).
     !> If F exceeds G in magnitude, C will be positive.
     !> Below, wp=>dp stands for double precision from LA_CONSTANTS module.

     pure subroutine la_dlartg(f,g,c,s,r)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(dp),intent(out) :: c,r,s
        real(dp),intent(in) :: f,g
        ! Local Scalars
        real(dp) :: d,f1,fs,g1,gs,p,u,uu
        ! Intrinsic Functions
        intrinsic :: abs,sign,sqrt
        ! Executable Statements
        f1 = abs(f)
        g1 = abs(g)
        if (g == zero) then
           c = one
           s = zero
           r = f
        else if (f == zero) then
           c = zero
           s = sign(one,g)
           r = g1
     else if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) &
               then
           d = sqrt(f*f + g*g)
           p = one/d
           c = f1*p
           s = g*sign(p,f)
           r = sign(d,f)
        else
           u = min(safmax,max(safmin,f1,g1))
           uu = one/u
           fs = f*uu
           gs = g*uu
           d = sqrt(fs*fs + gs*gs)
           p = one/d
           c = abs(fs)*p
           s = gs*sign(p,f)
           r = sign(d,f)*u
        end if
        return
     end subroutine la_dlartg
#ifdef LA_WITH_XDP
     !> !
     !>
     !> XLARTG: generates a plane rotation so that
     !> [  C  S  ]  .  [ F ]  =  [ R ]
     !> [ -S  C  ]     [ G ]     [ 0 ]
     !> where C**2 + S**2 = 1.
     !> The mathematical formulas used for C and S are
     !> R = sign(F) * sqrt(F**2 + G**2)
     !> C = F / R
     !> S = G / R
     !> Hence C >= 0. The algorithm used to compute these quantities
     !> incorporates scaling to avoid overflow or underflow in computing the
     !> square root of the sum of squares.
     !> This version is discontinuous in R at F = 0 but it returns the same
     !> C and S as YLARTG for complex inputs (F,0) and (G,0).
     !> This is a more accurate version of the BLAS1 routine XROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0 and (G .ne. 0), then C=0 and S=sign(1,G) without doing any
     !> floating point operations (saves work in XBDSQR when
     !> there are zeros on the diagonal).
     !> If F exceeds G in magnitude, C will be positive.
     !> Below, wp=>xdp stands for extended precision from LA_CONSTANTS module.

     pure subroutine la_xlartg(f,g,c,s,r)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(xdp),intent(out) :: c,r,s
        real(xdp),intent(in) :: f,g
        ! Local Scalars
        real(xdp) :: d,f1,fs,g1,gs,p,u,uu
        ! Intrinsic Functions
        intrinsic :: abs,sign,sqrt
        ! Executable Statements
        f1 = abs(f)
        g1 = abs(g)
        if (g == zero) then
           c = one
           s = zero
           r = f
        else if (f == zero) then
           c = zero
           s = sign(one,g)
           r = g1
     else if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) &
               then
           d = sqrt(f*f + g*g)
           p = one/d
           c = f1*p
           s = g*sign(p,f)
           r = sign(d,f)
        else
           u = min(safmax,max(safmin,f1,g1))
           uu = one/u
           fs = f*uu
           gs = g*uu
           d = sqrt(fs*fs + gs*gs)
           p = one/d
           c = abs(fs)*p
           s = gs*sign(p,f)
           r = sign(d,f)*u
        end if
        return
     end subroutine la_xlartg
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> QLARTG: generates a plane rotation so that
     !> [  C  S  ]  .  [ F ]  =  [ R ]
     !> [ -S  C  ]     [ G ]     [ 0 ]
     !> where C**2 + S**2 = 1.
     !> The mathematical formulas used for C and S are
     !> R = sign(F) * sqrt(F**2 + G**2)
     !> C = F / R
     !> S = G / R
     !> Hence C >= 0. The algorithm used to compute these quantities
     !> incorporates scaling to avoid overflow or underflow in computing the
     !> square root of the sum of squares.
     !> This version is discontinuous in R at F = 0 but it returns the same
     !> C and S as WLARTG for complex inputs (F,0) and (G,0).
     !> This is a more accurate version of the BLAS1 routine QROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0 and (G .ne. 0), then C=0 and S=sign(1,G) without doing any
     !> floating point operations (saves work in QBDSQR when
     !> there are zeros on the diagonal).
     !> If F exceeds G in magnitude, C will be positive.
     !> Below, wp=>qp stands for quad precision from LA_CONSTANTS module.

     pure subroutine la_qlartg(f,g,c,s,r)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(qp),intent(out) :: c,r,s
        real(qp),intent(in) :: f,g
        ! Local Scalars
        real(qp) :: d,f1,fs,g1,gs,p,u,uu
        ! Intrinsic Functions
        intrinsic :: abs,sign,sqrt
        ! Executable Statements
        f1 = abs(f)
        g1 = abs(g)
        if (g == zero) then
           c = one
           s = zero
           r = f
        else if (f == zero) then
           c = zero
           s = sign(one,g)
           r = g1
     else if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) &
               then
           d = sqrt(f*f + g*g)
           p = one/d
           c = f1*p
           s = g*sign(p,f)
           r = sign(d,f)
        else
           u = min(safmax,max(safmin,f1,g1))
           uu = one/u
           fs = f*uu
           gs = g*uu
           d = sqrt(fs*fs + gs*gs)
           p = one/d
           c = abs(fs)*p
           s = gs*sign(p,f)
           r = sign(d,f)*u
        end if
        return
     end subroutine la_qlartg
#endif

     !> SLARTGP: generates a plane rotation so that
     !> [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
     !> [ -SN  CS  ]     [ G ]     [ 0 ]
     !> This is a slower, more accurate version of the Level 1 BLAS routine SROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then CS=(+/-)1 and SN=0.
     !> If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.
     !> The sign is chosen so that R >= 0.

     pure subroutine la_slartgp(f,g,cs,sn,r)
        use la_constants_sp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(out) :: cs,r,sn
           real(sp),intent(in) :: f,g
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i
           real(sp) :: eps,f1,g1,safmin,safmn2,safmx2,scale
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,sign,sqrt
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Executable Statements
           ! if( first ) then
              safmin = la_slamch('S')
              eps = la_slamch('E')
              safmn2 = la_slamch('B')**int(log(safmin/eps)/log(la_slamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
              ! first = .false.
           ! end if
           if (g == zero) then
              cs = sign(one,f)
              sn = zero
              r = abs(f)
           else if (f == zero) then
              cs = zero
              sn = sign(one,g)
              r = abs(g)
           else
              f1 = f
              g1 = g
              scale = max(abs(f1),abs(g1))
              if (scale >= safmx2) then
                 count = 0
                 10 continue
                 count = count + 1
                 f1 = f1*safmn2
                 g1 = g1*safmn2
                 scale = max(abs(f1),abs(g1))
                 if (scale >= safmx2 .and. count < 20) go to 10
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmx2
                 end do
              else if (scale <= safmn2) then
                 count = 0
                 30 continue
                 count = count + 1
                 f1 = f1*safmx2
                 g1 = g1*safmx2
                 scale = max(abs(f1),abs(g1))
                 if (scale <= safmn2) go to 30
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmn2
                 end do
              else
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
              end if
              if (r < zero) then
                 cs = -cs
                 sn = -sn
                 r = -r
              end if
           end if
           return
     end subroutine la_slartgp
     !> DLARTGP: generates a plane rotation so that
     !> [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
     !> [ -SN  CS  ]     [ G ]     [ 0 ]
     !> This is a slower, more accurate version of the Level 1 BLAS routine DROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then CS=(+/-)1 and SN=0.
     !> If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.
     !> The sign is chosen so that R >= 0.

     pure subroutine la_dlartgp(f,g,cs,sn,r)
        use la_constants_dp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(out) :: cs,r,sn
           real(dp),intent(in) :: f,g
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i
           real(dp) :: eps,f1,g1,safmin,safmn2,safmx2,scale
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,sign,sqrt
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Executable Statements
           ! if( first ) then
              safmin = la_dlamch('S')
              eps = la_dlamch('E')
              safmn2 = la_dlamch('B')**int(log(safmin/eps)/log(la_dlamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
              ! first = .false.
           ! end if
           if (g == zero) then
              cs = sign(one,f)
              sn = zero
              r = abs(f)
           else if (f == zero) then
              cs = zero
              sn = sign(one,g)
              r = abs(g)
           else
              f1 = f
              g1 = g
              scale = max(abs(f1),abs(g1))
              if (scale >= safmx2) then
                 count = 0
                 10 continue
                 count = count + 1
                 f1 = f1*safmn2
                 g1 = g1*safmn2
                 scale = max(abs(f1),abs(g1))
                 if (scale >= safmx2 .and. count < 20) go to 10
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmx2
                 end do
              else if (scale <= safmn2) then
                 count = 0
                 30 continue
                 count = count + 1
                 f1 = f1*safmx2
                 g1 = g1*safmx2
                 scale = max(abs(f1),abs(g1))
                 if (scale <= safmn2) go to 30
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmn2
                 end do
              else
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
              end if
              if (r < zero) then
                 cs = -cs
                 sn = -sn
                 r = -r
              end if
           end if
           return
     end subroutine la_dlartgp
#ifdef LA_WITH_XDP
     !> XLARTGP: generates a plane rotation so that
     !> [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
     !> [ -SN  CS  ]     [ G ]     [ 0 ]
     !> This is a slower, more accurate version of the Level 1 BLAS routine XROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then CS=(+/-)1 and SN=0.
     !> If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.
     !> The sign is chosen so that R >= 0.

     pure subroutine la_xlartgp(f,g,cs,sn,r)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(out) :: cs,r,sn
           real(xdp),intent(in) :: f,g
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i
           real(xdp) :: eps,f1,g1,safmin,safmn2,safmx2,scale
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,sign,sqrt
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Executable Statements
           ! if( first ) then
              safmin = la_xlamch('S')
              eps = la_xlamch('E')
              safmn2 = la_xlamch('B')**int(log(safmin/eps)/log(la_xlamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
              ! first = .false.
           ! end if
           if (g == zero) then
              cs = sign(one,f)
              sn = zero
              r = abs(f)
           else if (f == zero) then
              cs = zero
              sn = sign(one,g)
              r = abs(g)
           else
              f1 = f
              g1 = g
              scale = max(abs(f1),abs(g1))
              if (scale >= safmx2) then
                 count = 0
                 10 continue
                 count = count + 1
                 f1 = f1*safmn2
                 g1 = g1*safmn2
                 scale = max(abs(f1),abs(g1))
                 if (scale >= safmx2 .and. count < 20) go to 10
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmx2
                 end do
              else if (scale <= safmn2) then
                 count = 0
                 30 continue
                 count = count + 1
                 f1 = f1*safmx2
                 g1 = g1*safmx2
                 scale = max(abs(f1),abs(g1))
                 if (scale <= safmn2) go to 30
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmn2
                 end do
              else
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
              end if
              if (r < zero) then
                 cs = -cs
                 sn = -sn
                 r = -r
              end if
           end if
           return
     end subroutine la_xlartgp
#endif
#ifdef LA_WITH_QP
     !> QLARTGP: generates a plane rotation so that
     !> [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
     !> [ -SN  CS  ]     [ G ]     [ 0 ]
     !> This is a slower, more accurate version of the Level 1 BLAS routine QROTG,
     !> with the following other differences:
     !> F and G are unchanged on return.
     !> If G=0, then CS=(+/-)1 and SN=0.
     !> If F=0 and (G .ne. 0), then CS=0 and SN=(+/-)1.
     !> The sign is chosen so that R >= 0.

     pure subroutine la_qlartgp(f,g,cs,sn,r)
        use la_constants_qp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(out) :: cs,r,sn
           real(qp),intent(in) :: f,g
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i
           real(qp) :: eps,f1,g1,safmin,safmn2,safmx2,scale
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,sign,sqrt
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Executable Statements
           ! if( first ) then
              safmin = la_qlamch('S')
              eps = la_qlamch('E')
              safmn2 = la_qlamch('B')**int(log(safmin/eps)/log(la_qlamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
              ! first = .false.
           ! end if
           if (g == zero) then
              cs = sign(one,f)
              sn = zero
              r = abs(f)
           else if (f == zero) then
              cs = zero
              sn = sign(one,g)
              r = abs(g)
           else
              f1 = f
              g1 = g
              scale = max(abs(f1),abs(g1))
              if (scale >= safmx2) then
                 count = 0
                 10 continue
                 count = count + 1
                 f1 = f1*safmn2
                 g1 = g1*safmn2
                 scale = max(abs(f1),abs(g1))
                 if (scale >= safmx2 .and. count < 20) go to 10
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmx2
                 end do
              else if (scale <= safmn2) then
                 count = 0
                 30 continue
                 count = count + 1
                 f1 = f1*safmx2
                 g1 = g1*safmx2
                 scale = max(abs(f1),abs(g1))
                 if (scale <= safmn2) go to 30
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
                 do i = 1,count
                    r = r*safmn2
                 end do
              else
                 r = sqrt(f1**2 + g1**2)
                 cs = f1/r
                 sn = g1/r
              end if
              if (r < zero) then
                 cs = -cs
                 sn = -sn
                 r = -r
              end if
           end if
           return
     end subroutine la_qlartgp
#endif

     !> SLARTV: applies a vector of real plane rotations to elements of the
     !> real vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (  c(i)  s(i) ) ( x(i) )
     !> ( y(i) )    ( -s(i)  c(i) ) ( y(i) )

     pure subroutine la_slartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: c(*),s(*)
           real(sp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(sp) :: xi,yi
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - s(ic)*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_slartv
     !> DLARTV: applies a vector of real plane rotations to elements of the
     !> real vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (  c(i)  s(i) ) ( x(i) )
     !> ( y(i) )    ( -s(i)  c(i) ) ( y(i) )

     pure subroutine la_dlartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: c(*),s(*)
           real(dp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(dp) :: xi,yi
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - s(ic)*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_dlartv
#ifdef LA_WITH_XDP
     !> XLARTV: applies a vector of real plane rotations to elements of the
     !> real vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (  c(i)  s(i) ) ( x(i) )
     !> ( y(i) )    ( -s(i)  c(i) ) ( y(i) )

     pure subroutine la_xlartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(xdp),intent(in) :: c(*),s(*)
           real(xdp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(xdp) :: xi,yi
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - s(ic)*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_xlartv
#endif
#ifdef LA_WITH_QP
     !> QLARTV: applies a vector of real plane rotations to elements of the
     !> real vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (  c(i)  s(i) ) ( x(i) )
     !> ( y(i) )    ( -s(i)  c(i) ) ( y(i) )

     pure subroutine la_qlartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(qp),intent(in) :: c(*),s(*)
           real(qp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           real(qp) :: xi,yi
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - s(ic)*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_qlartv
#endif

     !> SLASR: applies a sequence of plane rotations to a real matrix A,
     !> from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_slasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: c(*),s(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(sp) :: ctemp,stemp,temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('SLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_slasr
     !> DLASR: applies a sequence of plane rotations to a real matrix A,
     !> from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_dlasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: c(*),s(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(dp) :: ctemp,stemp,temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('DLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_dlasr
#ifdef LA_WITH_XDP
     !> XLASR: applies a sequence of plane rotations to a real matrix A,
     !> from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_xlasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: c(*),s(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(xdp) :: ctemp,stemp,temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('XLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_xlasr
#endif
#ifdef LA_WITH_QP
     !> QLASR: applies a sequence of plane rotations to a real matrix A,
     !> from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_qlasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: c(*),s(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(qp) :: ctemp,stemp,temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('QLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_qlasr
#endif

     !> CLACRT: performs the operation
     !> (  c  s )( x )  ==> ( x )
     !> ( -s  c )( y )      ( y )
     !> where c and s are complex and the vectors x and y are complex.

     pure subroutine la_clacrt(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           complex(sp),intent(in) :: c,s
           ! Array Arguments
           complex(sp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(sp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
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
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              ctemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - s*cx(i)
              cx(i) = ctemp
           end do
           return
     end subroutine la_clacrt
     !> ZLACRT: performs the operation
     !> (  c  s )( x )  ==> ( x )
     !> ( -s  c )( y )      ( y )
     !> where c and s are complex and the vectors x and y are complex.

     pure subroutine la_zlacrt(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           complex(dp),intent(in) :: c,s
           ! Array Arguments
           complex(dp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(dp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
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
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              ctemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - s*cx(i)
              cx(i) = ctemp
           end do
           return
     end subroutine la_zlacrt
#ifdef LA_WITH_XDP
     !> YLACRT: performs the operation
     !> (  c  s )( x )  ==> ( x )
     !> ( -s  c )( y )      ( y )
     !> where c and s are complex and the vectors x and y are complex.

     pure subroutine la_ylacrt(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           complex(xdp),intent(in) :: c,s
           ! Array Arguments
           complex(xdp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(xdp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
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
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              ctemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - s*cx(i)
              cx(i) = ctemp
           end do
           return
     end subroutine la_ylacrt
#endif
#ifdef LA_WITH_QP
     !> WLACRT: performs the operation
     !> (  c  s )( x )  ==> ( x )
     !> ( -s  c )( y )      ( y )
     !> where c and s are complex and the vectors x and y are complex.

     pure subroutine la_wlacrt(n,cx,incx,cy,incy,c,s)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,incy,n
           complex(qp),intent(in) :: c,s
           ! Array Arguments
           complex(qp),intent(inout) :: cx(*),cy(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ix,iy
           complex(qp) :: ctemp
           ! Executable Statements
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) go to 20
           ! code for unequal increments or equal increments not equal to 1
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
           return
           ! code for both increments equal to 1
           20 continue
           do i = 1,n
              ctemp = c*cx(i) + s*cy(i)
              cy(i) = c*cy(i) - s*cx(i)
              cx(i) = ctemp
           end do
           return
     end subroutine la_wlacrt
#endif

     !> CLAR2V: applies a vector of complex plane rotations with real cosines
     !> from both sides to a sequence of 2-by-2 complex Hermitian matrices,
     !> defined by the elements of the vectors x, y and z. For i = 1,2,...,n
     !> (       x(i)  z(i) ) :=
     !> ( conjg(z(i)) y(i) )
     !> (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
     !> ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )

     pure subroutine la_clar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(sp),intent(in) :: c(*)
           complex(sp),intent(in) :: s(*)
           complex(sp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(sp) :: ci,sii,sir,t1i,t1r,t5,t6,xi,yi,zii,zir
           complex(sp) :: si,t2,t3,t4,zi
           ! Intrinsic Functions
           intrinsic :: aimag,cmplx,conjg,real
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = real(x(ix),KIND=sp)
              yi = real(y(ix),KIND=sp)
              zi = z(ix)
              zir = real(zi,KIND=sp)
              zii = aimag(zi)
              ci = c(ic)
              si = s(ic)
              sir = real(si,KIND=sp)
              sii = aimag(si)
              t1r = sir*zir - sii*zii
              t1i = sir*zii + sii*zir
              t2 = ci*zi
              t3 = t2 - conjg(si)*xi
              t4 = conjg(t2) + si*yi
              t5 = ci*xi + t1r
              t6 = ci*yi - t1r
              x(ix) = ci*t5 + (sir*real(t4,KIND=sp) + sii*aimag(t4))
              y(ix) = ci*t6 - (sir*real(t3,KIND=sp) - sii*aimag(t3))
              z(ix) = ci*t3 + conjg(si)*cmplx(t6,t1i,KIND=sp)
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_clar2v
     !> ZLAR2V: applies a vector of complex plane rotations with real cosines
     !> from both sides to a sequence of 2-by-2 complex Hermitian matrices,
     !> defined by the elements of the vectors x, y and z. For i = 1,2,...,n
     !> (       x(i)  z(i) ) :=
     !> ( conjg(z(i)) y(i) )
     !> (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
     !> ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )

     pure subroutine la_zlar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(dp),intent(in) :: c(*)
           complex(dp),intent(in) :: s(*)
           complex(dp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(dp) :: ci,sii,sir,t1i,t1r,t5,t6,xi,yi,zii,zir
           complex(dp) :: si,t2,t3,t4,zi
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = real(x(ix),KIND=dp)
              yi = real(y(ix),KIND=dp)
              zi = z(ix)
              zir = real(zi,KIND=dp)
              zii = aimag(zi)
              ci = c(ic)
              si = s(ic)
              sir = real(si,KIND=dp)
              sii = aimag(si)
              t1r = sir*zir - sii*zii
              t1i = sir*zii + sii*zir
              t2 = ci*zi
              t3 = t2 - conjg(si)*xi
              t4 = conjg(t2) + si*yi
              t5 = ci*xi + t1r
              t6 = ci*yi - t1r
              x(ix) = ci*t5 + (sir*real(t4,KIND=dp) + sii*aimag(t4))
              y(ix) = ci*t6 - (sir*real(t3,KIND=dp) - sii*aimag(t3))
              z(ix) = ci*t3 + conjg(si)*cmplx(t6,t1i,KIND=dp)
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_zlar2v
#ifdef LA_WITH_XDP
     !> YLAR2V: applies a vector of complex plane rotations with real cosines
     !> from both sides to a sequence of 2-by-2 complex Hermitian matrices,
     !> defined by the elements of the vectors x, y and z. For i = 1,2,...,n
     !> (       x(i)  z(i) ) :=
     !> ( conjg(z(i)) y(i) )
     !> (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
     !> ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )

     pure subroutine la_ylar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(xdp),intent(in) :: c(*)
           complex(xdp),intent(in) :: s(*)
           complex(xdp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(xdp) :: ci,sii,sir,t1i,t1r,t5,t6,xi,yi,zii,zir
           complex(xdp) :: si,t2,t3,t4,zi
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = real(x(ix),KIND=xdp)
              yi = real(y(ix),KIND=xdp)
              zi = z(ix)
              zir = real(zi,KIND=xdp)
              zii = aimag(zi)
              ci = c(ic)
              si = s(ic)
              sir = real(si,KIND=xdp)
              sii = aimag(si)
              t1r = sir*zir - sii*zii
              t1i = sir*zii + sii*zir
              t2 = ci*zi
              t3 = t2 - conjg(si)*xi
              t4 = conjg(t2) + si*yi
              t5 = ci*xi + t1r
              t6 = ci*yi - t1r
              x(ix) = ci*t5 + (sir*real(t4,KIND=xdp) + sii*aimag(t4))
              y(ix) = ci*t6 - (sir*real(t3,KIND=xdp) - sii*aimag(t3))
              z(ix) = ci*t3 + conjg(si)*cmplx(t6,t1i,KIND=xdp)
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_ylar2v
#endif
#ifdef LA_WITH_QP
     !> WLAR2V: applies a vector of complex plane rotations with real cosines
     !> from both sides to a sequence of 2-by-2 complex Hermitian matrices,
     !> defined by the elements of the vectors x, y and z. For i = 1,2,...,n
     !> (       x(i)  z(i) ) :=
     !> ( conjg(z(i)) y(i) )
     !> (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
     !> ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )

     pure subroutine la_wlar2v(n,x,y,z,incx,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,n
           ! Array Arguments
           real(qp),intent(in) :: c(*)
           complex(qp),intent(in) :: s(*)
           complex(qp),intent(inout) :: x(*),y(*),z(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix
           real(qp) :: ci,sii,sir,t1i,t1r,t5,t6,xi,yi,zii,zir
           complex(qp) :: si,t2,t3,t4,zi
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,aimag
           ! Executable Statements
           ix = 1
           ic = 1
           do i = 1,n
              xi = real(x(ix),KIND=qp)
              yi = real(y(ix),KIND=qp)
              zi = z(ix)
              zir = real(zi,KIND=qp)
              zii = aimag(zi)
              ci = c(ic)
              si = s(ic)
              sir = real(si,KIND=qp)
              sii = aimag(si)
              t1r = sir*zir - sii*zii
              t1i = sir*zii + sii*zir
              t2 = ci*zi
              t3 = t2 - conjg(si)*xi
              t4 = conjg(t2) + si*yi
              t5 = ci*xi + t1r
              t6 = ci*yi - t1r
              x(ix) = ci*t5 + (sir*real(t4,KIND=qp) + sii*aimag(t4))
              y(ix) = ci*t6 - (sir*real(t3,KIND=qp) - sii*aimag(t3))
              z(ix) = ci*t3 + conjg(si)*cmplx(t6,t1i,KIND=qp)
              ix = ix + incx
              ic = ic + incc
           end do
           return
     end subroutine la_wlar2v
#endif

     !> !
     !>
     !> CLARTG: generates a plane rotation so that
     !> [  C         S  ] . [ F ]  =  [ R ]
     !> [ -conjg(S)  C  ]   [ G ]     [ 0 ]
     !> where C is real and C**2 + |S|**2 = 1.
     !> The mathematical formulas used for C and S are
     !> sgn(x) = {  x / |x|,   x != 0
     !> {  1,         x = 0
     !> R = sgn(F) * sqrt(|F|**2 + |G|**2)
     !> C = |F| / sqrt(|F|**2 + |G|**2)
     !> S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
     !> When F and G are real, the formulas simplify to C = F/R and
     !> S = G/R, and the returned values of C, S, and R should be
     !> identical to those returned by CLARTG.
     !> The algorithm used to compute these quantities incorporates scaling
     !> to avoid overflow or underflow in computing the square root of the
     !> sum of squares.
     !> This is a faster version of the BLAS1 routine CROTG, except for
     !> the following differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0, then C=0 and S is chosen so that R is real.
     !> Below, wp=>sp stands for single precision from LA_CONSTANTS module.

     pure subroutine la_clartg(f,g,c,s,r)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(sp),intent(out) :: c
        complex(sp),intent(in) :: f,g
        complex(sp),intent(out) :: r,s
        ! Local Scalars
        real(sp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(sp) :: fs,gs,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(sp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=sp)**2 + aimag(t)**2
        ! Executable Statements
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
        return
     end subroutine la_clartg
     !> !
     !>
     !> ZLARTG: generates a plane rotation so that
     !> [  C         S  ] . [ F ]  =  [ R ]
     !> [ -conjg(S)  C  ]   [ G ]     [ 0 ]
     !> where C is real and C**2 + |S|**2 = 1.
     !> The mathematical formulas used for C and S are
     !> sgn(x) = {  x / |x|,   x != 0
     !> {  1,         x = 0
     !> R = sgn(F) * sqrt(|F|**2 + |G|**2)
     !> C = |F| / sqrt(|F|**2 + |G|**2)
     !> S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
     !> When F and G are real, the formulas simplify to C = F/R and
     !> S = G/R, and the returned values of C, S, and R should be
     !> identical to those returned by DLARTG.
     !> The algorithm used to compute these quantities incorporates scaling
     !> to avoid overflow or underflow in computing the square root of the
     !> sum of squares.
     !> This is a faster version of the BLAS1 routine ZROTG, except for
     !> the following differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0, then C=0 and S is chosen so that R is real.
     !> Below, wp=>dp stands for double precision from LA_CONSTANTS module.

     pure subroutine la_zlartg(f,g,c,s,r)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(dp),intent(out) :: c
        complex(dp),intent(in) :: f,g
        complex(dp),intent(out) :: r,s
        ! Local Scalars
        real(dp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(dp) :: fs,gs,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(dp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=dp)**2 + aimag(t)**2
        ! Executable Statements
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
        return
     end subroutine la_zlartg
#ifdef LA_WITH_XDP
     !> !
     !>
     !> YLARTG: generates a plane rotation so that
     !> [  C         S  ] . [ F ]  =  [ R ]
     !> [ -conjg(S)  C  ]   [ G ]     [ 0 ]
     !> where C is real and C**2 + |S|**2 = 1.
     !> The mathematical formulas used for C and S are
     !> sgn(x) = {  x / |x|,   x != 0
     !> {  1,         x = 0
     !> R = sgn(F) * sqrt(|F|**2 + |G|**2)
     !> C = |F| / sqrt(|F|**2 + |G|**2)
     !> S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
     !> When F and G are real, the formulas simplify to C = F/R and
     !> S = G/R, and the returned values of C, S, and R should be
     !> identical to those returned by XLARTG.
     !> The algorithm used to compute these quantities incorporates scaling
     !> to avoid overflow or underflow in computing the square root of the
     !> sum of squares.
     !> This is a faster version of the BLAS1 routine YROTG, except for
     !> the following differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0, then C=0 and S is chosen so that R is real.
     !> Below, wp=>xdp stands for extended precision from LA_CONSTANTS module.

     pure subroutine la_ylartg(f,g,c,s,r)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(xdp),intent(out) :: c
        complex(xdp),intent(in) :: f,g
        complex(xdp),intent(out) :: r,s
        ! Local Scalars
        real(xdp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(xdp) :: fs,gs,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(xdp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=xdp)**2 + aimag(t)**2
        ! Executable Statements
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
        return
     end subroutine la_ylartg
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> WLARTG: generates a plane rotation so that
     !> [  C         S  ] . [ F ]  =  [ R ]
     !> [ -conjg(S)  C  ]   [ G ]     [ 0 ]
     !> where C is real and C**2 + |S|**2 = 1.
     !> The mathematical formulas used for C and S are
     !> sgn(x) = {  x / |x|,   x != 0
     !> {  1,         x = 0
     !> R = sgn(F) * sqrt(|F|**2 + |G|**2)
     !> C = |F| / sqrt(|F|**2 + |G|**2)
     !> S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
     !> When F and G are real, the formulas simplify to C = F/R and
     !> S = G/R, and the returned values of C, S, and R should be
     !> identical to those returned by QLARTG.
     !> The algorithm used to compute these quantities incorporates scaling
     !> to avoid overflow or underflow in computing the square root of the
     !> sum of squares.
     !> This is a faster version of the BLAS1 routine WROTG, except for
     !> the following differences:
     !> F and G are unchanged on return.
     !> If G=0, then C=1 and S=0.
     !> If F=0, then C=0 and S is chosen so that R is real.
     !> Below, wp=>qp stands for quad precision from LA_CONSTANTS module.

     pure subroutine la_wlartg(f,g,c,s,r)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! february 2021
        ! Scalar Arguments
        real(qp),intent(out) :: c
        complex(qp),intent(in) :: f,g
        complex(qp),intent(out) :: r,s
        ! Local Scalars
        real(qp) :: d,f1,f2,g1,g2,h2,p,u,uu,v,vv,w
        complex(qp) :: fs,gs,t
        ! Intrinsic Functions
        intrinsic :: abs,aimag,conjg,max,min,real,sqrt
        ! Statement Functions
        real(qp) :: abssq
        ! Statement Function Definitions
        abssq(t) = real(t,KIND=qp)**2 + aimag(t)**2
        ! Executable Statements
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
        return
     end subroutine la_wlartg
#endif

     !> CLARTV: applies a vector of complex plane rotations with real cosines
     !> to elements of the complex vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
     !> ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )

     pure subroutine la_clartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(sp),intent(in) :: c(*)
           complex(sp),intent(in) :: s(*)
           complex(sp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           complex(sp) :: xi,yi
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - conjg(s(ic))*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_clartv
     !> ZLARTV: applies a vector of complex plane rotations with real cosines
     !> to elements of the complex vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
     !> ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )

     pure subroutine la_zlartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(dp),intent(in) :: c(*)
           complex(dp),intent(in) :: s(*)
           complex(dp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           complex(dp) :: xi,yi
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - conjg(s(ic))*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_zlartv
#ifdef LA_WITH_XDP
     !> YLARTV: applies a vector of complex plane rotations with real cosines
     !> to elements of the complex vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
     !> ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )

     pure subroutine la_ylartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(xdp),intent(in) :: c(*)
           complex(xdp),intent(in) :: s(*)
           complex(xdp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           complex(xdp) :: xi,yi
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - conjg(s(ic))*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_ylartv
#endif
#ifdef LA_WITH_QP
     !> WLARTV: applies a vector of complex plane rotations with real cosines
     !> to elements of the complex vectors x and y. For i = 1,2,...,n
     !> ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
     !> ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )

     pure subroutine la_wlartv(n,x,incx,y,incy,c,s,incc)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(qp),intent(in) :: c(*)
           complex(qp),intent(in) :: s(*)
           complex(qp),intent(inout) :: x(*),y(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ic,ix,iy
           complex(qp) :: xi,yi
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           ix = 1
           iy = 1
           ic = 1
           do i = 1,n
              xi = x(ix)
              yi = y(iy)
              x(ix) = c(ic)*xi + s(ic)*yi
              y(iy) = c(ic)*yi - conjg(s(ic))*xi
              ix = ix + incx
              iy = iy + incy
              ic = ic + incc
           end do
           return
     end subroutine la_wlartv
#endif

     !> CLASR: applies a sequence of real plane rotations to a complex matrix
     !> A, from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_clasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(in) :: c(*),s(*)
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(sp) :: ctemp,stemp
           complex(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('CLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_clasr
     !> ZLASR: applies a sequence of real plane rotations to a complex matrix
     !> A, from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_zlasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(in) :: c(*),s(*)
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(dp) :: ctemp,stemp
           complex(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('ZLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_zlasr
#ifdef LA_WITH_XDP
     !> YLASR: applies a sequence of real plane rotations to a complex matrix
     !> A, from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_ylasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(in) :: c(*),s(*)
           complex(xdp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(xdp) :: ctemp,stemp
           complex(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('YLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_ylasr
#endif
#ifdef LA_WITH_QP
     !> WLASR: applies a sequence of real plane rotations to a complex matrix
     !> A, from either the left or the right.
     !> When SIDE = 'L', the transformation takes the form
     !> A := P*A
     !> and when SIDE = 'R', the transformation takes the form
     !> A := A*P**T
     !> where P is an orthogonal matrix consisting of a sequence of z plane
     !> rotations, with z = M when SIDE = 'L' and z = N when SIDE = 'R',
     !> and P**T is the transpose of P.
     !> When DIRECT = 'F' (Forward sequence), then
     !> P = P(z-1) * ... * P(2) * P(1)
     !> and when DIRECT = 'B' (Backward sequence), then
     !> P = P(1) * P(2) * ... * P(z-1)
     !> where P(k) is a plane rotation matrix defined by the 2-by-2 rotation
     !> R(k) = (  c(k)  s(k) )
     !> = ( -s(k)  c(k) ).
     !> When PIVOT = 'V' (Variable pivot), the rotation is performed
     !> for the plane (k,k+1), i.e., P(k) has the form
     !> P(k) = (  1                                            )
     !> (       ...                                     )
     !> (              1                                )
     !> (                   c(k)  s(k)                  )
     !> (                  -s(k)  c(k)                  )
     !> (                                1              )
     !> (                                     ...       )
     !> (                                            1  )
     !> where R(k) appears as a rank-2 modification to the identity matrix in
     !> rows and columns k and k+1.
     !> When PIVOT = 'T' (Top pivot), the rotation is performed for the
     !> plane (1,k+1), so P(k) has the form
     !> P(k) = (  c(k)                    s(k)                 )
     !> (         1                                     )
     !> (              ...                              )
     !> (                     1                         )
     !> ( -s(k)                    c(k)                 )
     !> (                                 1             )
     !> (                                      ...      )
     !> (                                             1 )
     !> where R(k) appears in rows and columns 1 and k+1.
     !> Similarly, when PIVOT = 'B' (Bottom pivot), the rotation is
     !> performed for the plane (k,z), giving P(k) the form
     !> P(k) = ( 1                                             )
     !> (      ...                                      )
     !> (             1                                 )
     !> (                  c(k)                    s(k) )
     !> (                         1                     )
     !> (                              ...              )
     !> (                                     1         )
     !> (                 -s(k)                    c(k) )
     !> where R(k) appears in rows and columns k and z.  The rotations are
     !> performed without ever forming P(k) explicitly.

     pure subroutine la_wlasr(side,pivot,direct,m,n,c,s,a,lda)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: direct,pivot,side
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(in) :: c(*),s(*)
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,j
           real(qp) :: ctemp,stemp
           complex(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (.not. (la_lsame(side,'L') .or. la_lsame(side,'R'))) then
              info = 1
           else if (.not. (la_lsame(pivot,'V') .or. la_lsame(pivot,'T') .or. &
                     la_lsame(pivot,'B'))) then
              info = 2
           else if (.not. (la_lsame(direct,'F') .or. la_lsame(direct,'B'))) &
                     then
              info = 3
           else if (m < 0) then
              info = 4
           else if (n < 0) then
              info = 5
           else if (lda < max(1,m)) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('WLASR ',info)
              return
           end if
           ! quick return if possible
           if ((m == 0) .or. (n == 0)) return
           if (la_lsame(side,'L')) then
              ! form  p * a
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j + 1,i)
                             a(j + 1,i) = ctemp*temp - stemp*a(j,i)
                             a(j,i) = stemp*temp + ctemp*a(j,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,m
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = ctemp*temp - stemp*a(1,i)
                             a(1,i) = stemp*temp + ctemp*a(1,i)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,m - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = m - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,n
                             temp = a(j,i)
                             a(j,i) = stemp*a(m,i) + ctemp*temp
                             a(m,i) = ctemp*a(m,i) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           else if (la_lsame(side,'R')) then
              ! form a * p**t
              if (la_lsame(pivot,'V')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j + 1)
                             a(i,j + 1) = ctemp*temp - stemp*a(i,j)
                             a(i,j) = stemp*temp + ctemp*a(i,j)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'T')) then
                 if (la_lsame(direct,'F')) then
                    do j = 2,n
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n,2,-1
                       ctemp = c(j - 1)
                       stemp = s(j - 1)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = ctemp*temp - stemp*a(i,1)
                             a(i,1) = stemp*temp + ctemp*a(i,1)
                          end do
                       end if
                    end do
                 end if
              else if (la_lsame(pivot,'B')) then
                 if (la_lsame(direct,'F')) then
                    do j = 1,n - 1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 else if (la_lsame(direct,'B')) then
                    do j = n - 1,1,-1
                       ctemp = c(j)
                       stemp = s(j)
                       if ((ctemp /= one) .or. (stemp /= zero)) then
                          do i = 1,m
                             temp = a(i,j)
                             a(i,j) = stemp*a(i,n) + ctemp*temp
                             a(i,n) = ctemp*a(i,n) - stemp*temp
                          end do
                       end if
                    end do
                 end if
              end if
           end if
           return
     end subroutine la_wlasr
#endif

     !> CLARGV: generates a vector of complex plane rotations with real
     !> cosines, determined by elements of the complex vectors x and y.
     !> For i = 1,2,...,n
     !> (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
     !> ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
     !> where c(i)**2 + ABS(s(i))**2 = 1
     !> The following conventions are used (these are the same as in CLARTG,
     !> but differ from the BLAS1 routine CROTG):
     !> If y(i)=0, then c(i)=1 and s(i)=0.
     !> If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.

     pure subroutine la_clargv(n,x,incx,y,incy,c,incc)
        use la_constants_sp,only:zero,one,two,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(sp),intent(out) :: c(*)
           complex(sp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i,ic,ix,iy,j
           real(sp) :: cs,d,di,dr,eps,f2,f2s,g2,g2s,safmin,safmn2,safmx2,scale
           complex(sp) :: f,ff,fs,g,gs,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,int,log,max,real,sqrt
           ! Statement Functions
           real(sp) :: abs1,abssq
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Statement Function Definitions
           abs1(ff) = max(abs(real(ff,KIND=sp)),abs(aimag(ff)))
           abssq(ff) = real(ff,KIND=sp)**2 + aimag(ff)**2
           ! Executable Statements
           ! if( first ) then
              ! first = .false.
              safmin = la_slamch('S')
              eps = la_slamch('E')
              safmn2 = la_slamch('B')**int(log(safmin/eps)/log(la_slamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
           ! end if
           ix = 1
           iy = 1
           ic = 1
           loop_60: do i = 1,n
              f = x(ix)
              g = y(iy)
              ! use identical algorithm as in la_clartg
              scale = max(abs1(f),abs1(g))
              fs = f
              gs = g
              count = 0
              if (scale >= safmx2) then
              10 continue
                 count = count + 1
                 fs = fs*safmn2
                 gs = gs*safmn2
                 scale = scale*safmn2
                 if (scale >= safmx2 .and. count < 20) go to 10
              else if (scale <= safmn2) then
                 if (g == czero) then
                    cs = one
                    sn = czero
                    r = f
                    go to 50
                 end if
                 20 continue
                 count = count - 1
                 fs = fs*safmx2
                 gs = gs*safmx2
                 scale = scale*safmx2
                 if (scale <= safmn2) go to 20
              end if
              f2 = abssq(fs)
              g2 = abssq(gs)
              if (f2 <= max(g2,one)*safmin) then
                 ! this is a rare case: f is very small.
                 if (f == czero) then
                    cs = zero
                    r = la_slapy2(real(g,KIND=sp),aimag(g))
                    ! do complex/real division explicitly with two real
                    ! divisions
                    d = la_slapy2(real(gs,KIND=sp),aimag(gs))
                    sn = cmplx(real(gs,KIND=sp)/d,-aimag(gs)/d,KIND=sp)
                    go to 50
                 end if
                 f2s = la_slapy2(real(fs,KIND=sp),aimag(fs))
                 ! g2 and g2s are accurate
                 ! g2 is at least safmin, and g2s is at least safmn2
                 g2s = sqrt(g2)
                 ! error in cs from underflow in f2s is at most
                 ! unfl / safmn2 .lt. sqrt(unfl*eps) .lt. eps
                 ! if max(g2,one)=g2, then f2 .lt. g2*safmin,
                 ! and so cs .lt. sqrt(safmin)
                 ! if max(g2,one)=one, then f2 .lt. safmin
                 ! and so cs .lt. sqrt(safmin)/safmn2 = sqrt(eps)
                 ! therefore, cs = f2s/g2s / sqrt( 1 + (f2s/g2s)**2 ) = f2s/g2s
                 cs = f2s/g2s
                 ! make sure abs(ff) = 1
                 ! do complex/real division explicitly with 2 real divisions
                 if (abs1(f) > one) then
                    d = la_slapy2(real(f,KIND=sp),aimag(f))
                    ff = cmplx(real(f,KIND=sp)/d,aimag(f)/d,KIND=sp)
                 else
                    dr = safmx2*real(f,KIND=sp)
                    di = safmx2*aimag(f)
                    d = la_slapy2(dr,di)
                    ff = cmplx(dr/d,di/d,KIND=sp)
                 end if
                 sn = ff*cmplx(real(gs,KIND=sp)/g2s,-aimag(gs)/g2s,KIND=sp)
                 r = cs*f + sn*g
              else
                 ! this is the most common case.
                 ! neither f2 nor f2/g2 are less than safmin
                 ! f2s cannot overflow, and it is accurate
                 f2s = sqrt(one + g2/f2)
                 ! do the f2s(real)*fs(complex) multiply with two real
                 ! multiplies
                 r = cmplx(f2s*real(fs,KIND=sp),f2s*aimag(fs),KIND=sp)
                 cs = one/f2s
                 d = f2 + g2
                 ! do complex/real division explicitly with two real divisions
                 sn = cmplx(real(r,KIND=sp)/d,aimag(r)/d,KIND=sp)
                 sn = sn*conjg(gs)
                 if (count /= 0) then
                    if (count > 0) then
                       do j = 1,count
                          r = r*safmx2
                       end do
                    else
                       do j = 1,-count
                          r = r*safmn2
                       end do
                    end if
                 end if
              end if
              50 continue
              c(ic) = cs
              y(iy) = sn
              x(ix) = r
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_60
           return
     end subroutine la_clargv
     !> ZLARGV: generates a vector of complex plane rotations with real
     !> cosines, determined by elements of the complex vectors x and y.
     !> For i = 1,2,...,n
     !> (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
     !> ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
     !> where c(i)**2 + ABS(s(i))**2 = 1
     !> The following conventions are used (these are the same as in ZLARTG,
     !> but differ from the BLAS1 routine ZROTG):
     !> If y(i)=0, then c(i)=1 and s(i)=0.
     !> If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.

     pure subroutine la_zlargv(n,x,incx,y,incy,c,incc)
        use la_constants_dp,only:zero,one,two,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(dp),intent(out) :: c(*)
           complex(dp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i,ic,ix,iy,j
           real(dp) :: cs,d,di,dr,eps,f2,f2s,g2,g2s,safmin,safmn2,safmx2,scale
           complex(dp) :: f,ff,fs,g,gs,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,int,log,max,sqrt
           ! Statement Functions
           real(dp) :: abs1,abssq
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Statement Function Definitions
           abs1(ff) = max(abs(real(ff,KIND=dp)),abs(aimag(ff)))
           abssq(ff) = real(ff,KIND=dp)**2 + aimag(ff)**2
           ! Executable Statements
           ! if( first ) then
              ! first = .false.
              safmin = la_dlamch('S')
              eps = la_dlamch('E')
              safmn2 = la_dlamch('B')**int(log(safmin/eps)/log(la_dlamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
           ! end if
           ix = 1
           iy = 1
           ic = 1
           loop_60: do i = 1,n
              f = x(ix)
              g = y(iy)
              ! use identical algorithm as in la_zlartg
              scale = max(abs1(f),abs1(g))
              fs = f
              gs = g
              count = 0
              if (scale >= safmx2) then
              10 continue
                 count = count + 1
                 fs = fs*safmn2
                 gs = gs*safmn2
                 scale = scale*safmn2
                 if (scale >= safmx2 .and. count < 20) go to 10
              else if (scale <= safmn2) then
                 if (g == czero) then
                    cs = one
                    sn = czero
                    r = f
                    go to 50
                 end if
                 20 continue
                 count = count - 1
                 fs = fs*safmx2
                 gs = gs*safmx2
                 scale = scale*safmx2
                 if (scale <= safmn2) go to 20
              end if
              f2 = abssq(fs)
              g2 = abssq(gs)
              if (f2 <= max(g2,one)*safmin) then
                 ! this is a rare case: f is very small.
                 if (f == czero) then
                    cs = zero
                    r = la_dlapy2(real(g,KIND=dp),aimag(g))
                    ! do complex/real division explicitly with two real
                    ! divisions
                    d = la_dlapy2(real(gs,KIND=dp),aimag(gs))
                    sn = cmplx(real(gs,KIND=dp)/d,-aimag(gs)/d,KIND=dp)
                    go to 50
                 end if
                 f2s = la_dlapy2(real(fs,KIND=dp),aimag(fs))
                 ! g2 and g2s are accurate
                 ! g2 is at least safmin, and g2s is at least safmn2
                 g2s = sqrt(g2)
                 ! error in cs from underflow in f2s is at most
                 ! unfl / safmn2 .lt. sqrt(unfl*eps) .lt. eps
                 ! if max(g2,one)=g2, then f2 .lt. g2*safmin,
                 ! and so cs .lt. sqrt(safmin)
                 ! if max(g2,one)=one, then f2 .lt. safmin
                 ! and so cs .lt. sqrt(safmin)/safmn2 = sqrt(eps)
                 ! therefore, cs = f2s/g2s / sqrt( 1 + (f2s/g2s)**2 ) = f2s/g2s
                 cs = f2s/g2s
                 ! make sure abs(ff) = 1
                 ! do complex/real division explicitly with 2 real divisions
                 if (abs1(f) > one) then
                    d = la_dlapy2(real(f,KIND=dp),aimag(f))
                    ff = cmplx(real(f,KIND=dp)/d,aimag(f)/d,KIND=dp)
                 else
                    dr = safmx2*real(f,KIND=dp)
                    di = safmx2*aimag(f)
                    d = la_dlapy2(dr,di)
                    ff = cmplx(dr/d,di/d,KIND=dp)
                 end if
                 sn = ff*cmplx(real(gs,KIND=dp)/g2s,-aimag(gs)/g2s,KIND=dp)
                 r = cs*f + sn*g
              else
                 ! this is the most common case.
                 ! neither f2 nor f2/g2 are less than safmin
                 ! f2s cannot overflow, and it is accurate
                 f2s = sqrt(one + g2/f2)
                 ! do the f2s(real)*fs(complex) multiply with two real
                 ! multiplies
                 r = cmplx(f2s*real(fs,KIND=dp),f2s*aimag(fs),KIND=dp)
                 cs = one/f2s
                 d = f2 + g2
                 ! do complex/real division explicitly with two real divisions
                 sn = cmplx(real(r,KIND=dp)/d,aimag(r)/d,KIND=dp)
                 sn = sn*conjg(gs)
                 if (count /= 0) then
                    if (count > 0) then
                       do j = 1,count
                          r = r*safmx2
                       end do
                    else
                       do j = 1,-count
                          r = r*safmn2
                       end do
                    end if
                 end if
              end if
              50 continue
              c(ic) = cs
              y(iy) = sn
              x(ix) = r
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_60
           return
     end subroutine la_zlargv
#ifdef LA_WITH_XDP
     !> YLARGV: generates a vector of complex plane rotations with real
     !> cosines, determined by elements of the complex vectors x and y.
     !> For i = 1,2,...,n
     !> (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
     !> ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
     !> where c(i)**2 + ABS(s(i))**2 = 1
     !> The following conventions are used (these are the same as in YLARTG,
     !> but differ from the BLAS1 routine YROTG):
     !> If y(i)=0, then c(i)=1 and s(i)=0.
     !> If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.

     pure subroutine la_ylargv(n,x,incx,y,incy,c,incc)
        use la_constants_xdp,only:zero,one,two,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(xdp),intent(out) :: c(*)
           complex(xdp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i,ic,ix,iy,j
           real(xdp) :: cs,d,di,dr,eps,f2,f2s,g2,g2s,safmin,safmn2,safmx2,scale
           complex(xdp) :: f,ff,fs,g,gs,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,int,log,max,sqrt
           ! Statement Functions
           real(xdp) :: abs1,abssq
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Statement Function Definitions
           abs1(ff) = max(abs(real(ff,KIND=xdp)),abs(aimag(ff)))
           abssq(ff) = real(ff,KIND=xdp)**2 + aimag(ff)**2
           ! Executable Statements
           ! if( first ) then
              ! first = .false.
              safmin = la_xlamch('S')
              eps = la_xlamch('E')
              safmn2 = la_xlamch('B')**int(log(safmin/eps)/log(la_xlamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
           ! end if
           ix = 1
           iy = 1
           ic = 1
           loop_60: do i = 1,n
              f = x(ix)
              g = y(iy)
              ! use identical algorithm as in la_ylartg
              scale = max(abs1(f),abs1(g))
              fs = f
              gs = g
              count = 0
              if (scale >= safmx2) then
              10 continue
                 count = count + 1
                 fs = fs*safmn2
                 gs = gs*safmn2
                 scale = scale*safmn2
                 if (scale >= safmx2 .and. count < 20) go to 10
              else if (scale <= safmn2) then
                 if (g == czero) then
                    cs = one
                    sn = czero
                    r = f
                    go to 50
                 end if
                 20 continue
                 count = count - 1
                 fs = fs*safmx2
                 gs = gs*safmx2
                 scale = scale*safmx2
                 if (scale <= safmn2) go to 20
              end if
              f2 = abssq(fs)
              g2 = abssq(gs)
              if (f2 <= max(g2,one)*safmin) then
                 ! this is a rare case: f is very small.
                 if (f == czero) then
                    cs = zero
                    r = la_xlapy2(real(g,KIND=xdp),aimag(g))
                    ! do complex/real division explicitly with two real
                    ! divisions
                    d = la_xlapy2(real(gs,KIND=xdp),aimag(gs))
                    sn = cmplx(real(gs,KIND=xdp)/d,-aimag(gs)/d,KIND=xdp)
                    go to 50
                 end if
                 f2s = la_xlapy2(real(fs,KIND=xdp),aimag(fs))
                 ! g2 and g2s are accurate
                 ! g2 is at least safmin, and g2s is at least safmn2
                 g2s = sqrt(g2)
                 ! error in cs from underflow in f2s is at most
                 ! unfl / safmn2 .lt. sqrt(unfl*eps) .lt. eps
                 ! if max(g2,one)=g2, then f2 .lt. g2*safmin,
                 ! and so cs .lt. sqrt(safmin)
                 ! if max(g2,one)=one, then f2 .lt. safmin
                 ! and so cs .lt. sqrt(safmin)/safmn2 = sqrt(eps)
                 ! therefore, cs = f2s/g2s / sqrt( 1 + (f2s/g2s)**2 ) = f2s/g2s
                 cs = f2s/g2s
                 ! make sure abs(ff) = 1
                 ! do complex/real division explicitly with 2 real divisions
                 if (abs1(f) > one) then
                    d = la_xlapy2(real(f,KIND=xdp),aimag(f))
                    ff = cmplx(real(f,KIND=xdp)/d,aimag(f)/d,KIND=xdp)
                 else
                    dr = safmx2*real(f,KIND=xdp)
                    di = safmx2*aimag(f)
                    d = la_xlapy2(dr,di)
                    ff = cmplx(dr/d,di/d,KIND=xdp)
                 end if
                 sn = ff*cmplx(real(gs,KIND=xdp)/g2s,-aimag(gs)/g2s,KIND=xdp)
                 r = cs*f + sn*g
              else
                 ! this is the most common case.
                 ! neither f2 nor f2/g2 are less than safmin
                 ! f2s cannot overflow, and it is accurate
                 f2s = sqrt(one + g2/f2)
                 ! do the f2s(real)*fs(complex) multiply with two real
                 ! multiplies
                 r = cmplx(f2s*real(fs,KIND=xdp),f2s*aimag(fs),KIND=xdp)
                 cs = one/f2s
                 d = f2 + g2
                 ! do complex/real division explicitly with two real divisions
                 sn = cmplx(real(r,KIND=xdp)/d,aimag(r)/d,KIND=xdp)
                 sn = sn*conjg(gs)
                 if (count /= 0) then
                    if (count > 0) then
                       do j = 1,count
                          r = r*safmx2
                       end do
                    else
                       do j = 1,-count
                          r = r*safmn2
                       end do
                    end if
                 end if
              end if
              50 continue
              c(ic) = cs
              y(iy) = sn
              x(ix) = r
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_60
           return
     end subroutine la_ylargv
#endif
#ifdef LA_WITH_QP
     !> WLARGV: generates a vector of complex plane rotations with real
     !> cosines, determined by elements of the complex vectors x and y.
     !> For i = 1,2,...,n
     !> (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
     !> ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
     !> where c(i)**2 + ABS(s(i))**2 = 1
     !> The following conventions are used (these are the same as in WLARTG,
     !> but differ from the BLAS1 routine WROTG):
     !> If y(i)=0, then c(i)=1 and s(i)=0.
     !> If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.

     pure subroutine la_wlargv(n,x,incx,y,incy,c,incc)
        use la_constants_qp,only:zero,one,two,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incc,incx,incy,n
           ! Array Arguments
           real(qp),intent(out) :: c(*)
           complex(qp),intent(inout) :: x(*),y(*)
        ! =====================================================================

           ! Local Scalars
           ! logical            first
           integer(ilp) :: count,i,ic,ix,iy,j
           real(qp) :: cs,d,di,dr,eps,f2,f2s,g2,g2s,safmin,safmn2,safmx2,scale
           complex(qp) :: f,ff,fs,g,gs,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,int,log,max,sqrt
           ! Statement Functions
           real(qp) :: abs1,abssq
           ! Save Statement
           ! save               first, safmx2, safmin, safmn2
           ! Data Statements
           ! data               first / .true. /
           ! Statement Function Definitions
           abs1(ff) = max(abs(real(ff,KIND=qp)),abs(aimag(ff)))
           abssq(ff) = real(ff,KIND=qp)**2 + aimag(ff)**2
           ! Executable Statements
           ! if( first ) then
              ! first = .false.
              safmin = la_qlamch('S')
              eps = la_qlamch('E')
              safmn2 = la_qlamch('B')**int(log(safmin/eps)/log(la_qlamch('B')) &
                         /two,KIND=ilp)
              safmx2 = one/safmn2
           ! end if
           ix = 1
           iy = 1
           ic = 1
           loop_60: do i = 1,n
              f = x(ix)
              g = y(iy)
              ! use identical algorithm as in la_wlartg
              scale = max(abs1(f),abs1(g))
              fs = f
              gs = g
              count = 0
              if (scale >= safmx2) then
              10 continue
                 count = count + 1
                 fs = fs*safmn2
                 gs = gs*safmn2
                 scale = scale*safmn2
                 if (scale >= safmx2 .and. count < 20) go to 10
              else if (scale <= safmn2) then
                 if (g == czero) then
                    cs = one
                    sn = czero
                    r = f
                    go to 50
                 end if
                 20 continue
                 count = count - 1
                 fs = fs*safmx2
                 gs = gs*safmx2
                 scale = scale*safmx2
                 if (scale <= safmn2) go to 20
              end if
              f2 = abssq(fs)
              g2 = abssq(gs)
              if (f2 <= max(g2,one)*safmin) then
                 ! this is a rare case: f is very small.
                 if (f == czero) then
                    cs = zero
                    r = la_qlapy2(real(g,KIND=qp),aimag(g))
                    ! do complex/real division explicitly with two real
                    ! divisions
                    d = la_qlapy2(real(gs,KIND=qp),aimag(gs))
                    sn = cmplx(real(gs,KIND=qp)/d,-aimag(gs)/d,KIND=qp)
                    go to 50
                 end if
                 f2s = la_qlapy2(real(fs,KIND=qp),aimag(fs))
                 ! g2 and g2s are accurate
                 ! g2 is at least safmin, and g2s is at least safmn2
                 g2s = sqrt(g2)
                 ! error in cs from underflow in f2s is at most
                 ! unfl / safmn2 .lt. sqrt(unfl*eps) .lt. eps
                 ! if max(g2,one)=g2, then f2 .lt. g2*safmin,
                 ! and so cs .lt. sqrt(safmin)
                 ! if max(g2,one)=one, then f2 .lt. safmin
                 ! and so cs .lt. sqrt(safmin)/safmn2 = sqrt(eps)
                 ! therefore, cs = f2s/g2s / sqrt( 1 + (f2s/g2s)**2 ) = f2s/g2s
                 cs = f2s/g2s
                 ! make sure abs(ff) = 1
                 ! do complex/real division explicitly with 2 real divisions
                 if (abs1(f) > one) then
                    d = la_qlapy2(real(f,KIND=qp),aimag(f))
                    ff = cmplx(real(f,KIND=qp)/d,aimag(f)/d,KIND=qp)
                 else
                    dr = safmx2*real(f,KIND=qp)
                    di = safmx2*aimag(f)
                    d = la_qlapy2(dr,di)
                    ff = cmplx(dr/d,di/d,KIND=qp)
                 end if
                 sn = ff*cmplx(real(gs,KIND=qp)/g2s,-aimag(gs)/g2s,KIND=qp)
                 r = cs*f + sn*g
              else
                 ! this is the most common case.
                 ! neither f2 nor f2/g2 are less than safmin
                 ! f2s cannot overflow, and it is accurate
                 f2s = sqrt(one + g2/f2)
                 ! do the f2s(real)*fs(complex) multiply with two real
                 ! multiplies
                 r = cmplx(f2s*real(fs,KIND=qp),f2s*aimag(fs),KIND=qp)
                 cs = one/f2s
                 d = f2 + g2
                 ! do complex/real division explicitly with two real divisions
                 sn = cmplx(real(r,KIND=qp)/d,aimag(r)/d,KIND=qp)
                 sn = sn*conjg(gs)
                 if (count /= 0) then
                    if (count > 0) then
                       do j = 1,count
                          r = r*safmx2
                       end do
                    else
                       do j = 1,-count
                          r = r*safmn2
                       end do
                    end if
                 end if
              end if
              50 continue
              c(ic) = cs
              y(iy) = sn
              x(ix) = r
              ic = ic + incc
              iy = iy + incy
              ix = ix + incx
           end do loop_60
           return
     end subroutine la_wlargv
#endif

end module la_lapack_givens_jacobi_rot
