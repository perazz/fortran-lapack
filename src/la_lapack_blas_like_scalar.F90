!> BLAS-like scalar: complex division, Pythagorean sums, NaN tests
module la_lapack_blas_like_scalar
     use la_constants
     use la_lapack_auxiliary
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slaisnan
     public :: la_slapy3
     public :: la_sisnan
     public :: la_slapy2
     public :: la_sladiv
     public :: la_dlaisnan
     public :: la_dlapy3
     public :: la_disnan
     public :: la_dlapy2
     public :: la_dladiv
#ifdef LA_WITH_XDP
     public :: la_xlaisnan
     public :: la_xlapy3
     public :: la_xisnan
     public :: la_xlapy2
     public :: la_xladiv
#endif
#ifdef LA_WITH_QP
     public :: la_qlaisnan
     public :: la_qlapy3
     public :: la_qisnan
     public :: la_qlapy2
     public :: la_qladiv
#endif
     public :: la_cladiv
     public :: la_zladiv
#ifdef LA_WITH_XDP
     public :: la_yladiv
#endif
#ifdef LA_WITH_QP
     public :: la_wladiv
#endif

     contains

     !> This routine is not for general use.  It exists solely to avoid
     !> over-optimization in SISNAN.
     !> SLAISNAN: checks for NaNs by comparing its two arguments for
     !> inequality.  NaN is the only floating-point value where NaN != NaN
     !> returns .TRUE.  To check for NaNs, pass the same variable as both
     !> arguments.
     !> A compiler must assume that the two arguments are
     !> not the same variable, and the test will not be optimized away.
     !> Interprocedural or whole-program optimization may delete this
     !> test.  The ISNAN functions will be replaced by the correct
     !> Fortran 03 intrinsic once the intrinsic is widely available.

     pure logical(lk) function la_slaisnan(sin1,sin2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: sin1,sin2
        ! =====================================================================
        ! Executable Statements
           la_slaisnan = (sin1 /= sin2)
           return
     end function la_slaisnan
     !> This routine is not for general use.  It exists solely to avoid
     !> over-optimization in DISNAN.
     !> DLAISNAN: checks for NaNs by comparing its two arguments for
     !> inequality.  NaN is the only floating-point value where NaN != NaN
     !> returns .TRUE.  To check for NaNs, pass the same variable as both
     !> arguments.
     !> A compiler must assume that the two arguments are
     !> not the same variable, and the test will not be optimized away.
     !> Interprocedural or whole-program optimization may delete this
     !> test.  The ISNAN functions will be replaced by the correct
     !> Fortran 03 intrinsic once the intrinsic is widely available.

     pure logical(lk) function la_dlaisnan(din1,din2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: din1,din2
        ! =====================================================================
        ! Executable Statements
           la_dlaisnan = (din1 /= din2)
           return
     end function la_dlaisnan
#ifdef LA_WITH_XDP
     !> This routine is not for general use.  It exists solely to avoid
     !> over-optimization in XISNAN.
     !> XLAISNAN: checks for NaNs by comparing its two arguments for
     !> inequality.  NaN is the only floating-point value where NaN != NaN
     !> returns .TRUE.  To check for NaNs, pass the same variable as both
     !> arguments.
     !> A compiler must assume that the two arguments are
     !> not the same variable, and the test will not be optimized away.
     !> Interprocedural or whole-program optimization may delete this
     !> test.  The ISNAN functions will be replaced by the correct
     !> Fortran 03 intrinsic once the intrinsic is widely available.

     pure logical(lk) function la_xlaisnan(xin1,xin2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: xin1,xin2
        ! =====================================================================
        ! Executable Statements
           la_xlaisnan = (xin1 /= xin2)
           return
     end function la_xlaisnan
#endif
#ifdef LA_WITH_QP
     !> This routine is not for general use.  It exists solely to avoid
     !> over-optimization in QISNAN.
     !> QLAISNAN: checks for NaNs by comparing its two arguments for
     !> inequality.  NaN is the only floating-point value where NaN != NaN
     !> returns .TRUE.  To check for NaNs, pass the same variable as both
     !> arguments.
     !> A compiler must assume that the two arguments are
     !> not the same variable, and the test will not be optimized away.
     !> Interprocedural or whole-program optimization may delete this
     !> test.  The ISNAN functions will be replaced by the correct
     !> Fortran 03 intrinsic once the intrinsic is widely available.

     pure logical(lk) function la_qlaisnan(qin1,qin2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: qin1,qin2
        ! =====================================================================
        ! Executable Statements
           la_qlaisnan = (qin1 /= qin2)
           return
     end function la_qlaisnan
#endif

     !> SLAPY3: returns sqrt(x**2+y**2+z**2), taking care not to cause
     !> unnecessary overflow and unnecessary underflow.

     pure real(sp) function la_slapy3(x,y,z)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: x,y,z
        ! =====================================================================

           ! Local Scalars
           real(sp) :: w,xabs,yabs,zabs,hugeval
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           hugeval = la_slamch('OVERFLOW')
           xabs = abs(x)
           yabs = abs(y)
           zabs = abs(z)
           w = max(xabs,yabs,zabs)
           if (w == zero .or. w > hugeval) then
           ! w can be zero for max(0,nan,0)
           ! adding all three entries together will make sure
           ! nan will not disappear.
              la_slapy3 = xabs + yabs + zabs
           else
              la_slapy3 = w*sqrt((xabs/w)**2 + (yabs/w)**2 + (zabs/w)**2)
           end if
           return
     end function la_slapy3
     !> DLAPY3: returns sqrt(x**2+y**2+z**2), taking care not to cause
     !> unnecessary overflow and unnecessary underflow.

     pure real(dp) function la_dlapy3(x,y,z)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: x,y,z
        ! =====================================================================

           ! Local Scalars
           real(dp) :: w,xabs,yabs,zabs,hugeval
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           hugeval = la_dlamch('OVERFLOW')
           xabs = abs(x)
           yabs = abs(y)
           zabs = abs(z)
           w = max(xabs,yabs,zabs)
           if (w == zero .or. w > hugeval) then
           ! w can be zero for max(0,nan,0)
           ! adding all three entries together will make sure
           ! nan will not disappear.
              la_dlapy3 = xabs + yabs + zabs
           else
              la_dlapy3 = w*sqrt((xabs/w)**2 + (yabs/w)**2 + (zabs/w)**2)
           end if
           return
     end function la_dlapy3
#ifdef LA_WITH_XDP
     !> XLAPY3: returns sqrt(x**2+y**2+z**2), taking care not to cause
     !> unnecessary overflow and unnecessary underflow.

     pure real(xdp) function la_xlapy3(x,y,z)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: x,y,z
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: w,xabs,yabs,zabs,hugeval
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           hugeval = la_xlamch('OVERFLOW')
           xabs = abs(x)
           yabs = abs(y)
           zabs = abs(z)
           w = max(xabs,yabs,zabs)
           if (w == zero .or. w > hugeval) then
           ! w can be zero for max(0,nan,0)
           ! adding all three entries together will make sure
           ! nan will not disappear.
              la_xlapy3 = xabs + yabs + zabs
           else
              la_xlapy3 = w*sqrt((xabs/w)**2 + (yabs/w)**2 + (zabs/w)**2)
           end if
           return
     end function la_xlapy3
#endif
#ifdef LA_WITH_QP
     !> QLAPY3: returns sqrt(x**2+y**2+z**2), taking care not to cause
     !> unnecessary overflow and unnecessary underflow.

     pure real(qp) function la_qlapy3(x,y,z)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: x,y,z
        ! =====================================================================

           ! Local Scalars
           real(qp) :: w,xabs,yabs,zabs,hugeval
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           hugeval = la_qlamch('OVERFLOW')
           xabs = abs(x)
           yabs = abs(y)
           zabs = abs(z)
           w = max(xabs,yabs,zabs)
           if (w == zero .or. w > hugeval) then
           ! w can be zero for max(0,nan,0)
           ! adding all three entries together will make sure
           ! nan will not disappear.
              la_qlapy3 = xabs + yabs + zabs
           else
              la_qlapy3 = w*sqrt((xabs/w)**2 + (yabs/w)**2 + (zabs/w)**2)
           end if
           return
     end function la_qlapy3
#endif

     !> SISNAN: returns .TRUE. if its argument is NaN, and .FALSE.
     !> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
     !> future.

     pure logical(lk) function la_sisnan(sin)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: sin
        ! =====================================================================
        ! Executable Statements
           la_sisnan = la_slaisnan(sin,sin)
           return
     end function la_sisnan
     !> DISNAN: returns .TRUE. if its argument is NaN, and .FALSE.
     !> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
     !> future.

     pure logical(lk) function la_disnan(din)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: din
        ! =====================================================================
        ! Executable Statements
           la_disnan = la_dlaisnan(din,din)
           return
     end function la_disnan
#ifdef LA_WITH_XDP
     !> XISNAN: returns .TRUE. if its argument is NaN, and .FALSE.
     !> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
     !> future.

     pure logical(lk) function la_xisnan(xin)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: xin
        ! =====================================================================
        ! Executable Statements
           la_xisnan = la_xlaisnan(xin,xin)
           return
     end function la_xisnan
#endif
#ifdef LA_WITH_QP
     !> QISNAN: returns .TRUE. if its argument is NaN, and .FALSE.
     !> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
     !> future.

     pure logical(lk) function la_qisnan(qin)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: qin
        ! =====================================================================
        ! Executable Statements
           la_qisnan = la_qlaisnan(qin,qin)
           return
     end function la_qisnan
#endif

     !> SLAPY2: returns sqrt(x**2+y**2), taking care not to cause unnecessary
     !> overflow and unnecessary underflow.

     pure real(sp) function la_slapy2(x,y)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: x,y
        ! =====================================================================

           ! Local Scalars
           real(sp) :: w,xabs,yabs,z,hugeval
           logical(lk) :: x_is_nan,y_is_nan
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           x_is_nan = la_sisnan(x)
           y_is_nan = la_sisnan(y)
           if (x_is_nan) la_slapy2 = x
           if (y_is_nan) la_slapy2 = y
           hugeval = la_slamch('OVERFLOW')
           if (.not. (x_is_nan .or. y_is_nan)) then
              xabs = abs(x)
              yabs = abs(y)
              w = max(xabs,yabs)
              z = min(xabs,yabs)
              if (z == zero .or. w > hugeval) then
                 la_slapy2 = w
              else
                 la_slapy2 = w*sqrt(one + (z/w)**2)
              end if
           end if
           return
     end function la_slapy2
     !> DLAPY2: returns sqrt(x**2+y**2), taking care not to cause unnecessary
     !> overflow and unnecessary underflow.

     pure real(dp) function la_dlapy2(x,y)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: x,y
        ! =====================================================================

           ! Local Scalars
           real(dp) :: w,xabs,yabs,z,hugeval
           logical(lk) :: x_is_nan,y_is_nan
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           x_is_nan = la_disnan(x)
           y_is_nan = la_disnan(y)
           if (x_is_nan) la_dlapy2 = x
           if (y_is_nan) la_dlapy2 = y
           hugeval = la_dlamch('OVERFLOW')
           if (.not. (x_is_nan .or. y_is_nan)) then
              xabs = abs(x)
              yabs = abs(y)
              w = max(xabs,yabs)
              z = min(xabs,yabs)
              if (z == zero .or. w > hugeval) then
                 la_dlapy2 = w
              else
                 la_dlapy2 = w*sqrt(one + (z/w)**2)
              end if
           end if
           return
     end function la_dlapy2
#ifdef LA_WITH_XDP
     !> XLAPY2: returns sqrt(x**2+y**2), taking care not to cause unnecessary
     !> overflow and unnecessary underflow.

     pure real(xdp) function la_xlapy2(x,y)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: x,y
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: w,xabs,yabs,z,hugeval
           logical(lk) :: x_is_nan,y_is_nan
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           x_is_nan = la_xisnan(x)
           y_is_nan = la_xisnan(y)
           if (x_is_nan) la_xlapy2 = x
           if (y_is_nan) la_xlapy2 = y
           hugeval = la_xlamch('OVERFLOW')
           if (.not. (x_is_nan .or. y_is_nan)) then
              xabs = abs(x)
              yabs = abs(y)
              w = max(xabs,yabs)
              z = min(xabs,yabs)
              if (z == zero .or. w > hugeval) then
                 la_xlapy2 = w
              else
                 la_xlapy2 = w*sqrt(one + (z/w)**2)
              end if
           end if
           return
     end function la_xlapy2
#endif
#ifdef LA_WITH_QP
     !> QLAPY2: returns sqrt(x**2+y**2), taking care not to cause unnecessary
     !> overflow and unnecessary underflow.

     pure real(qp) function la_qlapy2(x,y)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: x,y
        ! =====================================================================

           ! Local Scalars
           real(qp) :: w,xabs,yabs,z,hugeval
           logical(lk) :: x_is_nan,y_is_nan
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           x_is_nan = la_qisnan(x)
           y_is_nan = la_qisnan(y)
           if (x_is_nan) la_qlapy2 = x
           if (y_is_nan) la_qlapy2 = y
           hugeval = la_qlamch('OVERFLOW')
           if (.not. (x_is_nan .or. y_is_nan)) then
              xabs = abs(x)
              yabs = abs(y)
              w = max(xabs,yabs)
              z = min(xabs,yabs)
              if (z == zero .or. w > hugeval) then
                 la_qlapy2 = w
              else
                 la_qlapy2 = w*sqrt(one + (z/w)**2)
              end if
           end if
           return
     end function la_qlapy2
#endif

     !> SLADIV: performs complex division in  real arithmetic
     !> a + i*b
     !> p + i*q = ---------
     !> c + i*d
     !> The algorithm is due to Michael Baudin and Robert L. Smith
     !> and can be found in the paper
     !> "A Robust Complex Division in Scilab"

     pure subroutine la_sladiv(a,b,c,d,p,q)
        use la_constants_sp,only:half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: a,b,c,d
           real(sp),intent(out) :: p,q
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: bs = 2.0_sp

           ! Local Scalars
           real(sp) :: aa,bb,cc,dd,ab,cd,s,ov,un,be,eps
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           aa = a
           bb = b
           cc = c
           dd = d
           ab = max(abs(a),abs(b))
           cd = max(abs(c),abs(d))
           s = one
           ov = la_slamch('OVERFLOW THRESHOLD')
           un = la_slamch('SAFE MINIMUM')
           eps = la_slamch('EPSILON')
           be = bs/(eps*eps)
           if (ab >= half*ov) then
              aa = half*aa
              bb = half*bb
              s = two*s
           end if
           if (cd >= half*ov) then
              cc = half*cc
              dd = half*dd
              s = half*s
           end if
           if (ab <= un*bs/eps) then
              aa = aa*be
              bb = bb*be
              s = s/be
           end if
           if (cd <= un*bs/eps) then
              cc = cc*be
              dd = dd*be
              s = s*be
           end if
           if (abs(d) <= abs(c)) then
              call la_sladiv1(aa,bb,cc,dd,p,q)
           else
              call la_sladiv1(bb,aa,dd,cc,p,q)
              q = -q
           end if
           p = p*s
           q = q*s
           return
     end subroutine la_sladiv
     !> DLADIV: performs complex division in  real arithmetic
     !> a + i*b
     !> p + i*q = ---------
     !> c + i*d
     !> The algorithm is due to Michael Baudin and Robert L. Smith
     !> and can be found in the paper
     !> "A Robust Complex Division in Scilab"

     pure subroutine la_dladiv(a,b,c,d,p,q)
        use la_constants_dp,only:half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: a,b,c,d
           real(dp),intent(out) :: p,q
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: bs = 2.0_dp

           ! Local Scalars
           real(dp) :: aa,bb,cc,dd,ab,cd,s,ov,un,be,eps
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           aa = a
           bb = b
           cc = c
           dd = d
           ab = max(abs(a),abs(b))
           cd = max(abs(c),abs(d))
           s = one
           ov = la_dlamch('OVERFLOW THRESHOLD')
           un = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('EPSILON')
           be = bs/(eps*eps)
           if (ab >= half*ov) then
              aa = half*aa
              bb = half*bb
              s = two*s
           end if
           if (cd >= half*ov) then
              cc = half*cc
              dd = half*dd
              s = half*s
           end if
           if (ab <= un*bs/eps) then
              aa = aa*be
              bb = bb*be
              s = s/be
           end if
           if (cd <= un*bs/eps) then
              cc = cc*be
              dd = dd*be
              s = s*be
           end if
           if (abs(d) <= abs(c)) then
              call la_dladiv1(aa,bb,cc,dd,p,q)
           else
              call la_dladiv1(bb,aa,dd,cc,p,q)
              q = -q
           end if
           p = p*s
           q = q*s
           return
     end subroutine la_dladiv
#ifdef LA_WITH_XDP
     !> XLADIV: performs complex division in  real arithmetic
     !> a + i*b
     !> p + i*q = ---------
     !> c + i*d
     !> The algorithm is due to Michael Baudin and Robert L. Smith
     !> and can be found in the paper
     !> "A Robust Complex Division in Scilab"

     pure subroutine la_xladiv(a,b,c,d,p,q)
        use la_constants_xdp,only:half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: a,b,c,d
           real(xdp),intent(out) :: p,q
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: bs = 2.0_xdp

           ! Local Scalars
           real(xdp) :: aa,bb,cc,dd,ab,cd,s,ov,un,be,eps
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           aa = a
           bb = b
           cc = c
           dd = d
           ab = max(abs(a),abs(b))
           cd = max(abs(c),abs(d))
           s = one
           ov = la_xlamch('OVERFLOW THRESHOLD')
           un = la_xlamch('SAFE MINIMUM')
           eps = la_xlamch('EPSILON')
           be = bs/(eps*eps)
           if (ab >= half*ov) then
              aa = half*aa
              bb = half*bb
              s = two*s
           end if
           if (cd >= half*ov) then
              cc = half*cc
              dd = half*dd
              s = half*s
           end if
           if (ab <= un*bs/eps) then
              aa = aa*be
              bb = bb*be
              s = s/be
           end if
           if (cd <= un*bs/eps) then
              cc = cc*be
              dd = dd*be
              s = s*be
           end if
           if (abs(d) <= abs(c)) then
              call la_xladiv1(aa,bb,cc,dd,p,q)
           else
              call la_xladiv1(bb,aa,dd,cc,p,q)
              q = -q
           end if
           p = p*s
           q = q*s
           return
     end subroutine la_xladiv
#endif
#ifdef LA_WITH_QP
     !> QLADIV: performs complex division in  real arithmetic
     !> a + i*b
     !> p + i*q = ---------
     !> c + i*d
     !> The algorithm is due to Michael Baudin and Robert L. Smith
     !> and can be found in the paper
     !> "A Robust Complex Division in Scilab"

     pure subroutine la_qladiv(a,b,c,d,p,q)
        use la_constants_qp,only:half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: a,b,c,d
           real(qp),intent(out) :: p,q
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: bs = 2.0_qp

           ! Local Scalars
           real(qp) :: aa,bb,cc,dd,ab,cd,s,ov,un,be,eps
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           aa = a
           bb = b
           cc = c
           dd = d
           ab = max(abs(a),abs(b))
           cd = max(abs(c),abs(d))
           s = one
           ov = la_qlamch('OVERFLOW THRESHOLD')
           un = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('EPSILON')
           be = bs/(eps*eps)
           if (ab >= half*ov) then
              aa = half*aa
              bb = half*bb
              s = two*s
           end if
           if (cd >= half*ov) then
              cc = half*cc
              dd = half*dd
              s = half*s
           end if
           if (ab <= un*bs/eps) then
              aa = aa*be
              bb = bb*be
              s = s/be
           end if
           if (cd <= un*bs/eps) then
              cc = cc*be
              dd = dd*be
              s = s*be
           end if
           if (abs(d) <= abs(c)) then
              call la_qladiv1(aa,bb,cc,dd,p,q)
           else
              call la_qladiv1(bb,aa,dd,cc,p,q)
              q = -q
           end if
           p = p*s
           q = q*s
           return
     end subroutine la_qladiv
#endif

     !> CLADIV: := X / Y, where X and Y are complex.  The computation of X / Y
     !> will not overflow on an intermediary step unless the results
     !> overflows.

     pure complex(sp) function la_cladiv(x,y)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: x,y
        ! =====================================================================
           ! Local Scalars
           real(sp) :: zi,zr
           ! Intrinsic Functions
           intrinsic :: aimag,cmplx,real
           ! Executable Statements
           call la_sladiv(real(x,KIND=sp),aimag(x),real(y,KIND=sp),aimag(y),zr,zi)

           la_cladiv = cmplx(zr,zi,KIND=sp)
           return
     end function la_cladiv
     !> ZLADIV = X / Y, where X and Y are complex.  The computation of X / Y
     !> will not overflow on an intermediary step unless the results
     !> overflows.

     pure complex(dp) function la_zladiv(x,y)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: x,y
        ! =====================================================================
           ! Local Scalars
           real(dp) :: zi,zr
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           call la_dladiv(real(x,KIND=dp),aimag(x),real(y,KIND=dp),aimag(y),zr,zi)

           la_zladiv = cmplx(zr,zi,KIND=dp)
           return
     end function la_zladiv
#ifdef LA_WITH_XDP
     !> YLADIV = X / Y, where X and Y are complex.  The computation of X / Y
     !> will not overflow on an intermediary step unless the results
     !> overflows.

     pure complex(xdp) function la_yladiv(x,y)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: x,y
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: zi,zr
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           call la_xladiv(real(x,KIND=xdp),aimag(x),real(y,KIND=xdp),aimag(y),zr,zi)

           la_yladiv = cmplx(zr,zi,KIND=xdp)
           return
     end function la_yladiv
#endif
#ifdef LA_WITH_QP
     !> WLADIV = X / Y, where X and Y are complex.  The computation of X / Y
     !> will not overflow on an intermediary step unless the results
     !> overflows.

     pure complex(qp) function la_wladiv(x,y)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: x,y
        ! =====================================================================
           ! Local Scalars
           real(qp) :: zi,zr
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           call la_qladiv(real(x,KIND=qp),aimag(x),real(y,KIND=qp),aimag(y),zr,zi)

           la_wladiv = cmplx(zr,zi,KIND=qp)
           return
     end function la_wladiv
#endif

end module la_lapack_blas_like_scalar
