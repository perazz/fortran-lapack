!> BLAS-like level 2: matrix-vector products, scaling, rank updates
module la_lapack_blas_like_l2
     use la_constants
     use la_blas_aux
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_scalar
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sla_wwaddw
     public :: la_sla_gbamv
     public :: la_sla_geamv
     public :: la_slascl
     public :: la_dla_wwaddw
     public :: la_dla_gbamv
     public :: la_dla_geamv
     public :: la_dlascl
#ifdef LA_WITH_XDP
     public :: la_xla_wwaddw
     public :: la_xla_gbamv
     public :: la_xla_geamv
     public :: la_xlascl
#endif
#ifdef LA_WITH_QP
     public :: la_qla_wwaddw
     public :: la_qla_gbamv
     public :: la_qla_geamv
     public :: la_qlascl
#endif
     public :: la_cla_gbamv
     public :: la_cla_geamv
     public :: la_cla_heamv
     public :: la_cla_wwaddw
     public :: la_clascl
     public :: la_cspmv
     public :: la_cspr
     public :: la_csymv
     public :: la_csyr
     public :: la_zla_gbamv
     public :: la_zla_geamv
     public :: la_zla_heamv
     public :: la_zla_wwaddw
     public :: la_zlascl
     public :: la_zspmv
     public :: la_zspr
     public :: la_zsymv
     public :: la_zsyr
#ifdef LA_WITH_XDP
     public :: la_yla_gbamv
     public :: la_yla_geamv
     public :: la_yla_heamv
     public :: la_yla_wwaddw
     public :: la_ylascl
     public :: la_yspmv
     public :: la_yspr
     public :: la_ysymv
     public :: la_ysyr
#endif
#ifdef LA_WITH_QP
     public :: la_wla_gbamv
     public :: la_wla_geamv
     public :: la_wla_heamv
     public :: la_wla_wwaddw
     public :: la_wlascl
     public :: la_wspmv
     public :: la_wspr
     public :: la_wsymv
     public :: la_wsyr
#endif

     contains

     !> SLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_sla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: x(*),y(*)
           real(sp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_sla_wwaddw
     !> DLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_dla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: x(*),y(*)
           real(dp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_dla_wwaddw
#ifdef LA_WITH_XDP
     !> XLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_xla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(inout) :: x(*),y(*)
           real(xdp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           real(xdp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_xla_wwaddw
#endif
#ifdef LA_WITH_QP
     !> QLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_qla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: x(*),y(*)
           real(qp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_qla_wwaddw
#endif

     !> SLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_sla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*),x(*)
           real(sp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(sp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('SLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_sla_gbamv
     !> DLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_dla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*),x(*)
           real(dp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(dp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('DLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_dla_gbamv
#ifdef LA_WITH_XDP
     !> XLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_xla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           real(xdp),intent(in) :: ab(ldab,*),x(*)
           real(xdp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(xdp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('XLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_xlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_xla_gbamv
#endif
#ifdef LA_WITH_QP
     !> QLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_qla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*),x(*)
           real(qp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(qp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('QLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = abs(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_qla_gbamv
#endif

     !> SLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_sla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n,trans
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),x(*)
           real(sp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(sp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('SLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_sla_geamv
     !> DLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_dla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n,trans
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),x(*)
           real(dp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(dp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('DLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_dla_geamv
#ifdef LA_WITH_XDP
     !> XLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_xla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n,trans
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),x(*)
           real(xdp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(xdp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('XLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_xlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_xla_geamv
#endif
#ifdef LA_WITH_QP
     !> QLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_qla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n,trans
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),x(*)
           real(qp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(qp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('QLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       jx = kx
                       do j = 1,lenx
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_qla_geamv
#endif

     !> SLASCL: multiplies the M by N real matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_slascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(sp),intent(in) :: cfrom,cto
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(sp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_sisnan(cfrom)) then
              info = -4
           else if (la_sisnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_slamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower half of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper half of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_slascl
     !> DLASCL: multiplies the M by N real matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_dlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(dp),intent(in) :: cfrom,cto
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(dp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_disnan(cfrom)) then
              info = -4
           else if (la_disnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower half of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper half of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_dlascl
#ifdef LA_WITH_XDP
     !> XLASCL: multiplies the M by N real matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_xlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(xdp),intent(in) :: cfrom,cto
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(xdp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_xisnan(cfrom)) then
              info = -4
           else if (la_xisnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('XLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_xlamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower half of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper half of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_xlascl
#endif
#ifdef LA_WITH_QP
     !> QLASCL: multiplies the M by N real matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_qlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(qp),intent(in) :: cfrom,cto
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(qp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_qisnan(cfrom)) then
              info = -4
           else if (la_qisnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower half of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper half of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_qlascl
#endif

     !> CLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_cla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           complex(sp),intent(in) :: ab(ldab,*),x(*)
           real(sp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(sp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           complex(sp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('CLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_cla_gbamv
     !> ZLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_zla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           complex(dp),intent(in) :: ab(ldab,*),x(*)
           real(dp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(dp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           complex(dp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_zla_gbamv
#ifdef LA_WITH_XDP
     !> YLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_yla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           complex(xdp),intent(in) :: ab(ldab,*),x(*)
           real(xdp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(xdp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           complex(xdp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('YLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_xlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_yla_gbamv
#endif
#ifdef LA_WITH_QP
     !> WLA_GBAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_wla_gbamv(trans,m,n,kl,ku,alpha,ab,ldab,x,incx,beta,y,incy)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,ldab,m,n,kl,ku,trans
           ! Array Arguments
           complex(qp),intent(in) :: ab(ldab,*),x(*)
           real(qp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(qp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny,kd,ke
           complex(qp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (kl < 0 .or. kl > m - 1) then
              info = 4
           else if (ku < 0 .or. ku > n - 1) then
              info = 5
           else if (ldab < kl + ku + 1) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('WLA_GBAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           kd = ku + 1
           ke = kl + 1
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(kd + i - j,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = max(i - kl,1),min(i + ku,lenx)
                          temp = cabs1(ab(ke - i + j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_wla_gbamv
#endif

     !> CLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_cla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n
           integer(ilp),intent(in) :: trans
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),x(*)
           real(sp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(sp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           complex(sp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('CLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == 0.0_sp) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == 0.0_sp) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= 0.0_sp) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_cla_geamv
     !> ZLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_zla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n
           integer(ilp),intent(in) :: trans
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),x(*)
           real(dp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(dp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           complex(dp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_zla_geamv
#ifdef LA_WITH_XDP
     !> YLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_yla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n
           integer(ilp),intent(in) :: trans
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),x(*)
           real(xdp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(xdp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           complex(xdp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=xdp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('YLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_xlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_yla_geamv
#endif
#ifdef LA_WITH_QP
     !> WLA_GEAMV:  performs one of the matrix-vector operations
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> m by n matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_wla_geamv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,m,n
           integer(ilp),intent(in) :: trans
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),x(*)
           real(qp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(qp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky,lenx,leny
           complex(qp) :: cdum
           ! Intrinsic Functions
           intrinsic :: max,abs,real,aimag,sign
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. ((trans == la_ilatrans('N')) .or. (trans == la_ilatrans('T')) &
                     .or. (trans == la_ilatrans('C')))) then
              info = 1
           else if (m < 0) then
              info = 2
           else if (n < 0) then
              info = 3
           else if (lda < max(1,m)) then
              info = 6
           else if (incx == 0) then
              info = 8
           else if (incy == 0) then
              info = 11
           end if
           if (info /= 0) then
              call la_xerbla('WLA_GEAMV ',info)
              return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (trans == la_ilatrans('N')) then
              lenx = n
              leny = m
           else
              lenx = m
              leny = n
           end if
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (lenx - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (leny - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(m*n) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (trans == la_ilatrans('N')) then
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,leny
                    if (beta == czero) then
                       symb_zero = .true.
                       y(iy) = czero
                    else if (y(iy) == czero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= czero) then
                       jx = kx
                       do j = 1,lenx
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(jx) == czero .or. temp == czero)

                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_wla_geamv
#endif

     !> CLA_SYAMV  performs the matrix-vector operation
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> n by n symmetric matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_cla_heamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),x(*)
           real(sp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(sp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           complex(sp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,abs,sign,real,aimag
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (uplo /= la_ilauplo('U') .and. uplo /= la_ilauplo('L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('CHEMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(n^2) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_cla_heamv
     !> ZLA_SYAMV  performs the matrix-vector operation
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> n by n symmetric matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_zla_heamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),x(*)
           real(dp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(dp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,abs,sign,real,aimag
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (uplo /= la_ilauplo('U') .and. uplo /= la_ilauplo('L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('ZHEMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(n^2) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_zla_heamv
#ifdef LA_WITH_XDP
     !> YLA_SYAMV  performs the matrix-vector operation
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> n by n symmetric matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_yla_heamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),x(*)
           real(xdp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(xdp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           complex(xdp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,abs,sign,real,aimag
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (uplo /= la_ilauplo('U') .and. uplo /= la_ilauplo('L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('YHEMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_xlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(n^2) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_yla_heamv
#endif
#ifdef LA_WITH_QP
     !> WLA_SYAMV  performs the matrix-vector operation
     !> y := alpha*abs(A)*abs(x) + beta*abs(y),
     !> where alpha and beta are scalars, x and y are vectors and A is an
     !> n by n symmetric matrix.
     !> This function is primarily used in calculating error bounds.
     !> To protect against underflow during evaluation, components in
     !> the resulting vector are perturbed away from zero by (N+1)
     !> times the underflow threshold.  To prevent unnecessarily large
     !> errors for block-structure embedded in general matrices,
     !> "symbolically" zero components are not perturbed.  A zero
     !> entry is considered "symbolic" if all multiplications involved
     !> in computing that entry have at least one zero multiplicand.

     subroutine la_wla_heamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),x(*)
           real(qp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(qp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,abs,sign,real,aimag
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (uplo /= la_ilauplo('U') .and. uplo /= la_ilauplo('L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('WHEMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! set safe1 essentially to be the underflow threshold times the
           ! number of additions in each row.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (n + 1)*safe1
           ! form  y := alpha*abs(a)*abs(x) + beta*abs(y).
           ! the o(n^2) symb_zero tests could be replaced by o(n) queries to
           ! the inexact flag.  still doesn't help change the iteration order
           ! to per-column.
           iy = ky
           if (incx == 1) then
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(j))*temp
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           else
              if (uplo == la_ilauplo('U')) then
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              else
                 do i = 1,n
                    if (beta == zero) then
                       symb_zero = .true.
                       y(iy) = zero
                    else if (y(iy) == zero) then
                       symb_zero = .true.
                    else
                       symb_zero = .false.
                       y(iy) = beta*abs(y(iy))
                    end if
                    jx = kx
                    if (alpha /= zero) then
                       do j = 1,i
                          temp = cabs1(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = cabs1(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*cabs1(x(jx))*temp
                          jx = jx + incx
                       end do
                    end if
                    if (.not. symb_zero) y(iy) = y(iy) + sign(safe1,y(iy))
                    iy = iy + incy
                 end do
              end if
           end if
           return
     end subroutine la_wla_heamv
#endif

     !> CLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_cla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(sp),intent(inout) :: x(*),y(*)
           complex(sp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           complex(sp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_cla_wwaddw
     !> ZLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_zla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(dp),intent(inout) :: x(*),y(*)
           complex(dp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           complex(dp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_zla_wwaddw
#ifdef LA_WITH_XDP
     !> YLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_yla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(xdp),intent(inout) :: x(*),y(*)
           complex(xdp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           complex(xdp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_yla_wwaddw
#endif
#ifdef LA_WITH_QP
     !> WLA_WWADDW: adds a vector W into a doubled-single vector (X, Y).
     !> This works for all extant IBM's hex and binary floating point
     !> arithmetic, but not for decimal.

     pure subroutine la_wla_wwaddw(n,x,y,w)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(qp),intent(inout) :: x(*),y(*)
           complex(qp),intent(in) :: w(*)
        ! =====================================================================
           ! Local Scalars
           complex(qp) :: s
           integer(ilp) :: i
           ! Executable Statements
           do 10 i = 1,n
             s = x(i) + w(i)
             s = (s + s) - s
             y(i) = ((x(i) - s) + w(i)) + y(i)
             x(i) = s
             10 continue
           return
     end subroutine la_wla_wwaddw
#endif

     !> CLASCL: multiplies the M by N complex matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_clascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(sp),intent(in) :: cfrom,cto
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(sp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_sisnan(cfrom)) then
              info = -4
           else if (la_sisnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_slamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower chalf of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper chalf of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_clascl
     !> ZLASCL: multiplies the M by N complex matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_zlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(dp),intent(in) :: cfrom,cto
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(dp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_disnan(cfrom)) then
              info = -4
           else if (la_disnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower chalf of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper chalf of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_zlascl
#ifdef LA_WITH_XDP
     !> YLASCL: multiplies the M by N complex matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_ylascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(xdp),intent(in) :: cfrom,cto
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(xdp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_xisnan(cfrom)) then
              info = -4
           else if (la_xisnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('YLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_xlamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower chalf of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper chalf of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_ylascl
#endif
#ifdef LA_WITH_QP
     !> WLASCL: multiplies the M by N complex matrix A by the real scalar
     !> CTO/CFROM.  This is done without over/underflow as long as the final
     !> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
     !> A may be full, upper triangular, lower triangular, upper Hessenberg,
     !> or banded.

     pure subroutine la_wlascl(type,kl,ku,cfrom,cto,m,n,a,lda,info)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: type
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,lda,m,n
           real(qp),intent(in) :: cfrom,cto
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: i,itype,j,k1,k2,k3,k4
           real(qp) :: bignum,cfrom1,cfromc,cto1,ctoc,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           if (la_lsame(type,'G')) then
              itype = 0
           else if (la_lsame(type,'L')) then
              itype = 1
           else if (la_lsame(type,'U')) then
              itype = 2
           else if (la_lsame(type,'H')) then
              itype = 3
           else if (la_lsame(type,'B')) then
              itype = 4
           else if (la_lsame(type,'Q')) then
              itype = 5
           else if (la_lsame(type,'Z')) then
              itype = 6
           else
              itype = -1
           end if
           if (itype == -1) then
              info = -1
           else if (cfrom == zero .or. la_qisnan(cfrom)) then
              info = -4
           else if (la_qisnan(cto)) then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0 .or. (itype == 4 .and. n /= m) .or. (itype == 5 .and. n /= m)) then
              info = -7
           else if (itype <= 3 .and. lda < max(1,m)) then
              info = -9
           else if (itype >= 4) then
              if (kl < 0 .or. kl > max(m - 1,0)) then
                 info = -2
              else if (ku < 0 .or. ku > max(n - 1,0) .or. ((itype == 4 .or. itype == 5) .and. kl /= ku) &
                        ) then
                 info = -3
              else if ((itype == 4 .and. lda < kl + 1) .or. (itype == 5 .and. lda < ku + 1) .or. (itype == 6 &
                        .and. lda < 2*kl + ku + 1)) then
                 info = -9
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WLASCL',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) return
           ! get machine parameters
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           cfromc = cfrom
           ctoc = cto
           10 continue
           cfrom1 = cfromc*smlnum
           if (cfrom1 == cfromc) then
              ! cfromc is an inf.  multiply by a correctly signed zero for
              ! finite ctoc, or a nan if ctoc is infinite.
              mul = ctoc/cfromc
              done = .true.
              cto1 = ctoc
           else
              cto1 = ctoc/bignum
              if (cto1 == ctoc) then
                 ! ctoc is either 0 or an inf.  in both cases, ctoc itself
                 ! serves as the correct multiplication factor.
                 mul = ctoc
                 done = .true.
                 cfromc = one
              else if (abs(cfrom1) > abs(ctoc) .and. ctoc /= zero) then
                 mul = smlnum
                 done = .false.
                 cfromc = cfrom1
              else if (abs(cto1) > abs(cfromc)) then
                 mul = bignum
                 done = .false.
                 ctoc = cto1
              else
                 mul = ctoc/cfromc
                 done = .true.
              end if
           end if
           if (itype == 0) then
              ! full matrix
              do j = 1,n
                 do i = 1,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 1) then
              ! lower triangular matrix
              do j = 1,n
                 do i = j,m
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 2) then
              ! upper triangular matrix
              do j = 1,n
                 do i = 1,min(j,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 3) then
              ! upper hessenberg matrix
              do j = 1,n
                 do i = 1,min(j + 1,m)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 4) then
              ! lower chalf of a symmetric band matrix
              k3 = kl + 1
              k4 = n + 1
              do j = 1,n
                 do i = 1,min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 5) then
              ! upper chalf of a symmetric band matrix
              k1 = ku + 2
              k3 = ku + 1
              do j = 1,n
                 do i = max(k1 - j,1),k3
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           else if (itype == 6) then
              ! band matrix
              k1 = kl + ku + 2
              k2 = kl + 1
              k3 = 2*kl + ku + 1
              k4 = kl + ku + 1 + m
              do j = 1,n
                 do i = max(k1 - j,k2),min(k3,k4 - j)
                    a(i,j) = a(i,j)*mul
                 end do
              end do
           end if
           if (.not. done) go to 10
           return
     end subroutine la_wlascl
#endif

     !> CSPMV:  performs the matrix-vector operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix, supplied in packed form.

     pure subroutine la_cspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,n
           complex(sp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(sp),intent(in) :: ap(*),x(*)
           complex(sp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           complex(sp) :: temp1,temp2
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 6
           else if (incy == 0) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('CSPMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  y  when ap contains the upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    k = kk
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + temp1*ap(kk + j - 1) + alpha*temp2
                    kk = kk + j
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do k = kk,kk + j - 2
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*ap(kk + j - 1) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                 end do
              end if
           else
              ! form  y  when ap contains the lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*ap(kk)
                    k = kk + 1
                    do i = j + 1,n
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + alpha*temp2
                    kk = kk + (n - j + 1)
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*ap(kk)
                    ix = jx
                    iy = jy
                    do k = kk + 1,kk + n - j
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + (n - j + 1)
                 end do
              end if
           end if
           return
     end subroutine la_cspmv
     !> ZSPMV:  performs the matrix-vector operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix, supplied in packed form.

     pure subroutine la_zspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,n
           complex(dp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(dp),intent(in) :: ap(*),x(*)
           complex(dp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           complex(dp) :: temp1,temp2
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 6
           else if (incy == 0) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('ZSPMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  y  when ap contains the upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    k = kk
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + temp1*ap(kk + j - 1) + alpha*temp2
                    kk = kk + j
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do k = kk,kk + j - 2
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*ap(kk + j - 1) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                 end do
              end if
           else
              ! form  y  when ap contains the lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*ap(kk)
                    k = kk + 1
                    do i = j + 1,n
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + alpha*temp2
                    kk = kk + (n - j + 1)
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*ap(kk)
                    ix = jx
                    iy = jy
                    do k = kk + 1,kk + n - j
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + (n - j + 1)
                 end do
              end if
           end if
           return
     end subroutine la_zspmv
#ifdef LA_WITH_XDP
     !> YSPMV:  performs the matrix-vector operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix, supplied in packed form.

     pure subroutine la_yspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,n
           complex(xdp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(xdp),intent(in) :: ap(*),x(*)
           complex(xdp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           complex(xdp) :: temp1,temp2
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 6
           else if (incy == 0) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('YSPMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  y  when ap contains the upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    k = kk
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + temp1*ap(kk + j - 1) + alpha*temp2
                    kk = kk + j
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do k = kk,kk + j - 2
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*ap(kk + j - 1) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                 end do
              end if
           else
              ! form  y  when ap contains the lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*ap(kk)
                    k = kk + 1
                    do i = j + 1,n
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + alpha*temp2
                    kk = kk + (n - j + 1)
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*ap(kk)
                    ix = jx
                    iy = jy
                    do k = kk + 1,kk + n - j
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + (n - j + 1)
                 end do
              end if
           end if
           return
     end subroutine la_yspmv
#endif
#ifdef LA_WITH_QP
     !> WSPMV:  performs the matrix-vector operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix, supplied in packed form.

     pure subroutine la_wspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,n
           complex(qp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(qp),intent(in) :: ap(*),x(*)
           complex(qp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           complex(qp) :: temp1,temp2
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 6
           else if (incy == 0) then
              info = 9
           end if
           if (info /= 0) then
              call la_xerbla('WSPMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  y  when ap contains the upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    k = kk
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + temp1*ap(kk + j - 1) + alpha*temp2
                    kk = kk + j
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do k = kk,kk + j - 2
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*ap(kk + j - 1) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + j
                 end do
              end if
           else
              ! form  y  when ap contains the lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*ap(kk)
                    k = kk + 1
                    do i = j + 1,n
                       y(i) = y(i) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(i)
                       k = k + 1
                    end do
                    y(j) = y(j) + alpha*temp2
                    kk = kk + (n - j + 1)
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*ap(kk)
                    ix = jx
                    iy = jy
                    do k = kk + 1,kk + n - j
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*ap(k)
                       temp2 = temp2 + ap(k)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                    kk = kk + (n - j + 1)
                 end do
              end if
           end if
           return
     end subroutine la_wspmv
#endif

     !> CSPR:    performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix, supplied in packed form.

     pure subroutine la_cspr(uplo,n,alpha,x,incx,ap)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,n
           complex(sp),intent(in) :: alpha
           ! Array Arguments
           complex(sp),intent(inout) :: ap(*)
           complex(sp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           complex(sp) :: temp
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           end if
           if (info /= 0) then
              call la_xerbla('CSPR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  a  when upper triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       k = kk
                       do i = 1,j - 1
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(j)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    kk = kk + j
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do k = kk,kk + j - 2
                          ap(k) = ap(k) + x(ix)*temp
                          ix = ix + incx
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(jx)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    jx = jx + incx
                    kk = kk + j
                 end do
              end if
           else
              ! form  a  when lower triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       ap(kk) = ap(kk) + temp*x(j)
                       k = kk + 1
                       do i = j + 1,n
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    kk = kk + n - j + 1
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ap(kk) = ap(kk) + temp*x(jx)
                       ix = jx
                       do k = kk + 1,kk + n - j
                          ix = ix + incx
                          ap(k) = ap(k) + x(ix)*temp
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    jx = jx + incx
                    kk = kk + n - j + 1
                 end do
              end if
           end if
           return
     end subroutine la_cspr
     !> ZSPR:    performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix, supplied in packed form.

     pure subroutine la_zspr(uplo,n,alpha,x,incx,ap)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,n
           complex(dp),intent(in) :: alpha
           ! Array Arguments
           complex(dp),intent(inout) :: ap(*)
           complex(dp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           complex(dp) :: temp
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           end if
           if (info /= 0) then
              call la_xerbla('ZSPR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  a  when upper triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       k = kk
                       do i = 1,j - 1
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(j)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    kk = kk + j
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do k = kk,kk + j - 2
                          ap(k) = ap(k) + x(ix)*temp
                          ix = ix + incx
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(jx)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    jx = jx + incx
                    kk = kk + j
                 end do
              end if
           else
              ! form  a  when lower triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       ap(kk) = ap(kk) + temp*x(j)
                       k = kk + 1
                       do i = j + 1,n
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    kk = kk + n - j + 1
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ap(kk) = ap(kk) + temp*x(jx)
                       ix = jx
                       do k = kk + 1,kk + n - j
                          ix = ix + incx
                          ap(k) = ap(k) + x(ix)*temp
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    jx = jx + incx
                    kk = kk + n - j + 1
                 end do
              end if
           end if
           return
     end subroutine la_zspr
#ifdef LA_WITH_XDP
     !> YSPR:    performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix, supplied in packed form.

     pure subroutine la_yspr(uplo,n,alpha,x,incx,ap)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,n
           complex(xdp),intent(in) :: alpha
           ! Array Arguments
           complex(xdp),intent(inout) :: ap(*)
           complex(xdp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           complex(xdp) :: temp
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           end if
           if (info /= 0) then
              call la_xerbla('YSPR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  a  when upper triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       k = kk
                       do i = 1,j - 1
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(j)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    kk = kk + j
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do k = kk,kk + j - 2
                          ap(k) = ap(k) + x(ix)*temp
                          ix = ix + incx
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(jx)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    jx = jx + incx
                    kk = kk + j
                 end do
              end if
           else
              ! form  a  when lower triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       ap(kk) = ap(kk) + temp*x(j)
                       k = kk + 1
                       do i = j + 1,n
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    kk = kk + n - j + 1
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ap(kk) = ap(kk) + temp*x(jx)
                       ix = jx
                       do k = kk + 1,kk + n - j
                          ix = ix + incx
                          ap(k) = ap(k) + x(ix)*temp
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    jx = jx + incx
                    kk = kk + n - j + 1
                 end do
              end if
           end if
           return
     end subroutine la_yspr
#endif
#ifdef LA_WITH_QP
     !> WSPR:    performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix, supplied in packed form.

     pure subroutine la_wspr(uplo,n,alpha,x,incx,ap)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,n
           complex(qp),intent(in) :: alpha
           ! Array Arguments
           complex(qp),intent(inout) :: ap(*)
           complex(qp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           complex(qp) :: temp
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           end if
           if (info /= 0) then
              call la_xerbla('WSPR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           kk = 1
           if (la_lsame(uplo,'U')) then
              ! form  a  when upper triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       k = kk
                       do i = 1,j - 1
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(j)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    kk = kk + j
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do k = kk,kk + j - 2
                          ap(k) = ap(k) + x(ix)*temp
                          ix = ix + incx
                       end do
                       ap(kk + j - 1) = ap(kk + j - 1) + x(jx)*temp
                    else
                       ap(kk + j - 1) = ap(kk + j - 1)
                    end if
                    jx = jx + incx
                    kk = kk + j
                 end do
              end if
           else
              ! form  a  when lower triangle is stored in ap.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       ap(kk) = ap(kk) + temp*x(j)
                       k = kk + 1
                       do i = j + 1,n
                          ap(k) = ap(k) + x(i)*temp
                          k = k + 1
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    kk = kk + n - j + 1
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ap(kk) = ap(kk) + temp*x(jx)
                       ix = jx
                       do k = kk + 1,kk + n - j
                          ix = ix + incx
                          ap(k) = ap(k) + x(ix)*temp
                       end do
                    else
                       ap(kk) = ap(kk)
                    end if
                    jx = jx + incx
                    kk = kk + n - j + 1
                 end do
              end if
           end if
           return
     end subroutine la_wspr
#endif

     !> CSYMV:  performs the matrix-vector  operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix.

     pure subroutine la_csymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,lda,n
           complex(sp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),x(*)
           complex(sp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,kx,ky
           complex(sp) :: temp1,temp2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('CSYMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           if (la_lsame(uplo,'U')) then
              ! form  y  when a is stored in upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + temp1*a(j,j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do i = 1,j - 1
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           else
              ! form  y  when a is stored in lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*a(j,j)
                    do i = j + 1,n
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*a(j,j)
                    ix = jx
                    iy = jy
                    do i = j + 1,n
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           end if
           return
     end subroutine la_csymv
     !> ZSYMV:  performs the matrix-vector  operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix.

     pure subroutine la_zsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,lda,n
           complex(dp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),x(*)
           complex(dp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,kx,ky
           complex(dp) :: temp1,temp2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('ZSYMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           if (la_lsame(uplo,'U')) then
              ! form  y  when a is stored in upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + temp1*a(j,j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do i = 1,j - 1
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           else
              ! form  y  when a is stored in lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*a(j,j)
                    do i = j + 1,n
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*a(j,j)
                    ix = jx
                    iy = jy
                    do i = j + 1,n
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           end if
           return
     end subroutine la_zsymv
#ifdef LA_WITH_XDP
     !> YSYMV:  performs the matrix-vector  operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix.

     pure subroutine la_ysymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,lda,n
           complex(xdp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),x(*)
           complex(xdp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,kx,ky
           complex(xdp) :: temp1,temp2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('YSYMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           if (la_lsame(uplo,'U')) then
              ! form  y  when a is stored in upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + temp1*a(j,j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do i = 1,j - 1
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           else
              ! form  y  when a is stored in lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*a(j,j)
                    do i = j + 1,n
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*a(j,j)
                    ix = jx
                    iy = jy
                    do i = j + 1,n
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           end if
           return
     end subroutine la_ysymv
#endif
#ifdef LA_WITH_QP
     !> WSYMV:  performs the matrix-vector  operation
     !> y := alpha*A*x + beta*y,
     !> where alpha and beta are scalars, x and y are n element vectors and
     !> A is an n by n symmetric matrix.

     pure subroutine la_wsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,incy,lda,n
           complex(qp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),x(*)
           complex(qp),intent(inout) :: y(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,iy,j,jx,jy,kx,ky
           complex(qp) :: temp1,temp2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (lda < max(1,n)) then
              info = 5
           else if (incx == 0) then
              info = 7
           else if (incy == 0) then
              info = 10
           end if
           if (info /= 0) then
              call la_xerbla('WSYMV ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! set up the start points in  x  and  y.
           if (incx > 0) then
              kx = 1
           else
              kx = 1 - (n - 1)*incx
           end if
           if (incy > 0) then
              ky = 1
           else
              ky = 1 - (n - 1)*incy
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           ! first form  y := beta*y.
           if (beta /= cone) then
              if (incy == 1) then
                 if (beta == czero) then
                    do i = 1,n
                       y(i) = czero
                    end do
                 else
                    do i = 1,n
                       y(i) = beta*y(i)
                    end do
                 end if
              else
                 iy = ky
                 if (beta == czero) then
                    do i = 1,n
                       y(iy) = czero
                       iy = iy + incy
                    end do
                 else
                    do i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
                    end do
                 end if
              end if
           end if
           if (alpha == czero) return
           if (la_lsame(uplo,'U')) then
              ! form  y  when a is stored in upper triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    do i = 1,j - 1
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + temp1*a(j,j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    ix = kx
                    iy = ky
                    do i = 1,j - 1
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                       ix = ix + incx
                       iy = iy + incy
                    end do
                    y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           else
              ! form  y  when a is stored in lower triangle.
              if ((incx == 1) .and. (incy == 1)) then
                 do j = 1,n
                    temp1 = alpha*x(j)
                    temp2 = czero
                    y(j) = y(j) + temp1*a(j,j)
                    do i = j + 1,n
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(i)
                    end do
                    y(j) = y(j) + alpha*temp2
                 end do
              else
                 jx = kx
                 jy = ky
                 do j = 1,n
                    temp1 = alpha*x(jx)
                    temp2 = czero
                    y(jy) = y(jy) + temp1*a(j,j)
                    ix = jx
                    iy = jy
                    do i = j + 1,n
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + a(i,j)*x(ix)
                    end do
                    y(jy) = y(jy) + alpha*temp2
                    jx = jx + incx
                    jy = jy + incy
                 end do
              end if
           end if
           return
     end subroutine la_wsymv
#endif

     !> CSYR:   performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix.

     pure subroutine la_csyr(uplo,n,alpha,x,incx,a,lda)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,lda,n
           complex(sp),intent(in) :: alpha
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,kx
           complex(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           else if (lda < max(1,n)) then
              info = 7
           end if
           if (info /= 0) then
              call la_xerbla('CSYR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           if (la_lsame(uplo,'U')) then
              ! form  a  when a is stored in upper triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = 1,j
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do i = 1,j
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           else
              ! form  a  when a is stored in lower triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = j,n
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = jx
                       do i = j,n
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           end if
           return
     end subroutine la_csyr
     !> ZSYR:   performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix.

     pure subroutine la_zsyr(uplo,n,alpha,x,incx,a,lda)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,lda,n
           complex(dp),intent(in) :: alpha
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,kx
           complex(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           else if (lda < max(1,n)) then
              info = 7
           end if
           if (info /= 0) then
              call la_xerbla('ZSYR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           if (la_lsame(uplo,'U')) then
              ! form  a  when a is stored in upper triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = 1,j
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do i = 1,j
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           else
              ! form  a  when a is stored in lower triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = j,n
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = jx
                       do i = j,n
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           end if
           return
     end subroutine la_zsyr
#ifdef LA_WITH_XDP
     !> YSYR:   performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix.

     pure subroutine la_ysyr(uplo,n,alpha,x,incx,a,lda)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,lda,n
           complex(xdp),intent(in) :: alpha
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,kx
           complex(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           else if (lda < max(1,n)) then
              info = 7
           end if
           if (info /= 0) then
              call la_xerbla('YSYR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           if (la_lsame(uplo,'U')) then
              ! form  a  when a is stored in upper triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = 1,j
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do i = 1,j
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           else
              ! form  a  when a is stored in lower triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = j,n
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = jx
                       do i = j,n
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           end if
           return
     end subroutine la_ysyr
#endif
#ifdef LA_WITH_QP
     !> WSYR:   performs the symmetric rank 1 operation
     !> A := alpha*x*x**H + A,
     !> where alpha is a complex scalar, x is an n element vector and A is an
     !> n by n symmetric matrix.

     pure subroutine la_wsyr(uplo,n,alpha,x,incx,a,lda)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: incx,lda,n
           complex(qp),intent(in) :: alpha
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: x(*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,info,ix,j,jx,kx
           complex(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = 1
           else if (n < 0) then
              info = 2
           else if (incx == 0) then
              info = 5
           else if (lda < max(1,n)) then
              info = 7
           end if
           if (info /= 0) then
              call la_xerbla('WSYR  ',info)
              return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
              kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
              kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           if (la_lsame(uplo,'U')) then
              ! form  a  when a is stored in upper triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = 1,j
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = kx
                       do i = 1,j
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           else
              ! form  a  when a is stored in lower triangle.
              if (incx == 1) then
                 do j = 1,n
                    if (x(j) /= czero) then
                       temp = alpha*x(j)
                       do i = j,n
                          a(i,j) = a(i,j) + x(i)*temp
                       end do
                    end if
                 end do
              else
                 jx = kx
                 do j = 1,n
                    if (x(jx) /= czero) then
                       temp = alpha*x(jx)
                       ix = jx
                       do i = j,n
                          a(i,j) = a(i,j) + x(ix)*temp
                          ix = ix + incx
                       end do
                    end if
                    jx = jx + incx
                 end do
              end if
           end if
           return
     end subroutine la_wsyr
#endif

end module la_lapack_blas_like_l2
