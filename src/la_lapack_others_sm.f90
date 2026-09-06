!> Extra-precise refinement helpers: condition numbers and pivot growth
module la_lapack_others_sm
     use la_constants
     use la_blas_aux
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_solve_aux
     use la_lapack_solve_chol_comp
     use la_lapack_solve_ldl_comp
     use la_lapack_solve_ldl_comp3
     use la_lapack_solve_lu_comp
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sla_gerpvgrw
     public :: la_sla_syamv
     public :: la_sla_syrcond
     public :: la_sla_syrpvgrw
     public :: la_dla_gerpvgrw
     public :: la_dla_syamv
     public :: la_dla_syrcond
     public :: la_dla_syrpvgrw
     public :: la_qla_gerpvgrw
     public :: la_qla_syamv
     public :: la_qla_syrcond
     public :: la_qla_syrpvgrw
     public :: la_cla_gerpvgrw
     public :: la_cla_syamv
     public :: la_cla_gbrcond_c
     public :: la_cla_gercond_c
     public :: la_cla_hercond_c
     public :: la_cla_porcond_c
     public :: la_cla_syrcond_c
     public :: la_cla_syrpvgrw
     public :: la_zla_gerpvgrw
     public :: la_zla_syamv
     public :: la_zla_gbrcond_c
     public :: la_zla_gercond_c
     public :: la_zla_hercond_c
     public :: la_zla_porcond_c
     public :: la_zla_syrcond_c
     public :: la_zla_syrpvgrw
     public :: la_wla_gerpvgrw
     public :: la_wla_syamv
     public :: la_wla_gbrcond_c
     public :: la_wla_gercond_c
     public :: la_wla_hercond_c
     public :: la_wla_porcond_c
     public :: la_wla_syrcond_c
     public :: la_wla_syrpvgrw

     contains

     !> SLA_GERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     pure real(sp) function la_sla_gerpvgrw(n,ncols,a,lda,af,ldaf)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,ncols,lda,ldaf
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),af(ldaf,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: amax,umax,rpvgrw
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           rpvgrw = one
           do j = 1,ncols
              amax = zero
              umax = zero
              do i = 1,n
                 amax = max(abs(a(i,j)),amax)
              end do
              do i = 1,j
                 umax = max(abs(af(i,j)),umax)
              end do
              if (umax /= 0.0_sp) then
                 rpvgrw = min(amax/umax,rpvgrw)
              end if
           end do
           la_sla_gerpvgrw = rpvgrw
     end function la_sla_gerpvgrw
     !> DLA_GERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     pure real(dp) function la_dla_gerpvgrw(n,ncols,a,lda,af,ldaf)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,ncols,lda,ldaf
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),af(ldaf,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: amax,umax,rpvgrw
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           rpvgrw = one
           do j = 1,ncols
              amax = zero
              umax = zero
              do i = 1,n
                 amax = max(abs(a(i,j)),amax)
              end do
              do i = 1,j
                 umax = max(abs(af(i,j)),umax)
              end do
              if (umax /= zero) then
                 rpvgrw = min(amax/umax,rpvgrw)
              end if
           end do
           la_dla_gerpvgrw = rpvgrw
     end function la_dla_gerpvgrw
     !> QLA_GERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     pure real(qp) function la_qla_gerpvgrw(n,ncols,a,lda,af,ldaf)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,ncols,lda,ldaf
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),af(ldaf,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: amax,umax,rpvgrw
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           rpvgrw = one
           do j = 1,ncols
              amax = zero
              umax = zero
              do i = 1,n
                 amax = max(abs(a(i,j)),amax)
              end do
              do i = 1,j
                 umax = max(abs(af(i,j)),umax)
              end do
              if (umax /= zero) then
                 rpvgrw = min(amax/umax,rpvgrw)
              end if
           end do
           la_qla_gerpvgrw = rpvgrw
     end function la_qla_gerpvgrw

     !> SLA_SYAMV:  performs the matrix-vector operation
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

     subroutine la_sla_syamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),x(*)
           real(sp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(sp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
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
              call la_xerbla('SLA_SYAMV',info)
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
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
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
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                       do j = i + 1,n
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
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
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
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
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
     end subroutine la_sla_syamv
     !> DLA_SYAMV:  performs the matrix-vector operation
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

     subroutine la_dla_syamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),x(*)
           real(dp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(dp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
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
              call la_xerbla('DLA_SYAMV',info)
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
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
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
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                       do j = i + 1,n
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
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
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
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
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
     end subroutine la_dla_syamv
     !> QLA_SYAMV:  performs the matrix-vector operation
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

     subroutine la_qla_syamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),x(*)
           real(qp),intent(inout) :: y(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: symb_zero
           real(qp) :: temp,safe1
           integer(ilp) :: i,info,iy,j,jx,kx,ky
           ! Intrinsic Functions
           intrinsic :: max,abs,sign
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
              call la_xerbla('QLA_SYAMV',info)
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
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                       do j = i + 1,n
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
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
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(j))*temp
                       end do
                       do j = i + 1,n
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
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
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
                          temp = abs(a(i,j))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
                          y(iy) = y(iy) + alpha*abs(x(jx))*temp
                          jx = jx + incx
                       end do
                       do j = i + 1,n
                          temp = abs(a(j,i))
                          symb_zero = symb_zero .and. (x(j) == zero .or. temp == zero)
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
     end subroutine la_qla_syamv

     !> SLA_SYRCOND: estimates the Skeel condition number of  op(A) * op2(C)
     !> where op2 is determined by CMODE as follows
     !> CMODE =  1    op2(C) = C
     !> CMODE =  0    op2(C) = I
     !> CMODE = -1    op2(C) = inv(C)
     !> The Skeel condition number cond(A) = norminf( |inv(A)||A| )
     !> is computed by computing scaling factors R such that
     !> diag(R)*A*op2(C) is row equilibrated and computing the standard
     !> infinity-norm condition number.

     real(sp) function la_sla_syrcond(uplo,n,a,lda,af,ldaf,ipiv,cmode,c,info,work, &
               iwork)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,lda,ldaf,cmode
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           integer(ilp),intent(in) :: ipiv(*)
           real(sp),intent(in) :: a(lda,*),af(ldaf,*),c(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           character :: normin
           integer(ilp) :: kase,i,j
           real(sp) :: ainvnm,smlnum,tmp
           logical(lk) :: up
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           la_sla_syrcond = zero
           info = 0
           if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SLA_SYRCOND',-info)
              return
           end if
           if (n == 0) then
              la_sla_syrcond = one
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute the equilibration matrix r such that
           ! inv(r)*a*c has unit 1-norm.
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (cmode == 1) then
                    do j = 1,i
                       tmp = tmp + abs(a(j,i)*c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j)*c(j))
                    end do
                 else if (cmode == 0) then
                    do j = 1,i
                       tmp = tmp + abs(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j))
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + abs(a(j,i)/c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j)/c(j))
                    end do
                 end if
                 work(2*n + i) = tmp
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (cmode == 1) then
                    do j = 1,i
                       tmp = tmp + abs(a(i,j)*c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i)*c(j))
                    end do
                 else if (cmode == 0) then
                    do j = 1,i
                       tmp = tmp + abs(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i))
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + abs(a(i,j)/c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i)/c(j))
                    end do
                 end if
                 work(2*n + i) = tmp
              end do
           end if
           ! estimate the norm of inv(op(a)).
           smlnum = la_slamch('SAFE MINIMUM')
           ainvnm = zero
           normin = 'N'
           kase = 0
           10 continue
           call la_slacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*work(2*n + i)
                 end do
                 if (up) then
                    call la_ssytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_ssytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (cmode == 1) then
                    do i = 1,n
                       work(i) = work(i)/c(i)
                    end do
                 else if (cmode == -1) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**t).
                 if (cmode == 1) then
                    do i = 1,n
                       work(i) = work(i)/c(i)
                    end do
                 else if (cmode == -1) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_ssytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_ssytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*work(2*n + i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= 0.0_sp) la_sla_syrcond = (1.0_sp/ainvnm)
           return
     end function la_sla_syrcond
     !> DLA_SYRCOND: estimates the Skeel condition number of  op(A) * op2(C)
     !> where op2 is determined by CMODE as follows
     !> CMODE =  1    op2(C) = C
     !> CMODE =  0    op2(C) = I
     !> CMODE = -1    op2(C) = inv(C)
     !> The Skeel condition number cond(A) = norminf( |inv(A)||A| )
     !> is computed by computing scaling factors R such that
     !> diag(R)*A*op2(C) is row equilibrated and computing the standard
     !> infinity-norm condition number.

     real(dp) function la_dla_syrcond(uplo,n,a,lda,af,ldaf,ipiv,cmode,c,info,work, &
               iwork)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,lda,ldaf,cmode
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           integer(ilp),intent(in) :: ipiv(*)
           real(dp),intent(in) :: a(lda,*),af(ldaf,*),c(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           character :: normin
           integer(ilp) :: kase,i,j
           real(dp) :: ainvnm,smlnum,tmp
           logical(lk) :: up
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           la_dla_syrcond = zero
           info = 0
           if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DLA_SYRCOND',-info)
              return
           end if
           if (n == 0) then
              la_dla_syrcond = one
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute the equilibration matrix r such that
           ! inv(r)*a*c has unit 1-norm.
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (cmode == 1) then
                    do j = 1,i
                       tmp = tmp + abs(a(j,i)*c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j)*c(j))
                    end do
                 else if (cmode == 0) then
                    do j = 1,i
                       tmp = tmp + abs(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j))
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + abs(a(j,i)/c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j)/c(j))
                    end do
                 end if
                 work(2*n + i) = tmp
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (cmode == 1) then
                    do j = 1,i
                       tmp = tmp + abs(a(i,j)*c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i)*c(j))
                    end do
                 else if (cmode == 0) then
                    do j = 1,i
                       tmp = tmp + abs(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i))
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + abs(a(i,j)/c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i)/c(j))
                    end do
                 end if
                 work(2*n + i) = tmp
              end do
           end if
           ! estimate the norm of inv(op(a)).
           smlnum = la_dlamch('SAFE MINIMUM')
           ainvnm = zero
           normin = 'N'
           kase = 0
           10 continue
           call la_dlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*work(2*n + i)
                 end do
                 if (up) then
                    call la_dsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_dsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (cmode == 1) then
                    do i = 1,n
                       work(i) = work(i)/c(i)
                    end do
                 else if (cmode == -1) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**t).
                 if (cmode == 1) then
                    do i = 1,n
                       work(i) = work(i)/c(i)
                    end do
                 else if (cmode == -1) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_dsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_dsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*work(2*n + i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_dla_syrcond = (one/ainvnm)
           return
     end function la_dla_syrcond
     !> QLA_SYRCOND: estimates the Skeel condition number of  op(A) * op2(C)
     !> where op2 is determined by CMODE as follows
     !> CMODE =  1    op2(C) = C
     !> CMODE =  0    op2(C) = I
     !> CMODE = -1    op2(C) = inv(C)
     !> The Skeel condition number cond(A) = norminf( |inv(A)||A| )
     !> is computed by computing scaling factors R such that
     !> diag(R)*A*op2(C) is row equilibrated and computing the standard
     !> infinity-norm condition number.

     real(qp) function la_qla_syrcond(uplo,n,a,lda,af,ldaf,ipiv,cmode,c,info,work, &
               iwork)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,lda,ldaf,cmode
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           integer(ilp),intent(in) :: ipiv(*)
           real(qp),intent(in) :: a(lda,*),af(ldaf,*),c(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           character :: normin
           integer(ilp) :: kase,i,j
           real(qp) :: ainvnm,smlnum,tmp
           logical(lk) :: up
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           la_qla_syrcond = zero
           info = 0
           if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QLA_SYRCOND',-info)
              return
           end if
           if (n == 0) then
              la_qla_syrcond = one
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute the equilibration matrix r such that
           ! inv(r)*a*c has unit 1-norm.
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (cmode == 1) then
                    do j = 1,i
                       tmp = tmp + abs(a(j,i)*c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j)*c(j))
                    end do
                 else if (cmode == 0) then
                    do j = 1,i
                       tmp = tmp + abs(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j))
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + abs(a(j,i)/c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(i,j)/c(j))
                    end do
                 end if
                 work(2*n + i) = tmp
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (cmode == 1) then
                    do j = 1,i
                       tmp = tmp + abs(a(i,j)*c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i)*c(j))
                    end do
                 else if (cmode == 0) then
                    do j = 1,i
                       tmp = tmp + abs(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i))
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + abs(a(i,j)/c(j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + abs(a(j,i)/c(j))
                    end do
                 end if
                 work(2*n + i) = tmp
              end do
           end if
           ! estimate the norm of inv(op(a)).
           smlnum = la_qlamch('SAFE MINIMUM')
           ainvnm = zero
           normin = 'N'
           kase = 0
           10 continue
           call la_qlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*work(2*n + i)
                 end do
                 if (up) then
                    call la_qsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_qsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (cmode == 1) then
                    do i = 1,n
                       work(i) = work(i)/c(i)
                    end do
                 else if (cmode == -1) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**t).
                 if (cmode == 1) then
                    do i = 1,n
                       work(i) = work(i)/c(i)
                    end do
                 else if (cmode == -1) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_qsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_qsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*work(2*n + i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_qla_syrcond = (one/ainvnm)
           return
     end function la_qla_syrcond

     !> SLA_SYRPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(sp) function la_sla_syrpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(sp),intent(in) :: a(lda,*),af(ldaf,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(sp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(abs(a(i,j)),work(n + i))
                    work(n + j) = max(abs(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(abs(a(i,j)),work(n + i))
                    work(n + j) = max(abs(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_ssytrs.
           ! calls to la_sswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(abs(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(abs(af(i,k)),work(k))
                       work(k - 1) = max(abs(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(abs(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(abs(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(abs(af(i,k)),work(k))
                       work(k + 1) = max(abs(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(abs(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= 0.0_sp) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= 0.0_sp) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_sla_syrpvgrw = rpvgrw
     end function la_sla_syrpvgrw
     !> DLA_SYRPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(dp) function la_dla_syrpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(dp),intent(in) :: a(lda,*),af(ldaf,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(dp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(abs(a(i,j)),work(n + i))
                    work(n + j) = max(abs(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(abs(a(i,j)),work(n + i))
                    work(n + j) = max(abs(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_dsytrs.
           ! calls to la_sswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(abs(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(abs(af(i,k)),work(k))
                       work(k - 1) = max(abs(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(abs(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(abs(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(abs(af(i,k)),work(k))
                       work(k + 1) = max(abs(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(abs(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_dla_syrpvgrw = rpvgrw
     end function la_dla_syrpvgrw
     !> QLA_SYRPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(qp) function la_qla_syrpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(qp),intent(in) :: a(lda,*),af(ldaf,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(qp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(abs(a(i,j)),work(n + i))
                    work(n + j) = max(abs(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(abs(a(i,j)),work(n + i))
                    work(n + j) = max(abs(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_qsytrs.
           ! calls to la_dswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(abs(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(abs(af(i,k)),work(k))
                       work(k - 1) = max(abs(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(abs(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(abs(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(abs(af(i,k)),work(k))
                       work(k + 1) = max(abs(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(abs(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_qla_syrpvgrw = rpvgrw
     end function la_qla_syrpvgrw

     !> CLA_GERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     pure real(sp) function la_cla_gerpvgrw(n,ncols,a,lda,af,ldaf)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,ncols,lda,ldaf
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: amax,umax,rpvgrw
           complex(sp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,min,abs,real,aimag
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           rpvgrw = one
           do j = 1,ncols
              amax = zero
              umax = zero
              do i = 1,n
                 amax = max(cabs1(a(i,j)),amax)
              end do
              do i = 1,j
                 umax = max(cabs1(af(i,j)),umax)
              end do
              if (umax /= 0.0_sp) then
                 rpvgrw = min(amax/umax,rpvgrw)
              end if
           end do
           la_cla_gerpvgrw = rpvgrw
     end function la_cla_gerpvgrw
     !> ZLA_GERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     pure real(dp) function la_zla_gerpvgrw(n,ncols,a,lda,af,ldaf)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,ncols,lda,ldaf
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: amax,umax,rpvgrw
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,min,abs,real,aimag
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           rpvgrw = one
           do j = 1,ncols
              amax = zero
              umax = zero
              do i = 1,n
                 amax = max(cabs1(a(i,j)),amax)
              end do
              do i = 1,j
                 umax = max(cabs1(af(i,j)),umax)
              end do
              if (umax /= zero) then
                 rpvgrw = min(amax/umax,rpvgrw)
              end if
           end do
           la_zla_gerpvgrw = rpvgrw
     end function la_zla_gerpvgrw
     !> WLA_GERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     pure real(qp) function la_wla_gerpvgrw(n,ncols,a,lda,af,ldaf)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,ncols,lda,ldaf
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: amax,umax,rpvgrw
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: max,min,abs,real,aimag
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           rpvgrw = one
           do j = 1,ncols
              amax = zero
              umax = zero
              do i = 1,n
                 amax = max(cabs1(a(i,j)),amax)
              end do
              do i = 1,j
                 umax = max(cabs1(af(i,j)),umax)
              end do
              if (umax /= zero) then
                 rpvgrw = min(amax/umax,rpvgrw)
              end if
           end do
           la_wla_gerpvgrw = rpvgrw
     end function la_wla_gerpvgrw

     !> CLA_SYAMV:  performs the matrix-vector operation
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

     subroutine la_cla_syamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n
           integer(ilp),intent(in) :: uplo
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
              call la_xerbla('CLA_SYAMV',info)
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
     end subroutine la_cla_syamv
     !> ZLA_SYAMV:  performs the matrix-vector operation
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

     subroutine la_zla_syamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n
           integer(ilp),intent(in) :: uplo
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
              call la_xerbla('ZLA_SYAMV',info)
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
     end subroutine la_zla_syamv
     !> WLA_SYAMV:  performs the matrix-vector operation
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

     subroutine la_wla_syamv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: incx,incy,lda,n
           integer(ilp),intent(in) :: uplo
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
              call la_xerbla('WLA_SYAMV',info)
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
     end subroutine la_wla_syamv

     !> CLA_GBRCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a REAL vector.

     real(sp) function la_cla_gbrcond_c(trans,n,kl,ku,ab,ldab,afb,ldafb,ipiv,c, &
               capply,info,work,rwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,kl,ku,ldab,ldafb
           integer(ilp) :: kd,ke
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: ab(ldab,*),afb(ldafb,*)
           complex(sp),intent(out) :: work(*)
           real(sp),intent(in) :: c(*)
           real(sp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: notrans
           integer(ilp) :: kase,i,j
           real(sp) :: ainvnm,anorm,tmp
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           la_cla_gbrcond_c = zero
           info = 0
           notrans = la_lsame(trans,'N')
           if (.not. notrans .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kl < 0 .or. kl > n - 1) then
              info = -3
           else if (ku < 0 .or. ku > n - 1) then
              info = -4
           else if (ldab < kl + ku + 1) then
              info = -6
           else if (ldafb < 2*kl + ku + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CLA_GBRCOND_C',-info)
              return
           end if
           ! compute norm of op(a)*op2(c).
           anorm = zero
           kd = ku + 1
           ke = kl + 1
           if (notrans) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(kd + i - j,j))/c(j)
                    end do
                 else
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(kd + i - j,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(ke - i + j,i))/c(j)
                    end do
                 else
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(ke - i + j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_cla_gbrcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (notrans) then
                    call la_cgbtrs('NO TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv,work,n, &
                              info)
                 else
                    call la_cgbtrs('CONJUGATE TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv, &
                              work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (notrans) then
                    call la_cgbtrs('CONJUGATE TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv, &
                              work,n,info)
                 else
                    call la_cgbtrs('NO TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv,work,n, &
                              info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_cla_gbrcond_c = one/ainvnm
           return
     end function la_cla_gbrcond_c
     !> ZLA_GBRCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector.

     real(dp) function la_zla_gbrcond_c(trans,n,kl,ku,ab,ldab,afb,ldafb,ipiv,c, &
               capply,info,work,rwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,kl,ku,ldab,ldafb
           integer(ilp) :: kd,ke
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: ab(ldab,*),afb(ldafb,*)
           complex(dp),intent(out) :: work(*)
           real(dp),intent(in) :: c(*)
           real(dp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: notrans
           integer(ilp) :: kase,i,j
           real(dp) :: ainvnm,anorm,tmp
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           la_zla_gbrcond_c = zero
           info = 0
           notrans = la_lsame(trans,'N')
           if (.not. notrans .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kl < 0 .or. kl > n - 1) then
              info = -3
           else if (ku < 0 .or. ku > n - 1) then
              info = -4
           else if (ldab < kl + ku + 1) then
              info = -6
           else if (ldafb < 2*kl + ku + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_GBRCOND_C',-info)
              return
           end if
           ! compute norm of op(a)*op2(c).
           anorm = zero
           kd = ku + 1
           ke = kl + 1
           if (notrans) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(kd + i - j,j))/c(j)
                    end do
                 else
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(kd + i - j,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(ke - i + j,i))/c(j)
                    end do
                 else
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(ke - i + j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_zla_gbrcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (notrans) then
                    call la_zgbtrs('NO TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv,work,n, &
                              info)
                 else
                    call la_zgbtrs('CONJUGATE TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv, &
                              work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (notrans) then
                    call la_zgbtrs('CONJUGATE TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv, &
                              work,n,info)
                 else
                    call la_zgbtrs('NO TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv,work,n, &
                              info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_zla_gbrcond_c = one/ainvnm
           return
     end function la_zla_gbrcond_c
     !> WLA_GBRCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a QUAD PRECISION vector.

     real(qp) function la_wla_gbrcond_c(trans,n,kl,ku,ab,ldab,afb,ldafb,ipiv,c, &
               capply,info,work,rwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,kl,ku,ldab,ldafb
           integer(ilp) :: kd,ke
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: ab(ldab,*),afb(ldafb,*)
           complex(qp),intent(out) :: work(*)
           real(qp),intent(in) :: c(*)
           real(qp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: notrans
           integer(ilp) :: kase,i,j
           real(qp) :: ainvnm,anorm,tmp
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           la_wla_gbrcond_c = zero
           info = 0
           notrans = la_lsame(trans,'N')
           if (.not. notrans .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kl < 0 .or. kl > n - 1) then
              info = -3
           else if (ku < 0 .or. ku > n - 1) then
              info = -4
           else if (ldab < kl + ku + 1) then
              info = -6
           else if (ldafb < 2*kl + ku + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WLA_GBRCOND_C',-info)
              return
           end if
           ! compute norm of op(a)*op2(c).
           anorm = zero
           kd = ku + 1
           ke = kl + 1
           if (notrans) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(kd + i - j,j))/c(j)
                    end do
                 else
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(kd + i - j,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(ke - i + j,i))/c(j)
                    end do
                 else
                    do j = max(i - kl,1),min(i + ku,n)
                       tmp = tmp + cabs1(ab(ke - i + j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_wla_gbrcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (notrans) then
                    call la_wgbtrs('NO TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv,work,n, &
                              info)
                 else
                    call la_wgbtrs('CONJUGATE TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv, &
                              work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (notrans) then
                    call la_wgbtrs('CONJUGATE TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv, &
                              work,n,info)
                 else
                    call la_wgbtrs('NO TRANSPOSE',n,kl,ku,1,afb,ldafb,ipiv,work,n, &
                              info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_wla_gbrcond_c = one/ainvnm
           return
     end function la_wla_gbrcond_c

     !> CLA_GERCOND_C: computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a REAL vector.

     real(sp) function la_cla_gercond_c(trans,n,a,lda,af,ldaf,ipiv,c,capply,info, &
               work,rwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(sp),intent(out) :: work(*)
           real(sp),intent(in) :: c(*)
           real(sp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: notrans
           integer(ilp) :: kase,i,j
           real(sp) :: ainvnm,anorm,tmp
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,aimag
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           la_cla_gercond_c = zero
           info = 0
           notrans = la_lsame(trans,'N')
           if (.not. notrans .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CLA_GERCOND_C',-info)
              return
           end if
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (notrans) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_cla_gercond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (notrans) then
                    call la_cgetrs('NO TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info)

                 else
                    call la_cgetrs('CONJUGATE TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info &
                              )
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (notrans) then
                    call la_cgetrs('CONJUGATE TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info &
                              )
                 else
                    call la_cgetrs('NO TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info)

                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_cla_gercond_c = one/ainvnm
           return
     end function la_cla_gercond_c
     !> ZLA_GERCOND_C: computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector.

     real(dp) function la_zla_gercond_c(trans,n,a,lda,af,ldaf,ipiv,c,capply,info, &
               work,rwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(dp),intent(out) :: work(*)
           real(dp),intent(in) :: c(*)
           real(dp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: notrans
           integer(ilp) :: kase,i,j
           real(dp) :: ainvnm,anorm,tmp
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,aimag
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           la_zla_gercond_c = zero
           info = 0
           notrans = la_lsame(trans,'N')
           if (.not. notrans .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_GERCOND_C',-info)
              return
           end if
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (notrans) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_zla_gercond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (notrans) then
                    call la_zgetrs('NO TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info)

                 else
                    call la_zgetrs('CONJUGATE TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info &
                              )
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (notrans) then
                    call la_zgetrs('CONJUGATE TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info &
                              )
                 else
                    call la_zgetrs('NO TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info)

                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_zla_gercond_c = one/ainvnm
           return
     end function la_zla_gercond_c
     !> WLA_GERCOND_C: computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a QUAD PRECISION vector.

     real(qp) function la_wla_gercond_c(trans,n,a,lda,af,ldaf,ipiv,c,capply,info, &
               work,rwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(qp),intent(out) :: work(*)
           real(qp),intent(in) :: c(*)
           real(qp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: notrans
           integer(ilp) :: kase,i,j
           real(qp) :: ainvnm,anorm,tmp
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,aimag
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           la_wla_gercond_c = zero
           info = 0
           notrans = la_lsame(trans,'N')
           if (.not. notrans .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WLA_GERCOND_C',-info)
              return
           end if
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (notrans) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_wla_gercond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (notrans) then
                    call la_wgetrs('NO TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info)

                 else
                    call la_wgetrs('CONJUGATE TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info &
                              )
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (notrans) then
                    call la_wgetrs('CONJUGATE TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info &
                              )
                 else
                    call la_wgetrs('NO TRANSPOSE',n,1,af,ldaf,ipiv,work,n,info)

                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_wla_gercond_c = one/ainvnm
           return
     end function la_wla_gercond_c

     !> CLA_HERCOND_C: computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a REAL vector.

     real(sp) function la_cla_hercond_c(uplo,n,a,lda,af,ldaf,ipiv,c,capply,info, &
               work,rwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(sp),intent(out) :: work(*)
           real(sp),intent(in) :: c(*)
           real(sp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase,i,j
           real(sp) :: ainvnm,anorm,tmp
           logical(lk) :: up,upper
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           la_cla_hercond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CLA_HERCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_cla_hercond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_chetrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_chetrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_chetrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_chetrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_cla_hercond_c = one/ainvnm
           return
     end function la_cla_hercond_c
     !> ZLA_HERCOND_C: computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector.

     real(dp) function la_zla_hercond_c(uplo,n,a,lda,af,ldaf,ipiv,c,capply,info,work, &
                rwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(dp),intent(out) :: work(*)
           real(dp),intent(in) :: c(*)
           real(dp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase,i,j
           real(dp) :: ainvnm,anorm,tmp
           logical(lk) :: up,upper
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           la_zla_hercond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_HERCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_zla_hercond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_zhetrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_zhetrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_zhetrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_zhetrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_zla_hercond_c = one/ainvnm
           return
     end function la_zla_hercond_c
     !> WLA_HERCOND_C: computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a QUAD PRECISION vector.

     real(qp) function la_wla_hercond_c(uplo,n,a,lda,af,ldaf,ipiv,c,capply,info,work, &
                rwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(qp),intent(out) :: work(*)
           real(qp),intent(in) :: c(*)
           real(qp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase,i,j
           real(qp) :: ainvnm,anorm,tmp
           logical(lk) :: up,upper
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           la_wla_hercond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WLA_HERCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_wla_hercond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_whetrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_whetrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_whetrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_whetrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_wla_hercond_c = one/ainvnm
           return
     end function la_wla_hercond_c

     !> CLA_PORCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a REAL vector

     real(sp) function la_cla_porcond_c(uplo,n,a,lda,af,ldaf,c,capply,info,work, &
               rwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(sp),intent(out) :: work(*)
           real(sp),intent(in) :: c(*)
           real(sp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase
           real(sp) :: ainvnm,anorm,tmp
           integer(ilp) :: i,j
           logical(lk) :: up,upper
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,aimag
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           la_cla_porcond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CLA_PORCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_cla_porcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_cpotrs('U',n,1,af,ldaf,work,n,info)
                 else
                    call la_cpotrs('L',n,1,af,ldaf,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_cpotrs('U',n,1,af,ldaf,work,n,info)
                 else
                    call la_cpotrs('L',n,1,af,ldaf,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_cla_porcond_c = one/ainvnm
           return
     end function la_cla_porcond_c
     !> ZLA_PORCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector

     real(dp) function la_zla_porcond_c(uplo,n,a,lda,af,ldaf,c,capply,info,work, &
               rwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(dp),intent(out) :: work(*)
           real(dp),intent(in) :: c(*)
           real(dp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase
           real(dp) :: ainvnm,anorm,tmp
           integer(ilp) :: i,j
           logical(lk) :: up,upper
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,aimag
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           la_zla_porcond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_PORCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_zla_porcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_zpotrs('U',n,1,af,ldaf,work,n,info)
                 else
                    call la_zpotrs('L',n,1,af,ldaf,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_zpotrs('U',n,1,af,ldaf,work,n,info)
                 else
                    call la_zpotrs('L',n,1,af,ldaf,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_zla_porcond_c = one/ainvnm
           return
     end function la_zla_porcond_c
     !> WLA_PORCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a QUAD PRECISION vector

     real(qp) function la_wla_porcond_c(uplo,n,a,lda,af,ldaf,c,capply,info,work, &
               rwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(qp),intent(out) :: work(*)
           real(qp),intent(in) :: c(*)
           real(qp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase
           real(qp) :: ainvnm,anorm,tmp
           integer(ilp) :: i,j
           logical(lk) :: up,upper
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,aimag
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           la_wla_porcond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WLA_PORCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_wla_porcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_wpotrs('U',n,1,af,ldaf,work,n,info)
                 else
                    call la_wpotrs('L',n,1,af,ldaf,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**h).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_wpotrs('U',n,1,af,ldaf,work,n,info)
                 else
                    call la_wpotrs('L',n,1,af,ldaf,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_wla_porcond_c = one/ainvnm
           return
     end function la_wla_porcond_c

     !> CLA_SYRCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a REAL vector.

     real(sp) function la_cla_syrcond_c(uplo,n,a,lda,af,ldaf,ipiv,c,capply,info, &
               work,rwork)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(sp),intent(out) :: work(*)
           real(sp),intent(in) :: c(*)
           real(sp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase
           real(sp) :: ainvnm,anorm,tmp
           integer(ilp) :: i,j
           logical(lk) :: up,upper
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           la_cla_syrcond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CLA_SYRCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_cla_syrcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_csytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_csytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**t).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_csytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_csytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_cla_syrcond_c = one/ainvnm
           return
     end function la_cla_syrcond_c
     !> ZLA_SYRCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector.

     real(dp) function la_zla_syrcond_c(uplo,n,a,lda,af,ldaf,ipiv,c,capply,info,work, &
                rwork)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(dp),intent(out) :: work(*)
           real(dp),intent(in) :: c(*)
           real(dp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase
           real(dp) :: ainvnm,anorm,tmp
           integer(ilp) :: i,j
           logical(lk) :: up,upper
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           la_zla_syrcond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZLA_SYRCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_zla_syrcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_zsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_zsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**t).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_zsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_zsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_zla_syrcond_c = one/ainvnm
           return
     end function la_zla_syrcond_c
     !> WLA_SYRCOND_C: Computes the infinity norm condition number of
     !> op(A) * inv(diag(C)) where C is a QUAD PRECISION vector.

     real(qp) function la_wla_syrcond_c(uplo,n,a,lda,af,ldaf,ipiv,c,capply,info,work, &
                rwork)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           logical(lk),intent(in) :: capply
           integer(ilp),intent(in) :: n,lda,ldaf
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
           complex(qp),intent(out) :: work(*)
           real(qp),intent(in) :: c(*)
           real(qp),intent(out) :: rwork(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: kase
           real(qp) :: ainvnm,anorm,tmp
           integer(ilp) :: i,j
           logical(lk) :: up,upper
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           la_wla_syrcond_c = zero
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (ldaf < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WLA_SYRCOND_C',-info)
              return
           end if
           up = .false.
           if (la_lsame(uplo,'U')) up = .true.
           ! compute norm of op(a)*op2(c).
           anorm = zero
           if (up) then
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(j,i))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(i,j))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           else
              do i = 1,n
                 tmp = zero
                 if (capply) then
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))/c(j)
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))/c(j)
                    end do
                 else
                    do j = 1,i
                       tmp = tmp + cabs1(a(i,j))
                    end do
                    do j = i + 1,n
                       tmp = tmp + cabs1(a(j,i))
                    end do
                 end if
                 rwork(i) = tmp
                 anorm = max(anorm,tmp)
              end do
           end if
           ! quick return if possible.
           if (n == 0) then
              la_wla_syrcond_c = one
              return
           else if (anorm == zero) then
              return
           end if
           ! estimate the norm of inv(op(a)).
           ainvnm = zero
           kase = 0
           10 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              if (kase == 2) then
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
                 if (up) then
                    call la_wsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_wsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by inv(c).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
              else
                 ! multiply by inv(c**t).
                 if (capply) then
                    do i = 1,n
                       work(i) = work(i)*c(i)
                    end do
                 end if
                 if (up) then
                    call la_wsytrs('U',n,1,af,ldaf,ipiv,work,n,info)
                 else
                    call la_wsytrs('L',n,1,af,ldaf,ipiv,work,n,info)
                 end if
                 ! multiply by r.
                 do i = 1,n
                    work(i) = work(i)*rwork(i)
                 end do
              end if
              go to 10
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) la_wla_syrcond_c = one/ainvnm
           return
     end function la_wla_syrcond_c

     !> CLA_SYRPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(sp) function la_cla_syrpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
           real(sp),intent(out) :: work(*)
           integer(ilp),intent(in) :: ipiv(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(sp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper
           complex(sp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_csytrs.
           ! calls to la_sswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= 0.0_sp) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= 0.0_sp) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_cla_syrpvgrw = rpvgrw
     end function la_cla_syrpvgrw
     !> ZLA_SYRPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(dp) function la_zla_syrpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
           real(dp),intent(out) :: work(*)
           integer(ilp),intent(in) :: ipiv(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(dp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_zsytrs.
           ! calls to la_sswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_zla_syrpvgrw = rpvgrw
     end function la_zla_syrpvgrw
     !> WLA_SYRPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(qp) function la_wla_syrpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
           real(qp),intent(out) :: work(*)
           integer(ilp),intent(in) :: ipiv(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(qp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_wsytrs.
           ! calls to la_dswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_wla_syrpvgrw = rpvgrw
     end function la_wla_syrpvgrw

end module la_lapack_others_sm
