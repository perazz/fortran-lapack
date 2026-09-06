!> Linear solve helpers: condition estimation, componentwise backward error
module la_lapack_solve_aux
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_lapack_aux
     use la_lapack_auxiliary
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slacn2
     public :: la_slacon
     public :: la_sla_lin_berr
     public :: la_dlacn2
     public :: la_dlacon
     public :: la_dla_lin_berr
     public :: la_qlacn2
     public :: la_qlacon
     public :: la_qla_lin_berr
     public :: la_cla_lin_berr
     public :: la_clacn2
     public :: la_clacon
     public :: la_zla_lin_berr
     public :: la_zlacn2
     public :: la_zlacon
     public :: la_wla_lin_berr
     public :: la_wlacn2
     public :: la_wlacon

     contains

     !> SLACN2: estimates the 1-norm of a square, real matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     pure subroutine la_slacn2(n,v,x,isgn,est,kase,isave)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(sp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(out) :: isgn(*)
           integer(ilp),intent(inout) :: isave(3)
           real(sp),intent(out) :: v(*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,jlast
           real(sp) :: altsgn,estold,temp,xs
           ! Intrinsic Functions
           intrinsic :: abs,nint,real
           ! Executable Statements
           if (kase == 0) then
              do i = 1,n
                 x(i) = one/real(n,KIND=sp)
              end do
              kase = 1
              isave(1) = 1
              return
           end if
           go to(20,40,70,110,140) isave(1)
           ! ................ entry   (isave( 1 ) = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 150
           end if
           est = la_sasum(n,x,1)
           do i = 1,n
              if (x(i) >= zero) then
                 x(i) = one
              else
                 x(i) = -one
              end if
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           isave(1) = 2
           return
           ! ................ entry   (isave( 1 ) = 2)
           ! first iteration.  x has been overwritten by transpose(a)*x.
           40 continue
           isave(2) = la_isamax(n,x,1)
           isave(3) = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = zero
           end do
           x(isave(2)) = one
           kase = 1
           isave(1) = 3
           return
           ! ................ entry   (isave( 1 ) = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_scopy(n,x,1,v,1)
           estold = est
           est = la_sasum(n,v,1)
           do i = 1,n
              if (x(i) >= zero) then
                 xs = one
              else
                 xs = -one
              end if
              if (nint(xs,KIND=ilp) /= isgn(i)) go to 90
           end do
           ! repeated sign vector detected, hence algorithm has converged.
           go to 120
           90 continue
           ! test for cycling.
           if (est <= estold) go to 120
           do i = 1,n
              if (x(i) >= zero) then
                 x(i) = one
              else
                 x(i) = -one
              end if
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           isave(1) = 4
           return
           ! ................ entry   (isave( 1 ) = 4)
           ! x has been overwritten by transpose(a)*x.
           110 continue
           jlast = isave(2)
           isave(2) = la_isamax(n,x,1)
           if ((x(jlast) /= abs(x(isave(2)))) .and. (isave(3) < itmax)) then
              isave(3) = isave(3) + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           120 continue
           altsgn = one
           do i = 1,n
              x(i) = altsgn*(one + real(i - 1,KIND=sp)/real(n - 1,KIND=sp))
              altsgn = -altsgn
           end do
           kase = 1
           isave(1) = 5
           return
           ! ................ entry   (isave( 1 ) = 5)
           ! x has been overwritten by a*x.
           140 continue
           temp = two*(la_sasum(n,x,1)/real(3*n,KIND=sp))
           if (temp > est) then
              call la_scopy(n,x,1,v,1)
              est = temp
           end if
           150 continue
           kase = 0
           return
     end subroutine la_slacn2
     !> DLACN2: estimates the 1-norm of a square, real matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     pure subroutine la_dlacn2(n,v,x,isgn,est,kase,isave)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(dp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(out) :: isgn(*)
           integer(ilp),intent(inout) :: isave(3)
           real(dp),intent(out) :: v(*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,jlast
           real(dp) :: altsgn,estold,temp,xs
           ! Intrinsic Functions
           intrinsic :: abs,real,nint
           ! Executable Statements
           if (kase == 0) then
              do i = 1,n
                 x(i) = one/real(n,KIND=dp)
              end do
              kase = 1
              isave(1) = 1
              return
           end if
           go to(20,40,70,110,140) isave(1)
           ! ................ entry   (isave( 1 ) = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 150
           end if
           est = la_dasum(n,x,1)
           do i = 1,n
              if (x(i) >= zero) then
                 x(i) = one
              else
                 x(i) = -one
              end if
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           isave(1) = 2
           return
           ! ................ entry   (isave( 1 ) = 2)
           ! first iteration.  x has been overwritten by transpose(a)*x.
           40 continue
           isave(2) = la_idamax(n,x,1)
           isave(3) = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = zero
           end do
           x(isave(2)) = one
           kase = 1
           isave(1) = 3
           return
           ! ................ entry   (isave( 1 ) = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_dcopy(n,x,1,v,1)
           estold = est
           est = la_dasum(n,v,1)
           do i = 1,n
              if (x(i) >= zero) then
                 xs = one
              else
                 xs = -one
              end if
              if (nint(xs,KIND=ilp) /= isgn(i)) go to 90
           end do
           ! repeated sign vector detected, hence algorithm has converged.
           go to 120
           90 continue
           ! test for cycling.
           if (est <= estold) go to 120
           do i = 1,n
              if (x(i) >= zero) then
                 x(i) = one
              else
                 x(i) = -one
              end if
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           isave(1) = 4
           return
           ! ................ entry   (isave( 1 ) = 4)
           ! x has been overwritten by transpose(a)*x.
           110 continue
           jlast = isave(2)
           isave(2) = la_idamax(n,x,1)
           if ((x(jlast) /= abs(x(isave(2)))) .and. (isave(3) < itmax)) then
              isave(3) = isave(3) + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           120 continue
           altsgn = one
           do i = 1,n
              x(i) = altsgn*(one + real(i - 1,KIND=dp)/real(n - 1,KIND=dp))
              altsgn = -altsgn
           end do
           kase = 1
           isave(1) = 5
           return
           ! ................ entry   (isave( 1 ) = 5)
           ! x has been overwritten by a*x.
           140 continue
           temp = two*(la_dasum(n,x,1)/real(3*n,KIND=dp))
           if (temp > est) then
              call la_dcopy(n,x,1,v,1)
              est = temp
           end if
           150 continue
           kase = 0
           return
     end subroutine la_dlacn2
     !> QLACN2: estimates the 1-norm of a square, real matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     pure subroutine la_qlacn2(n,v,x,isgn,est,kase,isave)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(qp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(out) :: isgn(*)
           integer(ilp),intent(inout) :: isave(3)
           real(qp),intent(out) :: v(*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,jlast
           real(qp) :: altsgn,estold,temp,xs
           ! Intrinsic Functions
           intrinsic :: abs,real,nint
           ! Executable Statements
           if (kase == 0) then
              do i = 1,n
                 x(i) = one/real(n,KIND=qp)
              end do
              kase = 1
              isave(1) = 1
              return
           end if
           go to(20,40,70,110,140) isave(1)
           ! ................ entry   (isave( 1 ) = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 150
           end if
           est = la_qasum(n,x,1)
           do i = 1,n
              if (x(i) >= zero) then
                 x(i) = one
              else
                 x(i) = -one
              end if
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           isave(1) = 2
           return
           ! ................ entry   (isave( 1 ) = 2)
           ! first iteration.  x has been overwritten by transpose(a)*x.
           40 continue
           isave(2) = la_iqamax(n,x,1)
           isave(3) = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = zero
           end do
           x(isave(2)) = one
           kase = 1
           isave(1) = 3
           return
           ! ................ entry   (isave( 1 ) = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_qcopy(n,x,1,v,1)
           estold = est
           est = la_qasum(n,v,1)
           do i = 1,n
              if (x(i) >= zero) then
                 xs = one
              else
                 xs = -one
              end if
              if (nint(xs,KIND=ilp) /= isgn(i)) go to 90
           end do
           ! repeated sign vector detected, hence algorithm has converged.
           go to 120
           90 continue
           ! test for cycling.
           if (est <= estold) go to 120
           do i = 1,n
              if (x(i) >= zero) then
                 x(i) = one
              else
                 x(i) = -one
              end if
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           isave(1) = 4
           return
           ! ................ entry   (isave( 1 ) = 4)
           ! x has been overwritten by transpose(a)*x.
           110 continue
           jlast = isave(2)
           isave(2) = la_iqamax(n,x,1)
           if ((x(jlast) /= abs(x(isave(2)))) .and. (isave(3) < itmax)) then
              isave(3) = isave(3) + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           120 continue
           altsgn = one
           do i = 1,n
              x(i) = altsgn*(one + real(i - 1,KIND=qp)/real(n - 1,KIND=qp))
              altsgn = -altsgn
           end do
           kase = 1
           isave(1) = 5
           return
           ! ................ entry   (isave( 1 ) = 5)
           ! x has been overwritten by a*x.
           140 continue
           temp = two*(la_qasum(n,x,1)/real(3*n,KIND=qp))
           if (temp > est) then
              call la_qcopy(n,x,1,v,1)
              est = temp
           end if
           150 continue
           kase = 0
           return
     end subroutine la_qlacn2

     !> SLACON: estimates the 1-norm of a square, real matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     subroutine la_slacon(n,v,x,isgn,est,kase)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(sp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(out) :: isgn(*)
           real(sp),intent(out) :: v(*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,iter,j,jlast,jump
           real(sp) :: altsgn,estold,temp
           ! Intrinsic Functions
           intrinsic :: abs,nint,real,sign
           ! Save Statement
           save
           ! Executable Statements
           if (kase == 0) then
              do i = 1,n
                 x(i) = one/real(n,KIND=sp)
              end do
              kase = 1
              jump = 1
              return
           end if
           go to(20,40,70,110,140) jump
           ! ................ entry   (jump = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 150
           end if
           est = la_sasum(n,x,1)
           do i = 1,n
              x(i) = sign(one,x(i))
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           jump = 2
           return
           ! ................ entry   (jump = 2)
           ! first iteration.  x has been overwritten by transpose(a)*x.
           40 continue
           j = la_isamax(n,x,1)
           iter = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = zero
           end do
           x(j) = one
           kase = 1
           jump = 3
           return
           ! ................ entry   (jump = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_scopy(n,x,1,v,1)
           estold = est
           est = la_sasum(n,v,1)
           do i = 1,n
              if (nint(sign(one,x(i)),KIND=ilp) /= isgn(i)) go to 90
           end do
           ! repeated sign vector detected, hence algorithm has converged.
           go to 120
           90 continue
           ! test for cycling.
           if (est <= estold) go to 120
           do i = 1,n
              x(i) = sign(one,x(i))
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           jump = 4
           return
           ! ................ entry   (jump = 4)
           ! x has been overwritten by transpose(a)*x.
           110 continue
           jlast = j
           j = la_isamax(n,x,1)
           if ((x(jlast) /= abs(x(j))) .and. (iter < itmax)) then
              iter = iter + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           120 continue
           altsgn = one
           do i = 1,n
              x(i) = altsgn*(one + real(i - 1,KIND=sp)/real(n - 1,KIND=sp))
              altsgn = -altsgn
           end do
           kase = 1
           jump = 5
           return
           ! ................ entry   (jump = 5)
           ! x has been overwritten by a*x.
           140 continue
           temp = two*(la_sasum(n,x,1)/real(3*n,KIND=sp))
           if (temp > est) then
              call la_scopy(n,x,1,v,1)
              est = temp
           end if
           150 continue
           kase = 0
           return
     end subroutine la_slacon
     !> DLACON: estimates the 1-norm of a square, real matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     subroutine la_dlacon(n,v,x,isgn,est,kase)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(dp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(out) :: isgn(*)
           real(dp),intent(out) :: v(*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,iter,j,jlast,jump
           real(dp) :: altsgn,estold,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,nint,sign
           ! Save Statement
           save
           ! Executable Statements
           if (kase == 0) then
              do i = 1,n
                 x(i) = one/real(n,KIND=dp)
              end do
              kase = 1
              jump = 1
              return
           end if
           go to(20,40,70,110,140) jump
           ! ................ entry   (jump = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 150
           end if
           est = la_dasum(n,x,1)
           do i = 1,n
              x(i) = sign(one,x(i))
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           jump = 2
           return
           ! ................ entry   (jump = 2)
           ! first iteration.  x has been overwritten by transpose(a)*x.
           40 continue
           j = la_idamax(n,x,1)
           iter = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = zero
           end do
           x(j) = one
           kase = 1
           jump = 3
           return
           ! ................ entry   (jump = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_dcopy(n,x,1,v,1)
           estold = est
           est = la_dasum(n,v,1)
           do i = 1,n
              if (nint(sign(one,x(i)),KIND=ilp) /= isgn(i)) go to 90
           end do
           ! repeated sign vector detected, hence algorithm has converged.
           go to 120
           90 continue
           ! test for cycling.
           if (est <= estold) go to 120
           do i = 1,n
              x(i) = sign(one,x(i))
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           jump = 4
           return
           ! ................ entry   (jump = 4)
           ! x has been overwritten by transpose(a)*x.
           110 continue
           jlast = j
           j = la_idamax(n,x,1)
           if ((x(jlast) /= abs(x(j))) .and. (iter < itmax)) then
              iter = iter + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           120 continue
           altsgn = one
           do i = 1,n
              x(i) = altsgn*(one + real(i - 1,KIND=dp)/real(n - 1,KIND=dp))
              altsgn = -altsgn
           end do
           kase = 1
           jump = 5
           return
           ! ................ entry   (jump = 5)
           ! x has been overwritten by a*x.
           140 continue
           temp = two*(la_dasum(n,x,1)/real(3*n,KIND=dp))
           if (temp > est) then
              call la_dcopy(n,x,1,v,1)
              est = temp
           end if
           150 continue
           kase = 0
           return
     end subroutine la_dlacon
     !> QLACON: estimates the 1-norm of a square, real matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     subroutine la_qlacon(n,v,x,isgn,est,kase)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(qp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(out) :: isgn(*)
           real(qp),intent(out) :: v(*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,iter,j,jlast,jump
           real(qp) :: altsgn,estold,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,nint,sign
           ! Save Statement
           save
           ! Executable Statements
           if (kase == 0) then
              do i = 1,n
                 x(i) = one/real(n,KIND=qp)
              end do
              kase = 1
              jump = 1
              return
           end if
           go to(20,40,70,110,140) jump
           ! ................ entry   (jump = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 150
           end if
           est = la_qasum(n,x,1)
           do i = 1,n
              x(i) = sign(one,x(i))
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           jump = 2
           return
           ! ................ entry   (jump = 2)
           ! first iteration.  x has been overwritten by transpose(a)*x.
           40 continue
           j = la_iqamax(n,x,1)
           iter = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = zero
           end do
           x(j) = one
           kase = 1
           jump = 3
           return
           ! ................ entry   (jump = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_qcopy(n,x,1,v,1)
           estold = est
           est = la_qasum(n,v,1)
           do i = 1,n
              if (nint(sign(one,x(i)),KIND=ilp) /= isgn(i)) go to 90
           end do
           ! repeated sign vector detected, hence algorithm has converged.
           go to 120
           90 continue
           ! test for cycling.
           if (est <= estold) go to 120
           do i = 1,n
              x(i) = sign(one,x(i))
              isgn(i) = nint(x(i),KIND=ilp)
           end do
           kase = 2
           jump = 4
           return
           ! ................ entry   (jump = 4)
           ! x has been overwritten by transpose(a)*x.
           110 continue
           jlast = j
           j = la_iqamax(n,x,1)
           if ((x(jlast) /= abs(x(j))) .and. (iter < itmax)) then
              iter = iter + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           120 continue
           altsgn = one
           do i = 1,n
              x(i) = altsgn*(one + real(i - 1,KIND=qp)/real(n - 1,KIND=qp))
              altsgn = -altsgn
           end do
           kase = 1
           jump = 5
           return
           ! ................ entry   (jump = 5)
           ! x has been overwritten by a*x.
           140 continue
           temp = two*(la_qasum(n,x,1)/real(3*n,KIND=qp))
           if (temp > est) then
              call la_qcopy(n,x,1,v,1)
              est = temp
           end if
           150 continue
           kase = 0
           return
     end subroutine la_qlacon

     !> SLA_LIN_BERR: computes componentwise relative backward error from
     !> the formula
     !> max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
     !> where abs(Z) is the componentwise absolute value of the matrix
     !> or vector Z.

     pure subroutine la_sla_lin_berr(n,nz,nrhs,res,ayb,berr)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,nz,nrhs
           ! Array Arguments
           real(sp),intent(in) :: ayb(n,nrhs)
           real(sp),intent(out) :: berr(nrhs)
           real(sp),intent(in) :: res(n,nrhs)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: tmp,safe1
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! adding safe1 to the numerator guards against spuriously zero
           ! residuals.  a similar safeguard is in the sla_yyamv routine used
           ! to compute ayb.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (nz + 1)*safe1
           do j = 1,nrhs
              berr(j) = zero
              do i = 1,n
                 if (ayb(i,j) /= 0.0_sp) then
                    tmp = (safe1 + abs(res(i,j)))/ayb(i,j)
                    berr(j) = max(berr(j),tmp)
                 end if
           ! if ayb is exactly 0.0_sp (and if computed by sla_yyamv), then we know
           ! the true residual also must be exactly zero.
              end do
           end do
     end subroutine la_sla_lin_berr
     !> DLA_LIN_BERR: computes component-wise relative backward error from
     !> the formula
     !> max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
     !> where abs(Z) is the component-wise absolute value of the matrix
     !> or vector Z.

     pure subroutine la_dla_lin_berr(n,nz,nrhs,res,ayb,berr)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,nz,nrhs
           ! Array Arguments
           real(dp),intent(in) :: ayb(n,nrhs)
           real(dp),intent(out) :: berr(nrhs)
           real(dp),intent(in) :: res(n,nrhs)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: tmp,safe1
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! adding safe1 to the numerator guards against spuriously zero
           ! residuals.  a similar safeguard is in the sla_yyamv routine used
           ! to compute ayb.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (nz + 1)*safe1
           do j = 1,nrhs
              berr(j) = zero
              do i = 1,n
                 if (ayb(i,j) /= zero) then
                    tmp = (safe1 + abs(res(i,j)))/ayb(i,j)
                    berr(j) = max(berr(j),tmp)
                 end if
           ! if ayb is exactly 0.0_dp (and if computed by sla_yyamv), then we know
           ! the true residual also must be exactly zero.
              end do
           end do
     end subroutine la_dla_lin_berr
     !> QLA_LIN_BERR: computes component-wise relative backward error from
     !> the formula
     !> max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
     !> where abs(Z) is the component-wise absolute value of the matrix
     !> or vector Z.

     pure subroutine la_qla_lin_berr(n,nz,nrhs,res,ayb,berr)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,nz,nrhs
           ! Array Arguments
           real(qp),intent(in) :: ayb(n,nrhs)
           real(qp),intent(out) :: berr(nrhs)
           real(qp),intent(in) :: res(n,nrhs)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: tmp,safe1
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! adding safe1 to the numerator guards against spuriously zero
           ! residuals.  a similar safeguard is in the sla_yyamv routine used
           ! to compute ayb.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (nz + 1)*safe1
           do j = 1,nrhs
              berr(j) = zero
              do i = 1,n
                 if (ayb(i,j) /= zero) then
                    tmp = (safe1 + abs(res(i,j)))/ayb(i,j)
                    berr(j) = max(berr(j),tmp)
                 end if
           ! if ayb is exactly 0.0_qp (and if computed by sla_yyamv), then we know
           ! the true residual also must be exactly zero.
              end do
           end do
     end subroutine la_qla_lin_berr

     !> CLA_LIN_BERR: computes componentwise relative backward error from
     !> the formula
     !> max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
     !> where abs(Z) is the componentwise absolute value of the matrix
     !> or vector Z.

     pure subroutine la_cla_lin_berr(n,nz,nrhs,res,ayb,berr)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,nz,nrhs
           ! Array Arguments
           real(sp),intent(in) :: ayb(n,nrhs)
           real(sp),intent(out) :: berr(nrhs)
           complex(sp),intent(in) :: res(n,nrhs)
        ! =====================================================================
           ! Local Scalars
           real(sp) :: tmp,safe1
           integer(ilp) :: i,j
           complex(sp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           complex(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=sp)) + abs(aimag(cdum))
           ! Executable Statements
           ! adding safe1 to the numerator guards against spuriously zero
           ! residuals.  a similar safeguard is in the cla_yyamv routine used
           ! to compute ayb.
           safe1 = la_slamch('SAFE MINIMUM')
           safe1 = (nz + 1)*safe1
           do j = 1,nrhs
              berr(j) = zero
              do i = 1,n
                 if (ayb(i,j) /= 0.0_sp) then
                    tmp = (safe1 + cabs1(res(i,j)))/ayb(i,j)
                    berr(j) = max(berr(j),tmp)
                 end if
           ! if ayb is exactly 0.0_sp (and if computed by cla_yyamv), then we know
           ! the true residual also must be exactly zero.
              end do
           end do
     end subroutine la_cla_lin_berr
     !> ZLA_LIN_BERR: computes componentwise relative backward error from
     !> the formula
     !> max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
     !> where abs(Z) is the componentwise absolute value of the matrix
     !> or vector Z.

     pure subroutine la_zla_lin_berr(n,nz,nrhs,res,ayb,berr)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,nz,nrhs
           ! Array Arguments
           real(dp),intent(in) :: ayb(n,nrhs)
           real(dp),intent(out) :: berr(nrhs)
           complex(dp),intent(in) :: res(n,nrhs)
        ! =====================================================================
           ! Local Scalars
           real(dp) :: tmp,safe1
           integer(ilp) :: i,j
           complex(dp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           complex(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=dp)) + abs(aimag(cdum))
           ! Executable Statements
           ! adding safe1 to the numerator guards against spuriously zero
           ! residuals.  a similar safeguard is in the cla_yyamv routine used
           ! to compute ayb.
           safe1 = la_dlamch('SAFE MINIMUM')
           safe1 = (nz + 1)*safe1
           do j = 1,nrhs
              berr(j) = zero
              do i = 1,n
                 if (ayb(i,j) /= zero) then
                    tmp = (safe1 + cabs1(res(i,j)))/ayb(i,j)
                    berr(j) = max(berr(j),tmp)
                 end if
           ! if ayb is exactly 0.0_dp (and if computed by cla_yyamv), then we know
           ! the true residual also must be exactly zero.
              end do
           end do
     end subroutine la_zla_lin_berr
     !> WLA_LIN_BERR: computes componentwise relative backward error from
     !> the formula
     !> max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
     !> where abs(Z) is the componentwise absolute value of the matrix
     !> or vector Z.

     pure subroutine la_wla_lin_berr(n,nz,nrhs,res,ayb,berr)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,nz,nrhs
           ! Array Arguments
           real(qp),intent(in) :: ayb(n,nrhs)
           real(qp),intent(out) :: berr(nrhs)
           complex(qp),intent(in) :: res(n,nrhs)
        ! =====================================================================
           ! Local Scalars
           real(qp) :: tmp,safe1
           integer(ilp) :: i,j
           complex(qp) :: cdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           complex(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(cdum) = abs(real(cdum,KIND=qp)) + abs(aimag(cdum))
           ! Executable Statements
           ! adding safe1 to the numerator guards against spuriously zero
           ! residuals.  a similar safeguard is in the cla_yyamv routine used
           ! to compute ayb.
           safe1 = la_qlamch('SAFE MINIMUM')
           safe1 = (nz + 1)*safe1
           do j = 1,nrhs
              berr(j) = zero
              do i = 1,n
                 if (ayb(i,j) /= zero) then
                    tmp = (safe1 + cabs1(res(i,j)))/ayb(i,j)
                    berr(j) = max(berr(j),tmp)
                 end if
           ! if ayb is exactly 0.0_qp (and if computed by cla_yyamv), then we know
           ! the true residual also must be exactly zero.
              end do
           end do
     end subroutine la_wla_lin_berr

     !> CLACN2: estimates the 1-norm of a square, complex matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     pure subroutine la_clacn2(n,v,x,est,kase,isave)
        use la_constants_sp,only:one,two,czero,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(sp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(inout) :: isave(3)
           complex(sp),intent(out) :: v(*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,jlast
           real(sp) :: absxi,altsgn,estold,safmin,temp
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,real
           ! Executable Statements
           safmin = la_slamch('SAFE MINIMUM')
           if (kase == 0) then
              do i = 1,n
                 x(i) = cmplx(one/real(n,KIND=sp),KIND=sp)
              end do
              kase = 1
              isave(1) = 1
              return
           end if
           go to(20,40,70,90,120) isave(1)
           ! ................ entry   (isave( 1 ) = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 130
           end if
           est = la_scsum1(n,x,1)
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=sp)/absxi,aimag(x(i))/absxi,KIND=sp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           isave(1) = 2
           return
           ! ................ entry   (isave( 1 ) = 2)
           ! first iteration.  x has been overwritten by ctrans(a)*x.
           40 continue
           isave(2) = la_icmax1(n,x,1)
           isave(3) = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = czero
           end do
           x(isave(2)) = cone
           kase = 1
           isave(1) = 3
           return
           ! ................ entry   (isave( 1 ) = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_ccopy(n,x,1,v,1)
           estold = est
           est = la_scsum1(n,v,1)
           ! test for cycling.
           if (est <= estold) go to 100
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=sp)/absxi,aimag(x(i))/absxi,KIND=sp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           isave(1) = 4
           return
           ! ................ entry   (isave( 1 ) = 4)
           ! x has been overwritten by ctrans(a)*x.
           90 continue
           jlast = isave(2)
           isave(2) = la_icmax1(n,x,1)
           if ((abs(x(jlast)) /= abs(x(isave(2)))) .and. (isave(3) < itmax)) &
                     then
              isave(3) = isave(3) + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           100 continue
           altsgn = one
           do i = 1,n
              x(i) = cmplx(altsgn*(one + real(i - 1,KIND=sp)/real(n - 1,KIND=sp)),KIND=sp)

              altsgn = -altsgn
           end do
           kase = 1
           isave(1) = 5
           return
           ! ................ entry   (isave( 1 ) = 5)
           ! x has been overwritten by a*x.
           120 continue
           temp = two*(la_scsum1(n,x,1)/real(3*n,KIND=sp))
           if (temp > est) then
              call la_ccopy(n,x,1,v,1)
              est = temp
           end if
           130 continue
           kase = 0
           return
     end subroutine la_clacn2
     !> ZLACN2: estimates the 1-norm of a square, complex matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     pure subroutine la_zlacn2(n,v,x,est,kase,isave)
        use la_constants_dp,only:one,two,czero,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(dp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(inout) :: isave(3)
           complex(dp),intent(out) :: v(*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,jlast
           real(dp) :: absxi,altsgn,estold,safmin,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag
           ! Executable Statements
           safmin = la_dlamch('SAFE MINIMUM')
           if (kase == 0) then
              do i = 1,n
                 x(i) = cmplx(one/real(n,KIND=dp),KIND=dp)
              end do
              kase = 1
              isave(1) = 1
              return
           end if
           go to(20,40,70,90,120) isave(1)
           ! ................ entry   (isave( 1 ) = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 130
           end if
           est = la_dzsum1(n,x,1)
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=dp)/absxi,aimag(x(i))/absxi,KIND=dp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           isave(1) = 2
           return
           ! ................ entry   (isave( 1 ) = 2)
           ! first iteration.  x has been overwritten by ctrans(a)*x.
           40 continue
           isave(2) = la_izmax1(n,x,1)
           isave(3) = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = czero
           end do
           x(isave(2)) = cone
           kase = 1
           isave(1) = 3
           return
           ! ................ entry   (isave( 1 ) = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_zcopy(n,x,1,v,1)
           estold = est
           est = la_dzsum1(n,v,1)
           ! test for cycling.
           if (est <= estold) go to 100
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=dp)/absxi,aimag(x(i))/absxi,KIND=dp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           isave(1) = 4
           return
           ! ................ entry   (isave( 1 ) = 4)
           ! x has been overwritten by ctrans(a)*x.
           90 continue
           jlast = isave(2)
           isave(2) = la_izmax1(n,x,1)
           if ((abs(x(jlast)) /= abs(x(isave(2)))) .and. (isave(3) < itmax)) &
                     then
              isave(3) = isave(3) + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           100 continue
           altsgn = one
           do i = 1,n
              x(i) = cmplx(altsgn*(one + real(i - 1,KIND=dp)/real(n - 1,KIND=dp)),KIND=dp)

              altsgn = -altsgn
           end do
           kase = 1
           isave(1) = 5
           return
           ! ................ entry   (isave( 1 ) = 5)
           ! x has been overwritten by a*x.
           120 continue
           temp = two*(la_dzsum1(n,x,1)/real(3*n,KIND=dp))
           if (temp > est) then
              call la_zcopy(n,x,1,v,1)
              est = temp
           end if
           130 continue
           kase = 0
           return
     end subroutine la_zlacn2
     !> WLACN2: estimates the 1-norm of a square, complex matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     pure subroutine la_wlacn2(n,v,x,est,kase,isave)
        use la_constants_qp,only:one,two,czero,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(qp),intent(inout) :: est
           ! Array Arguments
           integer(ilp),intent(inout) :: isave(3)
           complex(qp),intent(out) :: v(*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,jlast
           real(qp) :: absxi,altsgn,estold,safmin,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag
           ! Executable Statements
           safmin = la_qlamch('SAFE MINIMUM')
           if (kase == 0) then
              do i = 1,n
                 x(i) = cmplx(one/real(n,KIND=qp),KIND=qp)
              end do
              kase = 1
              isave(1) = 1
              return
           end if
           go to(20,40,70,90,120) isave(1)
           ! ................ entry   (isave( 1 ) = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 130
           end if
           est = la_qwsum1(n,x,1)
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=qp)/absxi,aimag(x(i))/absxi,KIND=qp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           isave(1) = 2
           return
           ! ................ entry   (isave( 1 ) = 2)
           ! first iteration.  x has been overwritten by ctrans(a)*x.
           40 continue
           isave(2) = la_iwmax1(n,x,1)
           isave(3) = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = czero
           end do
           x(isave(2)) = cone
           kase = 1
           isave(1) = 3
           return
           ! ................ entry   (isave( 1 ) = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_wcopy(n,x,1,v,1)
           estold = est
           est = la_qwsum1(n,v,1)
           ! test for cycling.
           if (est <= estold) go to 100
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=qp)/absxi,aimag(x(i))/absxi,KIND=qp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           isave(1) = 4
           return
           ! ................ entry   (isave( 1 ) = 4)
           ! x has been overwritten by ctrans(a)*x.
           90 continue
           jlast = isave(2)
           isave(2) = la_iwmax1(n,x,1)
           if ((abs(x(jlast)) /= abs(x(isave(2)))) .and. (isave(3) < itmax)) &
                     then
              isave(3) = isave(3) + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           100 continue
           altsgn = one
           do i = 1,n
              x(i) = cmplx(altsgn*(one + real(i - 1,KIND=qp)/real(n - 1,KIND=qp)),KIND=qp)

              altsgn = -altsgn
           end do
           kase = 1
           isave(1) = 5
           return
           ! ................ entry   (isave( 1 ) = 5)
           ! x has been overwritten by a*x.
           120 continue
           temp = two*(la_qwsum1(n,x,1)/real(3*n,KIND=qp))
           if (temp > est) then
              call la_wcopy(n,x,1,v,1)
              est = temp
           end if
           130 continue
           kase = 0
           return
     end subroutine la_wlacn2

     !> CLACON: estimates the 1-norm of a square, complex matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     subroutine la_clacon(n,v,x,est,kase)
        use la_constants_sp,only:one,two,czero,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(sp),intent(inout) :: est
           ! Array Arguments
           complex(sp),intent(out) :: v(n)
           complex(sp),intent(inout) :: x(n)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,iter,j,jlast,jump
           real(sp) :: absxi,altsgn,estold,safmin,temp
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,real
           ! Save Statement
           save
           ! Executable Statements
           safmin = la_slamch('SAFE MINIMUM')
           if (kase == 0) then
              do i = 1,n
                 x(i) = cmplx(one/real(n,KIND=sp),KIND=sp)
              end do
              kase = 1
              jump = 1
              return
           end if
           go to(20,40,70,90,120) jump
           ! ................ entry   (jump = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 130
           end if
           est = la_scsum1(n,x,1)
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=sp)/absxi,aimag(x(i))/absxi,KIND=sp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           jump = 2
           return
           ! ................ entry   (jump = 2)
           ! first iteration.  x has been overwritten by ctrans(a)*x.
           40 continue
           j = la_icmax1(n,x,1)
           iter = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = czero
           end do
           x(j) = cone
           kase = 1
           jump = 3
           return
           ! ................ entry   (jump = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_ccopy(n,x,1,v,1)
           estold = est
           est = la_scsum1(n,v,1)
           ! test for cycling.
           if (est <= estold) go to 100
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=sp)/absxi,aimag(x(i))/absxi,KIND=sp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           jump = 4
           return
           ! ................ entry   (jump = 4)
           ! x has been overwritten by ctrans(a)*x.
           90 continue
           jlast = j
           j = la_icmax1(n,x,1)
           if ((abs(x(jlast)) /= abs(x(j))) .and. (iter < itmax)) then
              iter = iter + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           100 continue
           altsgn = one
           do i = 1,n
              x(i) = cmplx(altsgn*(one + real(i - 1,KIND=sp)/real(n - 1,KIND=sp)),KIND=sp)

              altsgn = -altsgn
           end do
           kase = 1
           jump = 5
           return
           ! ................ entry   (jump = 5)
           ! x has been overwritten by a*x.
           120 continue
           temp = two*(la_scsum1(n,x,1)/real(3*n,KIND=sp))
           if (temp > est) then
              call la_ccopy(n,x,1,v,1)
              est = temp
           end if
           130 continue
           kase = 0
           return
     end subroutine la_clacon
     !> ZLACON: estimates the 1-norm of a square, complex matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     subroutine la_zlacon(n,v,x,est,kase)
        use la_constants_dp,only:one,two,czero,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(dp),intent(inout) :: est
           ! Array Arguments
           complex(dp),intent(out) :: v(n)
           complex(dp),intent(inout) :: x(n)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,iter,j,jlast,jump
           real(dp) :: absxi,altsgn,estold,safmin,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag
           ! Save Statement
           save
           ! Executable Statements
           safmin = la_dlamch('SAFE MINIMUM')
           if (kase == 0) then
              do i = 1,n
                 x(i) = cmplx(one/real(n,KIND=dp),KIND=dp)
              end do
              kase = 1
              jump = 1
              return
           end if
           go to(20,40,70,90,120) jump
           ! ................ entry   (jump = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 130
           end if
           est = la_dzsum1(n,x,1)
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=dp)/absxi,aimag(x(i))/absxi,KIND=dp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           jump = 2
           return
           ! ................ entry   (jump = 2)
           ! first iteration.  x has been overwritten by ctrans(a)*x.
           40 continue
           j = la_izmax1(n,x,1)
           iter = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = czero
           end do
           x(j) = cone
           kase = 1
           jump = 3
           return
           ! ................ entry   (jump = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_zcopy(n,x,1,v,1)
           estold = est
           est = la_dzsum1(n,v,1)
           ! test for cycling.
           if (est <= estold) go to 100
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=dp)/absxi,aimag(x(i))/absxi,KIND=dp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           jump = 4
           return
           ! ................ entry   (jump = 4)
           ! x has been overwritten by ctrans(a)*x.
           90 continue
           jlast = j
           j = la_izmax1(n,x,1)
           if ((abs(x(jlast)) /= abs(x(j))) .and. (iter < itmax)) then
              iter = iter + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           100 continue
           altsgn = one
           do i = 1,n
              x(i) = cmplx(altsgn*(one + real(i - 1,KIND=dp)/real(n - 1,KIND=dp)),KIND=dp)

              altsgn = -altsgn
           end do
           kase = 1
           jump = 5
           return
           ! ................ entry   (jump = 5)
           ! x has been overwritten by a*x.
           120 continue
           temp = two*(la_dzsum1(n,x,1)/real(3*n,KIND=dp))
           if (temp > est) then
              call la_zcopy(n,x,1,v,1)
              est = temp
           end if
           130 continue
           kase = 0
           return
     end subroutine la_zlacon
     !> WLACON: estimates the 1-norm of a square, complex matrix A.
     !> Reverse communication is used for evaluating matrix-vector products.

     subroutine la_wlacon(n,v,x,est,kase)
        use la_constants_qp,only:one,two,czero,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(inout) :: kase
           integer(ilp),intent(in) :: n
           real(qp),intent(inout) :: est
           ! Array Arguments
           complex(qp),intent(out) :: v(n)
           complex(qp),intent(inout) :: x(n)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           integer(ilp) :: i,iter,j,jlast,jump
           real(qp) :: absxi,altsgn,estold,safmin,temp
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,aimag
           ! Save Statement
           save
           ! Executable Statements
           safmin = la_qlamch('SAFE MINIMUM')
           if (kase == 0) then
              do i = 1,n
                 x(i) = cmplx(one/real(n,KIND=qp),KIND=qp)
              end do
              kase = 1
              jump = 1
              return
           end if
           go to(20,40,70,90,120) jump
           ! ................ entry   (jump = 1)
           ! first iteration.  x has been overwritten by a*x.
           20 continue
           if (n == 1) then
              v(1) = x(1)
              est = abs(v(1))
              ! ... quit
              go to 130
           end if
           est = la_qwsum1(n,x,1)
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=qp)/absxi,aimag(x(i))/absxi,KIND=qp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           jump = 2
           return
           ! ................ entry   (jump = 2)
           ! first iteration.  x has been overwritten by ctrans(a)*x.
           40 continue
           j = la_iwmax1(n,x,1)
           iter = 2
           ! main loop - iterations 2,3,...,itmax.
           50 continue
           do i = 1,n
              x(i) = czero
           end do
           x(j) = cone
           kase = 1
           jump = 3
           return
           ! ................ entry   (jump = 3)
           ! x has been overwritten by a*x.
           70 continue
           call la_wcopy(n,x,1,v,1)
           estold = est
           est = la_qwsum1(n,v,1)
           ! test for cycling.
           if (est <= estold) go to 100
           do i = 1,n
              absxi = abs(x(i))
              if (absxi > safmin) then
                 x(i) = cmplx(real(x(i),KIND=qp)/absxi,aimag(x(i))/absxi,KIND=qp)

              else
                 x(i) = cone
              end if
           end do
           kase = 2
           jump = 4
           return
           ! ................ entry   (jump = 4)
           ! x has been overwritten by ctrans(a)*x.
           90 continue
           jlast = j
           j = la_iwmax1(n,x,1)
           if ((abs(x(jlast)) /= abs(x(j))) .and. (iter < itmax)) then
              iter = iter + 1
              go to 50
           end if
           ! iteration complete.  final stage.
           100 continue
           altsgn = one
           do i = 1,n
              x(i) = cmplx(altsgn*(one + real(i - 1,KIND=qp)/real(n - 1,KIND=qp)),KIND=qp)

              altsgn = -altsgn
           end do
           kase = 1
           jump = 5
           return
           ! ................ entry   (jump = 5)
           ! x has been overwritten by a*x.
           120 continue
           temp = two*(la_qwsum1(n,x,1)/real(3*n,KIND=qp))
           if (temp > est) then
              call la_wcopy(n,x,1,v,1)
              est = temp
           end if
           130 continue
           kase = 0
           return
     end subroutine la_wlacon

end module la_lapack_solve_aux
