!> Cholesky drivers: positive definite, packed, banded and tridiagonal systems
module la_lapack_solve_chol
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level3_gen
     use la_blas_level3_sym
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_mnorm
     use la_lapack_solve_chol_comp
     use la_lapack_solve_ldl_comp4
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sppsv
     public :: la_sppsvx
     public :: la_sptsv
     public :: la_sptsvx
     public :: la_spbsv
     public :: la_spbsvx
     public :: la_sposv
     public :: la_sposvx
     public :: la_dppsv
     public :: la_dppsvx
     public :: la_dptsv
     public :: la_dptsvx
     public :: la_dsposv
     public :: la_dpbsv
     public :: la_dpbsvx
     public :: la_dposv
     public :: la_dposvx
     public :: la_qppsv
     public :: la_qppsvx
     public :: la_qptsv
     public :: la_qptsvx
     public :: la_qdposv
     public :: la_qpbsv
     public :: la_qpbsvx
     public :: la_qposv
     public :: la_qposvx
     public :: la_cppsv
     public :: la_cppsvx
     public :: la_cpbsv
     public :: la_cpbsvx
     public :: la_cposv
     public :: la_cposvx
     public :: la_cptsv
     public :: la_cptsvx
     public :: la_zppsv
     public :: la_zppsvx
     public :: la_zcposv
     public :: la_zpbsv
     public :: la_zpbsvx
     public :: la_zposv
     public :: la_zposvx
     public :: la_zptsv
     public :: la_zptsvx
     public :: la_wppsv
     public :: la_wppsvx
     public :: la_wzposv
     public :: la_wpbsv
     public :: la_wpbsvx
     public :: la_wposv
     public :: la_wposvx
     public :: la_wptsv
     public :: la_wptsvx

     contains

     !> SPPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T* U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_sppsv(uplo,n,nrhs,ap,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: ap(*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SPPSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_spptrf(uplo,n,ap,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_spptrs(uplo,n,nrhs,ap,b,ldb,info)
           end if
           return
     end subroutine la_sppsv
     !> DPPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T* U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_dppsv(uplo,n,nrhs,ap,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: ap(*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DPPSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_dpptrf(uplo,n,ap,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dpptrs(uplo,n,nrhs,ap,b,ldb,info)
           end if
           return
     end subroutine la_dppsv
     !> QPPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T* U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_qppsv(uplo,n,nrhs,ap,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: ap(*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QPPSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_qpptrf(uplo,n,ap,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qpptrs(uplo,n,nrhs,ap,b,ldb,info)
           end if
           return
     end subroutine la_qppsv

     !> SPPSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_sppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,rcond,ferr, &
                berr,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: afp(*),ap(*),b(ldb,*),s(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(sp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -7
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -8
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -10
                 else if (ldx < max(1,n)) then
                    info = -12
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SPPSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_sppequ(uplo,n,ap,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_slaqsp(uplo,n,ap,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t * u or a = l * l**t.
              call la_scopy(n*(n + 1)/2,ap,1,afp,1)
              call la_spptrf(uplo,n,afp,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_slansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_sppcon(uplo,n,afp,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_spptrs(uplo,n,nrhs,afp,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_spprfs(uplo,n,nrhs,ap,afp,b,ldb,x,ldx,ferr,berr,work,iwork, &
                     info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_sppsvx
     !> DPPSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,rcond,ferr, &
                berr,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: afp(*),ap(*),b(ldb,*),s(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(dp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -7
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -8
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -10
                 else if (ldx < max(1,n)) then
                    info = -12
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DPPSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_dppequ(uplo,n,ap,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_dlaqsp(uplo,n,ap,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t * u or a = l * l**t.
              call la_dcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_dpptrf(uplo,n,afp,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_dlansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_dppcon(uplo,n,afp,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dpptrs(uplo,n,nrhs,afp,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_dpprfs(uplo,n,nrhs,ap,afp,b,ldb,x,ldx,ferr,berr,work,iwork, &
                     info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dppsvx
     !> QPPSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,rcond,ferr, &
                berr,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: afp(*),ap(*),b(ldb,*),s(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(qp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -7
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -8
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -10
                 else if (ldx < max(1,n)) then
                    info = -12
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QPPSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_qppequ(uplo,n,ap,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_qlaqsp(uplo,n,ap,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t * u or a = l * l**t.
              call la_qcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_qpptrf(uplo,n,afp,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_qlansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_qppcon(uplo,n,afp,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qpptrs(uplo,n,nrhs,afp,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_qpprfs(uplo,n,nrhs,ap,afp,b,ldb,x,ldx,ferr,berr,work,iwork, &
                     info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qppsvx

     !> SPTSV: computes the solution to a real system of linear equations
     !> A*X = B, where A is an N-by-N symmetric positive definite tridiagonal
     !> matrix, and X and B are N-by-NRHS matrices.
     !> A is factored as A = L*D*L**T, and the factored form of A is then
     !> used to solve the system of equations.

     pure subroutine la_sptsv(n,nrhs,d,e,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: b(ldb,*),d(*),e(*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SPTSV ',-info)
              return
           end if
           ! compute the l*d*l**t (or u**t*d*u) factorization of a.
           call la_spttrf(n,d,e,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_spttrs(n,nrhs,d,e,b,ldb,info)
           end if
           return
     end subroutine la_sptsv
     !> DPTSV: computes the solution to a real system of linear equations
     !> A*X = B, where A is an N-by-N symmetric positive definite tridiagonal
     !> matrix, and X and B are N-by-NRHS matrices.
     !> A is factored as A = L*D*L**T, and the factored form of A is then
     !> used to solve the system of equations.

     pure subroutine la_dptsv(n,nrhs,d,e,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: b(ldb,*),d(*),e(*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DPTSV ',-info)
              return
           end if
           ! compute the l*d*l**t (or u**t*d*u) factorization of a.
           call la_dpttrf(n,d,e,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dpttrs(n,nrhs,d,e,b,ldb,info)
           end if
           return
     end subroutine la_dptsv
     !> QPTSV: computes the solution to a real system of linear equations
     !> A*X = B, where A is an N-by-N symmetric positive definite tridiagonal
     !> matrix, and X and B are N-by-NRHS matrices.
     !> A is factored as A = L*D*L**T, and the factored form of A is then
     !> used to solve the system of equations.

     pure subroutine la_qptsv(n,nrhs,d,e,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: b(ldb,*),d(*),e(*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QPTSV ',-info)
              return
           end if
           ! compute the l*d*l**t (or u**t*d*u) factorization of a.
           call la_qpttrf(n,d,e,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qpttrs(n,nrhs,d,e,b,ldb,info)
           end if
           return
     end subroutine la_qptsv

     !> SPTSVX: uses the factorization A = L*D*L**T to compute the solution
     !> to a real system of linear equations A*X = B, where A is an N-by-N
     !> symmetric positive definite tridiagonal matrix and X and B are
     !> N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_sptsvx(fact,n,nrhs,d,e,df,ef,b,ldb,x,ldx,rcond,ferr,berr, &
                work,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(in) :: b(ldb,*),d(*),e(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(sp),intent(inout) :: df(*),ef(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('SPTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the l*d*l**t (or u**t*d*u) factorization of a.
              call la_scopy(n,d,1,df,1)
              if (n > 1) call la_scopy(n - 1,e,1,ef,1)
              call la_spttrf(n,df,ef,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_slanst('1',n,d,e)
           ! compute the reciprocal of the condition number of a.
           call la_sptcon(n,df,ef,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_spttrs(n,nrhs,df,ef,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_sptrfs(n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,work,info)

           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_sptsvx
     !> DPTSVX: uses the factorization A = L*D*L**T to compute the solution
     !> to a real system of linear equations A*X = B, where A is an N-by-N
     !> symmetric positive definite tridiagonal matrix and X and B are
     !> N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_dptsvx(fact,n,nrhs,d,e,df,ef,b,ldb,x,ldx,rcond,ferr,berr, &
                work,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(in) :: b(ldb,*),d(*),e(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(dp),intent(inout) :: df(*),ef(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('DPTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the l*d*l**t (or u**t*d*u) factorization of a.
              call la_dcopy(n,d,1,df,1)
              if (n > 1) call la_dcopy(n - 1,e,1,ef,1)
              call la_dpttrf(n,df,ef,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_dlanst('1',n,d,e)
           ! compute the reciprocal of the condition number of a.
           call la_dptcon(n,df,ef,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dpttrs(n,nrhs,df,ef,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_dptrfs(n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,work,info)

           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dptsvx
     !> QPTSVX: uses the factorization A = L*D*L**T to compute the solution
     !> to a real system of linear equations A*X = B, where A is an N-by-N
     !> symmetric positive definite tridiagonal matrix and X and B are
     !> N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_qptsvx(fact,n,nrhs,d,e,df,ef,b,ldb,x,ldx,rcond,ferr,berr, &
                work,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(in) :: b(ldb,*),d(*),e(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
           real(qp),intent(inout) :: df(*),ef(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('QPTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the l*d*l**t (or u**t*d*u) factorization of a.
              call la_qcopy(n,d,1,df,1)
              if (n > 1) call la_qcopy(n - 1,e,1,ef,1)
              call la_qpttrf(n,df,ef,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_qlanst('1',n,d,e)
           ! compute the reciprocal of the condition number of a.
           call la_qptcon(n,df,ef,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qpttrs(n,nrhs,df,ef,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_qptrfs(n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,work,info)

           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qptsvx

     !> DSPOSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> DSPOSV first attempts to factorize the matrix in SINGLE PRECISION
     !> and use this factorization within an iterative refinement procedure
     !> to produce a solution with DOUBLE PRECISION normwise backward error
     !> quality (see below). If the approach fails the method switches to a
     !> DOUBLE PRECISION factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio SINGLE PRECISION performance over DOUBLE PRECISION
     !> performance is too small. A reasonable strategy should take the
     !> number of right-hand sides and the size of the matrix into account.
     !> This might be done with a call to ILAENV in the future. Up to now, we
     !> always try iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by DLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_dsposv(uplo,n,nrhs,a,lda,b,ldb,x,ldx,work,swork,iter,info)
        use la_constants_dp,only:negone,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(sp),intent(out) :: swork(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: b(ldb,*)
           real(dp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(dp),parameter :: bwdmax = 1.0e+00_dp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(dp) :: anrm,cte,eps,rnrm,xnrm
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DSPOSV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip single precision iterative refinement if a priori slower
           ! than double precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_dlansy('I',uplo,n,a,lda,work)
           eps = la_dlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=dp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from double precision to single precision and store the
           ! result in sx.
           call la_dlag2s(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from double precision to single precision and store the
           ! result in sa.
           call la_dlat2s(uplo,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the cholesky factorization of sa.
           call la_spotrf(uplo,n,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_spotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
           ! convert sx back to double precision
           call la_slag2d(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_dlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_dsymm('LEFT',uplo,n,nrhs,negone,a,lda,x,ldx,one,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = abs(x(la_idamax(n,x(1,i),1),i))
              rnrm = abs(work(la_idamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from double precision to single precision
              ! and store the result in sx.
              call la_dlag2s(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_spotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
              ! convert sx back to double precision and update the current
              ! iterate.
              call la_slag2d(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_daxpy(n,one,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_dlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_dsymm('L',uplo,n,nrhs,negone,a,lda,x,ldx,one,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = abs(x(la_idamax(n,x(1,i),1),i))
                 rnrm = abs(work(la_idamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow
           ! up on double precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to double precision.
           call la_dpotrf(uplo,n,a,lda,info)
           if (info /= 0) return
           call la_dlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_dpotrs(uplo,n,nrhs,a,lda,x,ldx,info)
           return
     end subroutine la_dsposv
     !> QDPOSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> QDPOSV first attempts to factorize the matrix in SINGLE PRECISION
     !> and use this factorization within an iterative refinement procedure
     !> to produce a solution with QUAD PRECISION normwise backward error
     !> quality (see below). If the approach fails the method switches to a
     !> QUAD PRECISION factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio SINGLE PRECISION performance over QUAD PRECISION
     !> performance is too small. A reasonable strategy should take the
     !> number of right-hand sides and the size of the matrix into account.
     !> This might be done with a call to ILAENV in the future. Up to now, we
     !> always try iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by QLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_qdposv(uplo,n,nrhs,a,lda,b,ldb,x,ldx,work,swork,iter,info)
        use la_constants_qp,only:negone,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(dp),intent(out) :: swork(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: b(ldb,*)
           real(qp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(qp),parameter :: bwdmax = 1.0e+00_qp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(qp) :: anrm,cte,eps,rnrm,xnrm
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QDPOSV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip double precision iterative refinement if a priori slower
           ! than quad precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_qlansy('I',uplo,n,a,lda,work)
           eps = la_qlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=qp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from quad precision to double precision and store the
           ! result in sx.
           call la_qlag2d(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from quad precision to double precision and store the
           ! result in sa.
           call la_qlat2d(uplo,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the cholesky factorization of sa.
           call la_dpotrf(uplo,n,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_dpotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
           ! convert sx back to quad precision
           call la_dlag2q(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_qlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_qsymm('LEFT',uplo,n,nrhs,negone,a,lda,x,ldx,one,work,n)
           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = abs(x(la_iqamax(n,x(1,i),1),i))
              rnrm = abs(work(la_iqamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from quad precision to double precision
              ! and store the result in sx.
              call la_qlag2d(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_dpotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
              ! convert sx back to quad precision and update the current
              ! iterate.
              call la_dlag2q(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_qaxpy(n,one,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_qlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_qsymm('L',uplo,n,nrhs,negone,a,lda,x,ldx,one,work,n)
              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = abs(x(la_iqamax(n,x(1,i),1),i))
                 rnrm = abs(work(la_iqamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow
           ! up on quad precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to quad precision.
           call la_qpotrf(uplo,n,a,lda,info)
           if (info /= 0) return
           call la_qlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_qpotrs(uplo,n,nrhs,a,lda,x,ldx,info)
           return
     end subroutine la_qdposv

     !> SPBSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T * U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular band matrix, and L is a lower
     !> triangular band matrix, with the same number of superdiagonals or
     !> subdiagonals as A.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_spbsv(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kd < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < kd + 1) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('SPBSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_spbtrf(uplo,n,kd,ab,ldab,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_spbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
           end if
           return
     end subroutine la_spbsv
     !> DPBSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T * U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular band matrix, and L is a lower
     !> triangular band matrix, with the same number of superdiagonals or
     !> subdiagonals as A.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_dpbsv(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kd < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < kd + 1) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DPBSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_dpbtrf(uplo,n,kd,ab,ldab,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
           end if
           return
     end subroutine la_dpbsv
     !> QPBSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T * U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular band matrix, and L is a lower
     !> triangular band matrix, with the same number of superdiagonals or
     !> subdiagonals as A.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_qpbsv(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kd < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < kd + 1) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QPBSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_qpbtrf(uplo,n,kd,ab,ldab,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
           end if
           return
     end subroutine la_qpbsv

     !> SPBSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_spbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x, &
               ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldafb,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),s(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ,upper
           integer(ilp) :: i,infequ,j,j1,j2
           real(sp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           upper = la_lsame(uplo,'U')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kd < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           else if (ldafb < kd + 1) then
              info = -9
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -10
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -11
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -13
                 else if (ldx < max(1,n)) then
                    info = -15
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SPBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_spbequ(uplo,n,kd,ab,ldab,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_slaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t *u or a = l*l**t.
              if (upper) then
                 do j = 1,n
                    j1 = max(j - kd,1)
                    call la_scopy(j - j1 + 1,ab(kd + 1 - j + j1,j),1,afb(kd + 1 - j + j1,j),1)

                 end do
              else
                 do j = 1,n
                    j2 = min(j + kd,n)
                    call la_scopy(j2 - j + 1,ab(1,j),1,afb(1,j),1)
                 end do
              end if
              call la_spbtrf(uplo,n,kd,afb,ldafb,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_slansb('1',uplo,n,kd,ab,ldab,work)
           ! compute the reciprocal of the condition number of a.
           call la_spbcon(uplo,n,kd,afb,ldafb,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_spbtrs(uplo,n,kd,nrhs,afb,ldafb,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_spbrfs(uplo,n,kd,nrhs,ab,ldab,afb,ldafb,b,ldb,x,ldx,ferr,berr, &
                      work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_spbsvx
     !> DPBSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dpbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x, &
               ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldafb,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),s(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ,upper
           integer(ilp) :: i,infequ,j,j1,j2
           real(dp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           upper = la_lsame(uplo,'U')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kd < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           else if (ldafb < kd + 1) then
              info = -9
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -10
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -11
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -13
                 else if (ldx < max(1,n)) then
                    info = -15
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DPBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_dpbequ(uplo,n,kd,ab,ldab,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_dlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t *u or a = l*l**t.
              if (upper) then
                 do j = 1,n
                    j1 = max(j - kd,1)
                    call la_dcopy(j - j1 + 1,ab(kd + 1 - j + j1,j),1,afb(kd + 1 - j + j1,j),1)

                 end do
              else
                 do j = 1,n
                    j2 = min(j + kd,n)
                    call la_dcopy(j2 - j + 1,ab(1,j),1,afb(1,j),1)
                 end do
              end if
              call la_dpbtrf(uplo,n,kd,afb,ldafb,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_dlansb('1',uplo,n,kd,ab,ldab,work)
           ! compute the reciprocal of the condition number of a.
           call la_dpbcon(uplo,n,kd,afb,ldafb,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dpbtrs(uplo,n,kd,nrhs,afb,ldafb,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_dpbrfs(uplo,n,kd,nrhs,ab,ldab,afb,ldafb,b,ldb,x,ldx,ferr,berr, &
                      work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dpbsvx
     !> QPBSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qpbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x, &
               ldx,rcond,ferr,berr,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldafb,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*),s(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ,upper
           integer(ilp) :: i,infequ,j,j1,j2
           real(qp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           upper = la_lsame(uplo,'U')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kd < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           else if (ldafb < kd + 1) then
              info = -9
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -10
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -11
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -13
                 else if (ldx < max(1,n)) then
                    info = -15
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QPBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_qpbequ(uplo,n,kd,ab,ldab,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_qlaqsb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t *u or a = l*l**t.
              if (upper) then
                 do j = 1,n
                    j1 = max(j - kd,1)
                    call la_qcopy(j - j1 + 1,ab(kd + 1 - j + j1,j),1,afb(kd + 1 - j + j1,j),1)

                 end do
              else
                 do j = 1,n
                    j2 = min(j + kd,n)
                    call la_qcopy(j2 - j + 1,ab(1,j),1,afb(1,j),1)
                 end do
              end if
              call la_qpbtrf(uplo,n,kd,afb,ldafb,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_qlansb('1',uplo,n,kd,ab,ldab,work)
           ! compute the reciprocal of the condition number of a.
           call la_qpbcon(uplo,n,kd,afb,ldafb,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qpbtrs(uplo,n,kd,nrhs,afb,ldafb,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_qpbrfs(uplo,n,kd,nrhs,ab,ldab,afb,ldafb,b,ldb,x,ldx,ferr,berr, &
                      work,iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qpbsvx

     !> SPOSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T* U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_sposv(uplo,n,nrhs,a,lda,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SPOSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_spotrf(uplo,n,a,lda,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_spotrs(uplo,n,nrhs,a,lda,b,ldb,info)
           end if
           return
     end subroutine la_sposv
     !> DPOSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T* U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_dposv(uplo,n,nrhs,a,lda,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DPOSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_dpotrf(uplo,n,a,lda,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
           end if
           return
     end subroutine la_dposv
     !> QPOSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**T* U,  if UPLO = 'U', or
     !> A = L * L**T,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_qposv(uplo,n,nrhs,a,lda,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QPOSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**t*u or a = l*l**t.
           call la_qpotrf(uplo,n,a,lda,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
           end if
           return
     end subroutine la_qposv

     !> SPOSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_sposvx(fact,uplo,n,nrhs,a,lda,af,ldaf,equed,s,b,ldb,x,ldx, &
               rcond,ferr,berr,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),s(*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(sp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -9
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -10
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -12
                 else if (ldx < max(1,n)) then
                    info = -14
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SPOSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_spoequ(n,a,lda,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_slaqsy(uplo,n,a,lda,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t *u or a = l*l**t.
              call la_slacpy(uplo,n,n,a,lda,af,ldaf)
              call la_spotrf(uplo,n,af,ldaf,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_slansy('1',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_spocon(uplo,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_spotrs(uplo,n,nrhs,af,ldaf,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_sporfs(uplo,n,nrhs,a,lda,af,ldaf,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_sposvx
     !> DPOSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dposvx(fact,uplo,n,nrhs,a,lda,af,ldaf,equed,s,b,ldb,x,ldx, &
               rcond,ferr,berr,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),s(*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(dp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -9
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -10
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -12
                 else if (ldx < max(1,n)) then
                    info = -14
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DPOSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_dpoequ(n,a,lda,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_dlaqsy(uplo,n,a,lda,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t *u or a = l*l**t.
              call la_dlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_dpotrf(uplo,n,af,ldaf,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_dlansy('1',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_dpocon(uplo,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dpotrs(uplo,n,nrhs,af,ldaf,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_dporfs(uplo,n,nrhs,a,lda,af,ldaf,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dposvx
     !> QPOSVX: uses the Cholesky factorization A = U**T*U or A = L*L**T to
     !> compute the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qposvx(fact,uplo,n,nrhs,a,lda,af,ldaf,equed,s,b,ldb,x,ldx, &
               rcond,ferr,berr,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*),s(*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(qp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -9
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -10
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -12
                 else if (ldx < max(1,n)) then
                    info = -14
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QPOSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_qpoequ(n,a,lda,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_qlaqsy(uplo,n,a,lda,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**t *u or a = l*l**t.
              call la_qlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_qpotrf(uplo,n,af,ldaf,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_qlansy('1',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_qpocon(uplo,n,af,ldaf,anorm,rcond,work,iwork,info)
           ! compute the solution matrix x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qpotrs(uplo,n,nrhs,af,ldaf,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_qporfs(uplo,n,nrhs,a,lda,af,ldaf,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qposvx

     !> CPPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H * U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_cppsv(uplo,n,nrhs,ap,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(inout) :: ap(*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CPPSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_cpptrf(uplo,n,ap,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_cpptrs(uplo,n,nrhs,ap,b,ldb,info)
           end if
           return
     end subroutine la_cppsv
     !> ZPPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H * U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_zppsv(uplo,n,nrhs,ap,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(inout) :: ap(*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZPPSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_zpptrf(uplo,n,ap,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zpptrs(uplo,n,nrhs,ap,b,ldb,info)
           end if
           return
     end subroutine la_zppsv
     !> WPPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H * U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular matrix and L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_wppsv(uplo,n,nrhs,ap,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(inout) :: ap(*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WPPSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_wpptrf(uplo,n,ap,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wpptrs(uplo,n,nrhs,ap,b,ldb,info)
           end if
           return
     end subroutine la_wppsv

     !> CPPSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_cppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,rcond,ferr, &
                berr,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(sp),intent(inout) :: s(*)
           complex(sp),intent(inout) :: afp(*),ap(*),b(ldb,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(sp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -7
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -8
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -10
                 else if (ldx < max(1,n)) then
                    info = -12
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CPPSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_cppequ(uplo,n,ap,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_claqhp(uplo,n,ap,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h * u or a = l * l**h.
              call la_ccopy(n*(n + 1)/2,ap,1,afp,1)
              call la_cpptrf(uplo,n,afp,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_cppcon(uplo,n,afp,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cpptrs(uplo,n,nrhs,afp,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_cpprfs(uplo,n,nrhs,ap,afp,b,ldb,x,ldx,ferr,berr,work,rwork, &
                     info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_cppsvx
     !> ZPPSVX: uses the Cholesky factorization A = U**H * U or A = L * L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,rcond,ferr, &
                berr,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(dp),intent(inout) :: s(*)
           complex(dp),intent(inout) :: afp(*),ap(*),b(ldb,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(dp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -7
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -8
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -10
                 else if (ldx < max(1,n)) then
                    info = -12
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZPPSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_zppequ(uplo,n,ap,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_zlaqhp(uplo,n,ap,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h * u or a = l * l**h.
              call la_zcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_zpptrf(uplo,n,afp,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zppcon(uplo,n,afp,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zpptrs(uplo,n,nrhs,afp,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_zpprfs(uplo,n,nrhs,ap,afp,b,ldb,x,ldx,ferr,berr,work,rwork, &
                     info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zppsvx
     !> WPPSVX: uses the Cholesky factorization A = U**H * U or A = L * L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix stored in
     !> packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wppsvx(fact,uplo,n,nrhs,ap,afp,equed,s,b,ldb,x,ldx,rcond,ferr, &
                berr,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(qp),intent(inout) :: s(*)
           complex(qp),intent(inout) :: afp(*),ap(*),b(ldb,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(qp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -7
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -8
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -10
                 else if (ldx < max(1,n)) then
                    info = -12
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WPPSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_wppequ(uplo,n,ap,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_wlaqhp(uplo,n,ap,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h * u or a = l * l**h.
              call la_wcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_wpptrf(uplo,n,afp,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_wppcon(uplo,n,afp,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wpptrs(uplo,n,nrhs,afp,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_wpprfs(uplo,n,nrhs,ap,afp,b,ldb,x,ldx,ferr,berr,work,rwork, &
                     info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_wppsvx

     !> ZCPOSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> ZCPOSV first attempts to factorize the matrix in COMPLEX and use this
     !> factorization within an iterative refinement procedure to produce a
     !> solution with COMPLEX*16 normwise backward error quality (see below).
     !> If the approach fails the method switches to a COMPLEX*16
     !> factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio COMPLEX performance over COMPLEX*16 performance is too
     !> small. A reasonable strategy should take the number of right-hand
     !> sides and the size of the matrix into account. This might be done
     !> with a call to ILAENV in the future. Up to now, we always try
     !> iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by DLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_zcposv(uplo,n,nrhs,a,lda,b,ldb,x,ldx,work,swork,rwork,iter, &
               info)
        use la_constants_dp,only:cone,cnegone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(sp),intent(out) :: swork(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: b(ldb,*)
           complex(dp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(dp),parameter :: bwdmax = 1.0e+00_dp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(dp) :: anrm,cte,eps,rnrm,xnrm
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZCPOSV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip single precision iterative refinement if a priori slower
           ! than double precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_zlanhe('I',uplo,n,a,lda,rwork)
           eps = la_dlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=dp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from double precision to single precision and store the
           ! result in sx.
           call la_zlag2c(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from double precision to single precision and store the
           ! result in sa.
           call la_zlat2c(uplo,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the cholesky factorization of sa.
           call la_cpotrf(uplo,n,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_cpotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
           ! convert sx back to complex*16
           call la_clag2z(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_zlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_zhemm('LEFT',uplo,n,nrhs,cnegone,a,lda,x,ldx,cone,work,n)

           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = cabs1(x(la_izamax(n,x(1,i),1),i))
              rnrm = cabs1(work(la_izamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from double precision to single precision
              ! and store the result in sx.
              call la_zlag2c(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_cpotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
              ! convert sx back to double precision and update the current
              ! iterate.
              call la_clag2z(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_zaxpy(n,cone,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_zlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_zhemm('L',uplo,n,nrhs,cnegone,a,lda,x,ldx,cone,work,n)

              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = cabs1(x(la_izamax(n,x(1,i),1),i))
                 rnrm = cabs1(work(la_izamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow
           ! up on double precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to double precision.
           call la_zpotrf(uplo,n,a,lda,info)
           if (info /= 0) return
           call la_zlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_zpotrs(uplo,n,nrhs,a,lda,x,ldx,info)
           return
     end subroutine la_zcposv
     !> WZPOSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> WZPOSV first attempts to factorize the matrix in COMPLEX and use this
     !> factorization within an iterative refinement procedure to produce a
     !> solution with COMPLEX*16 normwise backward error quality (see below).
     !> If the approach fails the method switches to a COMPLEX*16
     !> factorization and solve.
     !> The iterative refinement is not going to be a winning strategy if
     !> the ratio COMPLEX performance over COMPLEX*16 performance is too
     !> small. A reasonable strategy should take the number of right-hand
     !> sides and the size of the matrix into account. This might be done
     !> with a call to ILAENV in the future. Up to now, we always try
     !> iterative refinement.
     !> The iterative refinement process is stopped if
     !> ITER > ITERMAX
     !> or for all the RHS we have:
     !> RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
     !> where
     !> o ITER is the number of the current iteration in the iterative
     !> refinement process
     !> o RNRM is the infinity-norm of the residual
     !> o XNRM is the infinity-norm of the solution
     !> o ANRM is the infinity-operator-norm of the matrix A
     !> o EPS is the machine epsilon returned by QLAMCH('Epsilon')
     !> The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
     !> respectively.

     subroutine la_wzposv(uplo,n,nrhs,a,lda,b,ldb,x,ldx,work,swork,rwork,iter, &
               info)
        use la_constants_qp,only:cone,cnegone
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,iter
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(dp),intent(out) :: swork(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: b(ldb,*)
           complex(qp),intent(out) :: work(n,*),x(ldx,*)
        ! =====================================================================
           ! Parameters
           logical(lk),parameter :: doitref = .true.
           integer(ilp),parameter :: itermax = 30
           real(qp),parameter :: bwdmax = 1.0e+00_qp

           ! Local Scalars
           integer(ilp) :: i,iiter,ptsa,ptsx
           real(qp) :: anrm,cte,eps,rnrm,xnrm
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,max,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           info = 0
           iter = 0
           ! test the input parameters.
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           else if (ldx < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WZPOSV',-info)
              return
           end if
           ! quick return if (n==0).
           if (n == 0) return
           ! skip double precision iterative refinement if a priori slower
           ! than quad precision factorization.
           if (.not. doitref) then
              iter = -1
              go to 40
           end if
           ! compute some constants.
           anrm = la_wlanhe('I',uplo,n,a,lda,rwork)
           eps = la_qlamch('EPSILON')
           cte = anrm*eps*sqrt(real(n,KIND=qp))*bwdmax
           ! set the indices ptsa, ptsx for referencing sa and sx in swork.
           ptsa = 1
           ptsx = ptsa + n*n
           ! convert b from quad precision to double precision and store the
           ! result in sx.
           call la_wlag2z(n,nrhs,b,ldb,swork(ptsx),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! convert a from quad precision to double precision and store the
           ! result in sa.
           call la_wlat2z(uplo,n,a,lda,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -2
              go to 40
           end if
           ! compute the cholesky factorization of sa.
           call la_zpotrf(uplo,n,swork(ptsa),n,info)
           if (info /= 0) then
              iter = -3
              go to 40
           end if
           ! solve the system sa*sx = sb.
           call la_zpotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
           ! convert sx back to complex*16
           call la_zlag2w(n,nrhs,swork(ptsx),n,x,ldx,info)
           ! compute r = b - ax (r is work).
           call la_wlacpy('ALL',n,nrhs,b,ldb,work,n)
           call la_whemm('LEFT',uplo,n,nrhs,cnegone,a,lda,x,ldx,cone,work,n)

           ! check whether the nrhs normwise backward errors satisfy the
           ! stopping criterion. if yes, set iter=0 and return.
           do i = 1,nrhs
              xnrm = cabs1(x(la_iwamax(n,x(1,i),1),i))
              rnrm = cabs1(work(la_iwamax(n,work(1,i),1),i))
              if (rnrm > xnrm*cte) go to 10
           end do
           ! if we are here, the nrhs normwise backward errors satisfy the
           ! stopping criterion. we are good to exit.
           iter = 0
           return
           10 continue
           loop_30: do iiter = 1,itermax
              ! convert r (in work) from quad precision to double precision
              ! and store the result in sx.
              call la_wlag2z(n,nrhs,work,n,swork(ptsx),n,info)
              if (info /= 0) then
                 iter = -2
                 go to 40
              end if
              ! solve the system sa*sx = sr.
              call la_zpotrs(uplo,n,nrhs,swork(ptsa),n,swork(ptsx),n,info)
              ! convert sx back to quad precision and update the current
              ! iterate.
              call la_zlag2w(n,nrhs,swork(ptsx),n,work,n,info)
              do i = 1,nrhs
                 call la_waxpy(n,cone,work(1,i),1,x(1,i),1)
              end do
              ! compute r = b - ax (r is work).
              call la_wlacpy('ALL',n,nrhs,b,ldb,work,n)
              call la_whemm('L',uplo,n,nrhs,cnegone,a,lda,x,ldx,cone,work,n)

              ! check whether the nrhs normwise backward errors satisfy the
              ! stopping criterion. if yes, set iter=iiter>0 and return.
              do i = 1,nrhs
                 xnrm = cabs1(x(la_iwamax(n,x(1,i),1),i))
                 rnrm = cabs1(work(la_iwamax(n,work(1,i),1),i))
                 if (rnrm > xnrm*cte) go to 20
              end do
              ! if we are here, the nrhs normwise backward errors satisfy the
              ! stopping criterion, we are good to exit.
              iter = iiter
              return
              20 continue
           end do loop_30
           ! if we are at this place of the code, this is because we have
           ! performed iter=itermax iterations and never satisfied the
           ! stopping criterion, set up the iter flag accordingly and follow
           ! up on quad precision routine.
           iter = -itermax - 1
           40 continue
           ! single-precision iterative refinement failed to converge to a
           ! satisfactory solution, so we resort to quad precision.
           call la_wpotrf(uplo,n,a,lda,info)
           if (info /= 0) return
           call la_wlacpy('ALL',n,nrhs,b,ldb,x,ldx)
           call la_wpotrs(uplo,n,nrhs,a,lda,x,ldx,info)
           return
     end subroutine la_wzposv

     !> CPBSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H * U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular band matrix, and L is a lower
     !> triangular band matrix, with the same number of superdiagonals or
     !> subdiagonals as A.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_cpbsv(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kd < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < kd + 1) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CPBSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h*u or a = l*l**h.
           call la_cpbtrf(uplo,n,kd,ab,ldab,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_cpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
           end if
           return
     end subroutine la_cpbsv
     !> ZPBSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H * U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular band matrix, and L is a lower
     !> triangular band matrix, with the same number of superdiagonals or
     !> subdiagonals as A.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_zpbsv(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kd < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < kd + 1) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZPBSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_zpbtrf(uplo,n,kd,ab,ldab,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
           end if
           return
     end subroutine la_zpbsv
     !> WPBSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H * U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular band matrix, and L is a lower
     !> triangular band matrix, with the same number of superdiagonals or
     !> subdiagonals as A.  The factored form of A is then used to solve the
     !> system of equations A * X = B.

     pure subroutine la_wpbsv(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(inout) :: ab(ldab,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (kd < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldab < kd + 1) then
              info = -6
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WPBSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_wpbtrf(uplo,n,kd,ab,ldab,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wpbtrs(uplo,n,kd,nrhs,ab,ldab,b,ldb,info)
           end if
           return
     end subroutine la_wpbsv

     !> CPBSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_cpbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x, &
               ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldafb,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(sp),intent(inout) :: s(*)
           complex(sp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ,upper
           integer(ilp) :: i,infequ,j,j1,j2
           real(sp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           upper = la_lsame(uplo,'U')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kd < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           else if (ldafb < kd + 1) then
              info = -9
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -10
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -11
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -13
                 else if (ldx < max(1,n)) then
                    info = -15
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CPBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_cpbequ(uplo,n,kd,ab,ldab,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_claqhb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h *u or a = l*l**h.
              if (upper) then
                 do j = 1,n
                    j1 = max(j - kd,1)
                    call la_ccopy(j - j1 + 1,ab(kd + 1 - j + j1,j),1,afb(kd + 1 - j + j1,j),1)

                 end do
              else
                 do j = 1,n
                    j2 = min(j + kd,n)
                    call la_ccopy(j2 - j + 1,ab(1,j),1,afb(1,j),1)
                 end do
              end if
              call la_cpbtrf(uplo,n,kd,afb,ldafb,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clanhb('1',uplo,n,kd,ab,ldab,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_cpbcon(uplo,n,kd,afb,ldafb,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cpbtrs(uplo,n,kd,nrhs,afb,ldafb,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_cpbrfs(uplo,n,kd,nrhs,ab,ldab,afb,ldafb,b,ldb,x,ldx,ferr,berr, &
                      work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_cpbsvx
     !> ZPBSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zpbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x, &
               ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldafb,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(dp),intent(inout) :: s(*)
           complex(dp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ,upper
           integer(ilp) :: i,infequ,j,j1,j2
           real(dp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           upper = la_lsame(uplo,'U')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kd < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           else if (ldafb < kd + 1) then
              info = -9
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -10
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -11
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -13
                 else if (ldx < max(1,n)) then
                    info = -15
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZPBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_zpbequ(uplo,n,kd,ab,ldab,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_zlaqhb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h *u or a = l*l**h.
              if (upper) then
                 do j = 1,n
                    j1 = max(j - kd,1)
                    call la_zcopy(j - j1 + 1,ab(kd + 1 - j + j1,j),1,afb(kd + 1 - j + j1,j),1)

                 end do
              else
                 do j = 1,n
                    j2 = min(j + kd,n)
                    call la_zcopy(j2 - j + 1,ab(1,j),1,afb(1,j),1)
                 end do
              end if
              call la_zpbtrf(uplo,n,kd,afb,ldafb,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlanhb('1',uplo,n,kd,ab,ldab,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zpbcon(uplo,n,kd,afb,ldafb,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zpbtrs(uplo,n,kd,nrhs,afb,ldafb,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_zpbrfs(uplo,n,kd,nrhs,ab,ldab,afb,ldafb,b,ldb,x,ldx,ferr,berr, &
                      work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zpbsvx
     !> WPBSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite band matrix and X
     !> and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wpbsvx(fact,uplo,n,kd,nrhs,ab,ldab,afb,ldafb,equed,s,b,ldb,x, &
               ldx,rcond,ferr,berr,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldafb,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(qp),intent(inout) :: s(*)
           complex(qp),intent(inout) :: ab(ldab,*),afb(ldafb,*),b(ldb,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ,upper
           integer(ilp) :: i,infequ,j,j1,j2
           real(qp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           upper = la_lsame(uplo,'U')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (kd < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           else if (ldafb < kd + 1) then
              info = -9
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -10
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -11
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -13
                 else if (ldx < max(1,n)) then
                    info = -15
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WPBSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_wpbequ(uplo,n,kd,ab,ldab,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_wlaqhb(uplo,n,kd,ab,ldab,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right-hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h *u or a = l*l**h.
              if (upper) then
                 do j = 1,n
                    j1 = max(j - kd,1)
                    call la_wcopy(j - j1 + 1,ab(kd + 1 - j + j1,j),1,afb(kd + 1 - j + j1,j),1)

                 end do
              else
                 do j = 1,n
                    j2 = min(j + kd,n)
                    call la_wcopy(j2 - j + 1,ab(1,j),1,afb(1,j),1)
                 end do
              end if
              call la_wpbtrf(uplo,n,kd,afb,ldafb,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlanhb('1',uplo,n,kd,ab,ldab,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_wpbcon(uplo,n,kd,afb,ldafb,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wpbtrs(uplo,n,kd,nrhs,afb,ldafb,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_wpbrfs(uplo,n,kd,nrhs,ab,ldab,afb,ldafb,b,ldb,x,ldx,ferr,berr, &
                      work,rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_wpbsvx

     !> CPOSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H* U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular matrix and  L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_cposv(uplo,n,nrhs,a,lda,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CPOSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h*u or a = l*l**h.
           call la_cpotrf(uplo,n,a,lda,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_cpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
           end if
           return
     end subroutine la_cposv
     !> ZPOSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H* U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular matrix and  L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_zposv(uplo,n,nrhs,a,lda,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZPOSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_zpotrf(uplo,n,a,lda,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
           end if
           return
     end subroutine la_zposv
     !> WPOSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> The Cholesky decomposition is used to factor A as
     !> A = U**H* U,  if UPLO = 'U', or
     !> A = L * L**H,  if UPLO = 'L',
     !> where U is an upper triangular matrix and  L is a lower triangular
     !> matrix.  The factored form of A is then used to solve the system of
     !> equations A * X = B.

     pure subroutine la_wposv(uplo,n,nrhs,a,lda,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WPOSV ',-info)
              return
           end if
           ! compute the cholesky factorization a = u**h *u or a = l*l**h.
           call la_wpotrf(uplo,n,a,lda,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wpotrs(uplo,n,nrhs,a,lda,b,ldb,info)
           end if
           return
     end subroutine la_wposv

     !> CPOSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_cposvx(fact,uplo,n,nrhs,a,lda,af,ldaf,equed,s,b,ldb,x,ldx, &
               rcond,ferr,berr,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(sp),intent(inout) :: s(*)
           complex(sp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(sp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_slamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -9
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -10
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -12
                 else if (ldx < max(1,n)) then
                    info = -14
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CPOSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_cpoequ(n,a,lda,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_claqhe(uplo,n,a,lda,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h *u or a = l*l**h.
              call la_clacpy(uplo,n,n,a,lda,af,ldaf)
              call la_cpotrf(uplo,n,af,ldaf,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clanhe('1',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_cpocon(uplo,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cpotrs(uplo,n,nrhs,af,ldaf,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_cporfs(uplo,n,nrhs,a,lda,af,ldaf,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_cposvx
     !> ZPOSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zposvx(fact,uplo,n,nrhs,a,lda,af,ldaf,equed,s,b,ldb,x,ldx, &
               rcond,ferr,berr,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(dp),intent(inout) :: s(*)
           complex(dp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(dp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_dlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -9
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -10
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -12
                 else if (ldx < max(1,n)) then
                    info = -14
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZPOSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_zpoequ(n,a,lda,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_zlaqhe(uplo,n,a,lda,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h *u or a = l*l**h.
              call la_zlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_zpotrf(uplo,n,af,ldaf,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlanhe('1',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zpocon(uplo,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zpotrs(uplo,n,nrhs,af,ldaf,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_zporfs(uplo,n,nrhs,a,lda,af,ldaf,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zposvx
     !> WPOSVX: uses the Cholesky factorization A = U**H*U or A = L*L**H to
     !> compute the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian positive definite matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wposvx(fact,uplo,n,nrhs,a,lda,af,ldaf,equed,s,b,ldb,x,ldx, &
               rcond,ferr,berr,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(inout) :: equed
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(qp),intent(inout) :: s(*)
           complex(qp),intent(inout) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: equil,nofact,rcequ
           integer(ilp) :: i,infequ,j
           real(qp) :: amax,anorm,bignum,scond,smax,smin,smlnum
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           info = 0
           nofact = la_lsame(fact,'N')
           equil = la_lsame(fact,'E')
           if (nofact .or. equil) then
              equed = 'N'
              rcequ = .false.
           else
              rcequ = la_lsame(equed,'Y')
              smlnum = la_qlamch('SAFE MINIMUM')
              bignum = one/smlnum
           end if
           ! test the input parameters.
           if (.not. nofact .and. .not. equil .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           else if (ldaf < max(1,n)) then
              info = -8
           else if (la_lsame(fact,'F') .and. .not. (rcequ .or. la_lsame(equed,'N')) &
                      ) then
              info = -9
           else
              if (rcequ) then
                 smin = bignum
                 smax = zero
                 do j = 1,n
                    smin = min(smin,s(j))
                    smax = max(smax,s(j))
                 end do
                 if (smin <= zero) then
                    info = -10
                 else if (n > 0) then
                    scond = max(smin,smlnum)/min(smax,bignum)
                 else
                    scond = one
                 end if
              end if
              if (info == 0) then
                 if (ldb < max(1,n)) then
                    info = -12
                 else if (ldx < max(1,n)) then
                    info = -14
                 end if
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WPOSVX',-info)
              return
           end if
           if (equil) then
              ! compute row and column scalings to equilibrate the matrix a.
              call la_wpoequ(n,a,lda,s,scond,amax,infequ)
              if (infequ == 0) then
                 ! equilibrate the matrix.
                 call la_wlaqhe(uplo,n,a,lda,s,scond,amax,equed)
                 rcequ = la_lsame(equed,'Y')
              end if
           end if
           ! scale the right hand side.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = s(i)*b(i,j)
                 end do
              end do
           end if
           if (nofact .or. equil) then
              ! compute the cholesky factorization a = u**h *u or a = l*l**h.
              call la_wlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_wpotrf(uplo,n,af,ldaf,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlanhe('1',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_wpocon(uplo,n,af,ldaf,anorm,rcond,work,rwork,info)
           ! compute the solution matrix x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wpotrs(uplo,n,nrhs,af,ldaf,x,ldx,info)
           ! use iterative refinement to improve the computed solution and
           ! compute error bounds and backward error estimates for it.
           call la_wporfs(uplo,n,nrhs,a,lda,af,ldaf,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! transform the solution matrix x to a solution of the original
           ! system.
           if (rcequ) then
              do j = 1,nrhs
                 do i = 1,n
                    x(i,j) = s(i)*x(i,j)
                 end do
              end do
              do j = 1,nrhs
                 ferr(j) = ferr(j)/scond
              end do
           end if
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_wposvx

     !> CPTSV: computes the solution to a complex system of linear equations
     !> A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal
     !> matrix, and X and B are N-by-NRHS matrices.
     !> A is factored as A = L*D*L**H, and the factored form of A is then
     !> used to solve the system of equations.

     pure subroutine la_cptsv(n,nrhs,d,e,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(inout) :: d(*)
           complex(sp),intent(inout) :: b(ldb,*),e(*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CPTSV ',-info)
              return
           end if
           ! compute the l*d*l**h (or u**h*d*u) factorization of a.
           call la_cpttrf(n,d,e,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_cpttrs('LOWER',n,nrhs,d,e,b,ldb,info)
           end if
           return
     end subroutine la_cptsv
     !> ZPTSV: computes the solution to a complex system of linear equations
     !> A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal
     !> matrix, and X and B are N-by-NRHS matrices.
     !> A is factored as A = L*D*L**H, and the factored form of A is then
     !> used to solve the system of equations.

     pure subroutine la_zptsv(n,nrhs,d,e,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(inout) :: d(*)
           complex(dp),intent(inout) :: b(ldb,*),e(*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZPTSV ',-info)
              return
           end if
           ! compute the l*d*l**h (or u**h*d*u) factorization of a.
           call la_zpttrf(n,d,e,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zpttrs('LOWER',n,nrhs,d,e,b,ldb,info)
           end if
           return
     end subroutine la_zptsv
     !> WPTSV: computes the solution to a complex system of linear equations
     !> A*X = B, where A is an N-by-N Hermitian positive definite tridiagonal
     !> matrix, and X and B are N-by-NRHS matrices.
     !> A is factored as A = L*D*L**H, and the factored form of A is then
     !> used to solve the system of equations.

     pure subroutine la_wptsv(n,nrhs,d,e,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(inout) :: d(*)
           complex(qp),intent(inout) :: b(ldb,*),e(*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (nrhs < 0) then
              info = -2
           else if (ldb < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WPTSV ',-info)
              return
           end if
           ! compute the l*d*l**h (or u**h*d*u) factorization of a.
           call la_wpttrf(n,d,e,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wpttrs('LOWER',n,nrhs,d,e,b,ldb,info)
           end if
           return
     end subroutine la_wptsv

     !> CPTSVX: uses the factorization A = L*D*L**H to compute the solution
     !> to a complex system of linear equations A*X = B, where A is an
     !> N-by-N Hermitian positive definite tridiagonal matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_cptsvx(fact,n,nrhs,d,e,df,ef,b,ldb,x,ldx,rcond,ferr,berr, &
                work,rwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(sp),intent(in) :: d(*)
           real(sp),intent(inout) :: df(*)
           complex(sp),intent(in) :: b(ldb,*),e(*)
           complex(sp),intent(inout) :: ef(*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CPTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the l*d*l**h (or u**h*d*u) factorization of a.
              call la_scopy(n,d,1,df,1)
              if (n > 1) call la_ccopy(n - 1,e,1,ef,1)
              call la_cpttrf(n,df,ef,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clanht('1',n,d,e)
           ! compute the reciprocal of the condition number of a.
           call la_cptcon(n,df,ef,anorm,rcond,rwork,info)
           ! compute the solution vectors x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_cpttrs('LOWER',n,nrhs,df,ef,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_cptrfs('LOWER',n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_cptsvx
     !> ZPTSVX: uses the factorization A = L*D*L**H to compute the solution
     !> to a complex system of linear equations A*X = B, where A is an
     !> N-by-N Hermitian positive definite tridiagonal matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_zptsvx(fact,n,nrhs,d,e,df,ef,b,ldb,x,ldx,rcond,ferr,berr, &
                work,rwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(dp),intent(in) :: d(*)
           real(dp),intent(inout) :: df(*)
           complex(dp),intent(in) :: b(ldb,*),e(*)
           complex(dp),intent(inout) :: ef(*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZPTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the l*d*l**h (or u**h*d*u) factorization of a.
              call la_dcopy(n,d,1,df,1)
              if (n > 1) call la_zcopy(n - 1,e,1,ef,1)
              call la_zpttrf(n,df,ef,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlanht('1',n,d,e)
           ! compute the reciprocal of the condition number of a.
           call la_zptcon(n,df,ef,anorm,rcond,rwork,info)
           ! compute the solution vectors x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zpttrs('LOWER',n,nrhs,df,ef,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_zptrfs('LOWER',n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zptsvx
     !> WPTSVX: uses the factorization A = L*D*L**H to compute the solution
     !> to a complex system of linear equations A*X = B, where A is an
     !> N-by-N Hermitian positive definite tridiagonal matrix and X and B
     !> are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     pure subroutine la_wptsvx(fact,n,nrhs,d,e,df,ef,b,ldb,x,ldx,rcond,ferr,berr, &
                work,rwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           real(qp),intent(in) :: d(*)
           real(qp),intent(inout) :: df(*)
           complex(qp),intent(in) :: b(ldb,*),e(*)
           complex(qp),intent(inout) :: ef(*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WPTSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the l*d*l**h (or u**h*d*u) factorization of a.
              call la_qcopy(n,d,1,df,1)
              if (n > 1) call la_wcopy(n - 1,e,1,ef,1)
              call la_wpttrf(n,df,ef,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlanht('1',n,d,e)
           ! compute the reciprocal of the condition number of a.
           call la_wptcon(n,df,ef,anorm,rcond,rwork,info)
           ! compute the solution vectors x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wpttrs('LOWER',n,nrhs,df,ef,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_wptrfs('LOWER',n,nrhs,d,e,df,ef,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_wptsvx

end module la_lapack_solve_chol
