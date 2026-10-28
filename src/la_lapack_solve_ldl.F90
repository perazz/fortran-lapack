!> Symmetric and Hermitian indefinite drivers
module la_lapack_solve_ldl
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_mnorm
     use la_lapack_solve_ldl_comp
     use la_lapack_solve_ldl_comp2
     use la_lapack_solve_ldl_comp3
     use la_lapack_solve_ldl_comp4
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sspsv
     public :: la_sspsvx
     public :: la_ssysv_rk
     public :: la_ssysv_rook
     public :: la_ssysv
     public :: la_ssysvx
     public :: la_ssysv_aa
     public :: la_dspsv
     public :: la_dspsvx
     public :: la_dsysv_rk
     public :: la_dsysv_rook
     public :: la_dsysv
     public :: la_dsysvx
     public :: la_dsysv_aa
#ifdef LA_WITH_XDP
     public :: la_xspsv
     public :: la_xspsvx
     public :: la_xsysv_rk
     public :: la_xsysv_rook
     public :: la_xsysv
     public :: la_xsysvx
     public :: la_xsysv_aa
#endif
#ifdef LA_WITH_QP
     public :: la_qspsv
     public :: la_qspsvx
     public :: la_qsysv_rk
     public :: la_qsysv_rook
     public :: la_qsysv
     public :: la_qsysvx
     public :: la_qsysv_aa
#endif
     public :: la_cspsv
     public :: la_cspsvx
     public :: la_csysv
     public :: la_csysv_rk
     public :: la_csysv_rook
     public :: la_csysvx
     public :: la_chesv
     public :: la_chesv_rk
     public :: la_chesv_rook
     public :: la_chesvx
     public :: la_chpsv
     public :: la_chpsvx
     public :: la_chesv_aa
     public :: la_csysv_aa
     public :: la_zspsv
     public :: la_zspsvx
     public :: la_zsysv
     public :: la_zsysv_rk
     public :: la_zsysv_rook
     public :: la_zsysvx
     public :: la_zhesv
     public :: la_zhesv_rk
     public :: la_zhesv_rook
     public :: la_zhesvx
     public :: la_zhpsv
     public :: la_zhpsvx
     public :: la_zhesv_aa
     public :: la_zsysv_aa
#ifdef LA_WITH_XDP
     public :: la_yspsv
     public :: la_yspsvx
     public :: la_ysysv
     public :: la_ysysv_rk
     public :: la_ysysv_rook
     public :: la_ysysvx
     public :: la_yhesv
     public :: la_yhesv_rk
     public :: la_yhesv_rook
     public :: la_yhesvx
     public :: la_yhpsv
     public :: la_yhpsvx
     public :: la_yhesv_aa
     public :: la_ysysv_aa
#endif
#ifdef LA_WITH_QP
     public :: la_wspsv
     public :: la_wspsvx
     public :: la_wsysv
     public :: la_wsysv_rk
     public :: la_wsysv_rook
     public :: la_wsysvx
     public :: la_whesv
     public :: la_whesv_rk
     public :: la_whesv_rook
     public :: la_whesvx
     public :: la_whpsv
     public :: la_whpsvx
     public :: la_whesv_aa
     public :: la_wsysv_aa
#endif

     contains

     !> SSPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_sspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_ssptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_ssptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_sspsv
     !> DSPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_dspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_dsptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dsptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_dspsv
#ifdef LA_WITH_XDP
     !> XSPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_xspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: ap(*),b(ldb,*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('XSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_xsptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_xsptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_xspsv
#endif
#ifdef LA_WITH_QP
     !> QSPSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_qspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_qsptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qsptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_qspsv
#endif

     !> SSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_sspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,iwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: afp(*)
           real(sp),intent(in) :: ap(*),b(ldb,*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('SSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_scopy(n*(n + 1)/2,ap,1,afp,1)
              call la_ssptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_slansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_sspcon(uplo,n,afp,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ssptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_ssprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_sspsvx
     !> DSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,iwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: afp(*)
           real(dp),intent(in) :: ap(*),b(ldb,*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('DSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_dcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_dsptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_dlansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_dspcon(uplo,n,afp,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dsptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_dsprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_dspsvx
#ifdef LA_WITH_XDP
     !> XSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_xspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,iwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: afp(*)
           real(xdp),intent(in) :: ap(*),b(ldb,*)
           real(xdp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('XSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_xcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_xsptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_xlansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_xspcon(uplo,n,afp,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_xlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_xsptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_xsprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           return
     end subroutine la_xspsvx
#endif
#ifdef LA_WITH_QP
     !> QSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,iwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: afp(*)
           real(qp),intent(in) :: ap(*),b(ldb,*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('QSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_qcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_qsptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_qlansp('I',uplo,n,ap,work)
           ! compute the reciprocal of the condition number of a.
           call la_qspcon(uplo,n,afp,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qsptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_qsprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_qspsvx
#endif

     !> SSYSV_RK: computes the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> SSYTRF_RK is called to compute the factorization of a real
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine SSYTRS_3.

     pure subroutine la_ssysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_ssytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_ssytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_ssytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_ssysv_rk
     !> DSYSV_RK: computes the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> DSYTRF_RK is called to compute the factorization of a real
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine DSYTRS_3.

     pure subroutine la_dsysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_dsytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_dsytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_dsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_dsysv_rk
#ifdef LA_WITH_XDP
     !> XSYSV_RK: computes the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> XSYTRF_RK is called to compute the factorization of a real
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine XSYTRS_3.

     pure subroutine la_xsysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_xsytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_xsytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_xsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_xsysv_rk
#endif
#ifdef LA_WITH_QP
     !> QSYSV_RK: computes the solution to a real system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> QSYTRF_RK is called to compute the factorization of a real
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine QSYTRS_3.

     pure subroutine la_qsysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_qsytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_qsytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_qsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_qsysv_rk
#endif

     !> SSYSV_ROOK: computes the solution to a real system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> SSYTRF_ROOK is called to compute the factorization of a real
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling SSYTRS_ROOK.

     pure subroutine la_ssysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_ssytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_ssytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_ssytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_ssysv_rook
     !> DSYSV_ROOK: computes the solution to a real system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> DSYTRF_ROOK is called to compute the factorization of a real
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling DSYTRS_ROOK.

     pure subroutine la_dsysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_dsytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_dsytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_dsytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_dsysv_rook
#ifdef LA_WITH_XDP
     !> XSYSV_ROOK: computes the solution to a real system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> XSYTRF_ROOK is called to compute the factorization of a real
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling XSYTRS_ROOK.

     pure subroutine la_xsysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_xsytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_xsytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_xsytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_xsysv_rook
#endif
#ifdef LA_WITH_QP
     !> QSYSV_ROOK: computes the solution to a real system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> QSYTRF_ROOK is called to compute the factorization of a real
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling QSYTRS_ROOK.

     pure subroutine la_qsysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_qsytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_qsytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_qsytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_qsysv_rook
#endif

     !> SSYSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_ssysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_ssytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_ssytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_ssytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_ssytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_ssysv
     !> DSYSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_dsysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_dsytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_dsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_dsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_dsysv
#ifdef LA_WITH_XDP
     !> XSYSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_xsysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_xsytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_xsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_xsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_xsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_xsysv
#endif
#ifdef LA_WITH_QP
     !> QSYSV: computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_qsysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_qsytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = work(1)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_qsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_qsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_qsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_qsysv
#endif

     !> SSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a real system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_ssysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,iwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*),b(ldb,*)
           real(sp),intent(inout) :: af(ldaf,*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,3*n)
              if (nofact) then
                 nb = la_ilaenv(1,'SSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_slacpy(uplo,n,n,a,lda,af,ldaf)
              call la_ssytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_slansy('I',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_ssycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_slacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ssytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_ssyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_ssysvx
     !> DSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a real system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_dsysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,iwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*),b(ldb,*)
           real(dp),intent(inout) :: af(ldaf,*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,3*n)
              if (nofact) then
                 nb = la_ilaenv(1,'DSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_dlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_dsytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_dlansy('I',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_dsycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_dlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_dsytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_dsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_dsysvx
#ifdef LA_WITH_XDP
     !> XSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a real system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_xsysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,iwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: a(lda,*),b(ldb,*)
           real(xdp),intent(inout) :: af(ldaf,*)
           real(xdp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,3*n)
              if (nofact) then
                 nb = la_ilaenv(1,'XSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_xlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_xsytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_xlansy('I',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_xsycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_xlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_xsytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_xsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_xsysvx
#endif
#ifdef LA_WITH_QP
     !> QSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a real system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_qsysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,iwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*),b(ldb,*)
           real(qp),intent(inout) :: af(ldaf,*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,3*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,3*n)
              if (nofact) then
                 nb = la_ilaenv(1,'QSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_qlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_qsytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_qlansy('I',uplo,n,a,lda,work)
           ! compute the reciprocal of the condition number of a.
           call la_qsycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,iwork,info)
           ! compute the solution vectors x.
           call la_qlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_qsytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_qsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,iwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_qsysvx
#endif

     !> SSYSV computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_ssysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_ssytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_ssytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SSYSV_AA',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_ssytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_ssytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_ssysv_aa
     !> DSYSV computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_dsysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_dsytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_dsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_dsytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_dsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_dsysv_aa
#ifdef LA_WITH_XDP
     !> XSYSV computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_xsysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_xsytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_xsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_xsytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_xsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_xsysv_aa
#endif
#ifdef LA_WITH_QP
     !> QSYSV computes the solution to a real system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_qsysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_qsytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_qsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_qsytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_qsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_qsysv_aa
#endif

     !> CSPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_cspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_csptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_csptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_cspsv
     !> ZSPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_zspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_zsptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zsptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_zspsv
#ifdef LA_WITH_XDP
     !> YSPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_yspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: ap(*),b(ldb,*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('YSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_ysptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_ysptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_yspsv
#endif
#ifdef LA_WITH_QP
     !> WSPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is symmetric and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_wspsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WSPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_wsptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wsptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_wspsv
#endif

     !> CSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_cspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(inout) :: afp(*)
           complex(sp),intent(in) :: ap(*),b(ldb,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_ccopy(n*(n + 1)/2,ap,1,afp,1)
              call la_csptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clansp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_cspcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_csptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_csprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_cspsvx
     !> ZSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(inout) :: afp(*)
           complex(dp),intent(in) :: ap(*),b(ldb,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_zcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_zsptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlansp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zspcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zsptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_zsprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zspsvx
#ifdef LA_WITH_XDP
     !> YSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_yspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(xdp),intent(inout) :: afp(*)
           complex(xdp),intent(in) :: ap(*),b(ldb,*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('YSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_ycopy(n*(n + 1)/2,ap,1,afp,1)
              call la_ysptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_ylansp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_yspcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ysptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_ysprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           return
     end subroutine la_yspsvx
#endif
#ifdef LA_WITH_QP
     !> WSPSVX: uses the diagonal pivoting factorization A = U*D*U**T or
     !> A = L*D*L**T to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wspsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(inout) :: afp(*)
           complex(qp),intent(in) :: ap(*),b(ldb,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WSPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_wcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_wsptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlansp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_wspcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wsptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_wsprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_wspsvx
#endif

     !> CSYSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_csysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_csytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=sp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_csytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_csytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_csytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_csysv
     !> ZSYSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_zsysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_zsytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=dp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_zsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_zsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_zsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_zsysv
#ifdef LA_WITH_XDP
     !> YSYSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_ysysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_ysytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=xdp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_ysytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_ysytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_ysytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_ysysv
#endif
#ifdef LA_WITH_QP
     !> WSYSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_wsysv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_wsytrf(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=qp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WSYSV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_wsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_wsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_wsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_wsysv
#endif

     !> CSYSV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> CSYTRF_RK is called to compute the factorization of a complex
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine CSYTRS_3.

     pure subroutine la_csysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_csytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=sp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_csytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_csytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_csysv_rk
     !> ZSYSV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> ZSYTRF_RK is called to compute the factorization of a complex
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine ZSYTRS_3.

     pure subroutine la_zsysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_zsytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=dp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_zsytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_zsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_zsysv_rk
#ifdef LA_WITH_XDP
     !> YSYSV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> YSYTRF_RK is called to compute the factorization of a complex
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine YSYTRS_3.

     pure subroutine la_ysysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_ysytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=xdp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_ysytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_ysytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_ysysv_rk
#endif
#ifdef LA_WITH_QP
     !> WSYSV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N symmetric matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**T)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**T)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> WSYTRF_RK is called to compute the factorization of a complex
     !> symmetric matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine WSYTRS_3.

     pure subroutine la_wsysv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_wsytrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=qp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WSYSV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**t)*(p**t) or
           ! a = p*u*d*(u**t)*(p**t).
           call la_wsytrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_wsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_wsysv_rk
#endif

     !> CSYSV_ROOK: computes the solution to a complex system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> CSYTRF_ROOK is called to compute the factorization of a complex
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling CSYTRS_ROOK.

     pure subroutine la_csysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_csytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=sp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_csytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_csytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_csysv_rook
     !> ZSYSV_ROOK: computes the solution to a complex system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> ZSYTRF_ROOK is called to compute the factorization of a complex
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling ZSYTRS_ROOK.

     pure subroutine la_zsysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_zsytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=dp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_zsytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_zsytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_zsysv_rook
#ifdef LA_WITH_XDP
     !> YSYSV_ROOK: computes the solution to a complex system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> YSYTRF_ROOK is called to compute the factorization of a complex
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling YSYTRS_ROOK.

     pure subroutine la_ysysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_ysytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=xdp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_ysytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_ysytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_ysysv_rook
#endif
#ifdef LA_WITH_QP
     !> WSYSV_ROOK: computes the solution to a complex system of linear
     !> equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> WSYTRF_ROOK is called to compute the factorization of a complex
     !> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling WSYTRS_ROOK.

     pure subroutine la_wsysv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_wsytrf_rook(uplo,n,a,lda,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=qp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WSYSV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_wsytrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs_rook ( use level 2 blas)
              call la_wsytrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_wsysv_rook
#endif

     !> CSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_csysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: af(ldaf,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'CSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_clacpy(uplo,n,n,a,lda,af,ldaf)
              call la_csytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clansy('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_csycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_csytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_csyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_csysvx
     !> ZSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zsysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: af(ldaf,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'ZSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_zlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_zsytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlansy('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zsycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zsytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_zsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_zsysvx
#ifdef LA_WITH_XDP
     !> YSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_ysysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: af(ldaf,*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'YSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_ylacpy(uplo,n,n,a,lda,af,ldaf)
              call la_ysytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_ylansy('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_ysycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_ysytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_ysyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_ysysvx
#endif
#ifdef LA_WITH_QP
     !> WSYSVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_wsysvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: af(ldaf,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'WSYTRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WSYSVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**t or a = l*d*l**t.
              call la_wlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_wsytrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlansy('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_wsycon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_wsytrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_wsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_wsysvx
#endif

     !> CHESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_chesv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'CHETRF',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHESV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_chetrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_chetrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_chetrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_chesv
     !> ZHESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_zhesv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'ZHETRF',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHESV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_zhetrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_zhetrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_zhetrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_zhesv
#ifdef LA_WITH_XDP
     !> YHESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_yhesv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'YHETRF',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YHESV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_yhetrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_yhetrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_yhetrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_yhesv
#endif
#ifdef LA_WITH_QP
     !> WHESV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then
     !> used to solve the system of equations A * X = B.

     pure subroutine la_whesv(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'WHETRF',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHESV ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_whetrf(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              if (lwork < n) then
              ! solve with trs ( use level blas 2)
                 call la_whetrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
              else
              ! solve with trs2 ( use level blas 3)
                 call la_whetrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_whesv
#endif

     !> CHESV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**H)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**H)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> CHETRF_RK is called to compute the factorization of a complex
     !> Hermitian matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine CHETRS_3.

     pure subroutine la_chesv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_chetrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=sp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHESV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**t or a = l*d*l**t.
           call la_chetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_chetrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_chesv_rk
     !> ZHESV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**H)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**H)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> ZHETRF_RK is called to compute the factorization of a complex
     !> Hermitian matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine ZHETRS_3.

     pure subroutine la_zhesv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_zhetrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=dp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHESV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**h)*(p**t) or
           ! a = p*u*d*(u**h)*(p**t).
           call la_zhetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_zhetrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_zhesv_rk
#ifdef LA_WITH_XDP
     !> YHESV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**H)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**H)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> YHETRF_RK is called to compute the factorization of a complex
     !> Hermitian matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine YHETRS_3.

     pure subroutine la_yhesv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_yhetrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=xdp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YHESV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**h)*(p**t) or
           ! a = p*u*d*(u**h)*(p**t).
           call la_yhetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_yhetrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_yhesv_rk
#endif
#ifdef LA_WITH_QP
     !> WHESV_RK: computes the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix
     !> and X and B are N-by-NRHS matrices.
     !> The bounded Bunch-Kaufman (rook) diagonal pivoting method is used
     !> to factor A as
     !> A = P*U*D*(U**H)*(P**T),  if UPLO = 'U', or
     !> A = P*L*D*(L**H)*(P**T),  if UPLO = 'L',
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> WHETRF_RK is called to compute the factorization of a complex
     !> Hermitian matrix.  The factored form of A is then used to solve
     !> the system of equations A * X = B by calling BLAS3 routine WHETRS_3.

     pure subroutine la_whesv_rk(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -9
           else if (lwork < 1 .and. .not. lquery) then
              info = -11
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 call la_whetrf_rk(uplo,n,a,lda,e,ipiv,work,-1,info)
                 lwkopt = real(work(1),KIND=qp)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHESV_RK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = p*u*d*(u**h)*(p**t) or
           ! a = p*u*d*(u**h)*(p**t).
           call la_whetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b with blas3 solver, overwriting b with x.
              call la_whetrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_whesv_rk
#endif

     !> CHESV_ROOK: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The bounded Bunch-Kaufman ("rook") diagonal pivoting method is used
     !> to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> CHETRF_ROOK is called to compute the factorization of a complex
     !> Hermition matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling CHETRS_ROOK (uses BLAS 2).

     pure subroutine la_chesv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'CHETRF_ROOK',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHESV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_chetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs ( use level blas 2)
              call la_chetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_chesv_rook
     !> ZHESV_ROOK: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The bounded Bunch-Kaufman ("rook") diagonal pivoting method is used
     !> to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> ZHETRF_ROOK is called to compute the factorization of a complex
     !> Hermition matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling ZHETRS_ROOK (uses BLAS 2).

     pure subroutine la_zhesv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'ZHETRF_ROOK',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHESV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_zhetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs ( use level blas 2)
              call la_zhetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_zhesv_rook
#ifdef LA_WITH_XDP
     !> YHESV_ROOK: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The bounded Bunch-Kaufman ("rook") diagonal pivoting method is used
     !> to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> YHETRF_ROOK is called to compute the factorization of a complex
     !> Hermition matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling YHETRS_ROOK (uses BLAS 2).

     pure subroutine la_yhesv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'YHETRF_ROOK',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YHESV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_yhetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs ( use level blas 2)
              call la_yhetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_yhesv_rook
#endif
#ifdef LA_WITH_QP
     !> WHESV_ROOK: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> The bounded Bunch-Kaufman ("rook") diagonal pivoting method is used
     !> to factor A as
     !> A = U * D * U**T,  if UPLO = 'U', or
     !> A = L * D * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> WHETRF_ROOK is called to compute the factorization of a complex
     !> Hermition matrix A using the bounded Bunch-Kaufman ("rook") diagonal
     !> pivoting method.
     !> The factored form of A is then used to solve the system
     !> of equations A * X = B by calling WHETRS_ROOK (uses BLAS 2).

     pure subroutine la_whesv_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,nb
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < 1 .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              if (n == 0) then
                 lwkopt = 1
              else
                 nb = la_ilaenv(1,'WHETRF_ROOK',uplo,n,-1,-1,-1)
                 lwkopt = n*nb
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHESV_ROOK ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_whetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              ! solve with trs ( use level blas 2)
              call la_whetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
           end if
           work(1) = lwkopt
           return
     end subroutine la_whesv_rook
#endif

     !> CHESVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_chesvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: af(ldaf,*)
           complex(sp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(sp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'CHETRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHESVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_clacpy(uplo,n,n,a,lda,af,ldaf)
              call la_chetrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clanhe('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_checon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_chetrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_cherfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_chesvx
     !> ZHESVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zhesvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: af(ldaf,*)
           complex(dp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(dp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'ZHETRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHESVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_zlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_zhetrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlanhe('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zhecon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zhetrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_zherfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_zhesvx
#ifdef LA_WITH_XDP
     !> YHESVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_yhesvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: af(ldaf,*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'YHETRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YHESVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_ylacpy(uplo,n,n,a,lda,af,ldaf)
              call la_yhetrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_ylanhe('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_yhecon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_yhetrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_yherfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_yhesvx
#endif
#ifdef LA_WITH_QP
     !> WHESVX: uses the diagonal pivoting factorization to compute the
     !> solution to a complex system of linear equations A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_whesvx(fact,uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,rcond, &
               ferr,berr,work,lwork,rwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,lwork,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: af(ldaf,*)
           complex(qp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,nofact
           integer(ilp) :: lwkopt,nb
           real(qp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           lquery = (lwork == -1)
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
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
           else if (ldb < max(1,n)) then
              info = -11
           else if (ldx < max(1,n)) then
              info = -13
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -18
           end if
           if (info == 0) then
              lwkopt = max(1,2*n)
              if (nofact) then
                 nb = la_ilaenv(1,'WHETRF',uplo,n,-1,-1,-1)
                 lwkopt = max(lwkopt,n*nb)
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHESVX',-info)
              return
           else if (lquery) then
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_wlacpy(uplo,n,n,a,lda,af,ldaf)
              call la_whetrf(uplo,n,af,ldaf,ipiv,work,lwork,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlanhe('I',uplo,n,a,lda,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_whecon(uplo,n,af,ldaf,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_whetrs(uplo,n,nrhs,af,ldaf,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_wherfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr,berr, &
                     work,rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           work(1) = lwkopt
           return
     end subroutine la_whesvx
#endif

     !> CHPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is Hermitian and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_chpsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CHPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_chptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_chptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_chpsv
     !> ZHPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is Hermitian and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_zhpsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZHPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_zhptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zhptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_zhpsv
#ifdef LA_WITH_XDP
     !> YHPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is Hermitian and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_yhpsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: ap(*),b(ldb,*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('YHPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_yhptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_yhptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_yhpsv
#endif
#ifdef LA_WITH_QP
     !> WHPSV: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix stored in packed format and X
     !> and B are N-by-NRHS matrices.
     !> The diagonal pivoting method is used to factor A as
     !> A = U * D * U**H,  if UPLO = 'U', or
     !> A = L * D * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, D is Hermitian and block diagonal with 1-by-1
     !> and 2-by-2 diagonal blocks.  The factored form of A is then used to
     !> solve the system of equations A * X = B.

     pure subroutine la_whpsv(uplo,n,nrhs,ap,ipiv,b,ldb,info)
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
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
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WHPSV ',-info)
              return
           end if
           ! compute the factorization a = u*d*u**h or a = l*d*l**h.
           call la_whptrf(uplo,n,ap,ipiv,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_whptrs(uplo,n,nrhs,ap,ipiv,b,ldb,info)
           end if
           return
     end subroutine la_whpsv
#endif

     !> CHPSVX: uses the diagonal pivoting factorization A = U*D*U**H or
     !> A = L*D*L**H to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_chpsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_sp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(inout) :: afp(*)
           complex(sp),intent(in) :: ap(*),b(ldb,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CHPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_ccopy(n*(n + 1)/2,ap,1,afp,1)
              call la_chptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_clanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_chpcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_clacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_chptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_chprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_slamch('EPSILON')) info = n + 1
           return
     end subroutine la_chpsvx
     !> ZHPSVX: uses the diagonal pivoting factorization A = U*D*U**H or
     !> A = L*D*L**H to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_zhpsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_dp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(inout) :: afp(*)
           complex(dp),intent(in) :: ap(*),b(ldb,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZHPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_zcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_zhptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_zlanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_zhpcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_zlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_zhptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_zhprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_dlamch('EPSILON')) info = n + 1
           return
     end subroutine la_zhpsvx
#ifdef LA_WITH_XDP
     !> YHPSVX: uses the diagonal pivoting factorization A = U*D*U**H or
     !> A = L*D*L**H to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_yhpsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_xdp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(xdp),intent(inout) :: afp(*)
           complex(xdp),intent(in) :: ap(*),b(ldb,*)
           complex(xdp),intent(out) :: work(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nofact
           real(xdp) :: anorm
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nofact = la_lsame(fact,'N')
           if (.not. nofact .and. .not. la_lsame(fact,'F')) then
              info = -1
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('YHPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_ycopy(n*(n + 1)/2,ap,1,afp,1)
              call la_yhptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_ylanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_yhpcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_ylacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_yhptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_yhprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_xlamch('EPSILON')) info = n + 1
           return
     end subroutine la_yhpsvx
#endif
#ifdef LA_WITH_QP
     !> WHPSVX: uses the diagonal pivoting factorization A = U*D*U**H or
     !> A = L*D*L**H to compute the solution to a complex system of linear
     !> equations A * X = B, where A is an N-by-N Hermitian matrix stored
     !> in packed format and X and B are N-by-NRHS matrices.
     !> Error bounds on the solution and a condition estimate are also
     !> provided.

     subroutine la_whpsvx(fact,uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,rcond,ferr, &
               berr,work,rwork,info)
        use la_constants_qp
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: fact,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(inout) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(inout) :: afp(*)
           complex(qp),intent(in) :: ap(*),b(ldb,*)
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
           else if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) &
                     then
              info = -2
           else if (n < 0) then
              info = -3
           else if (nrhs < 0) then
              info = -4
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WHPSVX',-info)
              return
           end if
           if (nofact) then
              ! compute the factorization a = u*d*u**h or a = l*d*l**h.
              call la_wcopy(n*(n + 1)/2,ap,1,afp,1)
              call la_whptrf(uplo,n,afp,ipiv,info)
              ! return if info is non-zero.
              if (info > 0) then
                 rcond = zero
                 return
              end if
           end if
           ! compute the norm of the matrix a.
           anorm = la_wlanhp('I',uplo,n,ap,rwork)
           ! compute the reciprocal of the condition number of a.
           call la_whpcon(uplo,n,afp,ipiv,anorm,rcond,work,info)
           ! compute the solution vectors x.
           call la_wlacpy('FULL',n,nrhs,b,ldb,x,ldx)
           call la_whptrs(uplo,n,nrhs,afp,ipiv,x,ldx,info)
           ! use iterative refinement to improve the computed solutions and
           ! compute error bounds and backward error estimates for them.
           call la_whprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                     rwork,info)
           ! set info = n+1 if the matrix is singular to working precision.
           if (rcond < la_qlamch('EPSILON')) info = n + 1
           return
     end subroutine la_whpsvx
#endif

     !> CHESV_AA: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**H * T * U,  if UPLO = 'U', or
     !> A = L * T * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is Hermitian and tridiagonal. The factored form
     !> of A is then used to solve the system of equations A * X = B.

     pure subroutine la_chesv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_hetrf,lwkopt_hetrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_chetrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_hetrf = int(work(1),KIND=ilp)
              call la_chetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_hetrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_hetrf,lwkopt_hetrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHESV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**h*t*u or a = l*t*l**h.
           call la_chetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_chetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_chesv_aa
     !> ZHESV_AA: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**H * T * U,  if UPLO = 'U', or
     !> A = L * T * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is Hermitian and tridiagonal. The factored form
     !> of A is then used to solve the system of equations A * X = B.

     pure subroutine la_zhesv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_hetrf,lwkopt_hetrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_zhetrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_hetrf = int(work(1),KIND=ilp)
              call la_zhetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_hetrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_hetrf,lwkopt_hetrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHESV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**h*t*u or a = l*t*l**h.
           call la_zhetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zhetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_zhesv_aa
#ifdef LA_WITH_XDP
     !> YHESV_AA: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**H * T * U,  if UPLO = 'U', or
     !> A = L * T * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is Hermitian and tridiagonal. The factored form
     !> of A is then used to solve the system of equations A * X = B.

     pure subroutine la_yhesv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_hetrf,lwkopt_hetrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_yhetrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_hetrf = int(work(1),KIND=ilp)
              call la_yhetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_hetrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_hetrf,lwkopt_hetrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YHESV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**h*t*u or a = l*t*l**h.
           call la_yhetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_yhetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_yhesv_aa
#endif
#ifdef LA_WITH_QP
     !> WHESV_AA: computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**H * T * U,  if UPLO = 'U', or
     !> A = L * T * L**H,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is Hermitian and tridiagonal. The factored form
     !> of A is then used to solve the system of equations A * X = B.

     pure subroutine la_whesv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_hetrf,lwkopt_hetrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_whetrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_hetrf = int(work(1),KIND=ilp)
              call la_whetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_hetrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_hetrf,lwkopt_hetrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHESV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**h*t*u or a = l*t*l**h.
           call la_whetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_whetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_whesv_aa
#endif

     !> CSYSV computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_csysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_csytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_csytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_csytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_csytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_csysv_aa
     !> ZSYSV computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_zsysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_zsytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_zsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_zsytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_zsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_zsysv_aa
#ifdef LA_WITH_XDP
     !> YSYSV computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_ysysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_ysytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_ysytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_ysytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_ysytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_ysysv_aa
#endif
#ifdef LA_WITH_QP
     !> WSYSV computes the solution to a complex system of linear equations
     !> A * X = B,
     !> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
     !> matrices.
     !> Aasen's algorithm is used to factor A as
     !> A = U**T * T * U,  if UPLO = 'U', or
     !> A = L * T * L**T,  if UPLO = 'L',
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is symmetric tridiagonal. The factored
     !> form of A is then used to solve the system of equations A * X = B.

     pure subroutine la_wsysv_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,lwork,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: lwkopt,lwkopt_sytrf,lwkopt_sytrs
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1)
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(2*n,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info == 0) then
              call la_wsytrf_aa(uplo,n,a,lda,ipiv,work,-1,info)
              lwkopt_sytrf = int(work(1),KIND=ilp)
              call la_wsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,-1,info)
              lwkopt_sytrs = int(work(1),KIND=ilp)
              lwkopt = max(lwkopt_sytrf,lwkopt_sytrs)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WSYSV_AA ',-info)
              return
           else if (lquery) then
              return
           end if
           ! compute the factorization a = u**t*t*u or a = l*t*l**t.
           call la_wsytrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
           if (info == 0) then
              ! solve the system a*x = b, overwriting b with x.
              call la_wsytrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)

           end if
           work(1) = lwkopt
           return
     end subroutine la_wsysv_aa
#endif

end module la_lapack_solve_ldl
