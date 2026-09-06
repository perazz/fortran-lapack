!> BLAS-like level 3: rank-k updates and solves in RFP storage
module la_lapack_blas_like_l3
     use la_constants
     use la_blas_aux
     use la_blas_level3_gen
     use la_blas_level3_sym
     use la_blas_level3_tri
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slagtm
     public :: la_ssfrk
     public :: la_stfsm
     public :: la_dlagtm
     public :: la_dsfrk
     public :: la_dtfsm
#ifdef LA_WITH_XDP
     public :: la_xlagtm
     public :: la_xsfrk
     public :: la_xtfsm
#endif
#ifdef LA_WITH_QP
     public :: la_qlagtm
     public :: la_qsfrk
     public :: la_qtfsm
#endif
     public :: la_chfrk
     public :: la_clacrm
     public :: la_clagtm
     public :: la_clarcm
     public :: la_ctfsm
     public :: la_zhfrk
     public :: la_zlacrm
     public :: la_zlagtm
     public :: la_zlarcm
     public :: la_ztfsm
#ifdef LA_WITH_XDP
     public :: la_yhfrk
     public :: la_ylacrm
     public :: la_ylagtm
     public :: la_ylarcm
     public :: la_ytfsm
#endif
#ifdef LA_WITH_QP
     public :: la_whfrk
     public :: la_wlacrm
     public :: la_wlagtm
     public :: la_wlarcm
     public :: la_wtfsm
#endif

     contains

     !> SLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_slagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_sp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(in) :: alpha,beta
           ! Array Arguments
           real(sp),intent(inout) :: b(ldb,*)
           real(sp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b + a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b - a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_slagtm
     !> DLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_dlagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_dp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(in) :: alpha,beta
           ! Array Arguments
           real(dp),intent(inout) :: b(ldb,*)
           real(dp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b + a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b - a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_dlagtm
#ifdef LA_WITH_XDP
     !> XLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_xlagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_xdp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(in) :: alpha,beta
           ! Array Arguments
           real(xdp),intent(inout) :: b(ldb,*)
           real(xdp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b + a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b - a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_xlagtm
#endif
#ifdef LA_WITH_QP
     !> QLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_qlagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_qp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(in) :: alpha,beta
           ! Array Arguments
           real(qp),intent(inout) :: b(ldb,*)
           real(qp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b + a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else
                 ! compute b := b - a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_qlagtm
#endif

     !> Level 3 BLAS like routine for C in RFP Format.
     !> SSFRK: performs one of the symmetric rank--k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n symmetric
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_ssfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('SSFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_ssyrk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = zero
              end do
              return
           end if
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_ssyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_ssyrk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_sgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 't'
                       call la_ssyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_ssyrk('U','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_sgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_ssyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_ssyrk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_sgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n2,1), &
                                 lda,beta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 't'
                       call la_ssyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_ssyrk('U','T',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_sgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n2), &
                                 lda,beta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'l', and trans = 'n'
                       call la_ssyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_ssyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_sgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n1 + 1,1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 't', uplo = 'l', and trans = 't'
                       call la_ssyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_ssyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_sgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n1 + 1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'u', and trans = 'n'
                       call la_ssyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_ssyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_sgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    else
                       ! n is odd, transr = 't', uplo = 'u', and trans = 't'
                       call la_ssyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_ssyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_sgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_ssyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_ssyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_sgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 't'
                       call la_ssyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_ssyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_sgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_ssyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_ssyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_sgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 't'
                       call la_ssyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_ssyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_sgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'l', and trans = 'n'
                       call la_ssyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_ssyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_sgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'l', and trans = 't'
                       call la_ssyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_ssyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_sgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'u', and trans = 'n'
                       call la_ssyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_ssyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_sgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'u', and trans = 't'
                       call la_ssyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_ssyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_sgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_ssfrk
     !> Level 3 BLAS like routine for C in RFP Format.
     !> DSFRK: performs one of the symmetric rank--k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n symmetric
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_dsfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DSFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_dsyrk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = zero
              end do
              return
           end if
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_dsyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_dsyrk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_dgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 't'
                       call la_dsyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_dsyrk('U','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_dgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_dsyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_dsyrk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_dgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n2,1), &
                                 lda,beta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 't'
                       call la_dsyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_dsyrk('U','T',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_dgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n2), &
                                 lda,beta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'l', and trans = 'n'
                       call la_dsyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_dsyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_dgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n1 + 1,1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 't', uplo = 'l', and trans = 't'
                       call la_dsyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_dsyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_dgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n1 + 1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'u', and trans = 'n'
                       call la_dsyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_dsyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_dgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    else
                       ! n is odd, transr = 't', uplo = 'u', and trans = 't'
                       call la_dsyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_dsyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_dgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_dsyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_dsyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_dgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 't'
                       call la_dsyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_dsyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_dgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_dsyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_dsyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_dgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 't'
                       call la_dsyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_dsyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_dgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'l', and trans = 'n'
                       call la_dsyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_dsyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_dgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'l', and trans = 't'
                       call la_dsyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_dsyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_dgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'u', and trans = 'n'
                       call la_dsyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_dsyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_dgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'u', and trans = 't'
                       call la_dsyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_dsyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_dgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_dsfrk
#ifdef LA_WITH_XDP
     !> Level 3 BLAS like routine for C in RFP Format.
     !> XSFRK: performs one of the symmetric rank--k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n symmetric
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_xsfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('XSFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_xsyrk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = zero
              end do
              return
           end if
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_xsyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_xsyrk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_xgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 't'
                       call la_xsyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_xsyrk('U','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_xgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_xsyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_xsyrk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_xgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n2,1), &
                                 lda,beta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 't'
                       call la_xsyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_xsyrk('U','T',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_xgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n2), &
                                 lda,beta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'l', and trans = 'n'
                       call la_xsyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_xsyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_xgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n1 + 1,1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 't', uplo = 'l', and trans = 't'
                       call la_xsyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_xsyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_xgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n1 + 1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'u', and trans = 'n'
                       call la_xsyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_xsyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_xgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    else
                       ! n is odd, transr = 't', uplo = 'u', and trans = 't'
                       call la_xsyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_xsyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_xgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_xsyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_xsyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_xgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 't'
                       call la_xsyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_xsyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_xgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_xsyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_xsyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_xgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 't'
                       call la_xsyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_xsyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_xgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'l', and trans = 'n'
                       call la_xsyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_xsyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_xgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'l', and trans = 't'
                       call la_xsyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_xsyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_xgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'u', and trans = 'n'
                       call la_xsyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_xsyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_xgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'u', and trans = 't'
                       call la_xsyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_xsyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_xgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_xsfrk
#endif
#ifdef LA_WITH_QP
     !> Level 3 BLAS like routine for C in RFP Format.
     !> QSFRK: performs one of the symmetric rank--k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n symmetric
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_qsfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QSFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_qsyrk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = zero
              end do
              return
           end if
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_qsyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_qsyrk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_qgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 't'
                       call la_qsyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_qsyrk('U','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_qgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_qsyrk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_qsyrk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_qgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n2,1), &
                                 lda,beta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 't'
                       call la_qsyrk('L','T',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_qsyrk('U','T',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_qgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n2), &
                                 lda,beta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 't'
                 if (lower) then
                    ! n is odd, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'l', and trans = 'n'
                       call la_qsyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_qsyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_qgemm('N','T',n1,n2,k,alpha,a(1,1),lda,a(n1 + 1,1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 't', uplo = 'l', and trans = 't'
                       call la_qsyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_qsyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_qgemm('T','N',n1,n2,k,alpha,a(1,1),lda,a(1,n1 + 1), &
                                  lda,beta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 't', uplo = 'u', and trans = 'n'
                       call la_qsyrk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_qsyrk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_qgemm('N','T',n2,n1,k,alpha,a(n1 + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    else
                       ! n is odd, transr = 't', uplo = 'u', and trans = 't'
                       call la_qsyrk('U','T',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_qsyrk('L','T',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_qgemm('T','N',n2,n1,k,alpha,a(1,n1 + 1),lda,a(1,1), &
                                  lda,beta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_qsyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_qsyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_qgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 't'
                       call la_qsyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_qsyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_qgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_qsyrk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_qsyrk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_qgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 't'
                       call la_qsyrk('L','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_qsyrk('U','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_qgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 't'
                 if (lower) then
                    ! n is even, transr = 't', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'l', and trans = 'n'
                       call la_qsyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_qsyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_qgemm('N','T',nk,nk,k,alpha,a(1,1),lda,a(nk + 1,1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'l', and trans = 't'
                       call la_qsyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_qsyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_qgemm('T','N',nk,nk,k,alpha,a(1,1),lda,a(1,nk + 1), &
                                  lda,beta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 't', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 't', uplo = 'u', and trans = 'n'
                       call la_qsyrk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_qsyrk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_qgemm('N','T',nk,nk,k,alpha,a(nk + 1,1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    else
                       ! n is even, transr = 't', uplo = 'u', and trans = 't'
                       call la_qsyrk('U','T',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_qsyrk('L','T',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_qgemm('T','N',nk,nk,k,alpha,a(1,nk + 1),lda,a(1,1), &
                                  lda,beta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_qsfrk
#endif

     !> Level 3 BLAS like routine for A in RFP Format.
     !> STFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**T.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_stfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           real(sp),intent(in) :: alpha
           ! Array Arguments
           real(sp),intent(in) :: a(0:*)
           real(sp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('STFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_sp)
           if (alpha == zero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = zero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_strsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_strsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_sgemm('N','N',m2,n,m1,-one,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_strsm('L','U','T',diag,m2,n,one,a(m),m,b(m1, &
                                       0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_strsm('L','L','T',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_strsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_sgemm('T','N',m1,n,m2,-one,a(m1),m,b(m1,0), &
                                       ldb,alpha,b,ldb)
                             call la_strsm('L','L','T',diag,m1,n,one,a(0),m,b,ldb &
                                       )
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_strsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_sgemm('T','N',m2,n,m1,-one,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_strsm('L','U','T',diag,m2,n,one,a(m1),m,b(m1,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_strsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_sgemm('N','N',m1,n,m2,-one,a(0),m,b(m1,0),ldb, &
                                     alpha,b,ldb)
                          call la_strsm('L','L','T',diag,m1,n,one,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_strsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_strsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_sgemm('T','N',m2,n,m1,-one,a(m1*m1),m1,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_strsm('L','L','N',diag,m2,n,one,a(1),m1,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_strsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_strsm('L','L','T',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_sgemm('N','N',m1,n,m2,-one,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_strsm('L','U','N',diag,m1,n,one,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_strsm('L','U','T',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_sgemm('N','N',m2,n,m1,-one,a(0),m2,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_strsm('L','L','N',diag,m2,n,one,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_strsm('L','L','T',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_sgemm('T','N',m1,n,m2,-one,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_strsm('L','U','N',diag,m1,n,one,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_strsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_sgemm('N','N',k,n,k,-one,a(k + 1),m + 1,b,ldb,alpha, &
                                     b(k,0),ldb)
                          call la_strsm('L','U','T',diag,k,n,one,a(0),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_strsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_sgemm('T','N',k,n,k,-one,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_strsm('L','L','T',diag,k,n,one,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_strsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_sgemm('T','N',k,n,k,-one,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_strsm('L','U','T',diag,k,n,one,a(k),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_strsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_sgemm('N','N',k,n,k,-one,a(0),m + 1,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_strsm('L','L','T',diag,k,n,one,a(k + 1),m + 1,b,ldb &
                                    )
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_strsm('L','U','T',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_sgemm('T','N',k,n,k,-one,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_strsm('L','L','N',diag,k,n,one,a(0),k,b(k,0), &
                                    ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_strsm('L','L','T',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_sgemm('N','N',k,n,k,-one,a(k*(k + 1)),k,b(k,0), &
                                     ldb,alpha,b,ldb)
                          call la_strsm('L','U','N',diag,k,n,one,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_strsm('L','U','T',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_sgemm('N','N',k,n,k,-one,a(0),k,b,ldb,alpha,b( &
                                    k,0),ldb)
                          call la_strsm('L','L','N',diag,k,n,one,a(k*k),k,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_strsm('L','L','T',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_sgemm('T','N',k,n,k,-one,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_strsm('L','U','N',diag,k,n,one,a(k*(k + 1)),k,b, &
                                    ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_strsm('R','U','T',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_sgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_strsm('R','L','N',diag,m,n1,one,a(0),n,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          call la_strsm('R','L','T',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_sgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(n1),n, &
                                     alpha,b(0,n1),ldb)
                          call la_strsm('R','U','N',diag,m,n2,one,a(n),n,b(0,n1) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_strsm('R','L','T',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_sgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(0),n, &
                                    alpha,b(0,n1),ldb)
                          call la_strsm('R','U','N',diag,m,n2,one,a(n1),n,b(0,n1 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_strsm('R','U','T',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_sgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(0),n, &
                                     alpha,b(0,0),ldb)
                          call la_strsm('R','L','N',diag,m,n1,one,a(n2),n,b(0,0) &
                                    ,ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          call la_strsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_sgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,0),ldb)
                          call la_strsm('R','U','T',diag,m,n1,one,a(0),n1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          call la_strsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_sgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(n1*n1), &
                                     n1,alpha,b(0,n1),ldb)
                          call la_strsm('R','L','T',diag,m,n2,one,a(1),n1,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_strsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_sgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(0),n2, &
                                     alpha,b(0,n1),ldb)
                          call la_strsm('R','L','T',diag,m,n2,one,a(n1*n2),n2,b(0, &
                                     n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_strsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_sgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_strsm('R','U','T',diag,m,n1,one,a(n2*n2),n2,b(0, &
                                     0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_strsm('R','U','T',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_sgemm('N','N',m,k,k,-one,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_strsm('R','L','N',diag,m,k,one,a(1),n + 1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_strsm('R','L','T',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_sgemm('N','T',m,k,k,-one,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_strsm('R','U','N',diag,m,k,one,a(0),n + 1,b(0,k) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_strsm('R','L','T',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_sgemm('N','N',m,k,k,-one,b(0,0),ldb,a(0),n + 1, &
                                    alpha,b(0,k),ldb)
                          call la_strsm('R','U','N',diag,m,k,one,a(k),n + 1,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_strsm('R','U','T',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_sgemm('N','T',m,k,k,-one,b(0,k),ldb,a(0),n + 1, &
                                    alpha,b(0,0),ldb)
                          call la_strsm('R','L','N',diag,m,k,one,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_strsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_sgemm('N','T',m,k,k,-one,b(0,k),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,0),ldb)
                          call la_strsm('R','U','T',diag,m,k,one,a(k),k,b(0,0), &
                                    ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_strsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_sgemm('N','N',m,k,k,-one,b(0,0),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,k),ldb)
                          call la_strsm('R','L','T',diag,m,k,one,a(0),k,b(0,k), &
                                    ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_strsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_sgemm('N','T',m,k,k,-one,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_strsm('R','L','T',diag,m,k,one,a(k*k),k,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_strsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_sgemm('N','N',m,k,k,-one,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_strsm('R','U','T',diag,m,k,one,a((k + 1)*k),k,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_stfsm
     !> Level 3 BLAS like routine for A in RFP Format.
     !> DTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**T.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_dtfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           real(dp),intent(in) :: alpha
           ! Array Arguments
           real(dp),intent(in) :: a(0:*)
           real(dp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('DTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_dp)
           if (alpha == zero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = zero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_dtrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_dtrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_dgemm('N','N',m2,n,m1,-one,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_dtrsm('L','U','T',diag,m2,n,one,a(m),m,b(m1, &
                                       0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_dtrsm('L','L','T',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_dtrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_dgemm('T','N',m1,n,m2,-one,a(m1),m,b(m1,0), &
                                       ldb,alpha,b,ldb)
                             call la_dtrsm('L','L','T',diag,m1,n,one,a(0),m,b,ldb &
                                       )
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_dtrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_dgemm('T','N',m2,n,m1,-one,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_dtrsm('L','U','T',diag,m2,n,one,a(m1),m,b(m1,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_dtrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_dgemm('N','N',m1,n,m2,-one,a(0),m,b(m1,0),ldb, &
                                     alpha,b,ldb)
                          call la_dtrsm('L','L','T',diag,m1,n,one,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_dtrsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_dtrsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_dgemm('T','N',m2,n,m1,-one,a(m1*m1),m1,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_dtrsm('L','L','N',diag,m2,n,one,a(1),m1,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_dtrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_dtrsm('L','L','T',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_dgemm('N','N',m1,n,m2,-one,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_dtrsm('L','U','N',diag,m1,n,one,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_dtrsm('L','U','T',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_dgemm('N','N',m2,n,m1,-one,a(0),m2,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_dtrsm('L','L','N',diag,m2,n,one,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_dtrsm('L','L','T',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_dgemm('T','N',m1,n,m2,-one,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_dtrsm('L','U','N',diag,m1,n,one,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_dtrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_dgemm('N','N',k,n,k,-one,a(k + 1),m + 1,b,ldb,alpha, &
                                     b(k,0),ldb)
                          call la_dtrsm('L','U','T',diag,k,n,one,a(0),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_dtrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_dgemm('T','N',k,n,k,-one,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_dtrsm('L','L','T',diag,k,n,one,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_dtrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_dgemm('T','N',k,n,k,-one,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_dtrsm('L','U','T',diag,k,n,one,a(k),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_dtrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_dgemm('N','N',k,n,k,-one,a(0),m + 1,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_dtrsm('L','L','T',diag,k,n,one,a(k + 1),m + 1,b,ldb &
                                    )
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_dtrsm('L','U','T',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_dgemm('T','N',k,n,k,-one,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_dtrsm('L','L','N',diag,k,n,one,a(0),k,b(k,0), &
                                    ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_dtrsm('L','L','T',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_dgemm('N','N',k,n,k,-one,a(k*(k + 1)),k,b(k,0), &
                                     ldb,alpha,b,ldb)
                          call la_dtrsm('L','U','N',diag,k,n,one,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_dtrsm('L','U','T',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_dgemm('N','N',k,n,k,-one,a(0),k,b,ldb,alpha,b( &
                                    k,0),ldb)
                          call la_dtrsm('L','L','N',diag,k,n,one,a(k*k),k,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_dtrsm('L','L','T',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_dgemm('T','N',k,n,k,-one,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_dtrsm('L','U','N',diag,k,n,one,a(k*(k + 1)),k,b, &
                                    ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_dtrsm('R','U','T',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_dgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_dtrsm('R','L','N',diag,m,n1,one,a(0),n,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          call la_dtrsm('R','L','T',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_dgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(n1),n, &
                                     alpha,b(0,n1),ldb)
                          call la_dtrsm('R','U','N',diag,m,n2,one,a(n),n,b(0,n1) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_dtrsm('R','L','T',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_dgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(0),n, &
                                    alpha,b(0,n1),ldb)
                          call la_dtrsm('R','U','N',diag,m,n2,one,a(n1),n,b(0,n1 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_dtrsm('R','U','T',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_dgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(0),n, &
                                     alpha,b(0,0),ldb)
                          call la_dtrsm('R','L','N',diag,m,n1,one,a(n2),n,b(0,0) &
                                    ,ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          call la_dtrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_dgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,0),ldb)
                          call la_dtrsm('R','U','T',diag,m,n1,one,a(0),n1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          call la_dtrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_dgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(n1*n1), &
                                     n1,alpha,b(0,n1),ldb)
                          call la_dtrsm('R','L','T',diag,m,n2,one,a(1),n1,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_dtrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_dgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(0),n2, &
                                     alpha,b(0,n1),ldb)
                          call la_dtrsm('R','L','T',diag,m,n2,one,a(n1*n2),n2,b(0, &
                                     n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_dtrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_dgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_dtrsm('R','U','T',diag,m,n1,one,a(n2*n2),n2,b(0, &
                                     0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_dtrsm('R','U','T',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_dgemm('N','N',m,k,k,-one,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_dtrsm('R','L','N',diag,m,k,one,a(1),n + 1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_dtrsm('R','L','T',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_dgemm('N','T',m,k,k,-one,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_dtrsm('R','U','N',diag,m,k,one,a(0),n + 1,b(0,k) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_dtrsm('R','L','T',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_dgemm('N','N',m,k,k,-one,b(0,0),ldb,a(0),n + 1, &
                                    alpha,b(0,k),ldb)
                          call la_dtrsm('R','U','N',diag,m,k,one,a(k),n + 1,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_dtrsm('R','U','T',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_dgemm('N','T',m,k,k,-one,b(0,k),ldb,a(0),n + 1, &
                                    alpha,b(0,0),ldb)
                          call la_dtrsm('R','L','N',diag,m,k,one,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_dtrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_dgemm('N','T',m,k,k,-one,b(0,k),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,0),ldb)
                          call la_dtrsm('R','U','T',diag,m,k,one,a(k),k,b(0,0), &
                                    ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_dtrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_dgemm('N','N',m,k,k,-one,b(0,0),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,k),ldb)
                          call la_dtrsm('R','L','T',diag,m,k,one,a(0),k,b(0,k), &
                                    ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_dtrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_dgemm('N','T',m,k,k,-one,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_dtrsm('R','L','T',diag,m,k,one,a(k*k),k,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_dtrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_dgemm('N','N',m,k,k,-one,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_dtrsm('R','U','T',diag,m,k,one,a((k + 1)*k),k,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_dtfsm
#ifdef LA_WITH_XDP
     !> Level 3 BLAS like routine for A in RFP Format.
     !> XTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**T.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_xtfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           real(xdp),intent(in) :: alpha
           ! Array Arguments
           real(xdp),intent(in) :: a(0:*)
           real(xdp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('XTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_xdp)
           if (alpha == zero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = zero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_xtrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_xtrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_xgemm('N','N',m2,n,m1,-one,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_xtrsm('L','U','T',diag,m2,n,one,a(m),m,b(m1, &
                                       0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_xtrsm('L','L','T',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_xtrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_xgemm('T','N',m1,n,m2,-one,a(m1),m,b(m1,0), &
                                       ldb,alpha,b,ldb)
                             call la_xtrsm('L','L','T',diag,m1,n,one,a(0),m,b,ldb &
                                       )
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_xtrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_xgemm('T','N',m2,n,m1,-one,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_xtrsm('L','U','T',diag,m2,n,one,a(m1),m,b(m1,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_xtrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_xgemm('N','N',m1,n,m2,-one,a(0),m,b(m1,0),ldb, &
                                     alpha,b,ldb)
                          call la_xtrsm('L','L','T',diag,m1,n,one,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_xtrsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_xtrsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_xgemm('T','N',m2,n,m1,-one,a(m1*m1),m1,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_xtrsm('L','L','N',diag,m2,n,one,a(1),m1,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_xtrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_xtrsm('L','L','T',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_xgemm('N','N',m1,n,m2,-one,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_xtrsm('L','U','N',diag,m1,n,one,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_xtrsm('L','U','T',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_xgemm('N','N',m2,n,m1,-one,a(0),m2,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_xtrsm('L','L','N',diag,m2,n,one,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_xtrsm('L','L','T',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_xgemm('T','N',m1,n,m2,-one,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_xtrsm('L','U','N',diag,m1,n,one,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_xtrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_xgemm('N','N',k,n,k,-one,a(k + 1),m + 1,b,ldb,alpha, &
                                     b(k,0),ldb)
                          call la_xtrsm('L','U','T',diag,k,n,one,a(0),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_xtrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_xgemm('T','N',k,n,k,-one,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_xtrsm('L','L','T',diag,k,n,one,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_xtrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_xgemm('T','N',k,n,k,-one,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_xtrsm('L','U','T',diag,k,n,one,a(k),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_xtrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_xgemm('N','N',k,n,k,-one,a(0),m + 1,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_xtrsm('L','L','T',diag,k,n,one,a(k + 1),m + 1,b,ldb &
                                    )
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_xtrsm('L','U','T',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_xgemm('T','N',k,n,k,-one,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_xtrsm('L','L','N',diag,k,n,one,a(0),k,b(k,0), &
                                    ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_xtrsm('L','L','T',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_xgemm('N','N',k,n,k,-one,a(k*(k + 1)),k,b(k,0), &
                                     ldb,alpha,b,ldb)
                          call la_xtrsm('L','U','N',diag,k,n,one,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_xtrsm('L','U','T',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_xgemm('N','N',k,n,k,-one,a(0),k,b,ldb,alpha,b( &
                                    k,0),ldb)
                          call la_xtrsm('L','L','N',diag,k,n,one,a(k*k),k,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_xtrsm('L','L','T',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_xgemm('T','N',k,n,k,-one,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_xtrsm('L','U','N',diag,k,n,one,a(k*(k + 1)),k,b, &
                                    ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_xtrsm('R','U','T',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_xgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_xtrsm('R','L','N',diag,m,n1,one,a(0),n,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          call la_xtrsm('R','L','T',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_xgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(n1),n, &
                                     alpha,b(0,n1),ldb)
                          call la_xtrsm('R','U','N',diag,m,n2,one,a(n),n,b(0,n1) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_xtrsm('R','L','T',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_xgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(0),n, &
                                    alpha,b(0,n1),ldb)
                          call la_xtrsm('R','U','N',diag,m,n2,one,a(n1),n,b(0,n1 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_xtrsm('R','U','T',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_xgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(0),n, &
                                     alpha,b(0,0),ldb)
                          call la_xtrsm('R','L','N',diag,m,n1,one,a(n2),n,b(0,0) &
                                    ,ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          call la_xtrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_xgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,0),ldb)
                          call la_xtrsm('R','U','T',diag,m,n1,one,a(0),n1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          call la_xtrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_xgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(n1*n1), &
                                     n1,alpha,b(0,n1),ldb)
                          call la_xtrsm('R','L','T',diag,m,n2,one,a(1),n1,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_xtrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_xgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(0),n2, &
                                     alpha,b(0,n1),ldb)
                          call la_xtrsm('R','L','T',diag,m,n2,one,a(n1*n2),n2,b(0, &
                                     n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_xtrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_xgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_xtrsm('R','U','T',diag,m,n1,one,a(n2*n2),n2,b(0, &
                                     0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_xtrsm('R','U','T',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_xgemm('N','N',m,k,k,-one,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_xtrsm('R','L','N',diag,m,k,one,a(1),n + 1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_xtrsm('R','L','T',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_xgemm('N','T',m,k,k,-one,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_xtrsm('R','U','N',diag,m,k,one,a(0),n + 1,b(0,k) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_xtrsm('R','L','T',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_xgemm('N','N',m,k,k,-one,b(0,0),ldb,a(0),n + 1, &
                                    alpha,b(0,k),ldb)
                          call la_xtrsm('R','U','N',diag,m,k,one,a(k),n + 1,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_xtrsm('R','U','T',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_xgemm('N','T',m,k,k,-one,b(0,k),ldb,a(0),n + 1, &
                                    alpha,b(0,0),ldb)
                          call la_xtrsm('R','L','N',diag,m,k,one,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_xtrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_xgemm('N','T',m,k,k,-one,b(0,k),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,0),ldb)
                          call la_xtrsm('R','U','T',diag,m,k,one,a(k),k,b(0,0), &
                                    ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_xtrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_xgemm('N','N',m,k,k,-one,b(0,0),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,k),ldb)
                          call la_xtrsm('R','L','T',diag,m,k,one,a(0),k,b(0,k), &
                                    ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_xtrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_xgemm('N','T',m,k,k,-one,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_xtrsm('R','L','T',diag,m,k,one,a(k*k),k,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_xtrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_xgemm('N','N',m,k,k,-one,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_xtrsm('R','U','T',diag,m,k,one,a((k + 1)*k),k,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_xtfsm
#endif
#ifdef LA_WITH_QP
     !> Level 3 BLAS like routine for A in RFP Format.
     !> QTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**T.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_qtfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           real(qp),intent(in) :: alpha
           ! Array Arguments
           real(qp),intent(in) :: a(0:*)
           real(qp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'T')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('QTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_qp)
           if (alpha == zero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = zero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_qtrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_qtrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_qgemm('N','N',m2,n,m1,-one,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_qtrsm('L','U','T',diag,m2,n,one,a(m),m,b(m1, &
                                       0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_qtrsm('L','L','T',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_qtrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_qgemm('T','N',m1,n,m2,-one,a(m1),m,b(m1,0), &
                                       ldb,alpha,b,ldb)
                             call la_qtrsm('L','L','T',diag,m1,n,one,a(0),m,b,ldb &
                                       )
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_qtrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_qgemm('T','N',m2,n,m1,-one,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_qtrsm('L','U','T',diag,m2,n,one,a(m1),m,b(m1,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_qtrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_qgemm('N','N',m1,n,m2,-one,a(0),m,b(m1,0),ldb, &
                                     alpha,b,ldb)
                          call la_qtrsm('L','L','T',diag,m1,n,one,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_qtrsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_qtrsm('L','U','T',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_qgemm('T','N',m2,n,m1,-one,a(m1*m1),m1,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_qtrsm('L','L','N',diag,m2,n,one,a(1),m1,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          if (m == 1) then
                             call la_qtrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_qtrsm('L','L','T',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_qgemm('N','N',m1,n,m2,-one,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_qtrsm('L','U','N',diag,m1,n,one,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_qtrsm('L','U','T',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_qgemm('N','N',m2,n,m1,-one,a(0),m2,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_qtrsm('L','L','N',diag,m2,n,one,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_qtrsm('L','L','T',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_qgemm('T','N',m1,n,m2,-one,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_qtrsm('L','U','N',diag,m1,n,one,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_qtrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_qgemm('N','N',k,n,k,-one,a(k + 1),m + 1,b,ldb,alpha, &
                                     b(k,0),ldb)
                          call la_qtrsm('L','U','T',diag,k,n,one,a(0),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_qtrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_qgemm('T','N',k,n,k,-one,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_qtrsm('L','L','T',diag,k,n,one,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_qtrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_qgemm('T','N',k,n,k,-one,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_qtrsm('L','U','T',diag,k,n,one,a(k),m + 1,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_qtrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_qgemm('N','N',k,n,k,-one,a(0),m + 1,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_qtrsm('L','L','T',diag,k,n,one,a(k + 1),m + 1,b,ldb &
                                    )
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='l', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_qtrsm('L','U','T',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_qgemm('T','N',k,n,k,-one,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_qtrsm('L','L','N',diag,k,n,one,a(0),k,b(k,0), &
                                    ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_qtrsm('L','L','T',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_qgemm('N','N',k,n,k,-one,a(k*(k + 1)),k,b(k,0), &
                                     ldb,alpha,b,ldb)
                          call la_qtrsm('L','U','N',diag,k,n,one,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 't', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_qtrsm('L','U','T',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_qgemm('N','N',k,n,k,-one,a(0),k,b,ldb,alpha,b( &
                                    k,0),ldb)
                          call la_qtrsm('L','L','N',diag,k,n,one,a(k*k),k,b(k,0) &
                                    ,ldb)
                       else
                          ! side  ='l', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_qtrsm('L','L','T',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_qgemm('T','N',k,n,k,-one,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_qtrsm('L','U','N',diag,k,n,one,a(k*(k + 1)),k,b, &
                                    ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_qtrsm('R','U','T',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_qgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_qtrsm('R','L','N',diag,m,n1,one,a(0),n,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 't'
                          call la_qtrsm('R','L','T',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_qgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(n1),n, &
                                     alpha,b(0,n1),ldb)
                          call la_qtrsm('R','U','N',diag,m,n2,one,a(n),n,b(0,n1) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_qtrsm('R','L','T',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_qgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(0),n, &
                                    alpha,b(0,n1),ldb)
                          call la_qtrsm('R','U','N',diag,m,n2,one,a(n1),n,b(0,n1 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 't'
                          call la_qtrsm('R','U','T',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_qgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(0),n, &
                                     alpha,b(0,0),ldb)
                          call la_qtrsm('R','L','N',diag,m,n1,one,a(n2),n,b(0,0) &
                                    ,ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 'n'
                          call la_qtrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_qgemm('N','T',m,n1,n2,-one,b(0,n1),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,0),ldb)
                          call la_qtrsm('R','U','T',diag,m,n1,one,a(0),n1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'l', and
                          ! trans = 't'
                          call la_qtrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_qgemm('N','N',m,n2,n1,-one,b(0,0),ldb,a(n1*n1), &
                                     n1,alpha,b(0,n1),ldb)
                          call la_qtrsm('R','L','T',diag,m,n2,one,a(1),n1,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 'n'
                          call la_qtrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_qgemm('N','T',m,n2,n1,-one,b(0,0),ldb,a(0),n2, &
                                     alpha,b(0,n1),ldb)
                          call la_qtrsm('R','L','T',diag,m,n2,one,a(n1*n2),n2,b(0, &
                                     n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 't', uplo = 'u', and
                          ! trans = 't'
                          call la_qtrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_qgemm('N','N',m,n1,n2,-one,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_qtrsm('R','U','T',diag,m,n1,one,a(n2*n2),n2,b(0, &
                                     0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_qtrsm('R','U','T',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_qgemm('N','N',m,k,k,-one,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_qtrsm('R','L','N',diag,m,k,one,a(1),n + 1,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 't'
                          call la_qtrsm('R','L','T',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_qgemm('N','T',m,k,k,-one,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_qtrsm('R','U','N',diag,m,k,one,a(0),n + 1,b(0,k) &
                                    ,ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_qtrsm('R','L','T',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_qgemm('N','N',m,k,k,-one,b(0,0),ldb,a(0),n + 1, &
                                    alpha,b(0,k),ldb)
                          call la_qtrsm('R','U','N',diag,m,k,one,a(k),n + 1,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 't'
                          call la_qtrsm('R','U','T',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_qgemm('N','T',m,k,k,-one,b(0,k),ldb,a(0),n + 1, &
                                    alpha,b(0,0),ldb)
                          call la_qtrsm('R','L','N',diag,m,k,one,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 't'
                    if (lower) then
                       ! side  ='r', n is even, transr = 't', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 'n'
                          call la_qtrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_qgemm('N','T',m,k,k,-one,b(0,k),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,0),ldb)
                          call la_qtrsm('R','U','T',diag,m,k,one,a(k),k,b(0,0), &
                                    ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'l',
                          ! and trans = 't'
                          call la_qtrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_qgemm('N','N',m,k,k,-one,b(0,0),ldb,a((k + 1)*k &
                                    ),k,alpha,b(0,k),ldb)
                          call la_qtrsm('R','L','T',diag,m,k,one,a(0),k,b(0,k), &
                                    ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 't', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 'n'
                          call la_qtrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_qgemm('N','T',m,k,k,-one,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_qtrsm('R','L','T',diag,m,k,one,a(k*k),k,b(0,k) &
                                    ,ldb)
                       else
                          ! side  ='r', n is even, transr = 't', uplo = 'u',
                          ! and trans = 't'
                          call la_qtrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_qgemm('N','N',m,k,k,-one,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_qtrsm('R','U','T',diag,m,k,one,a((k + 1)*k),k,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_qtfsm
#endif

     !> Level 3 BLAS like routine for C in RFP Format.
     !> CHFRK: performs one of the Hermitian rank--k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n Hermitian
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_chfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           complex(sp) :: calpha,cbeta
           ! Intrinsic Functions
           intrinsic :: max,cmplx
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CHFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_cherk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = czero
              end do
              return
           end if
           calpha = cmplx(alpha,zero,KIND=sp)
           cbeta = cmplx(beta,zero,KIND=sp)
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_cherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_cherk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_cgemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_cherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_cherk('U','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_cgemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_cherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_cherk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_cgemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n2,1), &
                                 lda,cbeta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_cherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_cherk('U','C',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_cgemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n2), &
                                 lda,cbeta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 'c'
                 if (lower) then
                    ! n is odd, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_cherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_cherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_cgemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n1 + 1,1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_cherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_cherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_cgemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n1 + 1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_cherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_cherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_cgemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    else
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_cherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_cherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_cgemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_cherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_cherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_cgemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_cherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_cherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_cgemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_cherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_cherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_cgemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_cherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_cherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_cgemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 'c'
                 if (lower) then
                    ! n is even, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_cherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_cherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_cgemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_cherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_cherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_cgemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_cherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_cherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_cgemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_cherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_cherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_cgemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_chfrk
     !> Level 3 BLAS like routine for C in RFP Format.
     !> ZHFRK: performs one of the Hermitian rank--k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n Hermitian
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_zhfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           complex(dp) :: calpha,cbeta
           ! Intrinsic Functions
           intrinsic :: max,cmplx
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZHFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_zherk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = czero
              end do
              return
           end if
           calpha = cmplx(alpha,zero,KIND=dp)
           cbeta = cmplx(beta,zero,KIND=dp)
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_zherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_zherk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_zgemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_zherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_zherk('U','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_zgemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_zherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_zherk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_zgemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n2,1), &
                                 lda,cbeta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_zherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_zherk('U','C',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_zgemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n2), &
                                 lda,cbeta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 'c'
                 if (lower) then
                    ! n is odd, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_zherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_zherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_zgemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n1 + 1,1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_zherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_zherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_zgemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n1 + 1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_zherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_zherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_zgemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    else
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_zherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_zherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_zgemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_zherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_zherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_zgemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_zherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_zherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_zgemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_zherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_zherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_zgemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_zherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_zherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_zgemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 'c'
                 if (lower) then
                    ! n is even, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_zherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_zherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_zgemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_zherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_zherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_zgemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_zherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_zherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_zgemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_zherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_zherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_zgemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_zhfrk
#ifdef LA_WITH_XDP
     !> Level 3 BLAS like routine for C in RFP Format.
     !> YHFRK: performs one of the Hermitian rank--k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n Hermitian
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_yhfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           complex(xdp) :: calpha,cbeta
           ! Intrinsic Functions
           intrinsic :: max,cmplx
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('YHFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_yherk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = czero
              end do
              return
           end if
           calpha = cmplx(alpha,zero,KIND=xdp)
           cbeta = cmplx(beta,zero,KIND=xdp)
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_yherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_yherk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_ygemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_yherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_yherk('U','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_ygemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_yherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_yherk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_ygemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n2,1), &
                                 lda,cbeta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_yherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_yherk('U','C',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_ygemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n2), &
                                 lda,cbeta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 'c'
                 if (lower) then
                    ! n is odd, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_yherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_yherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_ygemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n1 + 1,1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_yherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_yherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_ygemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n1 + 1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_yherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_yherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_ygemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    else
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_yherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_yherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_ygemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_yherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_yherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_ygemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_yherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_yherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_ygemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_yherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_yherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_ygemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_yherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_yherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_ygemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 'c'
                 if (lower) then
                    ! n is even, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_yherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_yherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_ygemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_yherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_yherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_ygemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_yherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_yherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_ygemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_yherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_yherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_ygemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_yhfrk
#endif
#ifdef LA_WITH_QP
     !> Level 3 BLAS like routine for C in RFP Format.
     !> WHFRK: performs one of the Hermitian rank--k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where alpha and beta are real scalars, C is an n--by--n Hermitian
     !> matrix and A is an n--by--k matrix in the first case and a k--by--n
     !> matrix in the second case.

     pure subroutine la_whfrk(transr,uplo,trans,n,k,alpha,a,lda,beta,c)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,n
           character,intent(in) :: trans,transr,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: c(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,normaltransr,nisodd,notrans
           integer(ilp) :: info,nrowa,j,nk,n1,n2
           complex(qp) :: calpha,cbeta
           ! Intrinsic Functions
           intrinsic :: max,cmplx
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (notrans) then
              nrowa = n
           else
              nrowa = k
           end if
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (k < 0) then
              info = -5
           else if (lda < max(1,nrowa)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WHFRK ',-info)
              return
           end if
           ! quick return if possible.
           ! the quick return case: ((alpha==0).and.(beta/=zero)) is not
           ! done (it is in la_wherk for example) and left in the general case.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           if ((alpha == zero) .and. (beta == zero)) then
              do j = 1, ((n*(n + 1))/2)
                 c(j) = czero
              end do
              return
           end if
           calpha = cmplx(alpha,zero,KIND=qp)
           cbeta = cmplx(beta,zero,KIND=qp)
           ! c is n-by-n.
           ! if n is odd, set nisodd = .true., and n1 and n2.
           ! if n is even, nisodd = .false., and nk.
           if (mod(n,2) == 0) then
              nisodd = .false.
              nk = n/2
           else
              nisodd = .true.
              if (lower) then
                 n2 = n/2
                 n1 = n - n2
              else
                 n1 = n/2
                 n2 = n - n1
              end if
           end if
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                    ! n is odd, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_wherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_wherk('U','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_wgemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_wherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(1),n)

                       call la_wherk('U','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(n + 1) &
                                 ,n)
                       call la_wgemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(n1 + 1),n)
                    end if
                 else
                    ! n is odd, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_wherk('L','N',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_wherk('U','N',n2,k,alpha,a(n2,1),lda,beta,c(n1 + 1), &
                                  n)
                       call la_wgemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n2,1), &
                                 lda,cbeta,c(1),n)
                    else
                       ! n is odd, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_wherk('L','C',n1,k,alpha,a(1,1),lda,beta,c(n2 + 1), &
                                 n)
                       call la_wherk('U','C',n2,k,alpha,a(1,n2),lda,beta,c(n1 + 1), &
                                  n)
                       call la_wgemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n2), &
                                 lda,cbeta,c(1),n)
                    end if
                 end if
              else
                 ! n is odd, and transr = 'c'
                 if (lower) then
                    ! n is odd, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_wherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_wherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c(2), &
                                 n1)
                       call la_wgemm('N','C',n1,n2,k,calpha,a(1,1),lda,a(n1 + 1,1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    else
                       ! n is odd, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_wherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(1),n1 &
                                 )
                       call la_wherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c(2), &
                                 n1)
                       call la_wgemm('C','N',n1,n2,k,calpha,a(1,1),lda,a(1,n1 + 1) &
                                 ,lda,cbeta,c(n1*n1 + 1),n1)
                    end if
                 else
                    ! n is odd, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_wherk('U','N',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_wherk('L','N',n2,k,alpha,a(n1 + 1,1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_wgemm('N','C',n2,n1,k,calpha,a(n1 + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    else
                       ! n is odd, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_wherk('U','C',n1,k,alpha,a(1,1),lda,beta,c(n2*n2 + 1 &
                                 ),n2)
                       call la_wherk('L','C',n2,k,alpha,a(1,n1 + 1),lda,beta,c( &
                                 n1*n2 + 1),n2)
                       call la_wgemm('C','N',n2,n1,k,calpha,a(1,n1 + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),n2)
                    end if
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! n is even, transr = 'n', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'n'
                       call la_wherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_wherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 n + 1)
                       call la_wgemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'l', and trans = 'c'
                       call la_wherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(2),n + &
                                 1)
                       call la_wherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 n + 1)
                       call la_wgemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(nk + 2),n + 1)
                    end if
                 else
                    ! n is even, transr = 'n', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'n'
                       call la_wherk('L','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_wherk('U','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_wgemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(1),n + 1)
                    else
                       ! n is even, transr = 'n', uplo = 'u', and trans = 'c'
                       call la_wherk('L','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 2), &
                                 n + 1)
                       call la_wherk('U','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(nk + 1 &
                                 ),n + 1)
                       call la_wgemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(1),n + 1)
                    end if
                 end if
              else
                 ! n is even, and transr = 'c'
                 if (lower) then
                    ! n is even, transr = 'c', and uplo = 'l'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'n'
                       call la_wherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_wherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c(1), &
                                 nk)
                       call la_wgemm('N','C',nk,nk,k,calpha,a(1,1),lda,a(nk + 1,1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'l', and trans = 'c'
                       call la_wherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk + 1), &
                                 nk)
                       call la_wherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c(1), &
                                 nk)
                       call la_wgemm('C','N',nk,nk,k,calpha,a(1,1),lda,a(1,nk + 1) &
                                 ,lda,cbeta,c(((nk + 1)*nk) + 1),nk)
                    end if
                 else
                    ! n is even, transr = 'c', and uplo = 'u'
                    if (notrans) then
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'n'
                       call la_wherk('U','N',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_wherk('L','N',nk,k,alpha,a(nk + 1,1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_wgemm('N','C',nk,nk,k,calpha,a(nk + 1,1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    else
                       ! n is even, transr = 'c', uplo = 'u', and trans = 'c'
                       call la_wherk('U','C',nk,k,alpha,a(1,1),lda,beta,c(nk*(nk + &
                                 1) + 1),nk)
                       call la_wherk('L','C',nk,k,alpha,a(1,nk + 1),lda,beta,c( &
                                 nk*nk + 1),nk)
                       call la_wgemm('C','N',nk,nk,k,calpha,a(1,nk + 1),lda,a(1,1) &
                                 ,lda,cbeta,c(1),nk)
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_whfrk
#endif

     !> CLACRM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by N and complex; B is N by N and real;
     !> C is M by N and complex.

     pure subroutine la_clacrm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(sp),intent(in) :: b(ldb,*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: aimag,cmplx,real
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(a(i,j),KIND=sp)
              end do
           end do
           l = m*n + 1
           call la_sgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(a(i,j))
              end do
           end do
           call la_sgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=sp),rwork(l + (j - 1)*m + i - 1),KIND=sp)

              end do
           end do
           return
     end subroutine la_clacrm
     !> ZLACRM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by N and complex; B is N by N and real;
     !> C is M by N and complex.

     pure subroutine la_zlacrm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(dp),intent(in) :: b(ldb,*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(a(i,j),KIND=dp)
              end do
           end do
           l = m*n + 1
           call la_dgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(a(i,j))
              end do
           end do
           call la_dgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=dp),rwork(l + (j - 1)*m + i - 1),KIND=dp)

              end do
           end do
           return
     end subroutine la_zlacrm
#ifdef LA_WITH_XDP
     !> YLACRM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by N and complex; B is N by N and real;
     !> C is M by N and complex.

     pure subroutine la_ylacrm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(xdp),intent(in) :: b(ldb,*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(a(i,j),KIND=xdp)
              end do
           end do
           l = m*n + 1
           call la_xgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(a(i,j))
              end do
           end do
           call la_xgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=xdp),rwork(l + (j - 1)*m + i - 1),KIND=xdp)

              end do
           end do
           return
     end subroutine la_ylacrm
#endif
#ifdef LA_WITH_QP
     !> WLACRM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by N and complex; B is N by N and real;
     !> C is M by N and complex.

     pure subroutine la_wlacrm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(qp),intent(in) :: b(ldb,*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(a(i,j),KIND=qp)
              end do
           end do
           l = m*n + 1
           call la_qgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(a(i,j))
              end do
           end do
           call la_qgemm('N','N',m,n,n,one,rwork,m,b,ldb,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=qp),rwork(l + (j - 1)*m + i - 1),KIND=qp)

              end do
           end do
           return
     end subroutine la_wlacrm
#endif

     !> CLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_clagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_sp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(sp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(sp),intent(inout) :: b(ldb,*)
           complex(sp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b + a**t * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b + a**h * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j) + conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) + conjg(du(n - 1))*x(n - 1,j) + conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + conjg(du(i - 1))*x(i - 1,j) + conjg(d(i)) &
                                    *x(i,j) + conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b - a**t*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b - a**h*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j) - conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) - conjg(du(n - 1))*x(n - 1,j) - conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - conjg(du(i - 1))*x(i - 1,j) - conjg(d(i)) &
                                    *x(i,j) - conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_clagtm
     !> ZLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_zlagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_dp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(dp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(dp),intent(inout) :: b(ldb,*)
           complex(dp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b + a**t * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b + a**h * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j) + conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) + conjg(du(n - 1))*x(n - 1,j) + conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + conjg(du(i - 1))*x(i - 1,j) + conjg(d(i)) &
                                    *x(i,j) + conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b - a**t *x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b - a**h *x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j) - conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) - conjg(du(n - 1))*x(n - 1,j) - conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - conjg(du(i - 1))*x(i - 1,j) - conjg(d(i)) &
                                    *x(i,j) - conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_zlagtm
#ifdef LA_WITH_XDP
     !> YLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_ylagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_xdp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(xdp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(xdp),intent(inout) :: b(ldb,*)
           complex(xdp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b + a**t * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b + a**h * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j) + conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) + conjg(du(n - 1))*x(n - 1,j) + conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + conjg(du(i - 1))*x(i - 1,j) + conjg(d(i)) &
                                    *x(i,j) + conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b - a**t *x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b - a**h *x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j) - conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) - conjg(du(n - 1))*x(n - 1,j) - conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - conjg(du(i - 1))*x(i - 1,j) - conjg(d(i)) &
                                    *x(i,j) - conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_ylagtm
#endif
#ifdef LA_WITH_QP
     !> WLAGTM: performs a matrix-vector product of the form
     !> B := alpha * A * X + beta * B
     !> where A is a tridiagonal matrix of order N, B and X are N by NRHS
     !> matrices, and alpha and beta are real scalars, each of which may be
     !> 0., 1., or -1.

     pure subroutine la_wlagtm(trans,n,nrhs,alpha,dl,d,du,x,ldx,beta,b,ldb)
        use la_constants_qp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: trans
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           real(qp),intent(in) :: alpha,beta
           ! Array Arguments
           complex(qp),intent(inout) :: b(ldb,*)
           complex(qp),intent(in) :: d(*),dl(*),du(*),x(ldx,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (n == 0) return
           ! multiply b by beta if beta/=1.
           if (beta == zero) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = zero
                 end do
              end do
           else if (beta == -one) then
              do j = 1,nrhs
                 do i = 1,n
                    b(i,j) = -b(i,j)
                 end do
              end do
           end if
           if (alpha == one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b + a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + du(1)*x(2,j)
                       b(n,j) = b(n,j) + dl(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + dl(i - 1)*x(i - 1,j) + d(i)*x(i,j) + du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b + a**t * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) + d(1)*x(1,j) + dl(1)*x(2,j)
                       b(n,j) = b(n,j) + du(n - 1)*x(n - 1,j) + d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + du(i - 1)*x(i - 1,j) + d(i)*x(i,j) + dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b + a**h * x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) + conjg(d(1))*x(1,j) + conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) + conjg(du(n - 1))*x(n - 1,j) + conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) + conjg(du(i - 1))*x(i - 1,j) + conjg(d(i)) &
                                    *x(i,j) + conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           else if (alpha == -one) then
              if (la_lsame(trans,'N')) then
                 ! compute b := b - a*x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - du(1)*x(2,j)
                       b(n,j) = b(n,j) - dl(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - dl(i - 1)*x(i - 1,j) - d(i)*x(i,j) - du(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'T')) then
                 ! compute b := b - a**t *x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - d(1)*x(1,j)
                    else
                       b(1,j) = b(1,j) - d(1)*x(1,j) - dl(1)*x(2,j)
                       b(n,j) = b(n,j) - du(n - 1)*x(n - 1,j) - d(n)*x(n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - du(i - 1)*x(i - 1,j) - d(i)*x(i,j) - dl(i &
                                    )*x(i + 1,j)
                       end do
                    end if
                 end do
              else if (la_lsame(trans,'C')) then
                 ! compute b := b - a**h *x
                 do j = 1,nrhs
                    if (n == 1) then
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j)
                    else
                       b(1,j) = b(1,j) - conjg(d(1))*x(1,j) - conjg(dl(1))*x(2, &
                                 j)
                       b(n,j) = b(n,j) - conjg(du(n - 1))*x(n - 1,j) - conjg(d(n))*x( &
                                  n,j)
                       do i = 2,n - 1
                          b(i,j) = b(i,j) - conjg(du(i - 1))*x(i - 1,j) - conjg(d(i)) &
                                    *x(i,j) - conjg(dl(i))*x(i + 1,j)
                       end do
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_wlagtm
#endif

     !> CLARCM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by M and real; B is M by N and complex;
     !> C is M by N and complex.

     pure subroutine la_clarcm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: b(ldb,*)
           complex(sp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: aimag,cmplx,real
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(b(i,j),KIND=sp)
              end do
           end do
           l = m*n + 1
           call la_sgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(b(i,j))
              end do
           end do
           call la_sgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=sp),rwork(l + (j - 1)*m + i - 1),KIND=sp)

              end do
           end do
           return
     end subroutine la_clarcm
     !> ZLARCM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by M and real; B is M by N and complex;
     !> C is M by N and complex.

     pure subroutine la_zlarcm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: b(ldb,*)
           complex(dp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(b(i,j),KIND=dp)
              end do
           end do
           l = m*n + 1
           call la_dgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(b(i,j))
              end do
           end do
           call la_dgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=dp),rwork(l + (j - 1)*m + i - 1),KIND=dp)

              end do
           end do
           return
     end subroutine la_zlarcm
#ifdef LA_WITH_XDP
     !> YLARCM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by M and real; B is M by N and complex;
     !> C is M by N and complex.

     pure subroutine la_ylarcm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(in) :: b(ldb,*)
           complex(xdp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(b(i,j),KIND=xdp)
              end do
           end do
           l = m*n + 1
           call la_xgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(b(i,j))
              end do
           end do
           call la_xgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=xdp),rwork(l + (j - 1)*m + i - 1),KIND=xdp)

              end do
           end do
           return
     end subroutine la_ylarcm
#endif
#ifdef LA_WITH_QP
     !> WLARCM: performs a very simple matrix-matrix multiplication:
     !> C := A * B,
     !> where A is M by M and real; B is M by N and complex;
     !> C is M by N and complex.

     pure subroutine la_wlarcm(m,n,a,lda,b,ldb,c,ldc,rwork)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: b(ldb,*)
           complex(qp),intent(out) :: c(ldc,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,j,l
           ! Intrinsic Functions
           intrinsic :: real,cmplx,aimag
           ! Executable Statements
           ! quick return if possible.
           if ((m == 0) .or. (n == 0)) return
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = real(b(i,j),KIND=qp)
              end do
           end do
           l = m*n + 1
           call la_qgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = rwork(l + (j - 1)*m + i - 1)
              end do
           end do
           do j = 1,n
              do i = 1,m
                 rwork((j - 1)*m + i) = aimag(b(i,j))
              end do
           end do
           call la_qgemm('N','N',m,n,m,one,a,lda,rwork,m,zero,rwork(l),m)

           do j = 1,n
              do i = 1,m
                 c(i,j) = cmplx(real(c(i,j),KIND=qp),rwork(l + (j - 1)*m + i - 1),KIND=qp)

              end do
           end do
           return
     end subroutine la_wlarcm
#endif

     !> Level 3 BLAS like routine for A in RFP Format.
     !> CTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**H.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_ctfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           complex(sp),intent(in) :: alpha
           ! Array Arguments
           complex(sp),intent(in) :: a(0:*)
           complex(sp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_sp,0e+0_sp)
           if (alpha == czero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = czero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_ctrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_ctrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_cgemm('N','N',m2,n,m1,-cone,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_ctrsm('L','U','C',diag,m2,n,cone,a(m),m,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_ctrsm('L','L','C',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_ctrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_cgemm('C','N',m1,n,m2,-cone,a(m1),m,b(m1,0), &
                                        ldb,alpha,b,ldb)
                             call la_ctrsm('L','L','C',diag,m1,n,cone,a(0),m,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_ctrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_cgemm('C','N',m2,n,m1,-cone,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_ctrsm('L','U','C',diag,m2,n,cone,a(m1),m,b(m1, &
                                    0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_ctrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_cgemm('N','N',m1,n,m2,-cone,a(0),m,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_ctrsm('L','L','C',diag,m1,n,cone,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_ctrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_ctrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_cgemm('C','N',m2,n,m1,-cone,a(m1*m1),m1,b,ldb, &
                                        alpha,b(m1,0),ldb)
                             call la_ctrsm('L','L','N',diag,m2,n,cone,a(1),m1,b( &
                                       m1,0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_ctrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_ctrsm('L','L','C',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_cgemm('N','N',m1,n,m2,-cone,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_ctrsm('L','U','N',diag,m1,n,cone,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_ctrsm('L','U','C',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_cgemm('N','N',m2,n,m1,-cone,a(0),m2,b,ldb,alpha, &
                                     b(m1,0),ldb)
                          call la_ctrsm('L','L','N',diag,m2,n,cone,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_ctrsm('L','L','C',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_cgemm('C','N',m1,n,m2,-cone,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_ctrsm('L','U','N',diag,m1,n,cone,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_ctrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_cgemm('N','N',k,n,k,-cone,a(k + 1),m + 1,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_ctrsm('L','U','C',diag,k,n,cone,a(0),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_ctrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_cgemm('C','N',k,n,k,-cone,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_ctrsm('L','L','C',diag,k,n,cone,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_ctrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_cgemm('C','N',k,n,k,-cone,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_ctrsm('L','U','C',diag,k,n,cone,a(k),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_ctrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_cgemm('N','N',k,n,k,-cone,a(0),m + 1,b(k,0),ldb, &
                                     alpha,b,ldb)
                          call la_ctrsm('L','L','C',diag,k,n,cone,a(k + 1),m + 1,b, &
                                    ldb)
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_ctrsm('L','U','C',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_cgemm('C','N',k,n,k,-cone,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_ctrsm('L','L','N',diag,k,n,cone,a(0),k,b(k,0), &
                                     ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_ctrsm('L','L','C',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_cgemm('N','N',k,n,k,-cone,a(k*(k + 1)),k,b(k,0) &
                                    ,ldb,alpha,b,ldb)
                          call la_ctrsm('L','U','N',diag,k,n,cone,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_ctrsm('L','U','C',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_cgemm('N','N',k,n,k,-cone,a(0),k,b,ldb,alpha,b( &
                                     k,0),ldb)
                          call la_ctrsm('L','L','N',diag,k,n,cone,a(k*k),k,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_ctrsm('L','L','C',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_cgemm('C','N',k,n,k,-cone,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_ctrsm('L','U','N',diag,k,n,cone,a(k*(k + 1)),k,b, &
                                     ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_ctrsm('R','U','C',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_cgemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_ctrsm('R','L','N',diag,m,n1,cone,a(0),n,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          call la_ctrsm('R','L','C',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_cgemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(n1), &
                                    n,alpha,b(0,n1),ldb)
                          call la_ctrsm('R','U','N',diag,m,n2,cone,a(n),n,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_ctrsm('R','L','C',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_cgemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(0),n, &
                                     alpha,b(0,n1),ldb)
                          call la_ctrsm('R','U','N',diag,m,n2,cone,a(n1),n,b(0, &
                                    n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_ctrsm('R','U','C',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_cgemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n,alpha,b(0,0),ldb)
                          call la_ctrsm('R','L','N',diag,m,n1,cone,a(n2),n,b(0,0 &
                                    ),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          call la_ctrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_cgemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(n1*n1 &
                                    ),n1,alpha,b(0,0),ldb)
                          call la_ctrsm('R','U','C',diag,m,n1,cone,a(0),n1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          call la_ctrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_cgemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,n1),ldb)
                          call la_ctrsm('R','L','C',diag,m,n2,cone,a(1),n1,b(0, &
                                    n1),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_ctrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_cgemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(0), &
                                    n2,alpha,b(0,n1),ldb)
                          call la_ctrsm('R','L','C',diag,m,n2,cone,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_ctrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_cgemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_ctrsm('R','U','C',diag,m,n1,cone,a(n2*n2),n2,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_ctrsm('R','U','C',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_cgemm('N','N',m,k,k,-cone,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_ctrsm('R','L','N',diag,m,k,cone,a(1),n + 1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_ctrsm('R','L','C',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_cgemm('N','C',m,k,k,-cone,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_ctrsm('R','U','N',diag,m,k,cone,a(0),n + 1,b(0,k &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_ctrsm('R','L','C',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_cgemm('N','N',m,k,k,-cone,b(0,0),ldb,a(0),n + 1, &
                                     alpha,b(0,k),ldb)
                          call la_ctrsm('R','U','N',diag,m,k,cone,a(k),n + 1,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_ctrsm('R','U','C',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_cgemm('N','C',m,k,k,-cone,b(0,k),ldb,a(0),n + 1, &
                                     alpha,b(0,0),ldb)
                          call la_ctrsm('R','L','N',diag,m,k,cone,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_ctrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_cgemm('N','C',m,k,k,-cone,b(0,k),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,0),ldb)
                          call la_ctrsm('R','U','C',diag,m,k,cone,a(k),k,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_ctrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_cgemm('N','N',m,k,k,-cone,b(0,0),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,k),ldb)
                          call la_ctrsm('R','L','C',diag,m,k,cone,a(0),k,b(0,k), &
                                     ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_ctrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_cgemm('N','C',m,k,k,-cone,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_ctrsm('R','L','C',diag,m,k,cone,a(k*k),k,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_ctrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_cgemm('N','N',m,k,k,-cone,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_ctrsm('R','U','C',diag,m,k,cone,a((k + 1)*k),k,b( &
                                     0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_ctfsm
     !> Level 3 BLAS like routine for A in RFP Format.
     !> ZTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**H.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_ztfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           complex(dp),intent(in) :: alpha
           ! Array Arguments
           complex(dp),intent(in) :: a(0:*)
           complex(dp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_dp,0e+0_dp)
           if (alpha == czero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = czero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_ztrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_ztrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_zgemm('N','N',m2,n,m1,-cone,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_ztrsm('L','U','C',diag,m2,n,cone,a(m),m,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_ztrsm('L','L','C',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_ztrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_zgemm('C','N',m1,n,m2,-cone,a(m1),m,b(m1,0), &
                                        ldb,alpha,b,ldb)
                             call la_ztrsm('L','L','C',diag,m1,n,cone,a(0),m,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_ztrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_zgemm('C','N',m2,n,m1,-cone,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_ztrsm('L','U','C',diag,m2,n,cone,a(m1),m,b(m1, &
                                    0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_ztrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_zgemm('N','N',m1,n,m2,-cone,a(0),m,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_ztrsm('L','L','C',diag,m1,n,cone,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_ztrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_ztrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_zgemm('C','N',m2,n,m1,-cone,a(m1*m1),m1,b,ldb, &
                                        alpha,b(m1,0),ldb)
                             call la_ztrsm('L','L','N',diag,m2,n,cone,a(1),m1,b( &
                                       m1,0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_ztrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_ztrsm('L','L','C',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_zgemm('N','N',m1,n,m2,-cone,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_ztrsm('L','U','N',diag,m1,n,cone,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_ztrsm('L','U','C',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_zgemm('N','N',m2,n,m1,-cone,a(0),m2,b,ldb,alpha, &
                                     b(m1,0),ldb)
                          call la_ztrsm('L','L','N',diag,m2,n,cone,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_ztrsm('L','L','C',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_zgemm('C','N',m1,n,m2,-cone,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_ztrsm('L','U','N',diag,m1,n,cone,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_ztrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_zgemm('N','N',k,n,k,-cone,a(k + 1),m + 1,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_ztrsm('L','U','C',diag,k,n,cone,a(0),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_ztrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_zgemm('C','N',k,n,k,-cone,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_ztrsm('L','L','C',diag,k,n,cone,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_ztrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_zgemm('C','N',k,n,k,-cone,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_ztrsm('L','U','C',diag,k,n,cone,a(k),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_ztrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_zgemm('N','N',k,n,k,-cone,a(0),m + 1,b(k,0),ldb, &
                                     alpha,b,ldb)
                          call la_ztrsm('L','L','C',diag,k,n,cone,a(k + 1),m + 1,b, &
                                    ldb)
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_ztrsm('L','U','C',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_zgemm('C','N',k,n,k,-cone,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_ztrsm('L','L','N',diag,k,n,cone,a(0),k,b(k,0), &
                                     ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_ztrsm('L','L','C',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_zgemm('N','N',k,n,k,-cone,a(k*(k + 1)),k,b(k,0) &
                                    ,ldb,alpha,b,ldb)
                          call la_ztrsm('L','U','N',diag,k,n,cone,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_ztrsm('L','U','C',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_zgemm('N','N',k,n,k,-cone,a(0),k,b,ldb,alpha,b( &
                                     k,0),ldb)
                          call la_ztrsm('L','L','N',diag,k,n,cone,a(k*k),k,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_ztrsm('L','L','C',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_zgemm('C','N',k,n,k,-cone,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_ztrsm('L','U','N',diag,k,n,cone,a(k*(k + 1)),k,b, &
                                     ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_ztrsm('R','U','C',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_zgemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_ztrsm('R','L','N',diag,m,n1,cone,a(0),n,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          call la_ztrsm('R','L','C',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_zgemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(n1), &
                                    n,alpha,b(0,n1),ldb)
                          call la_ztrsm('R','U','N',diag,m,n2,cone,a(n),n,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_ztrsm('R','L','C',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_zgemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(0),n, &
                                     alpha,b(0,n1),ldb)
                          call la_ztrsm('R','U','N',diag,m,n2,cone,a(n1),n,b(0, &
                                    n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_ztrsm('R','U','C',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_zgemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n,alpha,b(0,0),ldb)
                          call la_ztrsm('R','L','N',diag,m,n1,cone,a(n2),n,b(0,0 &
                                    ),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          call la_ztrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_zgemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(n1*n1 &
                                    ),n1,alpha,b(0,0),ldb)
                          call la_ztrsm('R','U','C',diag,m,n1,cone,a(0),n1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          call la_ztrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_zgemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,n1),ldb)
                          call la_ztrsm('R','L','C',diag,m,n2,cone,a(1),n1,b(0, &
                                    n1),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_ztrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_zgemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(0), &
                                    n2,alpha,b(0,n1),ldb)
                          call la_ztrsm('R','L','C',diag,m,n2,cone,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_ztrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_zgemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_ztrsm('R','U','C',diag,m,n1,cone,a(n2*n2),n2,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_ztrsm('R','U','C',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_zgemm('N','N',m,k,k,-cone,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_ztrsm('R','L','N',diag,m,k,cone,a(1),n + 1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_ztrsm('R','L','C',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_zgemm('N','C',m,k,k,-cone,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_ztrsm('R','U','N',diag,m,k,cone,a(0),n + 1,b(0,k &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_ztrsm('R','L','C',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_zgemm('N','N',m,k,k,-cone,b(0,0),ldb,a(0),n + 1, &
                                     alpha,b(0,k),ldb)
                          call la_ztrsm('R','U','N',diag,m,k,cone,a(k),n + 1,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_ztrsm('R','U','C',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_zgemm('N','C',m,k,k,-cone,b(0,k),ldb,a(0),n + 1, &
                                     alpha,b(0,0),ldb)
                          call la_ztrsm('R','L','N',diag,m,k,cone,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_ztrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_zgemm('N','C',m,k,k,-cone,b(0,k),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,0),ldb)
                          call la_ztrsm('R','U','C',diag,m,k,cone,a(k),k,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_ztrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_zgemm('N','N',m,k,k,-cone,b(0,0),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,k),ldb)
                          call la_ztrsm('R','L','C',diag,m,k,cone,a(0),k,b(0,k), &
                                     ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_ztrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_zgemm('N','C',m,k,k,-cone,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_ztrsm('R','L','C',diag,m,k,cone,a(k*k),k,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_ztrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_zgemm('N','N',m,k,k,-cone,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_ztrsm('R','U','C',diag,m,k,cone,a((k + 1)*k),k,b( &
                                     0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_ztfsm
#ifdef LA_WITH_XDP
     !> Level 3 BLAS like routine for A in RFP Format.
     !> YTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**H.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_ytfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           complex(xdp),intent(in) :: alpha
           ! Array Arguments
           complex(xdp),intent(in) :: a(0:*)
           complex(xdp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('YTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_xdp,0e+0_xdp)
           if (alpha == czero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = czero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_ytrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_ytrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_ygemm('N','N',m2,n,m1,-cone,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_ytrsm('L','U','C',diag,m2,n,cone,a(m),m,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_ytrsm('L','L','C',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_ytrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_ygemm('C','N',m1,n,m2,-cone,a(m1),m,b(m1,0), &
                                        ldb,alpha,b,ldb)
                             call la_ytrsm('L','L','C',diag,m1,n,cone,a(0),m,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_ytrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_ygemm('C','N',m2,n,m1,-cone,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_ytrsm('L','U','C',diag,m2,n,cone,a(m1),m,b(m1, &
                                    0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_ytrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_ygemm('N','N',m1,n,m2,-cone,a(0),m,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_ytrsm('L','L','C',diag,m1,n,cone,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_ytrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_ytrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_ygemm('C','N',m2,n,m1,-cone,a(m1*m1),m1,b,ldb, &
                                        alpha,b(m1,0),ldb)
                             call la_ytrsm('L','L','N',diag,m2,n,cone,a(1),m1,b( &
                                       m1,0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_ytrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_ytrsm('L','L','C',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_ygemm('N','N',m1,n,m2,-cone,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_ytrsm('L','U','N',diag,m1,n,cone,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_ytrsm('L','U','C',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_ygemm('N','N',m2,n,m1,-cone,a(0),m2,b,ldb,alpha, &
                                     b(m1,0),ldb)
                          call la_ytrsm('L','L','N',diag,m2,n,cone,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_ytrsm('L','L','C',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_ygemm('C','N',m1,n,m2,-cone,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_ytrsm('L','U','N',diag,m1,n,cone,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_ytrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_ygemm('N','N',k,n,k,-cone,a(k + 1),m + 1,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_ytrsm('L','U','C',diag,k,n,cone,a(0),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_ytrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_ygemm('C','N',k,n,k,-cone,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_ytrsm('L','L','C',diag,k,n,cone,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_ytrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_ygemm('C','N',k,n,k,-cone,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_ytrsm('L','U','C',diag,k,n,cone,a(k),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_ytrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_ygemm('N','N',k,n,k,-cone,a(0),m + 1,b(k,0),ldb, &
                                     alpha,b,ldb)
                          call la_ytrsm('L','L','C',diag,k,n,cone,a(k + 1),m + 1,b, &
                                    ldb)
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_ytrsm('L','U','C',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_ygemm('C','N',k,n,k,-cone,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_ytrsm('L','L','N',diag,k,n,cone,a(0),k,b(k,0), &
                                     ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_ytrsm('L','L','C',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_ygemm('N','N',k,n,k,-cone,a(k*(k + 1)),k,b(k,0) &
                                    ,ldb,alpha,b,ldb)
                          call la_ytrsm('L','U','N',diag,k,n,cone,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_ytrsm('L','U','C',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_ygemm('N','N',k,n,k,-cone,a(0),k,b,ldb,alpha,b( &
                                     k,0),ldb)
                          call la_ytrsm('L','L','N',diag,k,n,cone,a(k*k),k,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_ytrsm('L','L','C',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_ygemm('C','N',k,n,k,-cone,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_ytrsm('L','U','N',diag,k,n,cone,a(k*(k + 1)),k,b, &
                                     ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_ytrsm('R','U','C',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_ygemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_ytrsm('R','L','N',diag,m,n1,cone,a(0),n,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          call la_ytrsm('R','L','C',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_ygemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(n1), &
                                    n,alpha,b(0,n1),ldb)
                          call la_ytrsm('R','U','N',diag,m,n2,cone,a(n),n,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_ytrsm('R','L','C',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_ygemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(0),n, &
                                     alpha,b(0,n1),ldb)
                          call la_ytrsm('R','U','N',diag,m,n2,cone,a(n1),n,b(0, &
                                    n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_ytrsm('R','U','C',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_ygemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n,alpha,b(0,0),ldb)
                          call la_ytrsm('R','L','N',diag,m,n1,cone,a(n2),n,b(0,0 &
                                    ),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          call la_ytrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_ygemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(n1*n1 &
                                    ),n1,alpha,b(0,0),ldb)
                          call la_ytrsm('R','U','C',diag,m,n1,cone,a(0),n1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          call la_ytrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_ygemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,n1),ldb)
                          call la_ytrsm('R','L','C',diag,m,n2,cone,a(1),n1,b(0, &
                                    n1),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_ytrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_ygemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(0), &
                                    n2,alpha,b(0,n1),ldb)
                          call la_ytrsm('R','L','C',diag,m,n2,cone,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_ytrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_ygemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_ytrsm('R','U','C',diag,m,n1,cone,a(n2*n2),n2,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_ytrsm('R','U','C',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_ygemm('N','N',m,k,k,-cone,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_ytrsm('R','L','N',diag,m,k,cone,a(1),n + 1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_ytrsm('R','L','C',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_ygemm('N','C',m,k,k,-cone,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_ytrsm('R','U','N',diag,m,k,cone,a(0),n + 1,b(0,k &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_ytrsm('R','L','C',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_ygemm('N','N',m,k,k,-cone,b(0,0),ldb,a(0),n + 1, &
                                     alpha,b(0,k),ldb)
                          call la_ytrsm('R','U','N',diag,m,k,cone,a(k),n + 1,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_ytrsm('R','U','C',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_ygemm('N','C',m,k,k,-cone,b(0,k),ldb,a(0),n + 1, &
                                     alpha,b(0,0),ldb)
                          call la_ytrsm('R','L','N',diag,m,k,cone,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_ytrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_ygemm('N','C',m,k,k,-cone,b(0,k),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,0),ldb)
                          call la_ytrsm('R','U','C',diag,m,k,cone,a(k),k,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_ytrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_ygemm('N','N',m,k,k,-cone,b(0,0),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,k),ldb)
                          call la_ytrsm('R','L','C',diag,m,k,cone,a(0),k,b(0,k), &
                                     ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_ytrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_ygemm('N','C',m,k,k,-cone,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_ytrsm('R','L','C',diag,m,k,cone,a(k*k),k,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_ytrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_ygemm('N','N',m,k,k,-cone,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_ytrsm('R','U','C',diag,m,k,cone,a((k + 1)*k),k,b( &
                                     0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_ytfsm
#endif
#ifdef LA_WITH_QP
     !> Level 3 BLAS like routine for A in RFP Format.
     !> WTFSM:  solves the matrix equation
     !> op( A )*X = alpha*B  or  X*op( A ) = alpha*B
     !> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     !> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     !> op( A ) = A   or   op( A ) = A**H.
     !> A is in Rectangular Full Packed (RFP) Format.
     !> The matrix X is overwritten on B.

     pure subroutine la_wtfsm(transr,side,uplo,trans,diag,m,n,alpha,a,b,ldb)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,diag,side,trans,uplo
           integer(ilp),intent(in) :: ldb,m,n
           complex(qp),intent(in) :: alpha
           ! Array Arguments
           complex(qp),intent(in) :: a(0:*)
           complex(qp),intent(inout) :: b(0:ldb - 1,0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,lside,misodd,nisodd,normaltransr,notrans
           integer(ilp) :: m1,m2,n1,n2,k,info,i,j
           ! Intrinsic Functions
           intrinsic :: max,mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lside = la_lsame(side,'L')
           lower = la_lsame(uplo,'L')
           notrans = la_lsame(trans,'N')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lside .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -3
           else if (.not. notrans .and. .not. la_lsame(trans,'C')) then
              info = -4
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -5
           else if (m < 0) then
              info = -6
           else if (n < 0) then
              info = -7
           else if (ldb < max(1,m)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WTFSM ',-info)
              return
           end if
           ! quick return when ( (n==0).or.(m==0) )
           if ((m == 0) .or. (n == 0)) return
           ! quick return when alpha==(0e+0_qp,0e+0_qp)
           if (alpha == czero) then
              do j = 0,n - 1
                 do i = 0,m - 1
                    b(i,j) = czero
                 end do
              end do
              return
           end if
           if (lside) then
              ! side = 'l'
              ! a is m-by-m.
              ! if m is odd, set nisodd = .true., and m1 and m2.
              ! if m is even, nisodd = .false., and m.
              if (mod(m,2) == 0) then
                 misodd = .false.
                 k = m/2
              else
                 misodd = .true.
                 if (lower) then
                    m2 = m/2
                    m1 = m - m2
                 else
                    m1 = m/2
                    m2 = m - m1
                 end if
              end if
              if (misodd) then
                 ! side = 'l' and n is odd
                 if (normaltransr) then
                    ! side = 'l', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_wtrsm('L','L','N',diag,m1,n,alpha,a,m,b,ldb)

                          else
                             call la_wtrsm('L','L','N',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                             call la_wgemm('N','N',m2,n,m1,-cone,a(m1),m,b,ldb, &
                                       alpha,b(m1,0),ldb)
                             call la_wtrsm('L','U','C',diag,m2,n,cone,a(m),m,b(m1, &
                                        0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_wtrsm('L','L','C',diag,m1,n,alpha,a(0),m,b, &
                                       ldb)
                          else
                             call la_wtrsm('L','U','N',diag,m2,n,alpha,a(m),m,b( &
                                       m1,0),ldb)
                             call la_wgemm('C','N',m1,n,m2,-cone,a(m1),m,b(m1,0), &
                                        ldb,alpha,b,ldb)
                             call la_wtrsm('L','L','C',diag,m1,n,cone,a(0),m,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_wtrsm('L','L','N',diag,m1,n,alpha,a(m2),m,b,ldb &
                                    )
                          call la_wgemm('C','N',m2,n,m1,-cone,a(0),m,b,ldb,alpha, &
                                    b(m1,0),ldb)
                          call la_wtrsm('L','U','C',diag,m2,n,cone,a(m1),m,b(m1, &
                                    0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_wtrsm('L','U','N',diag,m2,n,alpha,a(m1),m,b(m1, &
                                    0),ldb)
                          call la_wgemm('N','N',m1,n,m2,-cone,a(0),m,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_wtrsm('L','L','C',diag,m1,n,cone,a(m2),m,b,ldb)

                       end if
                    end if
                 else
                    ! side = 'l', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          if (m == 1) then
                             call la_wtrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_wtrsm('L','U','C',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                             call la_wgemm('C','N',m2,n,m1,-cone,a(m1*m1),m1,b,ldb, &
                                        alpha,b(m1,0),ldb)
                             call la_wtrsm('L','L','N',diag,m2,n,cone,a(1),m1,b( &
                                       m1,0),ldb)
                          end if
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          if (m == 1) then
                             call la_wtrsm('L','U','N',diag,m1,n,alpha,a(0),m1,b, &
                                       ldb)
                          else
                             call la_wtrsm('L','L','C',diag,m2,n,alpha,a(1),m1,b( &
                                       m1,0),ldb)
                             call la_wgemm('N','N',m1,n,m2,-cone,a(m1*m1),m1,b(m1, &
                                       0),ldb,alpha,b,ldb)
                             call la_wtrsm('L','U','N',diag,m1,n,cone,a(0),m1,b, &
                                       ldb)
                          end if
                       end if
                    else
                       ! side  ='l', n is odd, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_wtrsm('L','U','C',diag,m1,n,alpha,a(m2*m2),m2,b, &
                                    ldb)
                          call la_wgemm('N','N',m2,n,m1,-cone,a(0),m2,b,ldb,alpha, &
                                     b(m1,0),ldb)
                          call la_wtrsm('L','L','N',diag,m2,n,cone,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                       else
                          ! side  ='l', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_wtrsm('L','L','C',diag,m2,n,alpha,a(m1*m2),m2,b( &
                                    m1,0),ldb)
                          call la_wgemm('C','N',m1,n,m2,-cone,a(0),m2,b(m1,0), &
                                    ldb,alpha,b,ldb)
                          call la_wtrsm('L','U','N',diag,m1,n,cone,a(m2*m2),m2,b, &
                                    ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'l' and n is even
                 if (normaltransr) then
                    ! side = 'l', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_wtrsm('L','L','N',diag,k,n,alpha,a(1),m + 1,b,ldb &
                                    )
                          call la_wgemm('N','N',k,n,k,-cone,a(k + 1),m + 1,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_wtrsm('L','U','C',diag,k,n,cone,a(0),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_wtrsm('L','U','N',diag,k,n,alpha,a(0),m + 1,b(k, &
                                    0),ldb)
                          call la_wgemm('C','N',k,n,k,-cone,a(k + 1),m + 1,b(k,0), &
                                    ldb,alpha,b,ldb)
                          call la_wtrsm('L','L','C',diag,k,n,cone,a(1),m + 1,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'n', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_wtrsm('L','L','N',diag,k,n,alpha,a(k + 1),m + 1,b, &
                                    ldb)
                          call la_wgemm('C','N',k,n,k,-cone,a(0),m + 1,b,ldb,alpha, &
                                    b(k,0),ldb)
                          call la_wtrsm('L','U','C',diag,k,n,cone,a(k),m + 1,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_wtrsm('L','U','N',diag,k,n,alpha,a(k),m + 1,b(k, &
                                    0),ldb)
                          call la_wgemm('N','N',k,n,k,-cone,a(0),m + 1,b(k,0),ldb, &
                                     alpha,b,ldb)
                          call la_wtrsm('L','L','C',diag,k,n,cone,a(k + 1),m + 1,b, &
                                    ldb)
                       end if
                    end if
                 else
                    ! side = 'l', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='l', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_wtrsm('L','U','C',diag,k,n,alpha,a(k),k,b,ldb)

                          call la_wgemm('C','N',k,n,k,-cone,a(k*(k + 1)),k,b,ldb, &
                                    alpha,b(k,0),ldb)
                          call la_wtrsm('L','L','N',diag,k,n,cone,a(0),k,b(k,0), &
                                     ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_wtrsm('L','L','C',diag,k,n,alpha,a(0),k,b(k,0) &
                                    ,ldb)
                          call la_wgemm('N','N',k,n,k,-cone,a(k*(k + 1)),k,b(k,0) &
                                    ,ldb,alpha,b,ldb)
                          call la_wtrsm('L','U','N',diag,k,n,cone,a(k),k,b,ldb)

                       end if
                    else
                       ! side  ='l', n is even, transr = 'c', and uplo = 'u'
                       if (.not. notrans) then
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_wtrsm('L','U','C',diag,k,n,alpha,a(k*(k + 1)),k, &
                                    b,ldb)
                          call la_wgemm('N','N',k,n,k,-cone,a(0),k,b,ldb,alpha,b( &
                                     k,0),ldb)
                          call la_wtrsm('L','L','N',diag,k,n,cone,a(k*k),k,b(k,0 &
                                    ),ldb)
                       else
                          ! side  ='l', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_wtrsm('L','L','C',diag,k,n,alpha,a(k*k),k,b(k, &
                                    0),ldb)
                          call la_wgemm('C','N',k,n,k,-cone,a(0),k,b(k,0),ldb, &
                                    alpha,b,ldb)
                          call la_wtrsm('L','U','N',diag,k,n,cone,a(k*(k + 1)),k,b, &
                                     ldb)
                       end if
                    end if
                 end if
              end if
           else
              ! side = 'r'
              ! a is n-by-n.
              ! if n is odd, set nisodd = .true., and n1 and n2.
              ! if n is even, nisodd = .false., and k.
              if (mod(n,2) == 0) then
                 nisodd = .false.
                 k = n/2
              else
                 nisodd = .true.
                 if (lower) then
                    n2 = n/2
                    n1 = n - n2
                 else
                    n1 = n/2
                    n2 = n - n1
                 end if
              end if
              if (nisodd) then
                 ! side = 'r' and n is odd
                 if (normaltransr) then
                    ! side = 'r', n is odd, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'n'
                          call la_wtrsm('R','U','C',diag,m,n2,alpha,a(n),n,b(0, &
                                    n1),ldb)
                          call la_wgemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(n1), &
                                    n,alpha,b(0,0),ldb)
                          call la_wtrsm('R','L','N',diag,m,n1,cone,a(0),n,b(0,0) &
                                    ,ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'l', and
                          ! trans = 'c'
                          call la_wtrsm('R','L','C',diag,m,n1,alpha,a(0),n,b(0,0 &
                                    ),ldb)
                          call la_wgemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(n1), &
                                    n,alpha,b(0,n1),ldb)
                          call la_wtrsm('R','U','N',diag,m,n2,cone,a(n),n,b(0,n1 &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'n'
                          call la_wtrsm('R','L','C',diag,m,n1,alpha,a(n2),n,b(0, &
                                    0),ldb)
                          call la_wgemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(0),n, &
                                     alpha,b(0,n1),ldb)
                          call la_wtrsm('R','U','N',diag,m,n2,cone,a(n1),n,b(0, &
                                    n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'n', uplo = 'u', and
                          ! trans = 'c'
                          call la_wtrsm('R','U','C',diag,m,n2,alpha,a(n1),n,b(0, &
                                    n1),ldb)
                          call la_wgemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n,alpha,b(0,0),ldb)
                          call la_wtrsm('R','L','N',diag,m,n1,cone,a(n2),n,b(0,0 &
                                    ),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is odd, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'n'
                          call la_wtrsm('R','L','N',diag,m,n2,alpha,a(1),n1,b(0, &
                                    n1),ldb)
                          call la_wgemm('N','C',m,n1,n2,-cone,b(0,n1),ldb,a(n1*n1 &
                                    ),n1,alpha,b(0,0),ldb)
                          call la_wtrsm('R','U','C',diag,m,n1,cone,a(0),n1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'l', and
                          ! trans = 'c'
                          call la_wtrsm('R','U','N',diag,m,n1,alpha,a(0),n1,b(0, &
                                    0),ldb)
                          call la_wgemm('N','N',m,n2,n1,-cone,b(0,0),ldb,a(n1*n1) &
                                    ,n1,alpha,b(0,n1),ldb)
                          call la_wtrsm('R','L','C',diag,m,n2,cone,a(1),n1,b(0, &
                                    n1),ldb)
                       end if
                    else
                       ! side  ='r', n is odd, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'n'
                          call la_wtrsm('R','U','N',diag,m,n1,alpha,a(n2*n2),n2,b( &
                                    0,0),ldb)
                          call la_wgemm('N','C',m,n2,n1,-cone,b(0,0),ldb,a(0), &
                                    n2,alpha,b(0,n1),ldb)
                          call la_wtrsm('R','L','C',diag,m,n2,cone,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                       else
                          ! side  ='r', n is odd, transr = 'c', uplo = 'u', and
                          ! trans = 'c'
                          call la_wtrsm('R','L','N',diag,m,n2,alpha,a(n1*n2),n2,b( &
                                    0,n1),ldb)
                          call la_wgemm('N','N',m,n1,n2,-cone,b(0,n1),ldb,a(0), &
                                    n2,alpha,b(0,0),ldb)
                          call la_wtrsm('R','U','C',diag,m,n1,cone,a(n2*n2),n2,b( &
                                    0,0),ldb)
                       end if
                    end if
                 end if
              else
                 ! side = 'r' and n is even
                 if (normaltransr) then
                    ! side = 'r', n is even, and transr = 'n'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'n', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'n'
                          call la_wtrsm('R','U','C',diag,m,k,alpha,a(0),n + 1,b(0, &
                                    k),ldb)
                          call la_wgemm('N','N',m,k,k,-cone,b(0,k),ldb,a(k + 1),n + &
                                    1,alpha,b(0,0),ldb)
                          call la_wtrsm('R','L','N',diag,m,k,cone,a(1),n + 1,b(0,0 &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'l',
                          ! and trans = 'c'
                          call la_wtrsm('R','L','C',diag,m,k,alpha,a(1),n + 1,b(0, &
                                    0),ldb)
                          call la_wgemm('N','C',m,k,k,-cone,b(0,0),ldb,a(k + 1),n + &
                                    1,alpha,b(0,k),ldb)
                          call la_wtrsm('R','U','N',diag,m,k,cone,a(0),n + 1,b(0,k &
                                    ),ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'n', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'n'
                          call la_wtrsm('R','L','C',diag,m,k,alpha,a(k + 1),n + 1,b(0, &
                                     0),ldb)
                          call la_wgemm('N','N',m,k,k,-cone,b(0,0),ldb,a(0),n + 1, &
                                     alpha,b(0,k),ldb)
                          call la_wtrsm('R','U','N',diag,m,k,cone,a(k),n + 1,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'n', uplo = 'u',
                          ! and trans = 'c'
                          call la_wtrsm('R','U','C',diag,m,k,alpha,a(k),n + 1,b(0, &
                                    k),ldb)
                          call la_wgemm('N','C',m,k,k,-cone,b(0,k),ldb,a(0),n + 1, &
                                     alpha,b(0,0),ldb)
                          call la_wtrsm('R','L','N',diag,m,k,cone,a(k + 1),n + 1,b(0, &
                                    0),ldb)
                       end if
                    end if
                 else
                    ! side = 'r', n is even, and transr = 'c'
                    if (lower) then
                       ! side  ='r', n is even, transr = 'c', and uplo = 'l'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'n'
                          call la_wtrsm('R','L','N',diag,m,k,alpha,a(0),k,b(0,k) &
                                    ,ldb)
                          call la_wgemm('N','C',m,k,k,-cone,b(0,k),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,0),ldb)
                          call la_wtrsm('R','U','C',diag,m,k,cone,a(k),k,b(0,0), &
                                     ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'l',
                          ! and trans = 'c'
                          call la_wtrsm('R','U','N',diag,m,k,alpha,a(k),k,b(0,0) &
                                    ,ldb)
                          call la_wgemm('N','N',m,k,k,-cone,b(0,0),ldb,a((k + 1) &
                                    *k),k,alpha,b(0,k),ldb)
                          call la_wtrsm('R','L','C',diag,m,k,cone,a(0),k,b(0,k), &
                                     ldb)
                       end if
                    else
                       ! side  ='r', n is even, transr = 'c', and uplo = 'u'
                       if (notrans) then
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'n'
                          call la_wtrsm('R','U','N',diag,m,k,alpha,a((k + 1)*k),k, &
                                    b(0,0),ldb)
                          call la_wgemm('N','C',m,k,k,-cone,b(0,0),ldb,a(0),k, &
                                    alpha,b(0,k),ldb)
                          call la_wtrsm('R','L','C',diag,m,k,cone,a(k*k),k,b(0,k &
                                    ),ldb)
                       else
                          ! side  ='r', n is even, transr = 'c', uplo = 'u',
                          ! and trans = 'c'
                          call la_wtrsm('R','L','N',diag,m,k,alpha,a(k*k),k,b(0, &
                                    k),ldb)
                          call la_wgemm('N','N',m,k,k,-cone,b(0,k),ldb,a(0),k, &
                                    alpha,b(0,0),ldb)
                          call la_wtrsm('R','U','C',diag,m,k,cone,a((k + 1)*k),k,b( &
                                     0,0),ldb)
                       end if
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_wtfsm
#endif

end module la_lapack_blas_like_l3
