!> BLAS level 3: symmetric matrix-matrix operations
module la_blas_level3_sym
     use la_constants
     use la_blas_aux
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_ssymm
     public :: la_ssyr2k
     public :: la_ssyrk
     public :: la_dsymm
     public :: la_dsyr2k
     public :: la_dsyrk
#ifdef LA_WITH_XDP
     public :: la_xsymm
     public :: la_xsyr2k
     public :: la_xsyrk
#endif
#ifdef LA_WITH_QP
     public :: la_qsymm
     public :: la_qsyr2k
     public :: la_qsyrk
#endif
     public :: la_csymm
     public :: la_csyr2k
     public :: la_csyrk
     public :: la_zsymm
     public :: la_zsyr2k
     public :: la_zsyrk
#ifdef LA_WITH_XDP
     public :: la_ysymm
     public :: la_ysyr2k
     public :: la_ysyrk
#endif
#ifdef LA_WITH_QP
     public :: la_wsymm
     public :: la_wsyr2k
     public :: la_wsyrk
#endif

     contains

     !> SSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
     !> C are  m by n matrices.

     pure subroutine la_ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),b(ldb,*)
           real(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(sp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('SSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (beta == zero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = zero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == zero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_ssymm
     !> DSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
     !> C are  m by n matrices.

     pure subroutine la_dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),b(ldb,*)
           real(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(dp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('DSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (beta == zero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = zero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == zero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_dsymm
#ifdef LA_WITH_XDP
     !> XSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
     !> C are  m by n matrices.

     pure subroutine la_xsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),b(ldb,*)
           real(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(xdp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('XSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (beta == zero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = zero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == zero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_xsymm
#endif
#ifdef LA_WITH_QP
     !> QSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
     !> C are  m by n matrices.

     pure subroutine la_qsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),b(ldb,*)
           real(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(qp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('QSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (beta == zero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = zero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == zero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_qsymm
#endif

     !> SSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_ssyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),b(ldb,*)
           real(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(sp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('SSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_ssyr2k
     !> DSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),b(ldb,*)
           real(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(dp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('DSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_dsyr2k
#ifdef LA_WITH_XDP
     !> XSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_xsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),b(ldb,*)
           real(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(xdp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('XSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_xsyr2k
#endif
#ifdef LA_WITH_QP
     !> QSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_qsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),b(ldb,*)
           real(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(qp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('QSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= zero) .or. (b(j,l) /= zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = zero
                           temp2 = zero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_qsyr2k
#endif

     !> SSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_ssyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('SSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_ssyrk
     !> DSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('DSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_dsyrk
#ifdef LA_WITH_XDP
     !> XSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_xsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(xdp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('XSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_xsyrk
#endif
#ifdef LA_WITH_QP
     !> QSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_qsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T')) .and. ( &
                     .not. la_lsame(trans,'C'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('QSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = zero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == zero) then
                           do i = j,n
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= zero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_qsyrk
#endif

     !> CSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where  alpha and beta are scalars, A is a symmetric matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_csymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(sp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('CSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = czero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == czero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_csymm
     !> ZSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where  alpha and beta are scalars, A is a symmetric matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(dp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('ZSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = czero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == czero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_zsymm
#ifdef LA_WITH_XDP
     !> YSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where  alpha and beta are scalars, A is a symmetric matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_ysymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(xdp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('YSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = czero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == czero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_ysymm
#endif
#ifdef LA_WITH_QP
     !> WSYMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where  alpha and beta are scalars, A is a symmetric matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_wsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: lda,ldb,ldc,m,n
           character,intent(in) :: side,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(qp) :: temp1,temp2
           integer(ilp) :: i,info,j,k,nrowa
           logical(lk) :: upper
           
           ! set nrowa as the number of rows of a.
           if (la_lsame(side,'L')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = la_lsame(uplo,'U')
           ! test the input parameters.
           info = 0
           if ((.not. la_lsame(side,'L')) .and. (.not. la_lsame(side,'R'))) then
               info = 1
           else if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,m)) then
               info = 9
           else if (ldc < max(1,m)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('WSYMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = czero
                       end do
                   end do
               else
                   do j = 1,n
                       do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(side,'L')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = czero
                           do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta == czero) then
                       do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do
                   else
                       do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do
                   end if
                   do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
               end do loop_170
           end if
           return
     end subroutine la_wsymm
#endif

     !> CSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_csyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(sp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('CSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_csyr2k
     !> ZSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(dp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('ZSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_zsyr2k
#ifdef LA_WITH_XDP
     !> YSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_ysyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(xdp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('YSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_ysyr2k
#endif
#ifdef LA_WITH_QP
     !> WSYR2K:  performs one of the symmetric rank 2k operations
     !> C := alpha*A*B**T + alpha*B*A**T + beta*C,
     !> or
     !> C := alpha*A**T*B + alpha*B**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     !> matrices in the second case.

     pure subroutine la_wsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(qp) :: temp1,temp2
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldb < max(1,nrowa)) then
               info = 9
           else if (ldc < max(1,n)) then
               info = 12
           end if
           if (info /= 0) then
               call la_xerbla('WSYR2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_wsyr2k
#endif

     !> CSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_csyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('CSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_csyrk
     !> ZSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('ZSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_zsyrk
#ifdef LA_WITH_XDP
     !> YSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_ysyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(xdp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('YSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_ysyrk
#endif
#ifdef LA_WITH_QP
     !> WSYRK:  performs one of the symmetric rank k operations
     !> C := alpha*A*A**T + beta*C,
     !> or
     !> C := alpha*A**T*A + beta*C,
     !> where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     !> and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     !> in the second case.

     pure subroutine la_wsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,j,l,nrowa
           logical(lk) :: upper
           
           ! test the input parameters.
           if (la_lsame(trans,'N')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = la_lsame(uplo,'U')
           info = 0
           if ((.not. upper) .and. (.not. la_lsame(uplo,'L'))) then
               info = 1
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'T'))) &
                     then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1,nrowa)) then
               info = 7
           else if (ldc < max(1,n)) then
               info = 10
           end if
           if (info /= 0) then
               call la_xerbla('WSYRK ',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               else
                   if (beta == czero) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == czero) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           if (a(j,l) /= czero) then
                               temp = alpha*a(j,l)
                               do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_wsyrk
#endif

end module la_blas_level3_sym
