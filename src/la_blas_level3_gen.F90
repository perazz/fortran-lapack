!> BLAS level 3: general and Hermitian matrix-matrix operations
module la_blas_level3_gen
     use la_constants
     use la_blas_aux
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sgemm
     public :: la_dgemm
#ifdef LA_WITH_XDP
     public :: la_xgemm
#endif
#ifdef LA_WITH_QP
     public :: la_qgemm
#endif
     public :: la_cgemm
     public :: la_chemm
     public :: la_cher2k
     public :: la_cherk
     public :: la_zgemm
     public :: la_zhemm
     public :: la_zher2k
     public :: la_zherk
#ifdef LA_WITH_XDP
     public :: la_ygemm
     public :: la_yhemm
     public :: la_yher2k
     public :: la_yherk
#endif
#ifdef LA_WITH_QP
     public :: la_wgemm
     public :: la_whemm
     public :: la_wher2k
     public :: la_wherk
#endif

     contains

     !> SGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),b(ldb,*)
           real(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. la_lsame(transa,'C')) .and. (.not. la_lsame(transa, &
                     'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. la_lsame(transb,'C')) .and. (.not. la_lsame( &
                     transb,'T'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('SGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           ! and if  alpha.eq.zero.
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (nota) then
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_sgemm
     !> DGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),b(ldb,*)
           real(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. la_lsame(transa,'C')) .and. (.not. la_lsame(transa, &
                     'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. la_lsame(transb,'C')) .and. (.not. la_lsame( &
                     transb,'T'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('DGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           ! and if  alpha.eq.zero.
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (nota) then
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_dgemm
#ifdef LA_WITH_XDP
     !> XGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_xgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),b(ldb,*)
           real(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(xdp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. la_lsame(transa,'C')) .and. (.not. la_lsame(transa, &
                     'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. la_lsame(transb,'C')) .and. (.not. la_lsame( &
                     transb,'T'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('XGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           ! and if  alpha.eq.zero.
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (nota) then
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_xgemm
#endif
#ifdef LA_WITH_QP
     !> QGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_qgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),b(ldb,*)
           real(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: max
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. la_lsame(transa,'C')) .and. (.not. la_lsame(transa, &
                     'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. la_lsame(transb,'C')) .and. (.not. la_lsame( &
                     transb,'T'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('QGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           ! and if  alpha.eq.zero.
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (nota) then
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,m
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = zero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_qgemm
#endif

     !> CGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: conja,conjb,nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and nrowb  as the number of rows  of  a  and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           conja = la_lsame(transa,'C')
           conjb = la_lsame(transb,'C')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. conja) .and. (.not. la_lsame(transa,'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. conjb) .and. (.not. la_lsame(transb,'T'))) &
                     then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('CGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) &
                     return
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else if (conja) then
                 ! form  c := alpha*a**h*b + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else if (nota) then
               if (conjb) then
                 ! form  c := alpha*a*b**h + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*conjg(b(j,l))
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               end if
           else if (conja) then
               if (conjb) then
                 ! form  c := alpha*a**h*b**h + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**h*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(j,l)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (conjb) then
                 ! form  c := alpha*a**t*b**h + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_cgemm
     !> ZGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: conja,conjb,nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and nrowb  as the number of rows  of  a  and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           conja = la_lsame(transa,'C')
           conjb = la_lsame(transb,'C')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. conja) .and. (.not. la_lsame(transa,'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. conjb) .and. (.not. la_lsame(transb,'T'))) &
                     then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('ZGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) &
                     return
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else if (conja) then
                 ! form  c := alpha*a**h*b + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else if (nota) then
               if (conjb) then
                 ! form  c := alpha*a*b**h + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*conjg(b(j,l))
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               end if
           else if (conja) then
               if (conjb) then
                 ! form  c := alpha*a**h*b**h + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**h*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(j,l)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (conjb) then
                 ! form  c := alpha*a**t*b**h + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_zgemm
#ifdef LA_WITH_XDP
     !> YGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_ygemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Local Scalars
           complex(xdp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: conja,conjb,nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and nrowb  as the number of rows  of  a  and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           conja = la_lsame(transa,'C')
           conjb = la_lsame(transb,'C')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. conja) .and. (.not. la_lsame(transa,'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. conjb) .and. (.not. la_lsame(transb,'T'))) &
                     then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('YGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) &
                     return
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else if (conja) then
                 ! form  c := alpha*a**h*b + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else if (nota) then
               if (conjb) then
                 ! form  c := alpha*a*b**h + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*conjg(b(j,l))
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               end if
           else if (conja) then
               if (conjb) then
                 ! form  c := alpha*a**h*b**h + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**h*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(j,l)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (conjb) then
                 ! form  c := alpha*a**t*b**h + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_ygemm
#endif
#ifdef LA_WITH_QP
     !> WGEMM:  performs one of the matrix-matrix operations
     !> C := alpha*op( A )*op( B ) + beta*C,
     !> where  op( X ) is one of
     !> op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
     !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
     !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     pure subroutine la_wgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,m,n
           character,intent(in) :: transa,transb
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,j,l,nrowa,nrowb
           logical(lk) :: conja,conjb,nota,notb
           
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and nrowb  as the number of rows  of  a  and  b  respectively.
           nota = la_lsame(transa,'N')
           notb = la_lsame(transb,'N')
           conja = la_lsame(transa,'C')
           conjb = la_lsame(transb,'C')
           if (nota) then
               nrowa = m
           else
               nrowa = k
           end if
           if (notb) then
               nrowb = k
           else
               nrowb = n
           end if
           ! test the input parameters.
           info = 0
           if ((.not. nota) .and. (.not. conja) .and. (.not. la_lsame(transa,'T'))) then
               info = 1
           else if ((.not. notb) .and. (.not. conjb) .and. (.not. la_lsame(transb,'T'))) &
                     then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1,nrowa)) then
               info = 8
           else if (ldb < max(1,nrowb)) then
               info = 10
           else if (ldc < max(1,m)) then
               info = 13
           end if
           if (info /= 0) then
               call la_xerbla('WGEMM ',info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) &
                     return
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
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(l,j)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else if (conja) then
                 ! form  c := alpha*a**h*b + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else if (nota) then
               if (conjb) then
                 ! form  c := alpha*a*b**h + beta*c.
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*conjg(b(j,l))
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1,n
                       if (beta == czero) then
                           do i = 1,m
                               c(i,j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do
                       end if
                       do l = 1,k
                           temp = alpha*b(j,l)
                           do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do
                       end do
                   end do
               end if
           else if (conja) then
               if (conjb) then
                 ! form  c := alpha*a**h*b**h + beta*c.
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**h*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*b(j,l)
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               end if
           else
               if (conjb) then
                 ! form  c := alpha*a**t*b**h + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*conjg(b(j,l))
                           end do
                           if (beta == czero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1,n
                       do i = 1,m
                           temp = czero
                           do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
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
     end subroutine la_wgemm
#endif

     !> CHEMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars, A is an hermitian matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_chemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
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
           intrinsic :: conjg,max,real
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
               call la_xerbla('CHEMM ',info)
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=sp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=sp) + &
                                         alpha*temp2
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=sp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=sp) + &
                                         alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*real(a(j,j),KIND=sp)
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
                           temp1 = alpha*conjg(a(j,k))
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*conjg(a(j,k))
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
     end subroutine la_chemm
     !> ZHEMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars, A is an hermitian matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
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
           intrinsic :: real,conjg,max
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
               call la_xerbla('ZHEMM ',info)
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=dp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=dp) + &
                                         alpha*temp2
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=dp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=dp) + &
                                         alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*real(a(j,j),KIND=dp)
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
                           temp1 = alpha*conjg(a(j,k))
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*conjg(a(j,k))
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
     end subroutine la_zhemm
#ifdef LA_WITH_XDP
     !> YHEMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars, A is an hermitian matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_yhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
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
           intrinsic :: real,conjg,max
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
               call la_xerbla('YHEMM ',info)
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=xdp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=xdp) + &
                                         alpha*temp2
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=xdp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=xdp) + &
                                         alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*real(a(j,j),KIND=xdp)
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
                           temp1 = alpha*conjg(a(j,k))
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*conjg(a(j,k))
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
     end subroutine la_yhemm
#endif
#ifdef LA_WITH_QP
     !> WHEMM:  performs one of the matrix-matrix operations
     !> C := alpha*A*B + beta*C,
     !> or
     !> C := alpha*B*A + beta*C,
     !> where alpha and beta are scalars, A is an hermitian matrix and  B and
     !> C are m by n matrices.

     pure subroutine la_whemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
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
           intrinsic :: real,conjg,max
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
               call la_xerbla('WHEMM ',info)
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=qp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=qp) + &
                                         alpha*temp2
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
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
                           end do
                           if (beta == czero) then
                               c(i,j) = temp1*real(a(i,i),KIND=qp) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i),KIND=qp) + &
                                         alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1,n
                   temp1 = alpha*real(a(j,j),KIND=qp)
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
                           temp1 = alpha*conjg(a(j,k))
                       end if
                       do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do
                   end do
                   do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*conjg(a(j,k))
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
     end subroutine la_whemm
#endif

     !> CHER2K:  performs one of the hermitian rank 2k operations
     !> C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
     !> or
     !> C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
     !> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
     !> hermitian matrix and  A and B  are  n by k matrices in the first case
     !> and  k by n  matrices in the second case.

     pure subroutine la_cher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(sp),intent(in) :: alpha
           real(sp),intent(in) :: beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*),b(ldb,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: conjg,max,real
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('CHER2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == real(czero,KIND=sp)) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                       end do
                   end if
               else
                   if (beta == real(czero,KIND=sp)) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
               if (upper) then
                   do j = 1,n
                       if (beta == real(czero,KIND=sp)) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                       else
                           c(j,j) = real(c(j,j),KIND=sp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=sp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=sp)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == real(czero,KIND=sp)) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                       else
                           c(j,j) = real(c(j,j),KIND=sp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=sp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=sp)
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=sp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=sp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=sp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=sp)
                               end if
                           else
                               if (beta == real(czero,KIND=sp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=sp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=sp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=sp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=sp)
                               end if
                           else
                               if (beta == real(czero,KIND=sp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_cher2k
     !> ZHER2K:  performs one of the hermitian rank 2k operations
     !> C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
     !> or
     !> C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
     !> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
     !> hermitian matrix and  A and B  are  n by k matrices in the first case
     !> and  k by n  matrices in the second case.

     pure subroutine la_zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(dp),intent(in) :: alpha
           real(dp),intent(in) :: beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*),b(ldb,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('ZHER2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == real(czero,KIND=dp)) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                       end do
                   end if
               else
                   if (beta == real(czero,KIND=dp)) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
               if (upper) then
                   do j = 1,n
                       if (beta == real(czero,KIND=dp)) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                       else
                           c(j,j) = real(c(j,j),KIND=dp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=dp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=dp)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == real(czero,KIND=dp)) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                       else
                           c(j,j) = real(c(j,j),KIND=dp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=dp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=dp)
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=dp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=dp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=dp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=dp)
                               end if
                           else
                               if (beta == real(czero,KIND=dp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=dp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=dp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=dp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=dp)
                               end if
                           else
                               if (beta == real(czero,KIND=dp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_zher2k
#ifdef LA_WITH_XDP
     !> YHER2K:  performs one of the hermitian rank 2k operations
     !> C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
     !> or
     !> C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
     !> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
     !> hermitian matrix and  A and B  are  n by k matrices in the first case
     !> and  k by n  matrices in the second case.

     pure subroutine la_yher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(xdp),intent(in) :: alpha
           real(xdp),intent(in) :: beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('YHER2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == real(czero,KIND=xdp)) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                       end do
                   end if
               else
                   if (beta == real(czero,KIND=xdp)) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
               if (upper) then
                   do j = 1,n
                       if (beta == real(czero,KIND=xdp)) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                       else
                           c(j,j) = real(c(j,j),KIND=xdp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=xdp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=xdp)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == real(czero,KIND=xdp)) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                       else
                           c(j,j) = real(c(j,j),KIND=xdp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=xdp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=xdp)
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=xdp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=xdp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=xdp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=xdp)
                               end if
                           else
                               if (beta == real(czero,KIND=xdp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=xdp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=xdp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=xdp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=xdp)
                               end if
                           else
                               if (beta == real(czero,KIND=xdp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_yher2k
#endif
#ifdef LA_WITH_QP
     !> WHER2K:  performs one of the hermitian rank 2k operations
     !> C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
     !> or
     !> C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
     !> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
     !> hermitian matrix and  A and B  are  n by k matrices in the first case
     !> and  k by n  matrices in the second case.

     pure subroutine la_wher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           complex(qp),intent(in) :: alpha
           real(qp),intent(in) :: beta
           integer(ilp),intent(in) :: k,lda,ldb,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*),b(ldb,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('WHER2K',info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == real(czero,KIND=qp)) then
                       do j = 1,n
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                       end do
                   end if
               else
                   if (beta == real(czero,KIND=qp)) then
                       do j = 1,n
                           do i = j,n
                               c(i,j) = czero
                           end do
                       end do
                   else
                       do j = 1,n
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
               if (upper) then
                   do j = 1,n
                       if (beta == real(czero,KIND=qp)) then
                           do i = 1,j
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                       else
                           c(j,j) = real(c(j,j),KIND=qp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=qp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=qp)
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       if (beta == real(czero,KIND=qp)) then
                           do i = j,n
                               c(i,j) = czero
                           end do
                       else if (beta /= one) then
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                       else
                           c(j,j) = real(c(j,j),KIND=qp)
                       end if
                       do l = 1,k
                           if ((a(j,l) /= czero) .or. (b(j,l) /= czero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 + b(i,l)*temp2
                               end do
                               c(j,j) = real(c(j,j),KIND=qp) + real(a(j,l)*temp1 + b(j,l)*temp2, &
                                         KIND=qp)
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=qp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=qp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=qp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=qp)
                               end if
                           else
                               if (beta == real(czero,KIND=qp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               else
                   do j = 1,n
                       do i = j,n
                           temp1 = czero
                           temp2 = czero
                           do l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
                           end do
                           if (i == j) then
                               if (beta == real(czero,KIND=qp)) then
                                   c(j,j) = real(alpha*temp1 + conjg(alpha)*temp2,KIND=qp)
                               else
                                   c(j,j) = beta*real(c(j,j),KIND=qp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2,KIND=qp)
                               end if
                           else
                               if (beta == real(czero,KIND=qp)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               end if
           end if
           return
     end subroutine la_wher2k
#endif

     !> CHERK:  performs one of the hermitian rank k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
     !> matrix and  A  is an  n by k  matrix in the  first case and a  k by n
     !> matrix in the second case.

     pure subroutine la_cherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_sp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: cmplx,conjg,max,real
           ! Local Scalars
           complex(sp) :: temp
           real(sp) :: rtemp
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('CHERK ',info)
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
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=sp)
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
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**h + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                       else
                           c(j,j) = real(c(j,j),KIND=sp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=sp)) then
                               temp = alpha*conjg(a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                               c(j,j) = real(c(j,j),KIND=sp) + real(temp*a(i,l),KIND=sp)
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
                           c(j,j) = beta*real(c(j,j),KIND=sp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       else
                           c(j,j) = real(c(j,j),KIND=sp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=sp)) then
                               temp = alpha*conjg(a(j,l))
                               c(j,j) = real(c(j,j),KIND=sp) + real(temp*a(j,l),KIND=sp)
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j - 1
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=sp)
                       end if
                   end do
               else
                   do j = 1,n
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=sp)
                       end if
                       do i = j + 1,n
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
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
     end subroutine la_cherk
     !> ZHERK:  performs one of the hermitian rank k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
     !> matrix and  A  is an  n by k  matrix in the  first case and a  k by n
     !> matrix in the second case.

     pure subroutine la_zherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_dp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,max
           ! Local Scalars
           complex(dp) :: temp
           real(dp) :: rtemp
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('ZHERK ',info)
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
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=dp)
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
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**h + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                       else
                           c(j,j) = real(c(j,j),KIND=dp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=dp)) then
                               temp = alpha*conjg(a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                               c(j,j) = real(c(j,j),KIND=dp) + real(temp*a(i,l),KIND=dp)
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
                           c(j,j) = beta*real(c(j,j),KIND=dp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       else
                           c(j,j) = real(c(j,j),KIND=dp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=dp)) then
                               temp = alpha*conjg(a(j,l))
                               c(j,j) = real(c(j,j),KIND=dp) + real(temp*a(j,l),KIND=dp)
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j - 1
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=dp)
                       end if
                   end do
               else
                   do j = 1,n
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=dp)
                       end if
                       do i = j + 1,n
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
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
     end subroutine la_zherk
#ifdef LA_WITH_XDP
     !> YHERK:  performs one of the hermitian rank k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
     !> matrix and  A  is an  n by k  matrix in the  first case and a  k by n
     !> matrix in the second case.

     pure subroutine la_yherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_xdp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,max
           ! Local Scalars
           complex(xdp) :: temp
           real(xdp) :: rtemp
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('YHERK ',info)
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
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
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
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**h + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                       else
                           c(j,j) = real(c(j,j),KIND=xdp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=xdp)) then
                               temp = alpha*conjg(a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                               c(j,j) = real(c(j,j),KIND=xdp) + real(temp*a(i,l),KIND=xdp)
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
                           c(j,j) = beta*real(c(j,j),KIND=xdp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       else
                           c(j,j) = real(c(j,j),KIND=xdp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=xdp)) then
                               temp = alpha*conjg(a(j,l))
                               c(j,j) = real(c(j,j),KIND=xdp) + real(temp*a(j,l),KIND=xdp)
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j - 1
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=xdp)
                       end if
                   end do
               else
                   do j = 1,n
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=xdp)
                       end if
                       do i = j + 1,n
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
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
     end subroutine la_yherk
#endif
#ifdef LA_WITH_QP
     !> WHERK:  performs one of the hermitian rank k operations
     !> C := alpha*A*A**H + beta*C,
     !> or
     !> C := alpha*A**H*A + beta*C,
     !> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
     !> matrix and  A  is an  n by k  matrix in the  first case and a  k by n
     !> matrix in the second case.

     pure subroutine la_wherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
        use la_constants_qp
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: alpha,beta
           integer(ilp),intent(in) :: k,lda,ldc,n
           character,intent(in) :: trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: c(ldc,*)
        ! =====================================================================
           ! Intrinsic Functions
           intrinsic :: real,cmplx,conjg,max
           ! Local Scalars
           complex(qp) :: temp
           real(qp) :: rtemp
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
           else if ((.not. la_lsame(trans,'N')) .and. (.not. la_lsame(trans,'C'))) &
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
               call la_xerbla('WHERK ',info)
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
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=qp)
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
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (la_lsame(trans,'N')) then
              ! form  c := alpha*a*a**h + beta*c.
               if (upper) then
                   do j = 1,n
                       if (beta == zero) then
                           do i = 1,j
                               c(i,j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                       else
                           c(j,j) = real(c(j,j),KIND=qp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=qp)) then
                               temp = alpha*conjg(a(j,l))
                               do i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                               c(j,j) = real(c(j,j),KIND=qp) + real(temp*a(i,l),KIND=qp)
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
                           c(j,j) = beta*real(c(j,j),KIND=qp)
                           do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do
                       else
                           c(j,j) = real(c(j,j),KIND=qp)
                       end if
                       do l = 1,k
                           if (a(j,l) /= cmplx(zero,KIND=qp)) then
                               temp = alpha*conjg(a(j,l))
                               c(j,j) = real(c(j,j),KIND=qp) + real(temp*a(j,l),KIND=qp)
                               do i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*a + beta*c.
               if (upper) then
                   do j = 1,n
                       do i = 1,j - 1
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
                           end do
                           if (beta == zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=qp)
                       end if
                   end do
               else
                   do j = 1,n
                       rtemp = zero
                       do l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
                       end do
                       if (beta == zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j),KIND=qp)
                       end if
                       do i = j + 1,n
                           temp = zero
                           do l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
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
     end subroutine la_wherk
#endif

end module la_blas_level3_gen
