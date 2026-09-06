!> BLAS level 2: triangular matrix-vector operations
module la_blas_level2_tri
     use la_constants
     use la_blas_aux
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_stbmv
     public :: la_stbsv
     public :: la_stpmv
     public :: la_stpsv
     public :: la_strmv
     public :: la_strsv
     public :: la_dtbmv
     public :: la_dtbsv
     public :: la_dtpmv
     public :: la_dtpsv
     public :: la_dtrmv
     public :: la_dtrsv
     public :: la_qtbmv
     public :: la_qtbsv
     public :: la_qtpmv
     public :: la_qtpsv
     public :: la_qtrmv
     public :: la_qtrsv
     public :: la_ctbmv
     public :: la_ctbsv
     public :: la_ctpmv
     public :: la_ctpsv
     public :: la_ctrmv
     public :: la_ctrsv
     public :: la_ztbmv
     public :: la_ztbsv
     public :: la_ztpmv
     public :: la_ztpsv
     public :: la_ztrmv
     public :: la_ztrsv
     public :: la_wtbmv
     public :: la_wtbsv
     public :: la_wtpmv
     public :: la_wtpsv
     public :: la_wtrmv
     public :: la_wtrsv

     contains

     !> STBMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     pure subroutine la_stbmv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('STBMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
               ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do i = j - 1,max(1,j - k),-1
                               temp = temp + a(l + i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do i = j - 1,max(1,j - k),-1
                               temp = temp + a(l + i,j)*x(ix)
                               ix = ix - incx
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do i = j + 1,min(n,j + k)
                               temp = temp + a(l + i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do i = j + 1,min(n,j + k)
                               temp = temp + a(l + i,j)*x(ix)
                               ix = ix + incx
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_stbmv
     !> DTBMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     pure subroutine la_dtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('DTBMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
               ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do i = j - 1,max(1,j - k),-1
                               temp = temp + a(l + i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do i = j - 1,max(1,j - k),-1
                               temp = temp + a(l + i,j)*x(ix)
                               ix = ix - incx
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do i = j + 1,min(n,j + k)
                               temp = temp + a(l + i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do i = j + 1,min(n,j + k)
                               temp = temp + a(l + i,j)*x(ix)
                               ix = ix + incx
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_dtbmv
     !> QTBMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     pure subroutine la_qtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('QTBMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
               ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do i = j - 1,max(1,j - k),-1
                               temp = temp + a(l + i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do i = j - 1,max(1,j - k),-1
                               temp = temp + a(l + i,j)*x(ix)
                               ix = ix - incx
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do i = j + 1,min(n,j + k)
                               temp = temp + a(l + i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do i = j + 1,min(n,j + k)
                               temp = temp + a(l + i,j)*x(ix)
                               ix = ix + incx
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_qtbmv

     !> STBSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
     !> diagonals.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_stbsv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('STBSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do i = j - 1,max(1,j - k),-1
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           kx = kx - incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do i = j - 1,max(1,j - k),-1
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do i = j + 1,min(n,j + k)
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           kx = kx + incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do i = j + 1,min(n,j + k)
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t)*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           do i = max(1,j - k),j - 1
                               temp = temp - a(l + i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(kplus1,j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           do i = max(1,j - k),j - 1
                               temp = temp - a(l + i,j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(kplus1,j)
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           do i = min(n,j + k),j + 1,-1
                               temp = temp - a(l + i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(1,j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           do i = min(n,j + k),j + 1,-1
                               temp = temp - a(l + i,j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(1,j)
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_stbsv
     !> DTBSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
     !> diagonals.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_dtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('DTBSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do i = j - 1,max(1,j - k),-1
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           kx = kx - incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do i = j - 1,max(1,j - k),-1
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do i = j + 1,min(n,j + k)
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           kx = kx + incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do i = j + 1,min(n,j + k)
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t)*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           do i = max(1,j - k),j - 1
                               temp = temp - a(l + i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(kplus1,j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           do i = max(1,j - k),j - 1
                               temp = temp - a(l + i,j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(kplus1,j)
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           do i = min(n,j + k),j + 1,-1
                               temp = temp - a(l + i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(1,j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           do i = min(n,j + k),j + 1,-1
                               temp = temp - a(l + i,j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(1,j)
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_dtbsv
     !> QTBSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
     !> diagonals.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_qtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('QTBSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do i = j - 1,max(1,j - k),-1
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           kx = kx - incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do i = j - 1,max(1,j - k),-1
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do i = j + 1,min(n,j + k)
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           kx = kx + incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do i = j + 1,min(n,j + k)
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t)*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           do i = max(1,j - k),j - 1
                               temp = temp - a(l + i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(kplus1,j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           do i = max(1,j - k),j - 1
                               temp = temp - a(l + i,j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(kplus1,j)
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           do i = min(n,j + k),j + 1,-1
                               temp = temp - a(l + i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(1,j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           do i = min(n,j + k),j + 1,-1
                               temp = temp - a(l + i,j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(1,j)
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_qtbsv

     !> STPMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix, supplied in packed form.

     pure subroutine la_stpmv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('STPMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x:= a*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk + j - 1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk - n + j)
                           end if
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk - (n - (j + 1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk - n + j)
                           end if
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk - 1
                           do i = j - 1,1,-1
                               temp = temp + ap(k)*x(i)
                               k = k - 1
                           end do
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk - 1,kk - j + 1,-1
                               ix = ix - incx
                               temp = temp + ap(k)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk + 1
                           do i = j + 1,n
                               temp = temp + ap(k)*x(i)
                               k = k + 1
                           end do
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk + 1,kk + n - j
                               ix = ix + incx
                               temp = temp + ap(k)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_stpmv
     !> DTPMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix, supplied in packed form.

     pure subroutine la_dtpmv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('DTPMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x:= a*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk + j - 1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk - n + j)
                           end if
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk - (n - (j + 1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk - n + j)
                           end if
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk - 1
                           do i = j - 1,1,-1
                               temp = temp + ap(k)*x(i)
                               k = k - 1
                           end do
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk - 1,kk - j + 1,-1
                               ix = ix - incx
                               temp = temp + ap(k)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk + 1
                           do i = j + 1,n
                               temp = temp + ap(k)*x(i)
                               k = k + 1
                           end do
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk + 1,kk + n - j
                               ix = ix + incx
                               temp = temp + ap(k)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_dtpmv
     !> QTPMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix, supplied in packed form.

     pure subroutine la_qtpmv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('QTPMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x:= a*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk + j - 1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk - n + j)
                           end if
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk - (n - (j + 1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk - n + j)
                           end if
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk - 1
                           do i = j - 1,1,-1
                               temp = temp + ap(k)*x(i)
                               k = k - 1
                           end do
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk - 1,kk - j + 1,-1
                               ix = ix - incx
                               temp = temp + ap(k)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk + 1
                           do i = j + 1,n
                               temp = temp + ap(k)*x(i)
                               k = k + 1
                           end do
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk + 1,kk + n - j
                               ix = ix + incx
                               temp = temp + ap(k)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_qtpmv

     !> STPSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix, supplied in packed form.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_stpsv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('STPSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk
                           do i = 1,j - 1
                               temp = temp - ap(k)*x(i)
                               k = k + 1
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           do k = kk,kk + j - 2
                               temp = temp - ap(k)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk
                           do i = n,j + 1,-1
                               temp = temp - ap(k)*x(i)
                               k = k - 1
                           end do
                           if (nounit) temp = temp/ap(kk - n + j)
                           x(j) = temp
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do k = kk,kk - (n - (j + 1)),-1
                               temp = temp - ap(k)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/ap(kk - n + j)
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_stpsv
     !> DTPSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix, supplied in packed form.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_dtpsv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('DTPSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk
                           do i = 1,j - 1
                               temp = temp - ap(k)*x(i)
                               k = k + 1
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           do k = kk,kk + j - 2
                               temp = temp - ap(k)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk
                           do i = n,j + 1,-1
                               temp = temp - ap(k)*x(i)
                               k = k - 1
                           end do
                           if (nounit) temp = temp/ap(kk - n + j)
                           x(j) = temp
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do k = kk,kk - (n - (j + 1)),-1
                               temp = temp - ap(k)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/ap(kk - n + j)
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_dtpsv
     !> QTPSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix, supplied in packed form.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_qtpsv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('QTPSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk
                           do i = 1,j - 1
                               temp = temp - ap(k)*x(i)
                               k = k + 1
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           do k = kk,kk + j - 2
                               temp = temp - ap(k)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk
                           do i = n,j + 1,-1
                               temp = temp - ap(k)*x(i)
                               k = k - 1
                           end do
                           if (nounit) temp = temp/ap(kk - n + j)
                           x(j) = temp
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do k = kk,kk - (n - (j + 1)),-1
                               temp = temp - ap(k)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/ap(kk - n + j)
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_qtpsv

     !> STRMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix.

     pure subroutine la_strmv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('STRMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do i = j + 1,n
                               temp = temp + a(i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_strmv
     !> DTRMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix.

     pure subroutine la_dtrmv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('DTRMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do i = j + 1,n
                               temp = temp + a(i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_dtrmv
     !> QTRMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix.

     pure subroutine la_qtrmv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('QTRMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do i = j + 1,n
                               temp = temp + a(i,j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_qtrmv

     !> STRSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_strsv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('STRSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           do i = 1,j - 1
                               temp = temp - a(i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           do i = 1,j - 1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           do i = n,j + 1,-1
                               temp = temp - a(i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do i = n,j + 1,-1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_strsv
     !> DTRSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_dtrsv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('DTRSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           do i = 1,j - 1
                               temp = temp - a(i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           do i = 1,j - 1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           do i = n,j + 1,-1
                               temp = temp - a(i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do i = n,j + 1,-1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_dtrsv
     !> QTRSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_qtrsv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           real(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('QTRSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           do i = 1,j - 1
                               temp = temp - a(i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           do i = 1,j - 1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           do i = n,j + 1,-1
                               temp = temp - a(i,j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do i = n,j + 1,-1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_qtrsv

     !> CTBMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     pure subroutine la_ctbmv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('CTBMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
               ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + a(l + i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + conjg(a(l + i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + a(l + i,j)*x(ix)
                                   ix = ix - incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + conjg(a(l + i,j))*x(ix)
                                   ix = ix - incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do i = j + 1,min(n,j + k)
                                   temp = temp + a(l + i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do i = j + 1,min(n,j + k)
                                   temp = temp + conjg(a(l + i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do i = j + 1,min(n,j + k)
                                   temp = temp + a(l + i,j)*x(ix)
                                   ix = ix + incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do i = j + 1,min(n,j + k)
                                   temp = temp + conjg(a(l + i,j))*x(ix)
                                   ix = ix + incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ctbmv
     !> ZTBMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     pure subroutine la_ztbmv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('ZTBMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
               ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + a(l + i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + conjg(a(l + i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + a(l + i,j)*x(ix)
                                   ix = ix - incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + conjg(a(l + i,j))*x(ix)
                                   ix = ix - incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do i = j + 1,min(n,j + k)
                                   temp = temp + a(l + i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do i = j + 1,min(n,j + k)
                                   temp = temp + conjg(a(l + i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do i = j + 1,min(n,j + k)
                                   temp = temp + a(l + i,j)*x(ix)
                                   ix = ix + incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do i = j + 1,min(n,j + k)
                                   temp = temp + conjg(a(l + i,j))*x(ix)
                                   ix = ix + incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ztbmv
     !> WTBMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     pure subroutine la_wtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('WTBMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
               ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1,j - k),j - 1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(i) = x(i) + temp*a(l + i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n,j + k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + a(l + i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + conjg(a(l + i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + a(l + i,j)*x(ix)
                                   ix = ix - incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do i = j - 1,max(1,j - k),-1
                                   temp = temp + conjg(a(l + i,j))*x(ix)
                                   ix = ix - incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do i = j + 1,min(n,j + k)
                                   temp = temp + a(l + i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do i = j + 1,min(n,j + k)
                                   temp = temp + conjg(a(l + i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do i = j + 1,min(n,j + k)
                                   temp = temp + a(l + i,j)*x(ix)
                                   ix = ix + incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do i = j + 1,min(n,j + k)
                                   temp = temp + conjg(a(l + i,j))*x(ix)
                                   ix = ix + incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_wtbmv

     !> CTBSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
     !> diagonals.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_ctbsv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('CTBSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do i = j - 1,max(1,j - k),-1
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           kx = kx - incx
                           if (x(jx) /= czero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do i = j - 1,max(1,j - k),-1
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do i = j + 1,min(n,j + k)
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           kx = kx + incx
                           if (x(jx) /= czero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do i = j + 1,min(n,j + k)
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1,j - k),j - 1
                                   temp = temp - a(l + i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do i = max(1,j - k),j - 1
                                   temp = temp - conjg(a(l + i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1,j - k),j - 1
                                   temp = temp - a(l + i,j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do i = max(1,j - k),j - 1
                                   temp = temp - conjg(a(l + i,j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - a(l + i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(1,j)
                           else
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - conjg(a(l + i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - a(l + i,j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(1,j)
                           else
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - conjg(a(l + i,j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ctbsv
     !> ZTBSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
     !> diagonals.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_ztbsv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('ZTBSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do i = j - 1,max(1,j - k),-1
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           kx = kx - incx
                           if (x(jx) /= czero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do i = j - 1,max(1,j - k),-1
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do i = j + 1,min(n,j + k)
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           kx = kx + incx
                           if (x(jx) /= czero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do i = j + 1,min(n,j + k)
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1,j - k),j - 1
                                   temp = temp - a(l + i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do i = max(1,j - k),j - 1
                                   temp = temp - conjg(a(l + i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1,j - k),j - 1
                                   temp = temp - a(l + i,j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do i = max(1,j - k),j - 1
                                   temp = temp - conjg(a(l + i,j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - a(l + i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(1,j)
                           else
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - conjg(a(l + i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - a(l + i,j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(1,j)
                           else
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - conjg(a(l + i,j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ztbsv
     !> WTBSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular band matrix, with ( k + 1 )
     !> diagonals.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_wtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,k,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kplus1,kx,l
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < (k + 1)) then
               info = 7
           else if (incx == 0) then
               info = 9
           end if
           if (info /= 0) then
               call la_xerbla('WTBSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do i = j - 1,max(1,j - k),-1
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           kx = kx - incx
                           if (x(jx) /= czero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do i = j - 1,max(1,j - k),-1
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do i = j + 1,min(n,j + k)
                                   x(i) = x(i) - temp*a(l + i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           kx = kx + incx
                           if (x(jx) /= czero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do i = j + 1,min(n,j + k)
                                   x(ix) = x(ix) - temp*a(l + i,j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1,j - k),j - 1
                                   temp = temp - a(l + i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do i = max(1,j - k),j - 1
                                   temp = temp - conjg(a(l + i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1,j - k),j - 1
                                   temp = temp - a(l + i,j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do i = max(1,j - k),j - 1
                                   temp = temp - conjg(a(l + i,j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - a(l + i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(1,j)
                           else
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - conjg(a(l + i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - a(l + i,j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(1,j)
                           else
                               do i = min(n,j + k),j + 1,-1
                                   temp = temp - conjg(a(l + i,j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_wtbsv

     !> CTPMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix, supplied in packed form.

     pure subroutine la_ctpmv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: ap(*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('CTPMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x:= a*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               k = kk
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk + j - 1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               k = kk
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk - n + j)
                           end if
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk - (n - (j + 1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk - n + j)
                           end if
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j - 1,1,-1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j - 1,1,-1
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(ap(k))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j + 1,n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j + 1,n
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + conjg(ap(k))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ctpmv
     !> ZTPMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix, supplied in packed form.

     pure subroutine la_ztpmv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: ap(*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('ZTPMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x:= a*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               k = kk
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk + j - 1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               k = kk
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk - n + j)
                           end if
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk - (n - (j + 1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk - n + j)
                           end if
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j - 1,1,-1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j - 1,1,-1
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(ap(k))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j + 1,n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j + 1,n
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + conjg(ap(k))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ztpmv
     !> WTPMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix, supplied in packed form.

     pure subroutine la_wtpmv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: ap(*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('WTPMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x:= a*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               k = kk
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk + j - 1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               k = kk
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk - n + j)
                           end if
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk,kk - (n - (j + 1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*ap(kk - n + j)
                           end if
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j - 1,1,-1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j - 1,1,-1
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(ap(k))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j + 1,n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j + 1,n
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + conjg(ap(k))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_wtpmv

     !> CTPSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix, supplied in packed form.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_ctpsv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: ap(*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('CTPSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/conjg(ap(kk + j - 1))
                           end if
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk,kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do k = kk,kk + j - 2
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(ap(kk + j - 1))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/conjg(ap(kk - n + j))
                           end if
                           x(j) = temp
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk,kk - (n - (j + 1)),-1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do k = kk,kk - (n - (j + 1)),-1
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(ap(kk - n + j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ctpsv
     !> ZTPSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix, supplied in packed form.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_ztpsv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: ap(*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('ZTPSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/conjg(ap(kk + j - 1))
                           end if
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk,kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do k = kk,kk + j - 2
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(ap(kk + j - 1))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/conjg(ap(kk - n + j))
                           end if
                           x(j) = temp
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk,kk - (n - (j + 1)),-1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do k = kk,kk - (n - (j + 1)),-1
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(ap(kk - n + j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ztpsv
     !> WTPSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix, supplied in packed form.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_wtpsv(uplo,trans,diag,n,ap,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: ap(*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,k,kk,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call la_xerbla('WTPSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do
                   end if
               else
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do
                           end if
                           jx = jx + incx
                           kk = kk + (n - j + 1)
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/conjg(ap(kk + j - 1))
                           end if
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk,kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do k = kk,kk + j - 2
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(ap(kk + j - 1))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do
                   end if
               else
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/conjg(ap(kk - n + j))
                           end if
                           x(j) = temp
                           kk = kk - (n - j + 1)
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk,kk - (n - (j + 1)),-1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do k = kk,kk - (n - (j + 1)),-1
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(ap(kk - n + j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n - j + 1)
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_wtpsv

     !> CTRMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix.

     pure subroutine la_ctrmv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('CTRMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j - 1,1,-1
                                   temp = temp + a(i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j - 1,1,-1
                                   temp = temp + conjg(a(i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + a(i,j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(a(i,j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j + 1,n
                                   temp = temp + a(i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j + 1,n
                                   temp = temp + conjg(a(i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + a(i,j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + conjg(a(i,j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ctrmv
     !> ZTRMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix.

     pure subroutine la_ztrmv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('ZTRMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j - 1,1,-1
                                   temp = temp + a(i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j - 1,1,-1
                                   temp = temp + conjg(a(i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + a(i,j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(a(i,j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j + 1,n
                                   temp = temp + a(i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j + 1,n
                                   temp = temp + conjg(a(i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + a(i,j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + conjg(a(i,j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ztrmv
     !> WTRMV:  performs one of the matrix-vector operations
     !> x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     !> where x is an n element vector and  A is an n by n unit, or non-unit,
     !> upper or lower triangular matrix.

     pure subroutine la_wtrmv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('WTRMV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := a*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               temp = x(j)
                               do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               temp = x(j)
                               do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               temp = x(jx)
                               ix = kx
                               do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x  or  x := a**h*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j - 1,1,-1
                                   temp = temp + a(i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j - 1,1,-1
                                   temp = temp + conjg(a(i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + a(i,j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(a(i,j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j + 1,n
                                   temp = temp + a(i,j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j + 1,n
                                   temp = temp + conjg(a(i,j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + a(i,j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + conjg(a(i,j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_wtrmv

     !> CTRSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_ctrsv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_sp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(sp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('CTRSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - a(i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ctrsv
     !> ZTRSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_ztrsv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_dp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(dp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('ZTRSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - a(i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_ztrsv
     !> WTRSV:  solves one of the systems of equations
     !> A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     !> where b and x are n element vectors and A is an n by n unit, or
     !> non-unit, upper or lower triangular matrix.
     !> No test for singularity or near-singularity is included in this
     !> routine. Such tests must be performed before calling this routine.

     pure subroutine la_wtrsv(uplo,trans,diag,n,a,lda,x,incx)
        use la_constants_qp
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,lda,n
           character,intent(in) :: diag,trans,uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================
           
           ! Local Scalars
           complex(qp) :: temp
           integer(ilp) :: i,info,ix,j,jx,kx
           logical(lk) :: noconj,nounit
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! test the input parameters.
           info = 0
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
               info = 1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') &
                     .and. .not. la_lsame(trans,'C')) then
               info = 2
           else if (.not. la_lsame(diag,'U') .and. .not. la_lsame(diag,'N')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1,n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call la_xerbla('WTRSV ',info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = la_lsame(trans,'T')
           nounit = la_lsame(diag,'N')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (la_lsame(trans,'N')) then
              ! form  x := inv( a )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = n,1,-1
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n,1,-1
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1,n
                           if (x(j) /= czero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           if (x(jx) /= czero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (la_lsame(uplo,'U')) then
                   if (incx == 1) then
                       do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - a(i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1,n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = 1,j - 1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n,1,-1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(j,j)
                           else
                               do i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
     end subroutine la_wtrsv

end module la_blas_level2_tri
