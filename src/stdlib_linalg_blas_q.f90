module stdlib_linalg_blas_q
     use stdlib_linalg_constants
     use stdlib_linalg_blas_aux
     use stdlib_linalg_blas_s
     implicit none(type, external)
     private

     public :: sp, dp, qp, lk, ilp
     public :: stdlib_qasum
     public :: stdlib_qaxpy
     public :: stdlib_qcopy
     public :: stdlib_qdot
     public :: stdlib_qgbmv
     public :: stdlib_qgemm
     public :: stdlib_qgemv
     public :: stdlib_qger
     public :: stdlib_qnrm2
     public :: stdlib_qrot
     public :: stdlib_qrotg
     public :: stdlib_qrotm
     public :: stdlib_qrotmg
     public :: stdlib_qsbmv
     public :: stdlib_qscal
     public :: stdlib_qsdot
     public :: stdlib_qspmv
     public :: stdlib_qspr
     public :: stdlib_qspr2
     public :: stdlib_qswap
     public :: stdlib_qsymm
     public :: stdlib_qsymv
     public :: stdlib_qsyr
     public :: stdlib_qsyr2
     public :: stdlib_qsyr2k
     public :: stdlib_qsyrk
     public :: stdlib_qtbmv
     public :: stdlib_qtbsv
     public :: stdlib_qtpmv
     public :: stdlib_qtpsv
     public :: stdlib_qtrmm
     public :: stdlib_qtrmv
     public :: stdlib_qtrsm
     public :: stdlib_qtrsv
     public :: stdlib_qzasum
     public :: stdlib_qznrm2

     ! 128-bit real constants
     real(qp), parameter, private :: zero = 0.00_qp
     real(qp), parameter, private :: half = 0.50_qp
     real(qp), parameter, private :: one = 1.00_qp
     real(qp), parameter, private :: two = 2.00_qp
     real(qp), parameter, private :: three = 3.00_qp
     real(qp), parameter, private :: four = 4.00_qp
     real(qp), parameter, private :: eight = 8.00_qp
     real(qp), parameter, private :: ten = 10.00_qp

     ! 128-bit complex constants
     complex(qp), parameter, private :: czero = (0.0_qp, 0.0_qp)
     complex(qp), parameter, private :: chalf = (0.5_qp, 0.0_qp)
     complex(qp), parameter, private :: cone = (1.0_qp, 0.0_qp)

     ! 128-bit scaling constants
     integer, parameter, private :: maxexp = maxexponent(zero)
     integer, parameter, private :: minexp = minexponent(zero)
     real(qp), parameter, private :: rradix = real(radix(zero), dp)
     real(qp), parameter, private :: ulp = epsilon(zero)
     real(qp), parameter, private :: eps = ulp*half
     real(qp), parameter, private :: safmin = rradix**max(minexp - 1, 1 - maxexp)
     real(qp), parameter, private :: safmax = one/safmin
     real(qp), parameter, private :: smlnum = safmin/ulp
     real(qp), parameter, private :: bignum = safmax*ulp
     real(qp), parameter, private :: rtmin = sqrt(smlnum)
     real(qp), parameter, private :: rtmax = sqrt(bignum)

     ! 128-bit Blue's scaling constants
     ! ssml>=1/s and sbig==1/S with s,S as defined in https://doi.org/10.1145/355769.355771
     real(qp), parameter, private :: tsml = rradix**ceiling((minexp - 1)*half)
     real(qp), parameter, private :: tbig = rradix**floor((maxexp - digits(zero) + 1)*half)
     real(qp), parameter, private :: ssml = rradix**(-floor((minexp - digits(zero))*half))
     real(qp), parameter, private :: sbig = rradix**(-ceiling((maxexp + digits(zero) - 1)*half))

     contains

     ! QASUM takes the sum of the absolute values.

     real(qp) function stdlib_qasum(n, dx, incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, n
           ! .. array arguments ..
           real(qp) :: dx(*)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: dtemp
           integer(ilp) :: i, m, mp1, nincx
           ! .. intrinsic functions ..
           intrinsic :: abs, mod
           stdlib_qasum = zero
           dtemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n, 6)
              if (m /= 0) then
                 do i = 1, m
                    dtemp = dtemp + abs(dx(i))
                 end do
                 if (n < 6) then
                    stdlib_qasum = dtemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1, n, 6
                 dtemp = dtemp + abs(dx(i)) + abs(dx(i + 1)) + abs(dx(i + 2)) + abs(dx(i + 3)) + abs(dx(i + &
                           4)) + abs(dx(i + 5))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1, nincx, incx
                 dtemp = dtemp + abs(dx(i))
              end do
           end if
           stdlib_qasum = dtemp
           return
           ! end of stdlib_qasum
     end function stdlib_qasum

     ! QAXPY constant times a vector plus a vector.
     ! uses unrolled loops for increments equal to one.

     subroutine stdlib_qaxpy(n, da, dx, incx, dy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: da
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(qp) :: dx(*), dy(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, ix, iy, m, mp1
           ! .. intrinsic functions ..
           intrinsic :: mod
           if (n <= 0) return
           if (da == zero) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n, 4)
              if (m /= 0) then
                 do i = 1, m
                    dy(i) = dy(i) + da*dx(i)
                 end do
              end if
              if (n < 4) return
              mp1 = m + 1
              do i = mp1, n, 4
                 dy(i) = dy(i) + da*dx(i)
                 dy(i + 1) = dy(i + 1) + da*dx(i + 1)
                 dy(i + 2) = dy(i + 2) + da*dx(i + 2)
                 dy(i + 3) = dy(i + 3) + da*dx(i + 3)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
               dy(iy) = dy(iy) + da*dx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_qaxpy
     end subroutine stdlib_qaxpy

     ! QCOPY copies a vector, x, to a vector, y.
     ! uses unrolled loops for increments equal to 1.

     subroutine stdlib_qcopy(n, dx, incx, dy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(qp) :: dx(*), dy(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, ix, iy, m, mp1
           ! .. intrinsic functions ..
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n, 7)
              if (m /= 0) then
                 do i = 1, m
                    dy(i) = dx(i)
                 end do
                 if (n < 7) return
              end if
              mp1 = m + 1
              do i = mp1, n, 7
                 dy(i) = dx(i)
                 dy(i + 1) = dx(i + 1)
                 dy(i + 2) = dx(i + 2)
                 dy(i + 3) = dx(i + 3)
                 dy(i + 4) = dx(i + 4)
                 dy(i + 5) = dx(i + 5)
                 dy(i + 6) = dx(i + 6)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 dy(iy) = dx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_qcopy
     end subroutine stdlib_qcopy

     ! QDOT forms the dot product of two vectors.
     ! uses unrolled loops for increments equal to one.

     real(qp) function stdlib_qdot(n, dx, incx, dy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(qp) :: dx(*), dy(*)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: dtemp
           integer(ilp) :: i, ix, iy, m, mp1
           ! .. intrinsic functions ..
           intrinsic :: mod
           stdlib_qdot = zero
           dtemp = zero
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              ! clean-up loop
              m = mod(n, 5)
              if (m /= 0) then
                 do i = 1, m
                    dtemp = dtemp + dx(i)*dy(i)
                 end do
                 if (n < 5) then
                    stdlib_qdot = dtemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1, n, 5
               dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + &
                         dx(i + 4)*dy(i + 4)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 dtemp = dtemp + dx(ix)*dy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_qdot = dtemp
           return
           ! end of stdlib_qdot
     end function stdlib_qdot

     ! QGBMV  performs one of the matrix-vector operations
     ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
     ! where alpha and beta are scalars, x and y are vectors and A is an
     ! m by n band matrix, with kl sub-diagonals and ku super-diagonals.

     subroutine stdlib_qgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: incx, incy, kl, ku, lda, m, n
           character :: trans
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, iy, j, jx, jy, k, kup1, kx, ky, lenx, leny
           ! .. intrinsic functions ..
           intrinsic :: max, min
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 1
           else if (m < 0) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (kl < 0) then
               info = 4
           else if (ku < 0) then
               info = 5
           else if (lda < (kl + ku + 1)) then
               info = 8
           else if (incx == 0) then
               info = 10
           else if (incy == 0) then
               info = 13
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qgbmv ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (stdlib_lsame(trans, 'n')) then
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
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the band part of a.
           ! first form  y := beta*y.
           if (beta /= one) then
               if (incy == 1) then
                   if (beta == zero) then
                       do i = 1, leny
                           y(i) = zero
                       end do
                   else
                       do i = 1, leny
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == zero) then
                       do i = 1, leny
                           y(iy) = zero
                           iy = iy + incy
                       end do
                   else
                       do i = 1, leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do
                   end if
               end if
           end if
           if (alpha == zero) return
           kup1 = ku + 1
           if (stdlib_lsame(trans, 'n')) then
              ! form  y := alpha*a*x + y.
               jx = kx
               if (incy == 1) then
                   do j = 1, n
                       temp = alpha*x(jx)
                       k = kup1 - j
                       do i = max(1, j - ku), min(m, j + kl)
                           y(i) = y(i) + temp*a(k + i, j)
                       end do
                       jx = jx + incx
                   end do
               else
                   do j = 1, n
                       temp = alpha*x(jx)
                       iy = ky
                       k = kup1 - j
                       do i = max(1, j - ku), min(m, j + kl)
                           y(iy) = y(iy) + temp*a(k + i, j)
                           iy = iy + incy
                       end do
                       jx = jx + incx
                       if (j > ku) ky = ky + incy
                   end do
               end if
           else
              ! form  y := alpha*a**t*x + y.
               jy = ky
               if (incx == 1) then
                   do j = 1, n
                       temp = zero
                       k = kup1 - j
                       do i = max(1, j - ku), min(m, j + kl)
                           temp = temp + a(k + i, j)*x(i)
                       end do
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do
               else
                   do j = 1, n
                       temp = zero
                       ix = kx
                       k = kup1 - j
                       do i = max(1, j - ku), min(m, j + kl)
                           temp = temp + a(k + i, j)*x(ix)
                           ix = ix + incx
                       end do
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j > ku) kx = kx + incx
                   end do
               end if
           end if
           return
           ! end of stdlib_qgbmv
     end subroutine stdlib_qgbmv

     ! QGEMM  performs one of the matrix-matrix operations
     ! C := alpha*op( A )*op( B ) + beta*C,
     ! where  op( X ) is one of
     ! op( X ) = X   or   op( X ) = X**T,
     ! alpha and beta are scalars, and A, B and C are matrices, with op( A )
     ! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     subroutine stdlib_qgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldb, ldc, m, n
           character :: transa, transb
           ! .. array arguments ..
           real(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, j, l, nrowa, nrowb
           logical(lk) :: nota, notb

           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
           nota = stdlib_lsame(transa, 'n')
           notb = stdlib_lsame(transb, 'n')
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
           if ((.not. nota) .and. (.not. stdlib_lsame(transa, 'c')) .and. (.not. stdlib_lsame(transa, &
                     't'))) then
               info = 1
           else if ((.not. notb) .and. (.not. stdlib_lsame(transb, 'c')) .and. (.not. stdlib_lsame( &
                     transb, 't'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (k < 0) then
               info = 5
           else if (lda < max(1, nrowa)) then
               info = 8
           else if (ldb < max(1, nrowb)) then
               info = 10
           else if (ldc < max(1, m)) then
               info = 13
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qgemm ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) &
                     return
           ! and if  alpha.eq.zero.
           if (alpha == zero) then
               if (beta == zero) then
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = zero
                       end do
                   end do
               else
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = beta*c(i, j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (notb) then
               if (nota) then
                 ! form  c := alpha*a*b + beta*c.
                   do j = 1, n
                       if (beta == zero) then
                           do i = 1, m
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1, m
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           temp = alpha*b(l, j)
                           do i = 1, m
                               c(i, j) = c(i, j) + temp*a(i, l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1, n
                       do i = 1, m
                           temp = zero
                           do l = 1, k
                               temp = temp + a(l, i)*b(l, j)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           else
               if (nota) then
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1, n
                       if (beta == zero) then
                           do i = 1, m
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1, m
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           temp = alpha*b(j, l)
                           do i = 1, m
                               c(i, j) = c(i, j) + temp*a(i, l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1, n
                       do i = 1, m
                           temp = zero
                           do l = 1, k
                               temp = temp + a(l, i)*b(j, l)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_qgemm
     end subroutine stdlib_qgemm

     ! QGEMV  performs one of the matrix-vector operations
     ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
     ! where alpha and beta are scalars, x and y are vectors and A is an
     ! m by n matrix.

     subroutine stdlib_qgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: incx, incy, lda, m, n
           character :: trans
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 1
           else if (m < 0) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (lda < max(1, m)) then
               info = 6
           else if (incx == 0) then
               info = 8
           else if (incy == 0) then
               info = 11
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qgemv ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
           if (stdlib_lsame(trans, 'n')) then
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
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           ! first form  y := beta*y.
           if (beta /= one) then
               if (incy == 1) then
                   if (beta == zero) then
                       do i = 1, leny
                           y(i) = zero
                       end do
                   else
                       do i = 1, leny
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == zero) then
                       do i = 1, leny
                           y(iy) = zero
                           iy = iy + incy
                       end do
                   else
                       do i = 1, leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do
                   end if
               end if
           end if
           if (alpha == zero) return
           if (stdlib_lsame(trans, 'n')) then
              ! form  y := alpha*a*x + y.
               jx = kx
               if (incy == 1) then
                   do j = 1, n
                       temp = alpha*x(jx)
                       do i = 1, m
                           y(i) = y(i) + temp*a(i, j)
                       end do
                       jx = jx + incx
                   end do
               else
                   do j = 1, n
                       temp = alpha*x(jx)
                       iy = ky
                       do i = 1, m
                           y(iy) = y(iy) + temp*a(i, j)
                           iy = iy + incy
                       end do
                       jx = jx + incx
                   end do
               end if
           else
              ! form  y := alpha*a**t*x + y.
               jy = ky
               if (incx == 1) then
                   do j = 1, n
                       temp = zero
                       do i = 1, m
                           temp = temp + a(i, j)*x(i)
                       end do
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do
               else
                   do j = 1, n
                       temp = zero
                       ix = kx
                       do i = 1, m
                           temp = temp + a(i, j)*x(ix)
                           ix = ix + incx
                       end do
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_qgemv
     end subroutine stdlib_qgemv

     ! QGER   performs the rank 1 operation
     ! A := alpha*x*y**T + A,
     ! where alpha is a scalar, x is an m element vector, y is an n element
     ! vector and A is an m by n matrix.

     subroutine stdlib_qger(m, n, alpha, x, incx, y, incy, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, incy, lda, m, n
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jy, kx
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (m < 0) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (incx == 0) then
               info = 5
           else if (incy == 0) then
               info = 7
           else if (lda < max(1, m)) then
               info = 9
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qger  ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (alpha == zero)) return
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (incy > 0) then
               jy = 1
           else
               jy = 1 - (n - 1)*incy
           end if
           if (incx == 1) then
               do j = 1, n
                   if (y(jy) /= zero) then
                       temp = alpha*y(jy)
                       do i = 1, m
                           a(i, j) = a(i, j) + x(i)*temp
                       end do
                   end if
                   jy = jy + incy
               end do
           else
               if (incx > 0) then
                   kx = 1
               else
                   kx = 1 - (m - 1)*incx
               end if
               do j = 1, n
                   if (y(jy) /= zero) then
                       temp = alpha*y(jy)
                       ix = kx
                       do i = 1, m
                           a(i, j) = a(i, j) + x(ix)*temp
                           ix = ix + incx
                       end do
                   end if
                   jy = jy + incy
               end do
           end if
           return
           ! end of stdlib_qger
     end subroutine stdlib_qger

     ! !
     ! QNRM2 returns the euclidean norm of a vector via the function
     ! name, so that
     ! QNRM2 := sqrt( x'*x )

     function stdlib_qnrm2(n, x, incx)
        integer, parameter :: wp = kind(1._qp)
        real(qp) :: stdlib_qnrm2
        ! -- reference blas level1 routine (version 3.9.1_qp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! .. constants ..
        real(qp), parameter :: zero = 0.0_qp
        real(qp), parameter :: one = 1.0_qp
        real(qp), parameter :: maxn = huge(0.0_qp)
        ! .. blue's scaling constants ..
     real(qp), parameter :: tsml = real(radix(0._qp), wp)**ceiling((minexponent(0._qp) - 1) &
                *0.5_qp)
     real(qp), parameter :: tbig = real(radix(0._qp), wp)**floor((maxexponent(0._qp) - &
               digits(0._qp) + 1)*0.5_qp)
     real(qp), parameter :: ssml = real(radix(0._qp), wp)**(-floor((minexponent(0._qp) - &
               digits(0._qp))*0.5_qp))
     real(qp), parameter :: sbig = real(radix(0._qp), wp)**(-ceiling((maxexponent(0._qp) &
               + digits(0._qp) - 1)*0.5_qp))
        ! .. scalar arguments ..
     integer(ilp) :: incx, n
        ! .. array arguments ..
        real(qp) :: x(*)
        ! .. local scalars ..
     integer(ilp) :: i, ix
     logical(lk) :: notbig
        real(qp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
        ! quick return if possible
        stdlib_qnrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1, n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        stdlib_qnrm2 = scl*sqrt(sumsq)
        return
     end function stdlib_qnrm2

     ! QROT applies a plane rotation.

     subroutine stdlib_qrot(n, dx, incx, dy, incy, c, s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: c, s
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(qp) :: dx(*), dy(*)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: dtemp
           integer(ilp) :: i, ix, iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1, n
                 dtemp = c*dx(i) + s*dy(i)
                 dy(i) = c*dy(i) - s*dx(i)
                 dx(i) = dtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 dtemp = c*dx(ix) + s*dy(iy)
                 dy(iy) = c*dy(iy) - s*dx(ix)
                 dx(ix) = dtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_qrot
     end subroutine stdlib_qrot

     ! !
     ! The computation uses the formulas
     ! sigma = sgn(a)    if |a| >  |b|
     ! = sgn(b)    if |b| >= |a|
     ! r = sigma*sqrt( a**2 + b**2 )
     ! c = 1; s = 0      if r = 0
     ! c = a/r; s = b/r  if r != 0
     ! The subroutine also computes
     ! z = s    if |a| > |b|,
     ! = 1/c  if |b| >= |a| and c != 0
     ! = 1    if c = 0
     ! This allows c and s to be reconstructed from z as follows:
     ! If z = 1, set c = 0, s = 1.
     ! If |z| < 1, set c = sqrt(1 - z**2) and s = z.
     ! If |z| > 1, set c = 1/z and s = sqrt( 1 - c**2).

     subroutine stdlib_qrotg(a, b, c, s)
        integer, parameter :: wp = kind(1._qp)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! .. constants ..
        real(qp), parameter :: zero = 0.0_qp
        real(qp), parameter :: one = 1.0_qp
        ! .. scaling constants ..
     real(qp), parameter :: safmin = real(radix(0._qp), wp)**max(minexponent(0._qp) - 1, 1 - &
               maxexponent(0._qp))
     real(qp), parameter :: safmax = real(radix(0._qp), wp)**max(1 - minexponent(0._qp), maxexponent( &
               0._qp) - 1)
        ! .. scalar arguments ..
        real(qp) :: a, b, c, s
        ! .. local scalars ..
        real(qp) :: anorm, bnorm, scl, sigma, r, z
        anorm = abs(a)
        bnorm = abs(b)
        if (bnorm == zero) then
           c = one
           s = zero
           b = zero
        else if (anorm == zero) then
           c = zero
           s = one
           a = b
           b = one
        else
           scl = min(safmax, max(safmin, anorm, bnorm))
           if (anorm > bnorm) then
              sigma = sign(one, a)
           else
              sigma = sign(one, b)
           end if
           r = sigma*(scl*sqrt((a/scl)**2 + (b/scl)**2))
           c = a/r
           s = b/r
           if (anorm > bnorm) then
              z = s
           else if (c /= zero) then
              z = one/c
           else
              z = one
           end if
           a = r
           b = z
        end if
        return
     end subroutine stdlib_qrotg

     ! APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
     ! (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
     ! (DY**T)
     ! QX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX >= 0, ELSE
     ! LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
     ! WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     ! QFLAG=-1._qp     DFLAG=0._qp        DFLAG=1._qp     DFLAG=-2.D0
     ! (DH11  DH12)    (1._qp  DH12)    (DH11  1._qp)    (1._qp  0._qp)
     ! H=(          )    (          )    (          )    (          )
     ! (DH21  DH22),   (DH21  1._qp),   (-1._qp DH22),   (0._qp  1._qp).
     ! SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.

     subroutine stdlib_qrotm(n, dx, incx, dy, incy, dparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(qp) :: dparam(5), dx(*), dy(*)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: dflag, dh11, dh12, dh21, dh22, two, w, z, zero
           integer(ilp) :: i, kx, ky, nsteps
           ! .. data statements ..
           data zero, two/0._qp, 2._qp/
           dflag = dparam(1)
           if (n <= 0 .or. (dflag + two == zero)) return
           if (incx == incy .and. incx > 0) then
              nsteps = n*incx
              if (dflag < zero) then
                 dh11 = dparam(2)
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 dh22 = dparam(5)
                 do i = 1, nsteps, incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w*dh11 + z*dh12
                    dy(i) = w*dh21 + z*dh22
                 end do
              else if (dflag == zero) then
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 do i = 1, nsteps, incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w + z*dh12
                    dy(i) = w*dh21 + z
                 end do
              else
                 dh11 = dparam(2)
                 dh22 = dparam(5)
                 do i = 1, nsteps, incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w*dh11 + z
                    dy(i) = -w + dh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              if (dflag < zero) then
                 dh11 = dparam(2)
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 dh22 = dparam(5)
                 do i = 1, n
                    w = dx(kx)
                    z = dy(ky)
                    dx(kx) = w*dh11 + z*dh12
                    dy(ky) = w*dh21 + z*dh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (dflag == zero) then
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 do i = 1, n
                    w = dx(kx)
                    z = dy(ky)
                    dx(kx) = w + z*dh12
                    dy(ky) = w*dh21 + z
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else
                  dh11 = dparam(2)
                  dh22 = dparam(5)
                  do i = 1, n
                     w = dx(kx)
                     z = dy(ky)
                     dx(kx) = w*dh11 + z
                     dy(ky) = -w + dh22*z
                     kx = kx + incx
                     ky = ky + incy
                 end do
              end if
           end if
           return
           ! end of stdlib_qrotm
     end subroutine stdlib_qrotm

     ! CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
     ! THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(DD1)*DX1,SQRT(DD2)    DY2)**T.
     ! WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
     ! QFLAG=-1._qp     DFLAG=0._qp        DFLAG=1._qp     DFLAG=-2.D0
     ! (DH11  DH12)    (1._qp  DH12)    (DH11  1._qp)    (1._qp  0._qp)
     ! H=(          )    (          )    (          )    (          )
     ! (DH21  DH22),   (DH21  1._qp),   (-1._qp DH22),   (0._qp  1._qp).
     ! LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
     ! RESPECTIVELY. (VALUES OF 1._qp, -1._qp, OR 0._qp IMPLIED BY THE
     ! VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
     ! THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
     ! INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
     ! OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.

     subroutine stdlib_qrotmg(dd1, dd2, dx1, dy1, dparam)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: dd1, dd2, dx1, dy1
           ! .. array arguments ..
           real(qp) :: dparam(5)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: dflag, dh11, dh12, dh21, dh22, dp1, dp2, dq1, dq2, dtemp, du, gam, gamsq, one, rgamsq, &
                     two, zero
           ! .. intrinsic functions ..
           intrinsic :: abs
           ! .. data statements ..
           data zero, one, two/0._qp, 1._qp, 2._qp/
           data gam, gamsq, rgamsq/4096._qp, 16777216._qp, 5.9604645e-8_qp/
           if (dd1 < zero) then
              ! go zero-h-d-and-dx1..
              dflag = -one
              dh11 = zero
              dh12 = zero
              dh21 = zero
              dh22 = zero
              dd1 = zero
              dd2 = zero
              dx1 = zero
           else
              ! case-dd1-nonnegative
              dp2 = dd2*dy1
              if (dp2 == zero) then
                 dflag = -two
                 dparam(1) = dflag
                 return
              end if
              ! regular-case..
              dp1 = dd1*dx1
              dq2 = dp2*dy1
              dq1 = dp1*dx1
              if (abs(dq1) > abs(dq2)) then
                 dh21 = -dy1/dx1
                 dh12 = dp2/dp1
                 du = one - dh12*dh21
                if (du > zero) then
                  dflag = zero
                  dd1 = dd1/du
                  dd2 = dd2/du
                  dx1 = dx1*du
                else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  dflag = -one
                  dh11 = zero
                  dh12 = zero
                  dh21 = zero
                  dh22 = zero
                  dd1 = zero
                  dd2 = zero
                  dx1 = zero
                end if
              else
                 if (dq2 < zero) then
                    ! go zero-h-d-and-dx1..
                    dflag = -one
                    dh11 = zero
                    dh12 = zero
                    dh21 = zero
                    dh22 = zero
                    dd1 = zero
                    dd2 = zero
                    dx1 = zero
                 else
                    dflag = one
                    dh11 = dp1/dp2
                    dh22 = dx1/dy1
                    du = one + dh11*dh22
                    dtemp = dd2/du
                    dd2 = dd1/du
                    dd1 = dtemp
                    dx1 = dy1*du
                 end if
              end if
           ! procedure..scale-check
              if (dd1 /= zero) then
                 do while ((dd1 <= rgamsq) .or. (dd1 >= gamsq))
                    if (dflag == zero) then
                       dh11 = one
                       dh22 = one
                       dflag = -one
                    else
                       dh21 = -one
                       dh12 = one
                       dflag = -one
                    end if
                    if (dd1 <= rgamsq) then
                       dd1 = dd1*gam**2
                       dx1 = dx1/gam
                       dh11 = dh11/gam
                       dh12 = dh12/gam
                    else
                       dd1 = dd1/gam**2
                       dx1 = dx1*gam
                       dh11 = dh11*gam
                       dh12 = dh12*gam
                    end if
                 end do
              end if
              if (dd2 /= zero) then
                 do while ((abs(dd2) <= rgamsq) .or. (abs(dd2) >= gamsq))
                    if (dflag == zero) then
                       dh11 = one
                       dh22 = one
                       dflag = -one
                    else
                       dh21 = -one
                       dh12 = one
                       dflag = -one
                    end if
                    if (abs(dd2) <= rgamsq) then
                       dd2 = dd2*gam**2
                       dh21 = dh21/gam
                       dh22 = dh22/gam
                    else
                       dd2 = dd2/gam**2
                       dh21 = dh21*gam
                       dh22 = dh22*gam
                    end if
                 end do
              end if
           end if
           if (dflag < zero) then
              dparam(2) = dh11
              dparam(3) = dh21
              dparam(4) = dh12
              dparam(5) = dh22
           else if (dflag == zero) then
              dparam(3) = dh21
              dparam(4) = dh12
           else
              dparam(2) = dh11
              dparam(5) = dh22
           end if
           dparam(1) = dflag
           return
           ! end of stdlib_qrotmg
     end subroutine stdlib_qrotmg

     ! QSBMV  performs the matrix-vector  operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n symmetric band matrix, with k super-diagonals.

     subroutine stdlib_qsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: incx, incy, k, lda, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
           ! .. intrinsic functions ..
           intrinsic :: max, min
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (k < 0) then
               info = 3
           else if (lda < (k + 1)) then
               info = 6
           else if (incx == 0) then
               info = 8
           else if (incy == 0) then
               info = 11
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsbmv ', info)
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
           ! start the operations. in this version the elements of the array a
           ! are accessed sequentially with one pass through a.
           ! first form  y := beta*y.
           if (beta /= one) then
               if (incy == 1) then
                   if (beta == zero) then
                       do i = 1, n
                           y(i) = zero
                       end do
                   else
                       do i = 1, n
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == zero) then
                       do i = 1, n
                           y(iy) = zero
                           iy = iy + incy
                       end do
                   else
                       do i = 1, n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do
                   end if
               end if
           end if
           if (alpha == zero) return
           if (stdlib_lsame(uplo, 'u')) then
              ! form  y  when upper triangle of a is stored.
               kplus1 = k + 1
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       l = kplus1 - j
                       do i = max(1, j - k), j - 1
                           y(i) = y(i) + temp1*a(l + i, j)
                           temp2 = temp2 + a(l + i, j)*x(i)
                       end do
                       y(j) = y(j) + temp1*a(kplus1, j) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       do i = max(1, j - k), j - 1
                           y(iy) = y(iy) + temp1*a(l + i, j)
                           temp2 = temp2 + a(l + i, j)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do
                       y(jy) = y(jy) + temp1*a(kplus1, j) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       if (j > k) then
                           kx = kx + incx
                           ky = ky + incy
                       end if
                   end do
               end if
           else
              ! form  y  when lower triangle of a is stored.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*a(1, j)
                       l = 1 - j
                       do i = j + 1, min(n, j + k)
                           y(i) = y(i) + temp1*a(l + i, j)
                           temp2 = temp2 + a(l + i, j)*x(i)
                       end do
                       y(j) = y(j) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*a(1, j)
                       l = 1 - j
                       ix = jx
                       iy = jy
                       do i = j + 1, min(n, j + k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l + i, j)
                           temp2 = temp2 + a(l + i, j)*x(ix)
                       end do
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_qsbmv
     end subroutine stdlib_qsbmv

     ! QSCAL scales a vector by a constant.
     ! uses unrolled loops for increment equal to 1.

     subroutine stdlib_qscal(n, da, dx, incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: da
           integer(ilp) :: incx, n
           ! .. array arguments ..
           real(qp) :: dx(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, m, mp1, nincx
           ! .. intrinsic functions ..
           intrinsic :: mod
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              ! clean-up loop
              m = mod(n, 5)
              if (m /= 0) then
                 do i = 1, m
                    dx(i) = da*dx(i)
                 end do
                 if (n < 5) return
              end if
              mp1 = m + 1
              do i = mp1, n, 5
                 dx(i) = da*dx(i)
                 dx(i + 1) = da*dx(i + 1)
                 dx(i + 2) = da*dx(i + 2)
                 dx(i + 3) = da*dx(i + 3)
                 dx(i + 4) = da*dx(i + 4)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1, nincx, incx
                 dx(i) = da*dx(i)
              end do
           end if
           return
           ! end of stdlib_qscal
     end subroutine stdlib_qscal

     ! Compute the inner product of two vectors with extended
     ! precision accumulation and result.
     ! Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
     ! QSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
     ! where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
     ! defined in a similar way using INCY.

     real(qp) function stdlib_qsdot(n, sx, incx, sy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(sp) :: sx(*), sy(*)
        ! authors:
        ! ========
        ! lawson, c. l., (jpl), hanson, r. j., (snla),
        ! kincaid, d. r., (u. of texas), krogh, f. t., (jpl)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, kx, ky, ns
           ! .. intrinsic functions ..
           intrinsic :: real
           stdlib_qsdot = zero
           if (n <= 0) return
           if (incx == incy .and. incx > 0) then
           ! code for equal, positive, non-unit increments.
              ns = n*incx
              do i = 1, ns, incx
                 stdlib_qsdot = stdlib_qsdot + real(sx(i), KIND=qp)*real(sy(i), KIND=qp)
              end do
           else
           ! code for unequal or nonpositive increments.
              kx = 1
              ky = 1
              if (incx < 0) kx = 1 + (1 - n)*incx
              if (incy < 0) ky = 1 + (1 - n)*incy
              do i = 1, n
                 stdlib_qsdot = stdlib_qsdot + real(sx(kx), KIND=qp)*real(sy(ky), KIND=qp)
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           return
           ! end of stdlib_qsdot
     end function stdlib_qsdot

     ! QSPMV  performs the matrix-vector operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n symmetric matrix, supplied in packed form.

     subroutine stdlib_qspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: incx, incy, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: ap(*), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, k, kk, kx, ky
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (incx == 0) then
               info = 6
           else if (incy == 0) then
               info = 9
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qspmv ', info)
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
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with one pass through ap.
           ! first form  y := beta*y.
           if (beta /= one) then
               if (incy == 1) then
                   if (beta == zero) then
                       do i = 1, n
                           y(i) = zero
                       end do
                   else
                       do i = 1, n
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == zero) then
                       do i = 1, n
                           y(iy) = zero
                           iy = iy + incy
                       end do
                   else
                       do i = 1, n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do
                   end if
               end if
           end if
           if (alpha == zero) return
           kk = 1
           if (stdlib_lsame(uplo, 'u')) then
              ! form  y  when ap contains the upper triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       k = kk
                       do i = 1, j - 1
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
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do k = kk, kk + j - 2
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
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*ap(kk)
                       k = kk + 1
                       do i = j + 1, n
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
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*ap(kk)
                       ix = jx
                       iy = jy
                       do k = kk + 1, kk + n - j
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
           ! end of stdlib_qspmv
     end subroutine stdlib_qspmv

     ! QSPR    performs the symmetric rank 1 operation
     ! A := alpha*x*x**T + A,
     ! where alpha is a real scalar, x is an n element vector and A is an
     ! n by n symmetric matrix, supplied in packed form.

     subroutine stdlib_qspr(uplo, n, alpha, x, incx, ap)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: ap(*), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, k, kk, kx
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (incx == 0) then
               info = 5
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qspr  ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == zero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with one pass through ap.
           kk = 1
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when upper triangle is stored in ap.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= zero) then
                           temp = alpha*x(j)
                           k = kk
                           do i = 1, j
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
                           end do
                       end if
                       kk = kk + j
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= zero) then
                           temp = alpha*x(jx)
                           ix = kx
                           do k = kk, kk + j - 1
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
                           end do
                       end if
                       jx = jx + incx
                       kk = kk + j
                   end do
               end if
           else
              ! form  a  when lower triangle is stored in ap.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= zero) then
                           temp = alpha*x(j)
                           k = kk
                           do i = j, n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
                           end do
                       end if
                       kk = kk + n - j + 1
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= zero) then
                           temp = alpha*x(jx)
                           ix = jx
                           do k = kk, kk + n - j
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
                           end do
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
                   end do
               end if
           end if
           return
           ! end of stdlib_qspr
     end subroutine stdlib_qspr

     ! QSPR2  performs the symmetric rank 2 operation
     ! A := alpha*x*y**T + alpha*y*x**T + A,
     ! where alpha is a scalar, x and y are n element vectors and A is an
     ! n by n symmetric matrix, supplied in packed form.

     subroutine stdlib_qspr2(uplo, n, alpha, x, incx, y, incy, ap)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, incy, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: ap(*), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, k, kk, kx, ky
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (incx == 0) then
               info = 5
           else if (incy == 0) then
               info = 7
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qspr2 ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == zero)) return
           ! set up the start points in x and y if the increments are not both
           ! unity.
           if ((incx /= 1) .or. (incy /= 1)) then
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
               jx = kx
               jy = ky
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with one pass through ap.
           kk = 1
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when upper triangle is stored in ap.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= zero) .or. (y(j) /= zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           k = kk
                           do i = 1, j
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
                           end do
                       end if
                       kk = kk + j
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = kx
                           iy = ky
                           do k = kk, kk + j - 1
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
                   end do
               end if
           else
              ! form  a  when lower triangle is stored in ap.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= zero) .or. (y(j) /= zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           k = kk
                           do i = j, n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
                           end do
                       end if
                       kk = kk + n - j + 1
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = jx
                           iy = jy
                           do k = kk, kk + n - j
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
                   end do
               end if
           end if
           return
           ! end of stdlib_qspr2
     end subroutine stdlib_qspr2

     ! QSWAP interchanges two vectors.
     ! uses unrolled loops for increments equal to 1.

     subroutine stdlib_qswap(n, dx, incx, dy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           real(qp) :: dx(*), dy(*)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: dtemp
           integer(ilp) :: i, ix, iy, m, mp1
           ! .. intrinsic functions ..
           intrinsic :: mod
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
             ! clean-up loop
              m = mod(n, 3)
              if (m /= 0) then
                 do i = 1, m
                    dtemp = dx(i)
                    dx(i) = dy(i)
                    dy(i) = dtemp
                 end do
                 if (n < 3) return
              end if
              mp1 = m + 1
              do i = mp1, n, 3
                 dtemp = dx(i)
                 dx(i) = dy(i)
                 dy(i) = dtemp
                 dtemp = dx(i + 1)
                 dx(i + 1) = dy(i + 1)
                 dy(i + 1) = dtemp
                 dtemp = dx(i + 2)
                 dx(i + 2) = dy(i + 2)
                 dy(i + 2) = dtemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 dtemp = dx(ix)
                 dx(ix) = dy(iy)
                 dy(iy) = dtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_qswap
     end subroutine stdlib_qswap

     ! QSYMM  performs one of the matrix-matrix operations
     ! C := alpha*A*B + beta*C,
     ! or
     ! C := alpha*B*A + beta*C,
     ! where alpha and beta are scalars,  A is a symmetric matrix and  B and
     ! C are  m by n matrices.

     subroutine stdlib_qsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: lda, ldb, ldc, m, n
           character :: side, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, j, k, nrowa
           logical(lk) :: upper

           ! set nrowa as the number of rows of a.
           if (stdlib_lsame(side, 'l')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = stdlib_lsame(uplo, 'u')
           ! test the input parameters.
           info = 0
           if ((.not. stdlib_lsame(side, 'l')) .and. (.not. stdlib_lsame(side, 'r'))) then
               info = 1
           else if ((.not. upper) .and. (.not. stdlib_lsame(uplo, 'l'))) then
               info = 2
           else if (m < 0) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1, nrowa)) then
               info = 7
           else if (ldb < max(1, m)) then
               info = 9
           else if (ldc < max(1, m)) then
               info = 12
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsymm ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == zero) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (beta == zero) then
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = zero
                       end do
                   end do
               else
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = beta*c(i, j)
                       end do
                   end do
               end if
               return
           end if
           ! start the operations.
           if (stdlib_lsame(side, 'l')) then
              ! form  c := alpha*a*b + beta*c.
               if (upper) then
                   do j = 1, n
                       do i = 1, m
                           temp1 = alpha*b(i, j)
                           temp2 = zero
                           do k = 1, i - 1
                               c(k, j) = c(k, j) + temp1*a(k, i)
                               temp2 = temp2 + b(k, j)*a(k, i)
                           end do
                           if (beta == zero) then
                               c(i, j) = temp1*a(i, i) + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + temp1*a(i, i) + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = m, 1, -1
                           temp1 = alpha*b(i, j)
                           temp2 = zero
                           do k = i + 1, m
                               c(k, j) = c(k, j) + temp1*a(k, i)
                               temp2 = temp2 + b(k, j)*a(k, i)
                           end do
                           if (beta == zero) then
                               c(i, j) = temp1*a(i, i) + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + temp1*a(i, i) + alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1, n
                   temp1 = alpha*a(j, j)
                   if (beta == zero) then
                       do i = 1, m
                           c(i, j) = temp1*b(i, j)
                       end do
                   else
                       do i = 1, m
                           c(i, j) = beta*c(i, j) + temp1*b(i, j)
                       end do
                   end if
                   do k = 1, j - 1
                       if (upper) then
                           temp1 = alpha*a(k, j)
                       else
                           temp1 = alpha*a(j, k)
                       end if
                       do i = 1, m
                           c(i, j) = c(i, j) + temp1*b(i, k)
                       end do
                   end do
                   do k = j + 1, n
                       if (upper) then
                           temp1 = alpha*a(j, k)
                       else
                           temp1 = alpha*a(k, j)
                       end if
                       do i = 1, m
                           c(i, j) = c(i, j) + temp1*b(i, k)
                       end do
                   end do
               end do loop_170
           end if
           return
           ! end of stdlib_qsymm
     end subroutine stdlib_qsymm

     ! QSYMV  performs the matrix-vector  operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n symmetric matrix.

     subroutine stdlib_qsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: incx, incy, lda, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kx, ky
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (lda < max(1, n)) then
               info = 5
           else if (incx == 0) then
               info = 7
           else if (incy == 0) then
               info = 10
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsymv ', info)
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
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the triangular part
           ! of a.
           ! first form  y := beta*y.
           if (beta /= one) then
               if (incy == 1) then
                   if (beta == zero) then
                       do i = 1, n
                           y(i) = zero
                       end do
                   else
                       do i = 1, n
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == zero) then
                       do i = 1, n
                           y(iy) = zero
                           iy = iy + incy
                       end do
                   else
                       do i = 1, n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do
                   end if
               end if
           end if
           if (alpha == zero) return
           if (stdlib_lsame(uplo, 'u')) then
              ! form  y  when a is stored in upper triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       do i = 1, j - 1
                           y(i) = y(i) + temp1*a(i, j)
                           temp2 = temp2 + a(i, j)*x(i)
                       end do
                       y(j) = y(j) + temp1*a(j, j) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do i = 1, j - 1
                           y(iy) = y(iy) + temp1*a(i, j)
                           temp2 = temp2 + a(i, j)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do
                       y(jy) = y(jy) + temp1*a(j, j) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           else
              ! form  y  when a is stored in lower triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*a(j, j)
                       do i = j + 1, n
                           y(i) = y(i) + temp1*a(i, j)
                           temp2 = temp2 + a(i, j)*x(i)
                       end do
                       y(j) = y(j) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*a(j, j)
                       ix = jx
                       iy = jy
                       do i = j + 1, n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i, j)
                           temp2 = temp2 + a(i, j)*x(ix)
                       end do
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_qsymv
     end subroutine stdlib_qsymv

     ! QSYR   performs the symmetric rank 1 operation
     ! A := alpha*x*x**T + A,
     ! where alpha is a real scalar, x is an n element vector and A is an
     ! n by n symmetric matrix.

     subroutine stdlib_qsyr(uplo, n, alpha, x, incx, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, lda, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kx
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (incx == 0) then
               info = 5
           else if (lda < max(1, n)) then
               info = 7
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsyr  ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == zero)) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the triangular part
           ! of a.
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when a is stored in upper triangle.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= zero) then
                           temp = alpha*x(j)
                           do i = 1, j
                               a(i, j) = a(i, j) + x(i)*temp
                           end do
                       end if
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= zero) then
                           temp = alpha*x(jx)
                           ix = kx
                           do i = 1, j
                               a(i, j) = a(i, j) + x(ix)*temp
                               ix = ix + incx
                           end do
                       end if
                       jx = jx + incx
                   end do
               end if
           else
              ! form  a  when a is stored in lower triangle.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= zero) then
                           temp = alpha*x(j)
                           do i = j, n
                               a(i, j) = a(i, j) + x(i)*temp
                           end do
                       end if
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= zero) then
                           temp = alpha*x(jx)
                           ix = jx
                           do i = j, n
                               a(i, j) = a(i, j) + x(ix)*temp
                               ix = ix + incx
                           end do
                       end if
                       jx = jx + incx
                   end do
               end if
           end if
           return
           ! end of stdlib_qsyr
     end subroutine stdlib_qsyr

     ! QSYR2  performs the symmetric rank 2 operation
     ! A := alpha*x*y**T + alpha*y*x**T + A,
     ! where alpha is a scalar, x and y are n element vectors and A is an n
     ! by n symmetric matrix.

     subroutine stdlib_qsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, incy, lda, n
           character :: uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kx, ky
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (n < 0) then
               info = 2
           else if (incx == 0) then
               info = 5
           else if (incy == 0) then
               info = 7
           else if (lda < max(1, n)) then
               info = 9
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsyr2 ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == zero)) return
           ! set up the start points in x and y if the increments are not both
           ! unity.
           if ((incx /= 1) .or. (incy /= 1)) then
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
               jx = kx
               jy = ky
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the triangular part
           ! of a.
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when a is stored in the upper triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= zero) .or. (y(j) /= zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           do i = 1, j
                               a(i, j) = a(i, j) + x(i)*temp1 + y(i)*temp2
                           end do
                       end if
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = kx
                           iy = ky
                           do i = 1, j
                               a(i, j) = a(i, j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do
                       end if
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           else
              ! form  a  when a is stored in the lower triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= zero) .or. (y(j) /= zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           do i = j, n
                               a(i, j) = a(i, j) + x(i)*temp1 + y(i)*temp2
                           end do
                       end if
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= zero) .or. (y(jy) /= zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = jx
                           iy = jy
                           do i = j, n
                               a(i, j) = a(i, j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do
                       end if
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_qsyr2
     end subroutine stdlib_qsyr2

     ! QSYR2K  performs one of the symmetric rank 2k operations
     ! C := alpha*A*B**T + alpha*B*A**T + beta*C,
     ! or
     ! C := alpha*A**T*B + alpha*B**T*A + beta*C,
     ! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     ! and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     ! matrices in the second case.

     subroutine stdlib_qsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldb, ldc, n
           character :: trans, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           real(qp) :: temp1, temp2
           integer(ilp) :: i, info, j, l, nrowa
           logical(lk) :: upper

           ! test the input parameters.
           if (stdlib_lsame(trans, 'n')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = stdlib_lsame(uplo, 'u')
           info = 0
           if ((.not. upper) .and. (.not. stdlib_lsame(uplo, 'l'))) then
               info = 1
           else if ((.not. stdlib_lsame(trans, 'n')) .and. (.not. stdlib_lsame(trans, 't')) .and. ( &
                     .not. stdlib_lsame(trans, 'c'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1, nrowa)) then
               info = 7
           else if (ldb < max(1, nrowa)) then
               info = 9
           else if (ldc < max(1, n)) then
               info = 12
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsyr2k', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = zero
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = beta*c(i, j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1, n
                           do i = j, n
                               c(i, j) = zero
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = j, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (stdlib_lsame(trans, 'n')) then
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
               if (upper) then
                   do j = 1, n
                       if (beta == zero) then
                           do i = 1, j
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1, j
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if ((a(j, l) /= zero) .or. (b(j, l) /= zero)) then
                               temp1 = alpha*b(j, l)
                               temp2 = alpha*a(j, l)
                               do i = 1, j
                                   c(i, j) = c(i, j) + a(i, l)*temp1 + b(i, l)*temp2
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       if (beta == zero) then
                           do i = j, n
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if ((a(j, l) /= zero) .or. (b(j, l) /= zero)) then
                               temp1 = alpha*b(j, l)
                               temp2 = alpha*a(j, l)
                               do i = j, n
                                   c(i, j) = c(i, j) + a(i, l)*temp1 + b(i, l)*temp2
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
               if (upper) then
                   do j = 1, n
                       do i = 1, j
                           temp1 = zero
                           temp2 = zero
                           do l = 1, k
                               temp1 = temp1 + a(l, i)*b(l, j)
                               temp2 = temp2 + b(l, i)*a(l, j)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp1 + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = j, n
                           temp1 = zero
                           temp2 = zero
                           do l = 1, k
                               temp1 = temp1 + a(l, i)*b(l, j)
                               temp2 = temp2 + b(l, i)*a(l, j)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp1 + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_qsyr2k
     end subroutine stdlib_qsyr2k

     ! QSYRK  performs one of the symmetric rank k operations
     ! C := alpha*A*A**T + beta*C,
     ! or
     ! C := alpha*A**T*A + beta*C,
     ! where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
     ! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     ! in the second case.

     subroutine stdlib_qsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldc, n
           character :: trans, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, j, l, nrowa
           logical(lk) :: upper

           ! test the input parameters.
           if (stdlib_lsame(trans, 'n')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = stdlib_lsame(uplo, 'u')
           info = 0
           if ((.not. upper) .and. (.not. stdlib_lsame(uplo, 'l'))) then
               info = 1
           else if ((.not. stdlib_lsame(trans, 'n')) .and. (.not. stdlib_lsame(trans, 't')) .and. ( &
                     .not. stdlib_lsame(trans, 'c'))) then
               info = 2
           else if (n < 0) then
               info = 3
           else if (k < 0) then
               info = 4
           else if (lda < max(1, nrowa)) then
               info = 7
           else if (ldc < max(1, n)) then
               info = 10
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qsyrk ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == zero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               if (upper) then
                   if (beta == zero) then
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = zero
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = beta*c(i, j)
                           end do
                       end do
                   end if
               else
                   if (beta == zero) then
                       do j = 1, n
                           do i = j, n
                               c(i, j) = zero
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = j, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (stdlib_lsame(trans, 'n')) then
              ! form  c := alpha*a*a**t + beta*c.
               if (upper) then
                   do j = 1, n
                       if (beta == zero) then
                           do i = 1, j
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1, j
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if (a(j, l) /= zero) then
                               temp = alpha*a(j, l)
                               do i = 1, j
                                   c(i, j) = c(i, j) + temp*a(i, l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       if (beta == zero) then
                           do i = j, n
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = j, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if (a(j, l) /= zero) then
                               temp = alpha*a(j, l)
                               do i = j, n
                                   c(i, j) = c(i, j) + temp*a(i, l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**t*a + beta*c.
               if (upper) then
                   do j = 1, n
                       do i = 1, j
                           temp = zero
                           do l = 1, k
                               temp = temp + a(l, i)*a(l, j)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = j, n
                           temp = zero
                           do l = 1, k
                               temp = temp + a(l, i)*a(l, j)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_qsyrk
     end subroutine stdlib_qsyrk

     ! QTBMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     subroutine stdlib_qtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, k, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kplus1, kx, l
           logical(lk) :: nounit
           ! .. intrinsic functions ..
           intrinsic :: max, min
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 2
           else if (.not. stdlib_lsame(diag, 'u') .and. .not. stdlib_lsame(diag, 'n')) then
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
               call stdlib_xerbla('stdlib_qtbmv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (stdlib_lsame(trans, 'n')) then
               ! form  x := a*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do i = max(1, j - k), j - 1
                                   x(i) = x(i) + temp*a(l + i, j)
                               end do
                               if (nounit) x(j) = x(j)*a(kplus1, j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do i = max(1, j - k), j - 1
                                   x(ix) = x(ix) + temp*a(l + i, j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(kplus1, j)
                           end if
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= zero) then
                               temp = x(j)
                               l = 1 - j
                               do i = min(n, j + k), j + 1, -1
                                   x(i) = x(i) + temp*a(l + i, j)
                               end do
                               if (nounit) x(j) = x(j)*a(1, j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do i = min(n, j + k), j + 1, -1
                                   x(ix) = x(ix) + temp*a(l + i, j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(1, j)
                           end if
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1, j)
                           do i = j - 1, max(1, j - k), -1
                               temp = temp + a(l + i, j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1, j)
                           do i = j - 1, max(1, j - k), -1
                               temp = temp + a(l + i, j)*x(ix)
                               ix = ix - incx
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           l = 1 - j
                           if (nounit) temp = temp*a(1, j)
                           do i = j + 1, min(n, j + k)
                               temp = temp + a(l + i, j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (nounit) temp = temp*a(1, j)
                           do i = j + 1, min(n, j + k)
                               temp = temp + a(l + i, j)*x(ix)
                               ix = ix + incx
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_qtbmv
     end subroutine stdlib_qtbmv

     ! QTBSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular band matrix, with ( k + 1 )
     ! diagonals.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.

     subroutine stdlib_qtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, k, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kplus1, kx, l
           logical(lk) :: nounit
           ! .. intrinsic functions ..
           intrinsic :: max, min
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 2
           else if (.not. stdlib_lsame(diag, 'u') .and. .not. stdlib_lsame(diag, 'n')) then
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
               call stdlib_xerbla('stdlib_qtbsv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with one pass through a.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := inv( a )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1, j)
                               temp = x(j)
                               do i = j - 1, max(1, j - k), -1
                                   x(i) = x(i) - temp*a(l + i, j)
                               end do
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           kx = kx - incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1, j)
                               temp = x(jx)
                               do i = j - 1, max(1, j - k), -1
                                   x(ix) = x(ix) - temp*a(l + i, j)
                                   ix = ix - incx
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1, j)
                               temp = x(j)
                               do i = j + 1, min(n, j + k)
                                   x(i) = x(i) - temp*a(l + i, j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           kx = kx + incx
                           if (x(jx) /= zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1, j)
                               temp = x(jx)
                               do i = j + 1, min(n, j + k)
                                   x(ix) = x(ix) - temp*a(l + i, j)
                                   ix = ix + incx
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t)*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           l = kplus1 - j
                           do i = max(1, j - k), j - 1
                               temp = temp - a(l + i, j)*x(i)
                           end do
                           if (nounit) temp = temp/a(kplus1, j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           do i = max(1, j - k), j - 1
                               temp = temp - a(l + i, j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(kplus1, j)
                           x(jx) = temp
                           jx = jx + incx
                           if (j > k) kx = kx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           l = 1 - j
                           do i = min(n, j + k), j + 1, -1
                               temp = temp - a(l + i, j)*x(i)
                           end do
                           if (nounit) temp = temp/a(1, j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           do i = min(n, j + k), j + 1, -1
                               temp = temp - a(l + i, j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(1, j)
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_qtbsv
     end subroutine stdlib_qtbsv

     ! QTPMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular matrix, supplied in packed form.

     subroutine stdlib_qtpmv(uplo, trans, diag, n, ap, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           real(qp) :: ap(*), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, k, kk, kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 2
           else if (.not. stdlib_lsame(diag, 'u') .and. .not. stdlib_lsame(diag, 'n')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qtpmv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x:= a*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = 1, j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do
                               if (nounit) x(j) = x(j)*ap(kk + j - 1)
                           end if
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk, kk + j - 2
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
                       do j = n, 1, -1
                           if (x(j) /= zero) then
                               temp = x(j)
                               k = kk
                               do i = n, j + 1, -1
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
                       do j = n, 1, -1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do k = kk, kk - (n - (j + 1)), -1
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
               if (stdlib_lsame(uplo, 'u')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk - 1
                           do i = j - 1, 1, -1
                               temp = temp + ap(k)*x(i)
                               k = k - 1
                           end do
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk - 1, kk - j + 1, -1
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
                       do j = 1, n
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk + 1
                           do i = j + 1, n
                               temp = temp + ap(k)*x(i)
                               k = k + 1
                           end do
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do k = kk + 1, kk + n - j
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
           ! end of stdlib_qtpmv
     end subroutine stdlib_qtpmv

     ! QTPSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular matrix, supplied in packed form.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.

     subroutine stdlib_qtpsv(uplo, trans, diag, n, ap, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           real(qp) :: ap(*), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, k, kk, kx
           logical(lk) :: nounit
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 2
           else if (.not. stdlib_lsame(diag, 'u') .and. .not. stdlib_lsame(diag, 'n')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (incx == 0) then
               info = 7
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qtpsv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := inv( a )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do i = j - 1, 1, -1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do
                           end if
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n, 1, -1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk - 1, kk - j + 1, -1
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
                       do j = 1, n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do i = j + 1, n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do
                           end if
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do k = kk + 1, kk + n - j
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
               if (stdlib_lsame(uplo, 'u')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           k = kk
                           do i = 1, j - 1
                               temp = temp - ap(k)*x(i)
                               k = k + 1
                           end do
                           if (nounit) temp = temp/ap(kk + j - 1)
                           x(j) = temp
                           kk = kk + j
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = kx
                           do k = kk, kk + j - 2
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
                       do j = n, 1, -1
                           temp = x(j)
                           k = kk
                           do i = n, j + 1, -1
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
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = kx
                           do k = kk, kk - (n - (j + 1)), -1
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
           ! end of stdlib_qtpsv
     end subroutine stdlib_qtpsv

     ! QTRMM  performs one of the matrix-matrix operations
     ! B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
     ! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     ! op( A ) = A   or   op( A ) = A**T.

     subroutine stdlib_qtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: lda, ldb, m, n
           character :: diag, side, transa, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), b(ldb, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, j, k, nrowa
           logical(lk) :: lside, nounit, upper

           ! test the input parameters.
           lside = stdlib_lsame(side, 'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
           nounit = stdlib_lsame(diag, 'n')
           upper = stdlib_lsame(uplo, 'u')
           info = 0
           if ((.not. lside) .and. (.not. stdlib_lsame(side, 'r'))) then
               info = 1
           else if ((.not. upper) .and. (.not. stdlib_lsame(uplo, 'l'))) then
               info = 2
           else if ((.not. stdlib_lsame(transa, 'n')) .and. (.not. stdlib_lsame(transa, 't')) .and. ( &
                     .not. stdlib_lsame(transa, 'c'))) then
               info = 3
           else if ((.not. stdlib_lsame(diag, 'u')) .and. (.not. stdlib_lsame(diag, 'n'))) &
                     then
               info = 4
           else if (m < 0) then
               info = 5
           else if (n < 0) then
               info = 6
           else if (lda < max(1, nrowa)) then
               info = 9
           else if (ldb < max(1, m)) then
               info = 11
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qtrmm ', info)
               return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               do j = 1, n
                   do i = 1, m
                       b(i, j) = zero
                   end do
               end do
               return
           end if
           ! start the operations.
           if (lside) then
               if (stdlib_lsame(transa, 'n')) then
                 ! form  b := alpha*a*b.
                   if (upper) then
                       do j = 1, n
                           do k = 1, m
                               if (b(k, j) /= zero) then
                                   temp = alpha*b(k, j)
                                   do i = 1, k - 1
                                       b(i, j) = b(i, j) + temp*a(i, k)
                                   end do
                                   if (nounit) temp = temp*a(k, k)
                                   b(k, j) = temp
                               end if
                           end do
                       end do
                   else
                       do j = 1, n
                           do k = m, 1, -1
                               if (b(k, j) /= zero) then
                                   temp = alpha*b(k, j)
                                   b(k, j) = temp
                                   if (nounit) b(k, j) = b(k, j)*a(k, k)
                                   do i = k + 1, m
                                       b(i, j) = b(i, j) + temp*a(i, k)
                                   end do
                               end if
                           end do
                       end do
                   end if
               else
                 ! form  b := alpha*a**t*b.
                   if (upper) then
                       do j = 1, n
                           do i = m, 1, -1
                               temp = b(i, j)
                               if (nounit) temp = temp*a(i, i)
                               do k = 1, i - 1
                                   temp = temp + a(k, i)*b(k, j)
                               end do
                               b(i, j) = alpha*temp
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = 1, m
                               temp = b(i, j)
                               if (nounit) temp = temp*a(i, i)
                               do k = i + 1, m
                                   temp = temp + a(k, i)*b(k, j)
                               end do
                               b(i, j) = alpha*temp
                           end do
                       end do
                   end if
               end if
           else
               if (stdlib_lsame(transa, 'n')) then
                 ! form  b := alpha*b*a.
                   if (upper) then
                       do j = n, 1, -1
                           temp = alpha
                           if (nounit) temp = temp*a(j, j)
                           do i = 1, m
                               b(i, j) = temp*b(i, j)
                           end do
                           do k = 1, j - 1
                               if (a(k, j) /= zero) then
                                   temp = alpha*a(k, j)
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                       end do
                   else
                       do j = 1, n
                           temp = alpha
                           if (nounit) temp = temp*a(j, j)
                           do i = 1, m
                               b(i, j) = temp*b(i, j)
                           end do
                           do k = j + 1, n
                               if (a(k, j) /= zero) then
                                   temp = alpha*a(k, j)
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                       end do
                   end if
               else
                 ! form  b := alpha*b*a**t.
                   if (upper) then
                       do k = 1, n
                           do j = 1, k - 1
                               if (a(j, k) /= zero) then
                                   temp = alpha*a(j, k)
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                           temp = alpha
                           if (nounit) temp = temp*a(k, k)
                           if (temp /= one) then
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                       end do
                   else
                       do k = n, 1, -1
                           do j = k + 1, n
                               if (a(j, k) /= zero) then
                                   temp = alpha*a(j, k)
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                           temp = alpha
                           if (nounit) temp = temp*a(k, k)
                           if (temp /= one) then
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_qtrmm
     end subroutine stdlib_qtrmm

     ! QTRMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular matrix.

     subroutine stdlib_qtrmv(uplo, trans, diag, n, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kx
           logical(lk) :: nounit
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 2
           else if (.not. stdlib_lsame(diag, 'u') .and. .not. stdlib_lsame(diag, 'n')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1, n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qtrmv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := a*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = 1, j - 1
                                   x(i) = x(i) + temp*a(i, j)
                               end do
                               if (nounit) x(j) = x(j)*a(j, j)
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = 1, j - 1
                                   x(ix) = x(ix) + temp*a(i, j)
                                   ix = ix + incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j, j)
                           end if
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= zero) then
                               temp = x(j)
                               do i = n, j + 1, -1
                                   x(i) = x(i) + temp*a(i, j)
                               end do
                               if (nounit) x(j) = x(j)*a(j, j)
                           end if
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           if (x(jx) /= zero) then
                               temp = x(jx)
                               ix = kx
                               do i = n, j + 1, -1
                                   x(ix) = x(ix) + temp*a(i, j)
                                   ix = ix - incx
                               end do
                               if (nounit) x(jx) = x(jx)*a(j, j)
                           end if
                           jx = jx - incx
                       end do
                   end if
               end if
           else
              ! form  x := a**t*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           if (nounit) temp = temp*a(j, j)
                           do i = j - 1, 1, -1
                               temp = temp + a(i, j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j, j)
                           do i = j - 1, 1, -1
                               ix = ix - incx
                               temp = temp + a(i, j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           if (nounit) temp = temp*a(j, j)
                           do i = j + 1, n
                               temp = temp + a(i, j)*x(i)
                           end do
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j, j)
                           do i = j + 1, n
                               ix = ix + incx
                               temp = temp + a(i, j)*x(ix)
                           end do
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_qtrmv
     end subroutine stdlib_qtrmv

     ! QTRSM  solves one of the matrix equations
     ! op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
     ! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     ! op( A ) = A   or   op( A ) = A**T.
     ! The matrix X is overwritten on B.

     subroutine stdlib_qtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: lda, ldb, m, n
           character :: diag, side, transa, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), b(ldb, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, j, k, nrowa
           logical(lk) :: lside, nounit, upper

           ! test the input parameters.
           lside = stdlib_lsame(side, 'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
           nounit = stdlib_lsame(diag, 'n')
           upper = stdlib_lsame(uplo, 'u')
           info = 0
           if ((.not. lside) .and. (.not. stdlib_lsame(side, 'r'))) then
               info = 1
           else if ((.not. upper) .and. (.not. stdlib_lsame(uplo, 'l'))) then
               info = 2
           else if ((.not. stdlib_lsame(transa, 'n')) .and. (.not. stdlib_lsame(transa, 't')) .and. ( &
                     .not. stdlib_lsame(transa, 'c'))) then
               info = 3
           else if ((.not. stdlib_lsame(diag, 'u')) .and. (.not. stdlib_lsame(diag, 'n'))) &
                     then
               info = 4
           else if (m < 0) then
               info = 5
           else if (n < 0) then
               info = 6
           else if (lda < max(1, nrowa)) then
               info = 9
           else if (ldb < max(1, m)) then
               info = 11
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qtrsm ', info)
               return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           ! and when  alpha.eq.zero.
           if (alpha == zero) then
               do j = 1, n
                   do i = 1, m
                       b(i, j) = zero
                   end do
               end do
               return
           end if
           ! start the operations.
           if (lside) then
               if (stdlib_lsame(transa, 'n')) then
                 ! form  b := alpha*inv( a )*b.
                   if (upper) then
                       do j = 1, n
                           if (alpha /= one) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = m, 1, -1
                               if (b(k, j) /= zero) then
                                   if (nounit) b(k, j) = b(k, j)/a(k, k)
                                   do i = 1, k - 1
                                       b(i, j) = b(i, j) - b(k, j)*a(i, k)
                                   end do
                               end if
                           end do
                       end do
                   else
                       do j = 1, n
                           if (alpha /= one) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = 1, m
                               if (b(k, j) /= zero) then
                                   if (nounit) b(k, j) = b(k, j)/a(k, k)
                                   do i = k + 1, m
                                       b(i, j) = b(i, j) - b(k, j)*a(i, k)
                                   end do
                               end if
                           end do
                       end do
                   end if
               else
                 ! form  b := alpha*inv( a**t )*b.
                   if (upper) then
                       do j = 1, n
                           do i = 1, m
                               temp = alpha*b(i, j)
                               do k = 1, i - 1
                                   temp = temp - a(k, i)*b(k, j)
                               end do
                               if (nounit) temp = temp/a(i, i)
                               b(i, j) = temp
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = m, 1, -1
                               temp = alpha*b(i, j)
                               do k = i + 1, m
                                   temp = temp - a(k, i)*b(k, j)
                               end do
                               if (nounit) temp = temp/a(i, i)
                               b(i, j) = temp
                           end do
                       end do
                   end if
               end if
           else
               if (stdlib_lsame(transa, 'n')) then
                 ! form  b := alpha*b*inv( a ).
                   if (upper) then
                       do j = 1, n
                           if (alpha /= one) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = 1, j - 1
                               if (a(k, j) /= zero) then
                                   do i = 1, m
                                       b(i, j) = b(i, j) - a(k, j)*b(i, k)
                                   end do
                               end if
                           end do
                           if (nounit) then
                               temp = one/a(j, j)
                               do i = 1, m
                                   b(i, j) = temp*b(i, j)
                               end do
                           end if
                       end do
                   else
                       do j = n, 1, -1
                           if (alpha /= one) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = j + 1, n
                               if (a(k, j) /= zero) then
                                   do i = 1, m
                                       b(i, j) = b(i, j) - a(k, j)*b(i, k)
                                   end do
                               end if
                           end do
                           if (nounit) then
                               temp = one/a(j, j)
                               do i = 1, m
                                   b(i, j) = temp*b(i, j)
                               end do
                           end if
                       end do
                   end if
               else
                 ! form  b := alpha*b*inv( a**t ).
                   if (upper) then
                       do k = n, 1, -1
                           if (nounit) then
                               temp = one/a(k, k)
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                           do j = 1, k - 1
                               if (a(j, k) /= zero) then
                                   temp = a(j, k)
                                   do i = 1, m
                                       b(i, j) = b(i, j) - temp*b(i, k)
                                   end do
                               end if
                           end do
                           if (alpha /= one) then
                               do i = 1, m
                                   b(i, k) = alpha*b(i, k)
                               end do
                           end if
                       end do
                   else
                       do k = 1, n
                           if (nounit) then
                               temp = one/a(k, k)
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                           do j = k + 1, n
                               if (a(j, k) /= zero) then
                                   temp = a(j, k)
                                   do i = 1, m
                                       b(i, j) = b(i, j) - temp*b(i, k)
                                   end do
                               end if
                           end do
                           if (alpha /= one) then
                               do i = 1, m
                                   b(i, k) = alpha*b(i, k)
                               end do
                           end if
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_qtrsm
     end subroutine stdlib_qtrsm

     ! QTRSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular matrix.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.

     subroutine stdlib_qtrsv(uplo, trans, diag, n, a, lda, x, incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           real(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           real(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kx
           logical(lk) :: nounit
           ! .. intrinsic functions ..
           intrinsic :: max
           ! test the input parameters.
           info = 0
           if (.not. stdlib_lsame(uplo, 'u') .and. .not. stdlib_lsame(uplo, 'l')) then
               info = 1
           else if (.not. stdlib_lsame(trans, 'n') .and. .not. stdlib_lsame(trans, 't') &
                     .and. .not. stdlib_lsame(trans, 'c')) then
               info = 2
           else if (.not. stdlib_lsame(diag, 'u') .and. .not. stdlib_lsame(diag, 'n')) then
               info = 3
           else if (n < 0) then
               info = 4
           else if (lda < max(1, n)) then
               info = 6
           else if (incx == 0) then
               info = 8
           end if
           if (info /= 0) then
               call stdlib_xerbla('stdlib_qtrsv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := inv( a )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j, j)
                               temp = x(j)
                               do i = j - 1, 1, -1
                                   x(i) = x(i) - temp*a(i, j)
                               end do
                           end if
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n, 1, -1
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j, j)
                               temp = x(jx)
                               ix = jx
                               do i = j - 1, 1, -1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i, j)
                               end do
                           end if
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= zero) then
                               if (nounit) x(j) = x(j)/a(j, j)
                               temp = x(j)
                               do i = j + 1, n
                                   x(i) = x(i) - temp*a(i, j)
                               end do
                           end if
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           if (x(jx) /= zero) then
                               if (nounit) x(jx) = x(jx)/a(j, j)
                               temp = x(jx)
                               ix = jx
                               do i = j + 1, n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i, j)
                               end do
                           end if
                           jx = jx + incx
                       end do
                   end if
               end if
           else
              ! form  x := inv( a**t )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           do i = 1, j - 1
                               temp = temp - a(i, j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j, j)
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = kx
                           do i = 1, j - 1
                               temp = temp - a(i, j)*x(ix)
                               ix = ix + incx
                           end do
                           if (nounit) temp = temp/a(j, j)
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           do i = n, j + 1, -1
                               temp = temp - a(i, j)*x(i)
                           end do
                           if (nounit) temp = temp/a(j, j)
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = kx
                           do i = n, j + 1, -1
                               temp = temp - a(i, j)*x(ix)
                               ix = ix - incx
                           end do
                           if (nounit) temp = temp/a(j, j)
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_qtrsv
     end subroutine stdlib_qtrsv

     ! QZASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
     ! returns a quad precision result.

     real(qp) function stdlib_qzasum(n, zx, incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, n
           ! .. array arguments ..
           complex(qp) :: zx(*)
        ! =====================================================================
           ! .. local scalars ..
           real(qp) :: stemp
           integer(ilp) :: i, nincx
           stdlib_qzasum = zero
           stemp = zero
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1, n
                 stemp = stemp + stdlib_qcabs1(zx(i))
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1, nincx, incx
                 stemp = stemp + stdlib_qcabs1(zx(i))
              end do
           end if
           stdlib_qzasum = stemp
           return
           ! end of stdlib_qzasum
     end function stdlib_qzasum

     ! !
     ! QZNRM2 returns the euclidean norm of a vector via the function
     ! name, so that
     ! QZNRM2 := sqrt( x**H*x )

     function stdlib_qznrm2(n, x, incx)
        integer, parameter :: wp = kind(1._qp)
        real(qp) :: stdlib_qznrm2
        ! -- reference blas level1 routine (version 3.9.1_qp) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
        ! .. constants ..
        real(qp), parameter :: zero = 0.0_qp
        real(qp), parameter :: one = 1.0_qp
        real(qp), parameter :: maxn = huge(0.0_qp)
        ! .. blue's scaling constants ..
     real(qp), parameter :: tsml = real(radix(0._qp), wp)**ceiling((minexponent(0._qp) - 1) &
                *0.5_qp)
     real(qp), parameter :: tbig = real(radix(0._qp), wp)**floor((maxexponent(0._qp) - &
               digits(0._qp) + 1)*0.5_qp)
     real(qp), parameter :: ssml = real(radix(0._qp), wp)**(-floor((minexponent(0._qp) - &
               digits(0._qp))*0.5_qp))
     real(qp), parameter :: sbig = real(radix(0._qp), wp)**(-ceiling((maxexponent(0._qp) &
               + digits(0._qp) - 1)*0.5_qp))
        ! .. scalar arguments ..
     integer(ilp) :: incx, n
        ! .. array arguments ..
        complex(qp) :: x(*)
        ! .. local scalars ..
     integer(ilp) :: i, ix
     logical(lk) :: notbig
        real(qp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
        ! quick return if possible
        stdlib_qznrm2 = zero
        if (n <= 0) return
        scl = one
        sumsq = zero
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1, n
           ax = abs(real(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if ((amed > zero) .or. (amed > maxn) .or. (amed /= amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range
           scl = one
           sumsq = amed
        end if
        stdlib_qznrm2 = scl*sqrt(sumsq)
        return
     end function stdlib_qznrm2

end module stdlib_linalg_blas_q