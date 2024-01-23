module stdlib_linalg_blas_w
     use stdlib_linalg_constants
     use stdlib_linalg_blas_aux
     use stdlib_linalg_blas_s
     use stdlib_linalg_blas_q
     use stdlib_linalg_blas_c
     implicit none(type, external)
     private

     public :: sp, dp, qp, lk, ilp
     public :: stdlib_waxpy
     public :: stdlib_wcopy
     public :: stdlib_wdotc
     public :: stdlib_wdotu
     public :: stdlib_wdrot
     public :: stdlib_wdscal
     public :: stdlib_wgbmv
     public :: stdlib_wgemm
     public :: stdlib_wgemv
     public :: stdlib_wgerc
     public :: stdlib_wgeru
     public :: stdlib_whbmv
     public :: stdlib_whemm
     public :: stdlib_whemv
     public :: stdlib_wher
     public :: stdlib_wher2
     public :: stdlib_wher2k
     public :: stdlib_wherk
     public :: stdlib_whpmv
     public :: stdlib_whpr
     public :: stdlib_whpr2
     public :: stdlib_wrotg
     public :: stdlib_wscal
     public :: stdlib_wswap
     public :: stdlib_wsymm
     public :: stdlib_wsyr2k
     public :: stdlib_wsyrk
     public :: stdlib_wtbmv
     public :: stdlib_wtbsv
     public :: stdlib_wtpmv
     public :: stdlib_wtpsv
     public :: stdlib_wtrmm
     public :: stdlib_wtrmv
     public :: stdlib_wtrsm
     public :: stdlib_wtrsv

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

     ! WAXPY constant times a vector plus a vector.

     subroutine stdlib_waxpy(n, za, zx, incx, zy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: za
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           complex(qp) :: zx(*), zy(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, ix, iy
           if (n <= 0) return
           if (stdlib_qcabs1(za) == zero) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1, n
                 zy(i) = zy(i) + za*zx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 zy(iy) = zy(iy) + za*zx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_waxpy
     end subroutine stdlib_waxpy

     ! WCOPY copies a vector, x, to a vector, y.

     subroutine stdlib_wcopy(n, zx, incx, zy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           complex(qp) :: zx(*), zy(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, ix, iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1, n
               zy(i) = zx(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 zy(iy) = zx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_wcopy
     end subroutine stdlib_wcopy

     ! WDOTC forms the dot product of two complex vectors
     ! WDOTC = X^H * Y

     complex(qp) function stdlib_wdotc(n, zx, incx, zy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           complex(qp) :: zx(*), zy(*)
        ! =====================================================================
           ! .. local scalars ..
           complex(qp) :: ztemp
           integer(ilp) :: i, ix, iy
           ! .. intrinsic functions ..
           intrinsic :: conjg
           ztemp = (zero, zero)
           stdlib_wdotc = (zero, zero)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1, n
                 ztemp = ztemp + conjg(zx(i))*zy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 ztemp = ztemp + conjg(zx(ix))*zy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_wdotc = ztemp
           return
           ! end of stdlib_wdotc
     end function stdlib_wdotc

     ! WDOTU forms the dot product of two complex vectors
     ! WDOTU = X^T * Y

     complex(qp) function stdlib_wdotu(n, zx, incx, zy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           complex(qp) :: zx(*), zy(*)
        ! =====================================================================
           ! .. local scalars ..
           complex(qp) :: ztemp
           integer(ilp) :: i, ix, iy
           ztemp = (zero, zero)
           stdlib_wdotu = (zero, zero)
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1, n
                 ztemp = ztemp + zx(i)*zy(i)
              end do
           else
              ! code for unequal increments or equal increments
                ! not equal to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 ztemp = ztemp + zx(ix)*zy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_wdotu = ztemp
           return
           ! end of stdlib_wdotu
     end function stdlib_wdotu

     ! Applies a plane rotation, where the cos and sin (c and s) are real
     ! and the vectors cx and cy are complex.
     ! jack dongarra, linpack, 3/11/78.

     subroutine stdlib_wdrot(n, zx, incx, zy, incy, c, s)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           real(qp) :: c, s
           ! .. array arguments ..
           complex(qp) :: zx(*), zy(*)
       ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, ix, iy
           complex(qp) :: ctemp
           ! .. executable statements ..
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
              ! code for both increments equal to 1
              do i = 1, n
                 ctemp = c*zx(i) + s*zy(i)
                 zy(i) = c*zy(i) - s*zx(i)
                 zx(i) = ctemp
              end do
           else
              ! code for unequal increments or equal increments not equal
                ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 ctemp = c*zx(ix) + s*zy(iy)
                 zy(iy) = c*zy(iy) - s*zx(ix)
                 zx(ix) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_wdrot
     end subroutine stdlib_wdrot

     ! WDSCAL scales a vector by a constant.

     subroutine stdlib_wdscal(n, da, zx, incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: da
           integer(ilp) :: incx, n
           ! .. array arguments ..
           complex(qp) :: zx(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, nincx
           ! .. intrinsic functions ..
           intrinsic :: dcmplx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1, n
                 zx(i) = cmplx(da, zero, KIND=qp)*zx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1, nincx, incx
                 zx(i) = cmplx(da, zero, KIND=qp)*zx(i)
              end do
           end if
           return
           ! end of stdlib_wdscal
     end subroutine stdlib_wdscal

     ! WGBMV  performs one of the matrix-vector operations
     ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
     ! y := alpha*A**H*x + beta*y,
     ! where alpha and beta are scalars, x and y are vectors and A is an
     ! m by n band matrix, with kl sub-diagonals and ku super-diagonals.

     subroutine stdlib_wgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: incx, incy, kl, ku, lda, m, n
           character :: trans
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, iy, j, jx, jy, k, kup1, kx, ky, lenx, leny
           logical(lk) :: noconj
           ! .. intrinsic functions ..
           intrinsic :: conjg, max, min
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
               call stdlib_xerbla('stdlib_wgbmv ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           noconj = stdlib_lsame(trans, 't')
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
           ! accessed sequentially with cone pass through the band part of a.
           ! first form  y := beta*y.
           if (beta /= cone) then
               if (incy == 1) then
                   if (beta == czero) then
                       do i = 1, leny
                           y(i) = czero
                       end do
                   else
                       do i = 1, leny
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == czero) then
                       do i = 1, leny
                           y(iy) = czero
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
           if (alpha == czero) return
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
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
               jy = ky
               if (incx == 1) then
                   do j = 1, n
                       temp = czero
                       k = kup1 - j
                       if (noconj) then
                           do i = max(1, j - ku), min(m, j + kl)
                               temp = temp + a(k + i, j)*x(i)
                           end do
                       else
                           do i = max(1, j - ku), min(m, j + kl)
                               temp = temp + conjg(a(k + i, j))*x(i)
                           end do
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do
               else
                   do j = 1, n
                       temp = czero
                       ix = kx
                       k = kup1 - j
                       if (noconj) then
                           do i = max(1, j - ku), min(m, j + kl)
                               temp = temp + a(k + i, j)*x(ix)
                               ix = ix + incx
                           end do
                       else
                           do i = max(1, j - ku), min(m, j + kl)
                               temp = temp + conjg(a(k + i, j))*x(ix)
                               ix = ix + incx
                           end do
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j > ku) kx = kx + incx
                   end do
               end if
           end if
           return
           ! end of stdlib_wgbmv
     end subroutine stdlib_wgbmv

     ! WGEMM  performs one of the matrix-matrix operations
     ! C := alpha*op( A )*op( B ) + beta*C,
     ! where  op( X ) is one of
     ! op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
     ! alpha and beta are scalars, and A, B and C are matrices, with op( A )
     ! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

     subroutine stdlib_wgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldb, ldc, m, n
           character :: transa, transb
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, j, l, nrowa, nrowb
           logical(lk) :: conja, conjb, nota, notb

           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and nrowb  as the number of rows  of  a  and  b  respectively.
           nota = stdlib_lsame(transa, 'n')
           notb = stdlib_lsame(transb, 'n')
           conja = stdlib_lsame(transa, 'c')
           conjb = stdlib_lsame(transb, 'c')
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
           if ((.not. nota) .and. (.not. conja) .and. (.not. stdlib_lsame(transa, 't'))) then
               info = 1
           else if ((.not. notb) .and. (.not. conjb) .and. (.not. stdlib_lsame(transb, 't'))) &
                     then
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
               call stdlib_xerbla('stdlib_wgemm ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) &
                     return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = czero
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
                       if (beta == czero) then
                           do i = 1, m
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
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
               else if (conja) then
                 ! form  c := alpha*a**h*b + beta*c.
                   do j = 1, n
                       do i = 1, m
                           temp = czero
                           do l = 1, k
                               temp = temp + conjg(a(l, i))*b(l, j)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b + beta*c
                   do j = 1, n
                       do i = 1, m
                           temp = czero
                           do l = 1, k
                               temp = temp + a(l, i)*b(l, j)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           else if (nota) then
               if (conjb) then
                 ! form  c := alpha*a*b**h + beta*c.
                   do j = 1, n
                       if (beta == czero) then
                           do i = 1, m
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1, m
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           temp = alpha*conjg(b(j, l))
                           do i = 1, m
                               c(i, j) = c(i, j) + temp*a(i, l)
                           end do
                       end do
                   end do
               else
                 ! form  c := alpha*a*b**t + beta*c
                   do j = 1, n
                       if (beta == czero) then
                           do i = 1, m
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
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
               end if
           else if (conja) then
               if (conjb) then
                 ! form  c := alpha*a**h*b**h + beta*c.
                   do j = 1, n
                       do i = 1, m
                           temp = czero
                           do l = 1, k
                               temp = temp + conjg(a(l, i))*conjg(b(j, l))
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**h*b**t + beta*c
                   do j = 1, n
                       do i = 1, m
                           temp = czero
                           do l = 1, k
                               temp = temp + conjg(a(l, i))*b(j, l)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           else
               if (conjb) then
                 ! form  c := alpha*a**t*b**h + beta*c
                   do j = 1, n
                       do i = 1, m
                           temp = czero
                           do l = 1, k
                               temp = temp + a(l, i)*conjg(b(j, l))
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               else
                 ! form  c := alpha*a**t*b**t + beta*c
                   do j = 1, n
                       do i = 1, m
                           temp = czero
                           do l = 1, k
                               temp = temp + a(l, i)*b(j, l)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_wgemm
     end subroutine stdlib_wgemm

     ! WGEMV  performs one of the matrix-vector operations
     ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
     ! y := alpha*A**H*x + beta*y,
     ! where alpha and beta are scalars, x and y are vectors and A is an
     ! m by n matrix.

     subroutine stdlib_wgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: incx, incy, lda, m, n
           character :: trans
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
           logical(lk) :: noconj
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
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
               call stdlib_xerbla('stdlib_wgemv ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           noconj = stdlib_lsame(trans, 't')
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
           ! accessed sequentially with cone pass through a.
           ! first form  y := beta*y.
           if (beta /= cone) then
               if (incy == 1) then
                   if (beta == czero) then
                       do i = 1, leny
                           y(i) = czero
                       end do
                   else
                       do i = 1, leny
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == czero) then
                       do i = 1, leny
                           y(iy) = czero
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
           if (alpha == czero) return
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
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
               jy = ky
               if (incx == 1) then
                   do j = 1, n
                       temp = czero
                       if (noconj) then
                           do i = 1, m
                               temp = temp + a(i, j)*x(i)
                           end do
                       else
                           do i = 1, m
                               temp = temp + conjg(a(i, j))*x(i)
                           end do
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do
               else
                   do j = 1, n
                       temp = czero
                       ix = kx
                       if (noconj) then
                           do i = 1, m
                               temp = temp + a(i, j)*x(ix)
                               ix = ix + incx
                           end do
                       else
                           do i = 1, m
                               temp = temp + conjg(a(i, j))*x(ix)
                               ix = ix + incx
                           end do
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_wgemv
     end subroutine stdlib_wgemv

     ! WGERC  performs the rank 1 operation
     ! A := alpha*x*y**H + A,
     ! where alpha is a scalar, x is an m element vector, y is an n element
     ! vector and A is an m by n matrix.

     subroutine stdlib_wgerc(m, n, alpha, x, incx, y, incy, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           integer(ilp) :: incx, incy, lda, m, n
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jy, kx
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
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
               call stdlib_xerbla('stdlib_wgerc ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (alpha == czero)) return
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (incy > 0) then
               jy = 1
           else
               jy = 1 - (n - 1)*incy
           end if
           if (incx == 1) then
               do j = 1, n
                   if (y(jy) /= czero) then
                       temp = alpha*conjg(y(jy))
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
                   if (y(jy) /= czero) then
                       temp = alpha*conjg(y(jy))
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
           ! end of stdlib_wgerc
     end subroutine stdlib_wgerc

     ! WGERU  performs the rank 1 operation
     ! A := alpha*x*y**T + A,
     ! where alpha is a scalar, x is an m element vector, y is an n element
     ! vector and A is an m by n matrix.

     subroutine stdlib_wgeru(m, n, alpha, x, incx, y, incy, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           integer(ilp) :: incx, incy, lda, m, n
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
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
               call stdlib_xerbla('stdlib_wgeru ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. (alpha == czero)) return
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (incy > 0) then
               jy = 1
           else
               jy = 1 - (n - 1)*incy
           end if
           if (incx == 1) then
               do j = 1, n
                   if (y(jy) /= czero) then
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
                   if (y(jy) /= czero) then
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
           ! end of stdlib_wgeru
     end subroutine stdlib_wgeru

     ! WHBMV  performs the matrix-vector  operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n hermitian band matrix, with k super-diagonals.

     subroutine stdlib_whbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: incx, incy, k, lda, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
           ! .. intrinsic functions ..
           intrinsic :: real, conjg, max, min
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
               call stdlib_xerbla('stdlib_whbmv ', info)
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
           ! start the operations. in this version the elements of the array a
           ! are accessed sequentially with cone pass through a.
           ! first form  y := beta*y.
           if (beta /= cone) then
               if (incy == 1) then
                   if (beta == czero) then
                       do i = 1, n
                           y(i) = czero
                       end do
                   else
                       do i = 1, n
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == czero) then
                       do i = 1, n
                           y(iy) = czero
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
           if (alpha == czero) return
           if (stdlib_lsame(uplo, 'u')) then
              ! form  y  when upper triangle of a is stored.
               kplus1 = k + 1
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = czero
                       l = kplus1 - j
                       do i = max(1, j - k), j - 1
                           y(i) = y(i) + temp1*a(l + i, j)
                           temp2 = temp2 + conjg(a(l + i, j))*x(i)
                       end do
                       y(j) = y(j) + temp1*real(a(kplus1, j), KIND=qp) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = czero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       do i = max(1, j - k), j - 1
                           y(iy) = y(iy) + temp1*a(l + i, j)
                           temp2 = temp2 + conjg(a(l + i, j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do
                       y(jy) = y(jy) + temp1*real(a(kplus1, j), KIND=qp) + alpha*temp2
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
                       temp2 = czero
                       y(j) = y(j) + temp1*real(a(1, j), KIND=qp)
                       l = 1 - j
                       do i = j + 1, min(n, j + k)
                           y(i) = y(i) + temp1*a(l + i, j)
                           temp2 = temp2 + conjg(a(l + i, j))*x(i)
                       end do
                       y(j) = y(j) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = czero
                       y(jy) = y(jy) + temp1*real(a(1, j), KIND=qp)
                       l = 1 - j
                       ix = jx
                       iy = jy
                       do i = j + 1, min(n, j + k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l + i, j)
                           temp2 = temp2 + conjg(a(l + i, j))*x(ix)
                       end do
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_whbmv
     end subroutine stdlib_whbmv

     ! WHEMM  performs one of the matrix-matrix operations
     ! C := alpha*A*B + beta*C,
     ! or
     ! C := alpha*B*A + beta*C,
     ! where alpha and beta are scalars, A is an hermitian matrix and  B and
     ! C are m by n matrices.

     subroutine stdlib_whemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: lda, ldb, ldc, m, n
           character :: side, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: real, conjg, max
           ! .. local scalars ..
           complex(qp) :: temp1, temp2
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
               call stdlib_xerbla('stdlib_whemm ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = czero
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
                           temp2 = czero
                           do k = 1, i - 1
                               c(k, j) = c(k, j) + temp1*a(k, i)
                               temp2 = temp2 + b(k, j)*conjg(a(k, i))
                           end do
                           if (beta == czero) then
                               c(i, j) = temp1*real(a(i, i), KIND=qp) + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + temp1*real(a(i, i), KIND=qp) + &
                                         alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = m, 1, -1
                           temp1 = alpha*b(i, j)
                           temp2 = czero
                           do k = i + 1, m
                               c(k, j) = c(k, j) + temp1*a(k, i)
                               temp2 = temp2 + b(k, j)*conjg(a(k, i))
                           end do
                           if (beta == czero) then
                               c(i, j) = temp1*real(a(i, i), KIND=qp) + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + temp1*real(a(i, i), KIND=qp) + &
                                         alpha*temp2
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*b*a + beta*c.
               loop_170: do j = 1, n
                   temp1 = alpha*real(a(j, j), KIND=qp)
                   if (beta == czero) then
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
                           temp1 = alpha*conjg(a(j, k))
                       end if
                       do i = 1, m
                           c(i, j) = c(i, j) + temp1*b(i, k)
                       end do
                   end do
                   do k = j + 1, n
                       if (upper) then
                           temp1 = alpha*conjg(a(j, k))
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
           ! end of stdlib_whemm
     end subroutine stdlib_whemm

     ! WHEMV  performs the matrix-vector  operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n hermitian matrix.

     subroutine stdlib_whemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: incx, incy, lda, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kx, ky
           ! .. intrinsic functions ..
           intrinsic :: real, conjg, max
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
               call stdlib_xerbla('stdlib_whemv ', info)
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
                       do i = 1, n
                           y(i) = czero
                       end do
                   else
                       do i = 1, n
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == czero) then
                       do i = 1, n
                           y(iy) = czero
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
           if (alpha == czero) return
           if (stdlib_lsame(uplo, 'u')) then
              ! form  y  when a is stored in upper triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = czero
                       do i = 1, j - 1
                           y(i) = y(i) + temp1*a(i, j)
                           temp2 = temp2 + conjg(a(i, j))*x(i)
                       end do
                       y(j) = y(j) + temp1*real(a(j, j), KIND=qp) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = czero
                       ix = kx
                       iy = ky
                       do i = 1, j - 1
                           y(iy) = y(iy) + temp1*a(i, j)
                           temp2 = temp2 + conjg(a(i, j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do
                       y(jy) = y(jy) + temp1*real(a(j, j), KIND=qp) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           else
              ! form  y  when a is stored in lower triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = czero
                       y(j) = y(j) + temp1*real(a(j, j), KIND=qp)
                       do i = j + 1, n
                           y(i) = y(i) + temp1*a(i, j)
                           temp2 = temp2 + conjg(a(i, j))*x(i)
                       end do
                       y(j) = y(j) + alpha*temp2
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = czero
                       y(jy) = y(jy) + temp1*real(a(j, j), KIND=qp)
                       ix = jx
                       iy = jy
                       do i = j + 1, n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i, j)
                           temp2 = temp2 + conjg(a(i, j))*x(ix)
                       end do
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_whemv
     end subroutine stdlib_whemv

     ! WHER   performs the hermitian rank 1 operation
     ! A := alpha*x*x**H + A,
     ! where alpha is a real scalar, x is an n element vector and A is an
     ! n by n hermitian matrix.

     subroutine stdlib_wher(uplo, n, alpha, x, incx, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, lda, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kx
           ! .. intrinsic functions ..
           intrinsic :: real, conjg, max
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
               call stdlib_xerbla('stdlib_wher  ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == real(czero, KIND=qp))) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when a is stored in upper triangle.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= czero) then
                           temp = alpha*conjg(x(j))
                           do i = 1, j - 1
                               a(i, j) = a(i, j) + x(i)*temp
                           end do
                           a(j, j) = real(a(j, j), KIND=qp) + real(x(j)*temp, KIND=qp)
                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= czero) then
                           temp = alpha*conjg(x(jx))
                           ix = kx
                           do i = 1, j - 1
                               a(i, j) = a(i, j) + x(ix)*temp
                               ix = ix + incx
                           end do
                           a(j, j) = real(a(j, j), KIND=qp) + real(x(jx)*temp, KIND=qp)
                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                       jx = jx + incx
                   end do
               end if
           else
              ! form  a  when a is stored in lower triangle.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= czero) then
                           temp = alpha*conjg(x(j))
                           a(j, j) = real(a(j, j), KIND=qp) + real(temp*x(j), KIND=qp)
                           do i = j + 1, n
                               a(i, j) = a(i, j) + x(i)*temp
                           end do
                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= czero) then
                           temp = alpha*conjg(x(jx))
                           a(j, j) = real(a(j, j), KIND=qp) + real(temp*x(jx), KIND=qp)
                           ix = jx
                           do i = j + 1, n
                               ix = ix + incx
                               a(i, j) = a(i, j) + x(ix)*temp
                           end do
                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                       jx = jx + incx
                   end do
               end if
           end if
           return
           ! end of stdlib_wher
     end subroutine stdlib_wher

     ! WHER2  performs the hermitian rank 2 operation
     ! A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
     ! where alpha is a scalar, x and y are n element vectors and A is an n
     ! by n hermitian matrix.

     subroutine stdlib_wher2(uplo, n, alpha, x, incx, y, incy, a, lda)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           integer(ilp) :: incx, incy, lda, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, kx, ky
           ! .. intrinsic functions ..
           intrinsic :: real, conjg, max
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
               call stdlib_xerbla('stdlib_wher2 ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
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
           ! accessed sequentially with cone pass through the triangular part
           ! of a.
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when a is stored in the upper triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= czero) .or. (y(j) /= czero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           do i = 1, j - 1
                               a(i, j) = a(i, j) + x(i)*temp1 + y(i)*temp2
                           end do
                           a(j, j) = real(a(j, j), KIND=qp) + real(x(j)*temp1 + y(j)*temp2, KIND=qp)

                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= czero) .or. (y(jy) /= czero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           do i = 1, j - 1
                               a(i, j) = a(i, j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do
                           a(j, j) = real(a(j, j), KIND=qp) + real(x(jx)*temp1 + y(jy)*temp2, KIND=qp)

                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           else
              ! form  a  when a is stored in the lower triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= czero) .or. (y(j) /= czero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           a(j, j) = real(a(j, j), KIND=qp) + real(x(j)*temp1 + y(j)*temp2, KIND=qp)

                           do i = j + 1, n
                               a(i, j) = a(i, j) + x(i)*temp1 + y(i)*temp2
                           end do
                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= czero) .or. (y(jy) /= czero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           a(j, j) = real(a(j, j), KIND=qp) + real(x(jx)*temp1 + y(jy)*temp2, KIND=qp)

                           ix = jx
                           iy = jy
                           do i = j + 1, n
                               ix = ix + incx
                               iy = iy + incy
                               a(i, j) = a(i, j) + x(ix)*temp1 + y(iy)*temp2
                           end do
                       else
                           a(j, j) = real(a(j, j), KIND=qp)
                       end if
                       jx = jx + incx
                       jy = jy + incy
                   end do
               end if
           end if
           return
           ! end of stdlib_wher2
     end subroutine stdlib_wher2

     ! WHER2K  performs one of the hermitian rank 2k operations
     ! C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
     ! or
     ! C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
     ! where  alpha and beta  are scalars with  beta  real,  C is an  n by n
     ! hermitian matrix and  A and B  are  n by k matrices in the first case
     ! and  k by n  matrices in the second case.

     subroutine stdlib_wher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           real(qp) :: beta
           integer(ilp) :: k, lda, ldb, ldc, n
           character :: trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: real, conjg, max
           ! .. local scalars ..
           complex(qp) :: temp1, temp2
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
           else if ((.not. stdlib_lsame(trans, 'n')) .and. (.not. stdlib_lsame(trans, 'c'))) &
                     then
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
               call stdlib_xerbla('stdlib_wher2k', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == one))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == real(czero, KIND=qp)) then
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = czero
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = 1, j - 1
                               c(i, j) = beta*c(i, j)
                           end do
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                       end do
                   end if
               else
                   if (beta == real(czero, KIND=qp)) then
                       do j = 1, n
                           do i = j, n
                               c(i, j) = czero
                           end do
                       end do
                   else
                       do j = 1, n
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                           do i = j + 1, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (stdlib_lsame(trans, 'n')) then
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
               if (upper) then
                   do j = 1, n
                       if (beta == real(czero, KIND=qp)) then
                           do i = 1, j
                               c(i, j) = czero
                           end do
                       else if (beta /= one) then
                           do i = 1, j - 1
                               c(i, j) = beta*c(i, j)
                           end do
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                       else
                           c(j, j) = real(c(j, j), KIND=qp)
                       end if
                       do l = 1, k
                           if ((a(j, l) /= czero) .or. (b(j, l) /= czero)) then
                               temp1 = alpha*conjg(b(j, l))
                               temp2 = conjg(alpha*a(j, l))
                               do i = 1, j - 1
                                   c(i, j) = c(i, j) + a(i, l)*temp1 + b(i, l)*temp2
                               end do
                               c(j, j) = real(c(j, j), KIND=qp) + real(a(j, l)*temp1 + b(j, l)*temp2, &
                                         KIND=qp)
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       if (beta == real(czero, KIND=qp)) then
                           do i = j, n
                               c(i, j) = czero
                           end do
                       else if (beta /= one) then
                           do i = j + 1, n
                               c(i, j) = beta*c(i, j)
                           end do
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                       else
                           c(j, j) = real(c(j, j), KIND=qp)
                       end if
                       do l = 1, k
                           if ((a(j, l) /= czero) .or. (b(j, l) /= czero)) then
                               temp1 = alpha*conjg(b(j, l))
                               temp2 = conjg(alpha*a(j, l))
                               do i = j + 1, n
                                   c(i, j) = c(i, j) + a(i, l)*temp1 + b(i, l)*temp2
                               end do
                               c(j, j) = real(c(j, j), KIND=qp) + real(a(j, l)*temp1 + b(j, l)*temp2, &
                                         KIND=qp)
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
               if (upper) then
                   do j = 1, n
                       do i = 1, j
                           temp1 = czero
                           temp2 = czero
                           do l = 1, k
                               temp1 = temp1 + conjg(a(l, i))*b(l, j)
                               temp2 = temp2 + conjg(b(l, i))*a(l, j)
                           end do
                           if (i == j) then
                               if (beta == real(czero, KIND=qp)) then
                                   c(j, j) = real(alpha*temp1 + conjg(alpha)*temp2, KIND=qp)
                               else
                                   c(j, j) = beta*real(c(j, j), KIND=qp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2, KIND=qp)
                               end if
                           else
                               if (beta == real(czero, KIND=qp)) then
                                   c(i, j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i, j) = beta*c(i, j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = j, n
                           temp1 = czero
                           temp2 = czero
                           do l = 1, k
                               temp1 = temp1 + conjg(a(l, i))*b(l, j)
                               temp2 = temp2 + conjg(b(l, i))*a(l, j)
                           end do
                           if (i == j) then
                               if (beta == real(czero, KIND=qp)) then
                                   c(j, j) = real(alpha*temp1 + conjg(alpha)*temp2, KIND=qp)
                               else
                                   c(j, j) = beta*real(c(j, j), KIND=qp) + real(alpha*temp1 + conjg( &
                                             alpha)*temp2, KIND=qp)
                               end if
                           else
                               if (beta == real(czero, KIND=qp)) then
                                   c(i, j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i, j) = beta*c(i, j) + alpha*temp1 + conjg(alpha)*temp2
                               end if
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_wher2k
     end subroutine stdlib_wher2k

     ! WHERK  performs one of the hermitian rank k operations
     ! C := alpha*A*A**H + beta*C,
     ! or
     ! C := alpha*A**H*A + beta*C,
     ! where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
     ! matrix and  A  is an  n by k  matrix in the  first case and a  k by n
     ! matrix in the second case.

     subroutine stdlib_wherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldc, n
           character :: trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: real, dcmplx, conjg, max
           ! .. local scalars ..
           complex(qp) :: temp
           real(qp) :: rtemp
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
           else if ((.not. stdlib_lsame(trans, 'n')) .and. (.not. stdlib_lsame(trans, 'c'))) &
                     then
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
               call stdlib_xerbla('stdlib_wherk ', info)
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
                           do i = 1, j - 1
                               c(i, j) = beta*c(i, j)
                           end do
                           c(j, j) = beta*real(c(j, j), KIND=qp)
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
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                           do i = j + 1, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end do
                   end if
               end if
               return
           end if
           ! start the operations.
           if (stdlib_lsame(trans, 'n')) then
              ! form  c := alpha*a*a**h + beta*c.
               if (upper) then
                   do j = 1, n
                       if (beta == zero) then
                           do i = 1, j
                               c(i, j) = zero
                           end do
                       else if (beta /= one) then
                           do i = 1, j - 1
                               c(i, j) = beta*c(i, j)
                           end do
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                       else
                           c(j, j) = real(c(j, j), KIND=qp)
                       end if
                       do l = 1, k
                           if (a(j, l) /= cmplx(zero, KIND=qp)) then
                               temp = alpha*conjg(a(j, l))
                               do i = 1, j - 1
                                   c(i, j) = c(i, j) + temp*a(i, l)
                               end do
                               c(j, j) = real(c(j, j), KIND=qp) + real(temp*a(i, l), KIND=qp)
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
                           c(j, j) = beta*real(c(j, j), KIND=qp)
                           do i = j + 1, n
                               c(i, j) = beta*c(i, j)
                           end do
                       else
                           c(j, j) = real(c(j, j), KIND=qp)
                       end if
                       do l = 1, k
                           if (a(j, l) /= cmplx(zero, KIND=qp)) then
                               temp = alpha*conjg(a(j, l))
                               c(j, j) = real(c(j, j), KIND=qp) + real(temp*a(j, l), KIND=qp)
                               do i = j + 1, n
                                   c(i, j) = c(i, j) + temp*a(i, l)
                               end do
                           end if
                       end do
                   end do
               end if
           else
              ! form  c := alpha*a**h*a + beta*c.
               if (upper) then
                   do j = 1, n
                       do i = 1, j - 1
                           temp = zero
                           do l = 1, k
                               temp = temp + conjg(a(l, i))*a(l, j)
                           end do
                           if (beta == zero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                       rtemp = zero
                       do l = 1, k
                           rtemp = rtemp + conjg(a(l, j))*a(l, j)
                       end do
                       if (beta == zero) then
                           c(j, j) = alpha*rtemp
                       else
                           c(j, j) = alpha*rtemp + beta*real(c(j, j), KIND=qp)
                       end if
                   end do
               else
                   do j = 1, n
                       rtemp = zero
                       do l = 1, k
                           rtemp = rtemp + conjg(a(l, j))*a(l, j)
                       end do
                       if (beta == zero) then
                           c(j, j) = alpha*rtemp
                       else
                           c(j, j) = alpha*rtemp + beta*real(c(j, j), KIND=qp)
                       end if
                       do i = j + 1, n
                           temp = zero
                           do l = 1, k
                               temp = temp + conjg(a(l, i))*a(l, j)
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
           ! end of stdlib_wherk
     end subroutine stdlib_wherk

     ! WHPMV  performs the matrix-vector operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n hermitian matrix, supplied in packed form.

     subroutine stdlib_whpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: incx, incy, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: ap(*), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, k, kk, kx, ky
           ! .. intrinsic functions ..
           intrinsic :: real, conjg
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
               call stdlib_xerbla('stdlib_whpmv ', info)
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
                       do i = 1, n
                           y(i) = czero
                       end do
                   else
                       do i = 1, n
                           y(i) = beta*y(i)
                       end do
                   end if
               else
                   iy = ky
                   if (beta == czero) then
                       do i = 1, n
                           y(iy) = czero
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
           if (alpha == czero) return
           kk = 1
           if (stdlib_lsame(uplo, 'u')) then
              ! form  y  when ap contains the upper triangle.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       temp1 = alpha*x(j)
                       temp2 = czero
                       k = kk
                       do i = 1, j - 1
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(i)
                           k = k + 1
                       end do
                       y(j) = y(j) + temp1*real(ap(kk + j - 1), KIND=qp) + alpha*temp2
                       kk = kk + j
                   end do
               else
                   jx = kx
                   jy = ky
                   do j = 1, n
                       temp1 = alpha*x(jx)
                       temp2 = czero
                       ix = kx
                       iy = ky
                       do k = kk, kk + j - 2
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do
                       y(jy) = y(jy) + temp1*real(ap(kk + j - 1), KIND=qp) + alpha*temp2
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
                       temp2 = czero
                       y(j) = y(j) + temp1*real(ap(kk), KIND=qp)
                       k = kk + 1
                       do i = j + 1, n
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(i)
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
                       temp2 = czero
                       y(jy) = y(jy) + temp1*real(ap(kk), KIND=qp)
                       ix = jx
                       iy = jy
                       do k = kk + 1, kk + n - j
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(ix)
                       end do
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + (n - j + 1)
                   end do
               end if
           end if
           return
           ! end of stdlib_whpmv
     end subroutine stdlib_whpmv

     ! WHPR    performs the hermitian rank 1 operation
     ! A := alpha*x*x**H + A,
     ! where alpha is a real scalar, x is an n element vector and A is an
     ! n by n hermitian matrix, supplied in packed form.

     subroutine stdlib_whpr(uplo, n, alpha, x, incx, ap)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           real(qp) :: alpha
           integer(ilp) :: incx, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: ap(*), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, k, kk, kx
           ! .. intrinsic functions ..
           intrinsic :: real, conjg
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
               call stdlib_xerbla('stdlib_whpr  ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == real(czero, KIND=qp))) return
           ! set the start point in x if the increment is not unity.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with cone pass through ap.
           kk = 1
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when upper triangle is stored in ap.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= czero) then
                           temp = alpha*conjg(x(j))
                           k = kk
                           do i = 1, j - 1
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
                           end do
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp) + real(x(j)*temp, KIND=qp)
                       else
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp)
                       end if
                       kk = kk + j
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= czero) then
                           temp = alpha*conjg(x(jx))
                           ix = kx
                           do k = kk, kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
                           end do
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp) + real(x(jx)*temp, KIND=qp)

                       else
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp)
                       end if
                       jx = jx + incx
                       kk = kk + j
                   end do
               end if
           else
              ! form  a  when lower triangle is stored in ap.
               if (incx == 1) then
                   do j = 1, n
                       if (x(j) /= czero) then
                           temp = alpha*conjg(x(j))
                           ap(kk) = real(ap(kk), KIND=qp) + real(temp*x(j), KIND=qp)
                           k = kk + 1
                           do i = j + 1, n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
                           end do
                       else
                           ap(kk) = real(ap(kk), KIND=qp)
                       end if
                       kk = kk + n - j + 1
                   end do
               else
                   jx = kx
                   do j = 1, n
                       if (x(jx) /= czero) then
                           temp = alpha*conjg(x(jx))
                           ap(kk) = real(ap(kk), KIND=qp) + real(temp*x(jx), KIND=qp)
                           ix = jx
                           do k = kk + 1, kk + n - j
                               ix = ix + incx
                               ap(k) = ap(k) + x(ix)*temp
                           end do
                       else
                           ap(kk) = real(ap(kk), KIND=qp)
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
                   end do
               end if
           end if
           return
           ! end of stdlib_whpr
     end subroutine stdlib_whpr

     ! WHPR2  performs the hermitian rank 2 operation
     ! A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
     ! where alpha is a scalar, x and y are n element vectors and A is an
     ! n by n hermitian matrix, supplied in packed form.

     subroutine stdlib_whpr2(uplo, n, alpha, x, incx, y, incy, ap)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           integer(ilp) :: incx, incy, n
           character :: uplo
           ! .. array arguments ..
           complex(qp) :: ap(*), x(*), y(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp1, temp2
           integer(ilp) :: i, info, ix, iy, j, jx, jy, k, kk, kx, ky
           ! .. intrinsic functions ..
           intrinsic :: real, conjg
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
               call stdlib_xerbla('stdlib_whpr2 ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (alpha == czero)) return
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
           ! are accessed sequentially with cone pass through ap.
           kk = 1
           if (stdlib_lsame(uplo, 'u')) then
              ! form  a  when upper triangle is stored in ap.
               if ((incx == 1) .and. (incy == 1)) then
                   do j = 1, n
                       if ((x(j) /= czero) .or. (y(j) /= czero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           k = kk
                           do i = 1, j - 1
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
                           end do
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp) + real(x(j)*temp1 + y(j)*temp2, &
                                     KIND=qp)
                       else
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp)
                       end if
                       kk = kk + j
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= czero) .or. (y(jy) /= czero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           do k = kk, kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp) + real(x(jx)*temp1 + y(jy)*temp2, &
                                     KIND=qp)
                       else
                           ap(kk + j - 1) = real(ap(kk + j - 1), KIND=qp)
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
                       if ((x(j) /= czero) .or. (y(j) /= czero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           ap(kk) = real(ap(kk), KIND=qp) + real(x(j)*temp1 + y(j)*temp2, KIND=qp)

                           k = kk + 1
                           do i = j + 1, n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
                           end do
                       else
                           ap(kk) = real(ap(kk), KIND=qp)
                       end if
                       kk = kk + n - j + 1
                   end do
               else
                   do j = 1, n
                       if ((x(jx) /= czero) .or. (y(jy) /= czero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           ap(kk) = real(ap(kk), KIND=qp) + real(x(jx)*temp1 + y(jy)*temp2, KIND=qp)

                           ix = jx
                           iy = jy
                           do k = kk + 1, kk + n - j
                               ix = ix + incx
                               iy = iy + incy
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                           end do
                       else
                           ap(kk) = real(ap(kk), KIND=qp)
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
                   end do
               end if
           end if
           return
           ! end of stdlib_whpr2
     end subroutine stdlib_whpr2

     ! !
     ! The computation uses the formulas
     ! |x| = sqrt( Re(x)**2 + Im(x)**2 )
     ! sgn(x) = x / |x|  if x /= 0
     ! = 1        if x  = 0
     ! c = |a| / sqrt(|a|**2 + |b|**2)
     ! s = sgn(a) * conjg(b) / sqrt(|a|**2 + |b|**2)
     ! When a and b are real and r /= 0, the formulas simplify to
     ! r = sgn(a)*sqrt(|a|**2 + |b|**2)
     ! c = a / r
     ! s = b / r
     ! the same as in DROTG when |a| > |b|.  When |b| >= |a|, the
     ! sign of c and s will be different from those computed by DROTG
     ! if the signs of a and b are not the same.

     subroutine stdlib_wrotg(a, b, c, s)
        integer, parameter :: wp = kind(1._qp)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! .. constants ..
        real(qp), parameter :: zero = 0.0_qp
        real(qp), parameter :: one = 1.0_qp
        complex(qp), parameter :: czero = 0.0_qp
        ! .. scaling constants ..
     real(qp), parameter :: safmin = real(radix(0._qp), wp)**max(minexponent(0._qp) - 1, 1 - &
               maxexponent(0._qp))
     real(qp), parameter :: safmax = real(radix(0._qp), wp)**max(1 - minexponent(0._qp), maxexponent( &
               0._qp) - 1)
     real(qp), parameter :: rtmin = sqrt(real(radix(0._qp), wp)**max(minexponent(0._qp) - 1, 1 - &
               maxexponent(0._qp))/epsilon(0._qp))
     real(qp), parameter :: rtmax = sqrt(real(radix(0._qp), wp)**max(1 - minexponent(0._qp), &
               maxexponent(0._qp) - 1)*epsilon(0._qp))
        ! .. scalar arguments ..
        real(qp) :: c
        complex(qp) :: a, b, s
        ! .. local scalars ..
        real(qp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
        complex(qp) :: f, fs, g, gs, r, t
        ! .. intrinsic functions ..
        intrinsic :: abs, aimag, conjg, max, min, real, sqrt
        ! .. statement functions ..
        real(qp) :: abssq
        ! .. statement function definitions ..
        abssq(t) = real(t)**2 + aimag(t)**2
        ! .. executable statements ..
        f = a
        g = b
        if (g == czero) then
           c = one
           s = czero
           r = f
        else if (f == czero) then
           c = zero
           g1 = max(abs(real(g)), abs(aimag(g)))
           if (g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              g2 = abssq(g)
              d = sqrt(g2)
              s = conjg(g)/d
              r = d
           else
              ! use scaled algorithm
              u = min(safmax, max(safmin, g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              d = sqrt(g2)
              s = conjg(gs)/d
              r = d*u
           end if
        else
           f1 = max(abs(real(f)), abs(aimag(f)))
           g1 = max(abs(real(g)), abs(aimag(g)))
     if (f1 > rtmin .and. f1 < rtmax .and. g1 > rtmin .and. g1 < rtmax) then
              ! use unscaled algorithm
              f2 = abssq(f)
              g2 = abssq(g)
              h2 = f2 + g2
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = f2*p
              s = conjg(g)*(f*p)
              r = f*(h2*p)
           else
              ! use scaled algorithm
              u = min(safmax, max(safmin, f1, g1))
              uu = one/u
              gs = g*uu
              g2 = abssq(gs)
              if (f1*uu < rtmin) then
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
                 v = min(safmax, max(safmin, f1))
                 vv = one/v
                 w = v*uu
                 fs = f*vv
                 f2 = abssq(fs)
                 h2 = f2*w**2 + g2
              else
                 ! otherwise use the same scaling for f and g.
                 w = one
                 fs = f*uu
                 f2 = abssq(fs)
                 h2 = f2 + g2
              end if
              if (f2 > rtmin .and. h2 < rtmax) then
                 d = sqrt(f2*h2)
              else
                 d = sqrt(f2)*sqrt(h2)
              end if
              p = 1/d
              c = (f2*p)*w
              s = conjg(gs)*(fs*p)
              r = (fs*(h2*p))*u
           end if
        end if
        a = r
        return
     end subroutine stdlib_wrotg

     ! WSCAL scales a vector by a constant.

     subroutine stdlib_wscal(n, za, zx, incx)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: za
           integer(ilp) :: incx, n
           ! .. array arguments ..
           complex(qp) :: zx(*)
        ! =====================================================================
           ! .. local scalars ..
           integer(ilp) :: i, nincx
           if (n <= 0 .or. incx <= 0) return
           if (incx == 1) then
              ! code for increment equal to 1
              do i = 1, n
                 zx(i) = za*zx(i)
              end do
           else
              ! code for increment not equal to 1
              nincx = n*incx
              do i = 1, nincx, incx
                 zx(i) = za*zx(i)
              end do
           end if
           return
           ! end of stdlib_wscal
     end subroutine stdlib_wscal

     ! WSWAP interchanges two vectors.

     subroutine stdlib_wswap(n, zx, incx, zy, incy)
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, incy, n
           ! .. array arguments ..
           complex(qp) :: zx(*), zy(*)
        ! =====================================================================
           ! .. local scalars ..
           complex(qp) :: ztemp
           integer(ilp) :: i, ix, iy
           if (n <= 0) return
           if (incx == 1 .and. incy == 1) then
             ! code for both increments equal to 1
              do i = 1, n
                 ztemp = zx(i)
                 zx(i) = zy(i)
                 zy(i) = ztemp
              end do
           else
             ! code for unequal increments or equal increments not equal
               ! to 1
              ix = 1
              iy = 1
              if (incx < 0) ix = (-n + 1)*incx + 1
              if (incy < 0) iy = (-n + 1)*incy + 1
              do i = 1, n
                 ztemp = zx(ix)
                 zx(ix) = zy(iy)
                 zy(iy) = ztemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
           ! end of stdlib_wswap
     end subroutine stdlib_wswap

     ! WSYMM  performs one of the matrix-matrix operations
     ! C := alpha*A*B + beta*C,
     ! or
     ! C := alpha*B*A + beta*C,
     ! where  alpha and beta are scalars, A is a symmetric matrix and  B and
     ! C are m by n matrices.

     subroutine stdlib_wsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: lda, ldb, ldc, m, n
           character :: side, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           complex(qp) :: temp1, temp2
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
               call stdlib_xerbla('stdlib_wsymm ', info)
               return
           end if
           ! quick return if possible.
           if ((m == 0) .or. (n == 0) .or. ((alpha == czero) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (beta == czero) then
                   do j = 1, n
                       do i = 1, m
                           c(i, j) = czero
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
                           temp2 = czero
                           do k = 1, i - 1
                               c(k, j) = c(k, j) + temp1*a(k, i)
                               temp2 = temp2 + b(k, j)*a(k, i)
                           end do
                           if (beta == czero) then
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
                           temp2 = czero
                           do k = i + 1, m
                               c(k, j) = c(k, j) + temp1*a(k, i)
                               temp2 = temp2 + b(k, j)*a(k, i)
                           end do
                           if (beta == czero) then
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
                   if (beta == czero) then
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
           ! end of stdlib_wsymm
     end subroutine stdlib_wsymm

     ! WSYR2K  performs one of the symmetric rank 2k operations
     ! C := alpha*A*B**T + alpha*B*A**T + beta*C,
     ! or
     ! C := alpha*A**T*B + alpha*B**T*A + beta*C,
     ! where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     ! and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     ! matrices in the second case.

     subroutine stdlib_wsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldb, ldc, n
           character :: trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           complex(qp) :: temp1, temp2
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
           else if ((.not. stdlib_lsame(trans, 'n')) .and. (.not. stdlib_lsame(trans, 't'))) &
                     then
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
               call stdlib_xerbla('stdlib_wsyr2k', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = czero
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
                   if (beta == czero) then
                       do j = 1, n
                           do i = j, n
                               c(i, j) = czero
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
                       if (beta == czero) then
                           do i = 1, j
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1, j
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if ((a(j, l) /= czero) .or. (b(j, l) /= czero)) then
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
                       if (beta == czero) then
                           do i = j, n
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if ((a(j, l) /= czero) .or. (b(j, l) /= czero)) then
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
                           temp1 = czero
                           temp2 = czero
                           do l = 1, k
                               temp1 = temp1 + a(l, i)*b(l, j)
                               temp2 = temp2 + b(l, i)*a(l, j)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp1 + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = j, n
                           temp1 = czero
                           temp2 = czero
                           do l = 1, k
                               temp1 = temp1 + a(l, i)*b(l, j)
                               temp2 = temp2 + b(l, i)*a(l, j)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp1 + alpha*temp2
                           else
                               c(i, j) = beta*c(i, j) + alpha*temp1 + alpha*temp2
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_wsyr2k
     end subroutine stdlib_wsyr2k

     ! WSYRK  performs one of the symmetric rank k operations
     ! C := alpha*A*A**T + beta*C,
     ! or
     ! C := alpha*A**T*A + beta*C,
     ! where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     ! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     ! in the second case.

     subroutine stdlib_wsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha, beta
           integer(ilp) :: k, lda, ldc, n
           character :: trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), c(ldc, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: max
           ! .. local scalars ..
           complex(qp) :: temp
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
           else if ((.not. stdlib_lsame(trans, 'n')) .and. (.not. stdlib_lsame(trans, 't'))) &
                     then
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
               call stdlib_xerbla('stdlib_wsyrk ', info)
               return
           end if
           ! quick return if possible.
           if ((n == 0) .or. (((alpha == czero) .or. (k == 0)) .and. (beta == cone))) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               if (upper) then
                   if (beta == czero) then
                       do j = 1, n
                           do i = 1, j
                               c(i, j) = czero
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
                   if (beta == czero) then
                       do j = 1, n
                           do i = j, n
                               c(i, j) = czero
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
                       if (beta == czero) then
                           do i = 1, j
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = 1, j
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if (a(j, l) /= czero) then
                               temp = alpha*a(j, l)
                               do i = 1, j
                                   c(i, j) = c(i, j) + temp*a(i, l)
                               end do
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       if (beta == czero) then
                           do i = j, n
                               c(i, j) = czero
                           end do
                       else if (beta /= cone) then
                           do i = j, n
                               c(i, j) = beta*c(i, j)
                           end do
                       end if
                       do l = 1, k
                           if (a(j, l) /= czero) then
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
                           temp = czero
                           do l = 1, k
                               temp = temp + a(l, i)*a(l, j)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               else
                   do j = 1, n
                       do i = j, n
                           temp = czero
                           do l = 1, k
                               temp = temp + a(l, i)*a(l, j)
                           end do
                           if (beta == czero) then
                               c(i, j) = alpha*temp
                           else
                               c(i, j) = alpha*temp + beta*c(i, j)
                           end if
                       end do
                   end do
               end if
           end if
           return
           ! end of stdlib_wsyrk
     end subroutine stdlib_wsyrk

     ! WTBMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular band matrix, with ( k + 1 ) diagonals.

     subroutine stdlib_wtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, k, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kplus1, kx, l
           logical(lk) :: noconj, nounit
           ! .. intrinsic functions ..
           intrinsic :: conjg, max, min
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
               call stdlib_xerbla('stdlib_wtbmv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = stdlib_lsame(trans, 't')
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (stdlib_lsame(trans, 'n')) then
               ! form  x := a*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
              ! form  x := a**t*x  or  x := a**h*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1, j)
                               do i = j - 1, max(1, j - k), -1
                                   temp = temp + a(l + i, j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1, j))
                               do i = j - 1, max(1, j - k), -1
                                   temp = temp + conjg(a(l + i, j))*x(i)
                               end do
                           end if
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
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1, j)
                               do i = j - 1, max(1, j - k), -1
                                   temp = temp + a(l + i, j)*x(ix)
                                   ix = ix - incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(kplus1, j))
                               do i = j - 1, max(1, j - k), -1
                                   temp = temp + conjg(a(l + i, j))*x(ix)
                                   ix = ix - incx
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1, j)
                               do i = j + 1, min(n, j + k)
                                   temp = temp + a(l + i, j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1, j))
                               do i = j + 1, min(n, j + k)
                                   temp = temp + conjg(a(l + i, j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1, j)
                               do i = j + 1, min(n, j + k)
                                   temp = temp + a(l + i, j)*x(ix)
                                   ix = ix + incx
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(1, j))
                               do i = j + 1, min(n, j + k)
                                   temp = temp + conjg(a(l + i, j))*x(ix)
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
           ! end of stdlib_wtbmv
     end subroutine stdlib_wtbmv

     ! WTBSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular band matrix, with ( k + 1 )
     ! diagonals.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.

     subroutine stdlib_wtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, k, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kplus1, kx, l
           logical(lk) :: noconj, nounit
           ! .. intrinsic functions ..
           intrinsic :: conjg, max, min
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
               call stdlib_xerbla('stdlib_wtbsv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = stdlib_lsame(trans, 't')
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with cone pass through a.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := inv( a )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kplus1 = k + 1
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1, j - k), j - 1
                                   temp = temp - a(l + i, j)*x(i)
                               end do
                               if (nounit) temp = temp/a(kplus1, j)
                           else
                               do i = max(1, j - k), j - 1
                                   temp = temp - conjg(a(l + i, j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1, j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               do i = max(1, j - k), j - 1
                                   temp = temp - a(l + i, j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(kplus1, j)
                           else
                               do i = max(1, j - k), j - 1
                                   temp = temp - conjg(a(l + i, j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(kplus1, j))
                           end if
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
                           if (noconj) then
                               do i = min(n, j + k), j + 1, -1
                                   temp = temp - a(l + i, j)*x(i)
                               end do
                               if (nounit) temp = temp/a(1, j)
                           else
                               do i = min(n, j + k), j + 1, -1
                                   temp = temp - conjg(a(l + i, j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(1, j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               do i = min(n, j + k), j + 1, -1
                                   temp = temp - a(l + i, j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(1, j)
                           else
                               do i = min(n, j + k), j + 1, -1
                                   temp = temp - conjg(a(l + i, j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(1, j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n - j) >= k) kx = kx - incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_wtbsv
     end subroutine stdlib_wtbsv

     ! WTPMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular matrix, supplied in packed form.

     subroutine stdlib_wtpmv(uplo, trans, diag, n, ap, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           complex(qp) :: ap(*), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, k, kk, kx
           logical(lk) :: noconj, nounit
           ! .. intrinsic functions ..
           intrinsic :: conjg
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
               call stdlib_xerbla('stdlib_wtpmv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = stdlib_lsame(trans, 't')
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x:= a*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
              ! form  x := a**t*x  or  x := a**h*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j - 1, 1, -1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j - 1, 1, -1
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k - 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk - j
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk - 1, kk - j + 1, -1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk - 1, kk - j + 1, -1
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
                       do j = 1, n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do i = j + 1, n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do i = j + 1, n
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k + 1
                               end do
                           end if
                           x(j) = temp
                           kk = kk + (n - j + 1)
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do k = kk + 1, kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do k = kk + 1, kk + n - j
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
           ! end of stdlib_wtpmv
     end subroutine stdlib_wtpmv

     ! WTPSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular matrix, supplied in packed form.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.

     subroutine stdlib_wtpsv(uplo, trans, diag, n, ap, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           complex(qp) :: ap(*), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, k, kk, kx
           logical(lk) :: noconj, nounit
           ! .. intrinsic functions ..
           intrinsic :: conjg
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
               call stdlib_xerbla('stdlib_wtpsv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = stdlib_lsame(trans, 't')
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with cone pass through ap.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := inv( a )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kk = (n*(n + 1))/2
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   kk = 1
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = 1, j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do i = 1, j - 1
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
                       do j = 1, n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk, kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/ap(kk + j - 1)
                           else
                               do k = kk, kk + j - 2
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
                       do j = n, 1, -1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do i = n, j + 1, -1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do i = n, j + 1, -1
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
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do k = kk, kk - (n - (j + 1)), -1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/ap(kk - n + j)
                           else
                               do k = kk, kk - (n - (j + 1)), -1
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
           ! end of stdlib_wtpsv
     end subroutine stdlib_wtpsv

     ! WTRMM  performs one of the matrix-matrix operations
     ! B := alpha*op( A )*B,   or   B := alpha*B*op( A )
     ! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     ! op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.

     subroutine stdlib_wtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           integer(ilp) :: lda, ldb, m, n
           character :: diag, side, transa, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, j, k, nrowa
           logical(lk) :: lside, noconj, nounit, upper

           ! test the input parameters.
           lside = stdlib_lsame(side, 'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
           noconj = stdlib_lsame(transa, 't')
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
               call stdlib_xerbla('stdlib_wtrmm ', info)
               return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               do j = 1, n
                   do i = 1, m
                       b(i, j) = czero
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
                               if (b(k, j) /= czero) then
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
                               if (b(k, j) /= czero) then
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
                 ! form  b := alpha*a**t*b   or   b := alpha*a**h*b.
                   if (upper) then
                       do j = 1, n
                           do i = m, 1, -1
                               temp = b(i, j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i, i)
                                   do k = 1, i - 1
                                       temp = temp + a(k, i)*b(k, j)
                                   end do
                               else
                                   if (nounit) temp = temp*conjg(a(i, i))
                                   do k = 1, i - 1
                                       temp = temp + conjg(a(k, i))*b(k, j)
                                   end do
                               end if
                               b(i, j) = alpha*temp
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = 1, m
                               temp = b(i, j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i, i)
                                   do k = i + 1, m
                                       temp = temp + a(k, i)*b(k, j)
                                   end do
                               else
                                   if (nounit) temp = temp*conjg(a(i, i))
                                   do k = i + 1, m
                                       temp = temp + conjg(a(k, i))*b(k, j)
                                   end do
                               end if
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
                               if (a(k, j) /= czero) then
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
                               if (a(k, j) /= czero) then
                                   temp = alpha*a(k, j)
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                       end do
                   end if
               else
                 ! form  b := alpha*b*a**t   or   b := alpha*b*a**h.
                   if (upper) then
                       loop_280: do k = 1, n
                           do j = 1, k - 1
                               if (a(j, k) /= czero) then
                                   if (noconj) then
                                       temp = alpha*a(j, k)
                                   else
                                       temp = alpha*conjg(a(j, k))
                                   end if
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k, k)
                               else
                                   temp = temp*conjg(a(k, k))
                               end if
                           end if
                           if (temp /= cone) then
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                       end do loop_280
                   else
                       loop_320: do k = n, 1, -1
                           do j = k + 1, n
                               if (a(j, k) /= czero) then
                                   if (noconj) then
                                       temp = alpha*a(j, k)
                                   else
                                       temp = alpha*conjg(a(j, k))
                                   end if
                                   do i = 1, m
                                       b(i, j) = b(i, j) + temp*b(i, k)
                                   end do
                               end if
                           end do
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k, k)
                               else
                                   temp = temp*conjg(a(k, k))
                               end if
                           end if
                           if (temp /= cone) then
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                       end do loop_320
                   end if
               end if
           end if
           return
           ! end of stdlib_wtrmm
     end subroutine stdlib_wtrmm

     ! WTRMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular matrix.

     subroutine stdlib_wtrmv(uplo, trans, diag, n, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kx
           logical(lk) :: noconj, nounit
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
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
               call stdlib_xerbla('stdlib_wtrmv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = stdlib_lsame(trans, 't')
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := a*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = 1, n
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
              ! form  x := a**t*x  or  x := a**h*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j, j)
                               do i = j - 1, 1, -1
                                   temp = temp + a(i, j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j, j))
                               do i = j - 1, 1, -1
                                   temp = temp + conjg(a(i, j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx + (n - 1)*incx
                       do j = n, 1, -1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j, j)
                               do i = j - 1, 1, -1
                                   ix = ix - incx
                                   temp = temp + a(i, j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j, j))
                               do i = j - 1, 1, -1
                                   ix = ix - incx
                                   temp = temp + conjg(a(i, j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j, j)
                               do i = j + 1, n
                                   temp = temp + a(i, j)*x(i)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j, j))
                               do i = j + 1, n
                                   temp = temp + conjg(a(i, j))*x(i)
                               end do
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j, j)
                               do i = j + 1, n
                                   ix = ix + incx
                                   temp = temp + a(i, j)*x(ix)
                               end do
                           else
                               if (nounit) temp = temp*conjg(a(j, j))
                               do i = j + 1, n
                                   ix = ix + incx
                                   temp = temp + conjg(a(i, j))*x(ix)
                               end do
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_wtrmv
     end subroutine stdlib_wtrmv

     ! WTRSM  solves one of the matrix equations
     ! op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
     ! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     ! op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
     ! The matrix X is overwritten on B.

     subroutine stdlib_wtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           complex(qp) :: alpha
           integer(ilp) :: lda, ldb, m, n
           character :: diag, side, transa, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), b(ldb, *)
        ! =====================================================================
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, j, k, nrowa
           logical(lk) :: lside, noconj, nounit, upper

           ! test the input parameters.
           lside = stdlib_lsame(side, 'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
           noconj = stdlib_lsame(transa, 't')
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
               call stdlib_xerbla('stdlib_wtrsm ', info)
               return
           end if
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           ! and when  alpha.eq.czero.
           if (alpha == czero) then
               do j = 1, n
                   do i = 1, m
                       b(i, j) = czero
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
                           if (alpha /= cone) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = m, 1, -1
                               if (b(k, j) /= czero) then
                                   if (nounit) b(k, j) = b(k, j)/a(k, k)
                                   do i = 1, k - 1
                                       b(i, j) = b(i, j) - b(k, j)*a(i, k)
                                   end do
                               end if
                           end do
                       end do
                   else
                       do j = 1, n
                           if (alpha /= cone) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = 1, m
                               if (b(k, j) /= czero) then
                                   if (nounit) b(k, j) = b(k, j)/a(k, k)
                                   do i = k + 1, m
                                       b(i, j) = b(i, j) - b(k, j)*a(i, k)
                                   end do
                               end if
                           end do
                       end do
                   end if
               else
                 ! form  b := alpha*inv( a**t )*b
                 ! or    b := alpha*inv( a**h )*b.
                   if (upper) then
                       do j = 1, n
                           do i = 1, m
                               temp = alpha*b(i, j)
                               if (noconj) then
                                   do k = 1, i - 1
                                       temp = temp - a(k, i)*b(k, j)
                                   end do
                                   if (nounit) temp = temp/a(i, i)
                               else
                                   do k = 1, i - 1
                                       temp = temp - conjg(a(k, i))*b(k, j)
                                   end do
                                   if (nounit) temp = temp/conjg(a(i, i))
                               end if
                               b(i, j) = temp
                           end do
                       end do
                   else
                       do j = 1, n
                           do i = m, 1, -1
                               temp = alpha*b(i, j)
                               if (noconj) then
                                   do k = i + 1, m
                                       temp = temp - a(k, i)*b(k, j)
                                   end do
                                   if (nounit) temp = temp/a(i, i)
                               else
                                   do k = i + 1, m
                                       temp = temp - conjg(a(k, i))*b(k, j)
                                   end do
                                   if (nounit) temp = temp/conjg(a(i, i))
                               end if
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
                           if (alpha /= cone) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = 1, j - 1
                               if (a(k, j) /= czero) then
                                   do i = 1, m
                                       b(i, j) = b(i, j) - a(k, j)*b(i, k)
                                   end do
                               end if
                           end do
                           if (nounit) then
                               temp = cone/a(j, j)
                               do i = 1, m
                                   b(i, j) = temp*b(i, j)
                               end do
                           end if
                       end do
                   else
                       do j = n, 1, -1
                           if (alpha /= cone) then
                               do i = 1, m
                                   b(i, j) = alpha*b(i, j)
                               end do
                           end if
                           do k = j + 1, n
                               if (a(k, j) /= czero) then
                                   do i = 1, m
                                       b(i, j) = b(i, j) - a(k, j)*b(i, k)
                                   end do
                               end if
                           end do
                           if (nounit) then
                               temp = cone/a(j, j)
                               do i = 1, m
                                   b(i, j) = temp*b(i, j)
                               end do
                           end if
                       end do
                   end if
               else
                 ! form  b := alpha*b*inv( a**t )
                 ! or    b := alpha*b*inv( a**h ).
                   if (upper) then
                       loop_330: do k = n, 1, -1
                           if (nounit) then
                               if (noconj) then
                                   temp = cone/a(k, k)
                               else
                                   temp = cone/conjg(a(k, k))
                               end if
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                           do j = 1, k - 1
                               if (a(j, k) /= czero) then
                                   if (noconj) then
                                       temp = a(j, k)
                                   else
                                       temp = conjg(a(j, k))
                                   end if
                                   do i = 1, m
                                       b(i, j) = b(i, j) - temp*b(i, k)
                                   end do
                               end if
                           end do
                           if (alpha /= cone) then
                               do i = 1, m
                                   b(i, k) = alpha*b(i, k)
                               end do
                           end if
                       end do loop_330
                   else
                       loop_380: do k = 1, n
                           if (nounit) then
                               if (noconj) then
                                   temp = cone/a(k, k)
                               else
                                   temp = cone/conjg(a(k, k))
                               end if
                               do i = 1, m
                                   b(i, k) = temp*b(i, k)
                               end do
                           end if
                           do j = k + 1, n
                               if (a(j, k) /= czero) then
                                   if (noconj) then
                                       temp = a(j, k)
                                   else
                                       temp = conjg(a(j, k))
                                   end if
                                   do i = 1, m
                                       b(i, j) = b(i, j) - temp*b(i, k)
                                   end do
                               end if
                           end do
                           if (alpha /= cone) then
                               do i = 1, m
                                   b(i, k) = alpha*b(i, k)
                               end do
                           end if
                       end do loop_380
                   end if
               end if
           end if
           return
           ! end of stdlib_wtrsm
     end subroutine stdlib_wtrsm

     ! WTRSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular matrix.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.

     subroutine stdlib_wtrsv(uplo, trans, diag, n, a, lda, x, incx)
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! .. scalar arguments ..
           integer(ilp) :: incx, lda, n
           character :: diag, trans, uplo
           ! .. array arguments ..
           complex(qp) :: a(lda, *), x(*)
        ! =====================================================================

           ! .. local scalars ..
           complex(qp) :: temp
           integer(ilp) :: i, info, ix, j, jx, kx
           logical(lk) :: noconj, nounit
           ! .. intrinsic functions ..
           intrinsic :: conjg, max
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
               call stdlib_xerbla('stdlib_wtrsv ', info)
               return
           end if
           ! quick return if possible.
           if (n == 0) return
           noconj = stdlib_lsame(trans, 't')
           nounit = stdlib_lsame(diag, 'n')
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
           if (incx <= 0) then
               kx = 1 - (n - 1)*incx
           else if (incx /= 1) then
               kx = 1
           end if
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with cone pass through a.
           if (stdlib_lsame(trans, 'n')) then
              ! form  x := inv( a )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = n, 1, -1
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
                           if (x(j) /= czero) then
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
                           if (x(jx) /= czero) then
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
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
               if (stdlib_lsame(uplo, 'u')) then
                   if (incx == 1) then
                       do j = 1, n
                           temp = x(j)
                           if (noconj) then
                               do i = 1, j - 1
                                   temp = temp - a(i, j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j, j)
                           else
                               do i = 1, j - 1
                                   temp = temp - conjg(a(i, j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j, j))
                           end if
                           x(j) = temp
                       end do
                   else
                       jx = kx
                       do j = 1, n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = 1, j - 1
                                   temp = temp - a(i, j)*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/a(j, j)
                           else
                               do i = 1, j - 1
                                   temp = temp - conjg(a(i, j))*x(ix)
                                   ix = ix + incx
                               end do
                               if (nounit) temp = temp/conjg(a(j, j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do
                   end if
               else
                   if (incx == 1) then
                       do j = n, 1, -1
                           temp = x(j)
                           if (noconj) then
                               do i = n, j + 1, -1
                                   temp = temp - a(i, j)*x(i)
                               end do
                               if (nounit) temp = temp/a(j, j)
                           else
                               do i = n, j + 1, -1
                                   temp = temp - conjg(a(i, j))*x(i)
                               end do
                               if (nounit) temp = temp/conjg(a(j, j))
                           end if
                           x(j) = temp
                       end do
                   else
                       kx = kx + (n - 1)*incx
                       jx = kx
                       do j = n, 1, -1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do i = n, j + 1, -1
                                   temp = temp - a(i, j)*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/a(j, j)
                           else
                               do i = n, j + 1, -1
                                   temp = temp - conjg(a(i, j))*x(ix)
                                   ix = ix - incx
                               end do
                               if (nounit) temp = temp/conjg(a(j, j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do
                   end if
               end if
           end if
           return
           ! end of stdlib_wtrsv
     end subroutine stdlib_wtrsv

end module stdlib_linalg_blas_w