module stdlib_linalg_blas_z
     use stdlib_linalg_constants
     use stdlib_linalg_blas_aux
     use stdlib_linalg_blas_s
     use stdlib_linalg_blas_d
     use stdlib_linalg_blas_c
     implicit none(type,external)
     private






     public :: sp,dp,lk,int32,int64
     public :: stdlib_zaxpy
     public :: stdlib_zcopy
     public :: stdlib_zdotc
     public :: stdlib_zdotu
     public :: stdlib_zdrot
     public :: stdlib_zdscal
     public :: stdlib_zgbmv
     public :: stdlib_zgemm
     public :: stdlib_zgemv
     public :: stdlib_zgerc
     public :: stdlib_zgeru
     public :: stdlib_zhbmv
     public :: stdlib_zhemm
     public :: stdlib_zhemv
     public :: stdlib_zher
     public :: stdlib_zher2
     public :: stdlib_zher2k
     public :: stdlib_zherk
     public :: stdlib_zhpmv
     public :: stdlib_zhpr
     public :: stdlib_zhpr2
     public :: stdlib_zrotg
     public :: stdlib_zscal
     public :: stdlib_zswap
     public :: stdlib_zsymm
     public :: stdlib_zsyr2k
     public :: stdlib_zsyrk
     public :: stdlib_ztbmv
     public :: stdlib_ztbsv
     public :: stdlib_ztpmv
     public :: stdlib_ztpsv
     public :: stdlib_ztrmm
     public :: stdlib_ztrmv
     public :: stdlib_ztrsm
     public :: stdlib_ztrsv


     contains
     
     
     ! ZAXPY constant times a vector plus a vector.
     subroutine stdlib_zaxpy(n,za,zx,incx,zy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) za
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*),zy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy
           ! ..
     
     
           if (n<=0) return
           if (stdlib_dcabs1(za)==0.0d0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 zy(i) = zy(i) + za*zx(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 zy(iy) = zy(iy) + za*zx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
     
           return
     
           ! end of stdlib_zaxpy
     
     end subroutine stdlib_zaxpy
     
     
     ! ZCOPY copies a vector, x, to a vector, y.
     subroutine stdlib_zcopy(n,zx,incx,zy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*),zy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
               zy(i) = zx(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 zy(iy) = zx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_zcopy
     
     end subroutine stdlib_zcopy
     
     
     ! ZDOTC forms the dot product of two complex vectors
     ! ZDOTC = X^H * Y
     complex(dp) function stdlib_zdotc(n,zx,incx,zy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*),zy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           complex(dp) ztemp
           integer(int32) i,ix,iy
           ! ..
           ! .. intrinsic functions ..
           intrinsic dconjg
           ! ..
           ztemp = (0.0d0,0.0d0)
           stdlib_zdotc = (0.0d0,0.0d0)
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 ztemp = ztemp + dconjg(zx(i))*zy(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 ztemp = ztemp + dconjg(zx(ix))*zy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_zdotc = ztemp
           return
     
           ! end of stdlib_zdotc
     
     end function stdlib_zdotc
     
     
     ! ZDOTU forms the dot product of two complex vectors
     ! ZDOTU = X^T * Y
     complex(dp) function stdlib_zdotu(n,zx,incx,zy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*),zy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           complex(dp) ztemp
           integer(int32) i,ix,iy
           ! ..
           ztemp = (0.0d0,0.0d0)
           stdlib_zdotu = (0.0d0,0.0d0)
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 ztemp = ztemp + zx(i)*zy(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 ztemp = ztemp + zx(ix)*zy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_zdotu = ztemp
           return
     
           ! end of stdlib_zdotu
     
     end function stdlib_zdotu
     
     
     ! Applies a plane rotation, where the cos and sin (c and s) are real
     ! and the vectors cx and cy are complex.
     ! jack dongarra, linpack, 3/11/78.
     subroutine stdlib_zdrot( n, zx, incx, zy, incy, c, s )
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32)            incx, incy, n
           real(dp)   c, s
           ! ..
           ! .. array arguments ..
           complex(dp)         zx( * ), zy( * )
           ! ..
     
       ! =====================================================================
     
           ! .. local scalars ..
           integer(int32)            i, ix, iy
           complex(dp)         ctemp
           ! ..
           ! .. executable statements ..
     
           if( n<=0 )return
           if( incx==1 .and. incy==1 ) then
     
              ! code for both increments equal to 1
     
              do i = 1, n
                 ctemp = c*zx( i ) + s*zy( i )
                 zy( i ) = c*zy( i ) - s*zx( i )
                 zx( i ) = ctemp
              end do
           else
     
              ! code for unequal increments or equal increments not equal
                ! to 1
     
              ix = 1
              iy = 1
              if( incx<0 )ix = ( -n+1 )*incx + 1
              if( incy<0 )iy = ( -n+1 )*incy + 1
              do i = 1, n
                 ctemp = c*zx( ix ) + s*zy( iy )
                 zy( iy ) = c*zy( iy ) - s*zx( ix )
                 zx( ix ) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_zdrot
     
     end subroutine stdlib_zdrot
     
     
     ! ZDSCAL scales a vector by a constant.
     subroutine stdlib_zdscal(n,da,zx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) da
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic dcmplx
           ! ..
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              do i = 1,n
                 zx(i) = dcmplx(da,0.0d0)*zx(i)
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 zx(i) = dcmplx(da,0.0d0)*zx(i)
              end do
           end if
           return
     
           ! end of stdlib_zdscal
     
     end subroutine stdlib_zdscal
     
     
     ! ZGBMV  performs one of the matrix-vector operations
     ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
     ! y := alpha*A**H*x + beta*y,
     ! where alpha and beta are scalars, x and y are vectors and A is an
     ! m by n band matrix, with kl sub-diagonals and ku super-diagonals.
     subroutine stdlib_zgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) incx,incy,kl,ku,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,k,kup1,kx,ky,lenx,leny
           logical(lk) noconj
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 1
           else if (m<0) then
               info = 2
           else if (n<0) then
               info = 3
           else if (kl<0) then
               info = 4
           else if (ku<0) then
               info = 5
           else if (lda< (kl+ku+1)) then
               info = 8
           else if (incx==0) then
               info = 10
           else if (incy==0) then
               info = 13
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zgbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           noconj = stdlib_lsame(trans,'t')
     
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
     
           if (stdlib_lsame(trans,'n')) then
               lenx = n
               leny = m
           else
               lenx = m
               leny = n
           end if
           if (incx>0) then
               kx = 1
           else
               kx = 1 - (lenx-1)*incx
           end if
           if (incy>0) then
               ky = 1
           else
               ky = 1 - (leny-1)*incy
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the band part of a.
     
           ! first form  y := beta*y.
     
           if (beta/=one) then
               if (incy==1) then
                   if (beta==zero) then
                       loop_10: do i = 1,leny
                           y(i) = zero
                       end do loop_10
                   else
                       loop_20: do i = 1,leny
                           y(i) = beta*y(i)
                       end do loop_20
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       loop_30: do i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
                       end do loop_30
                   else
                       loop_40: do i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do loop_40
                   end if
               end if
           end if
           if (alpha==zero) return
           kup1 = ku + 1
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   loop_60: do j = 1,n
                       temp = alpha*x(jx)
                       k = kup1 - j
                       loop_50: do i = max(1,j-ku),min(m,j+kl)
                           y(i) = y(i) + temp*a(k+i,j)
                       end do loop_50
                       jx = jx + incx
                   end do loop_60
               else
                   loop_80: do j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       k = kup1 - j
                       loop_70: do i = max(1,j-ku),min(m,j+kl)
                           y(iy) = y(iy) + temp*a(k+i,j)
                           iy = iy + incy
                       end do loop_70
                       jx = jx + incx
                       if (j>ku) ky = ky + incy
                   end do loop_80
               end if
           else
     
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
     
               jy = ky
               if (incx==1) then
                   loop_110: do j = 1,n
                       temp = zero
                       k = kup1 - j
                       if (noconj) then
                           loop_90: do i = max(1,j-ku),min(m,j+kl)
                               temp = temp + a(k+i,j)*x(i)
                           end do loop_90
                       else
                           loop_100: do i = max(1,j-ku),min(m,j+kl)
                               temp = temp + dconjg(a(k+i,j))*x(i)
                           end do loop_100
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do loop_110
               else
                   loop_140: do j = 1,n
                       temp = zero
                       ix = kx
                       k = kup1 - j
                       if (noconj) then
                           loop_120: do i = max(1,j-ku),min(m,j+kl)
                               temp = temp + a(k+i,j)*x(ix)
                               ix = ix + incx
                           end do loop_120
                       else
                           loop_130: do i = max(1,j-ku),min(m,j+kl)
                               temp = temp + dconjg(a(k+i,j))*x(ix)
                               ix = ix + incx
                           end do loop_130
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j>ku) kx = kx + incx
                   end do loop_140
               end if
           end if
     
           return
     
           ! end of stdlib_zgbmv
     
     end subroutine stdlib_zgbmv
     
     
     ! ZGEMM  performs one of the matrix-matrix operations
     ! C := alpha*op( A )*op( B ) + beta*C,
     ! where  op( X ) is one of
     ! op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
     ! alpha and beta are scalars, and A, B and C are matrices, with op( A )
     ! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
     subroutine stdlib_zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) k,lda,ldb,ldc,m,n
           character transa,transb
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,j,l,nrowa,nrowb
           logical(lk) conja,conjb,nota,notb
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and nrowb  as the number of rows  of  a  and  b  respectively.
     
           nota = stdlib_lsame(transa,'n')
           notb = stdlib_lsame(transb,'n')
           conja = stdlib_lsame(transa,'c')
           conjb = stdlib_lsame(transb,'c')
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
           if ((.not.nota) .and. (.not.conja) .and.(.not.stdlib_lsame(transa,'t'))) then
               info = 1
           else if ((.not.notb) .and. (.not.conjb) .and.(.not.stdlib_lsame(transb,'t'))) &
                     then
               info = 2
           else if (m<0) then
               info = 3
           else if (n<0) then
               info = 4
           else if (k<0) then
               info = 5
           else if (lda<max(1,nrowa)) then
               info = 8
           else if (ldb<max(1,nrowb)) then
               info = 10
           else if (ldc<max(1,m)) then
               info = 13
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zgemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.(((alpha==zero).or. (k==0)).and. (beta==one))) &
                     return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   loop_20: do j = 1,n
                       loop_10: do i = 1,m
                           c(i,j) = zero
                       end do loop_10
                   end do loop_20
               else
                   loop_40: do j = 1,n
                       loop_30: do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do loop_30
                   end do loop_40
               end if
               return
           end if
     
           ! start the operations.
     
           if (notb) then
               if (nota) then
     
                 ! form  c := alpha*a*b + beta*c.
     
                   loop_90: do j = 1,n
                       if (beta==zero) then
                           loop_50: do i = 1,m
                               c(i,j) = zero
                           end do loop_50
                       else if (beta/=one) then
                           loop_60: do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do loop_60
                       end if
                       loop_80: do l = 1,k
                           temp = alpha*b(l,j)
                           loop_70: do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do loop_70
                       end do loop_80
                   end do loop_90
               else if (conja) then
     
                 ! form  c := alpha*a**h*b + beta*c.
     
                   loop_120: do j = 1,n
                       loop_110: do i = 1,m
                           temp = zero
                           loop_100: do l = 1,k
                               temp = temp + dconjg(a(l,i))*b(l,j)
                           end do loop_100
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_110
                   end do loop_120
               else
     
                 ! form  c := alpha*a**t*b + beta*c
     
                   loop_150: do j = 1,n
                       loop_140: do i = 1,m
                           temp = zero
                           loop_130: do l = 1,k
                               temp = temp + a(l,i)*b(l,j)
                           end do loop_130
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_140
                   end do loop_150
               end if
           else if (nota) then
               if (conjb) then
     
                 ! form  c := alpha*a*b**h + beta*c.
     
                   loop_200: do j = 1,n
                       if (beta==zero) then
                           loop_160: do i = 1,m
                               c(i,j) = zero
                           end do loop_160
                       else if (beta/=one) then
                           loop_170: do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do loop_170
                       end if
                       loop_190: do l = 1,k
                           temp = alpha*dconjg(b(j,l))
                           loop_180: do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do loop_180
                       end do loop_190
                   end do loop_200
               else
     
                 ! form  c := alpha*a*b**t + beta*c
     
                   loop_250: do j = 1,n
                       if (beta==zero) then
                           loop_210: do i = 1,m
                               c(i,j) = zero
                           end do loop_210
                       else if (beta/=one) then
                           loop_220: do i = 1,m
                               c(i,j) = beta*c(i,j)
                           end do loop_220
                       end if
                       loop_240: do l = 1,k
                           temp = alpha*b(j,l)
                           loop_230: do i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
                           end do loop_230
                       end do loop_240
                   end do loop_250
               end if
           else if (conja) then
               if (conjb) then
     
                 ! form  c := alpha*a**h*b**h + beta*c.
     
                   loop_280: do j = 1,n
                       loop_270: do i = 1,m
                           temp = zero
                           loop_260: do l = 1,k
                               temp = temp + dconjg(a(l,i))*dconjg(b(j,l))
                           end do loop_260
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_270
                   end do loop_280
               else
     
                 ! form  c := alpha*a**h*b**t + beta*c
     
                   loop_310: do j = 1,n
                       loop_300: do i = 1,m
                           temp = zero
                           loop_290: do l = 1,k
                               temp = temp + dconjg(a(l,i))*b(j,l)
                           end do loop_290
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_300
                   end do loop_310
               end if
           else
               if (conjb) then
     
                 ! form  c := alpha*a**t*b**h + beta*c
     
                   loop_340: do j = 1,n
                       loop_330: do i = 1,m
                           temp = zero
                           loop_320: do l = 1,k
                               temp = temp + a(l,i)*dconjg(b(j,l))
                           end do loop_320
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_330
                   end do loop_340
               else
     
                 ! form  c := alpha*a**t*b**t + beta*c
     
                   loop_370: do j = 1,n
                       loop_360: do i = 1,m
                           temp = zero
                           loop_350: do l = 1,k
                               temp = temp + a(l,i)*b(j,l)
                           end do loop_350
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_360
                   end do loop_370
               end if
           end if
     
           return
     
           ! end of stdlib_zgemm
     
     end subroutine stdlib_zgemm
     
     
     ! ZGEMV  performs one of the matrix-vector operations
     ! y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
     ! y := alpha*A**H*x + beta*y,
     ! where alpha and beta are scalars, x and y are vectors and A is an
     ! m by n matrix.
     subroutine stdlib_zgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) incx,incy,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
           logical(lk) noconj
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 1
           else if (m<0) then
               info = 2
           else if (n<0) then
               info = 3
           else if (lda<max(1,m)) then
               info = 6
           else if (incx==0) then
               info = 8
           else if (incy==0) then
               info = 11
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zgemv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           noconj = stdlib_lsame(trans,'t')
     
           ! set  lenx  and  leny, the lengths of the vectors x and y, and set
           ! up the start points in  x  and  y.
     
           if (stdlib_lsame(trans,'n')) then
               lenx = n
               leny = m
           else
               lenx = m
               leny = n
           end if
           if (incx>0) then
               kx = 1
           else
               kx = 1 - (lenx-1)*incx
           end if
           if (incy>0) then
               ky = 1
           else
               ky = 1 - (leny-1)*incy
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
     
           ! first form  y := beta*y.
     
           if (beta/=one) then
               if (incy==1) then
                   if (beta==zero) then
                       loop_10: do i = 1,leny
                           y(i) = zero
                       end do loop_10
                   else
                       loop_20: do i = 1,leny
                           y(i) = beta*y(i)
                       end do loop_20
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       loop_30: do i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
                       end do loop_30
                   else
                       loop_40: do i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do loop_40
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   loop_60: do j = 1,n
                       temp = alpha*x(jx)
                       loop_50: do i = 1,m
                           y(i) = y(i) + temp*a(i,j)
                       end do loop_50
                       jx = jx + incx
                   end do loop_60
               else
                   loop_80: do j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       loop_70: do i = 1,m
                           y(iy) = y(iy) + temp*a(i,j)
                           iy = iy + incy
                       end do loop_70
                       jx = jx + incx
                   end do loop_80
               end if
           else
     
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
     
               jy = ky
               if (incx==1) then
                   loop_110: do j = 1,n
                       temp = zero
                       if (noconj) then
                           loop_90: do i = 1,m
                               temp = temp + a(i,j)*x(i)
                           end do loop_90
                       else
                           loop_100: do i = 1,m
                               temp = temp + dconjg(a(i,j))*x(i)
                           end do loop_100
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do loop_110
               else
                   loop_140: do j = 1,n
                       temp = zero
                       ix = kx
                       if (noconj) then
                           loop_120: do i = 1,m
                               temp = temp + a(i,j)*x(ix)
                               ix = ix + incx
                           end do loop_120
                       else
                           loop_130: do i = 1,m
                               temp = temp + dconjg(a(i,j))*x(ix)
                               ix = ix + incx
                           end do loop_130
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                   end do loop_140
               end if
           end if
     
           return
     
           ! end of stdlib_zgemv
     
     end subroutine stdlib_zgemv
     
     
     ! ZGERC  performs the rank 1 operation
     ! A := alpha*x*y**H + A,
     ! where alpha is a scalar, x is an m element vector, y is an n element
     ! vector and A is an m by n matrix.
     subroutine stdlib_zgerc(m,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           integer(int32) incx,incy,lda,m,n
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jy,kx
           ! ..
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (m<0) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 5
           else if (incy==0) then
               info = 7
           else if (lda<max(1,m)) then
               info = 9
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zgerc ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or. (alpha==zero)) return
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
     
           if (incy>0) then
               jy = 1
           else
               jy = 1 - (n-1)*incy
           end if
           if (incx==1) then
               loop_20: do j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*dconjg(y(jy))
                       loop_10: do i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
                       end do loop_10
                   end if
                   jy = jy + incy
               end do loop_20
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               loop_40: do j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*dconjg(y(jy))
                       ix = kx
                       loop_30: do i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
                       end do loop_30
                   end if
                   jy = jy + incy
               end do loop_40
           end if
     
           return
     
           ! end of stdlib_zgerc
     
     end subroutine stdlib_zgerc
     
     
     ! ZGERU  performs the rank 1 operation
     ! A := alpha*x*y**T + A,
     ! where alpha is a scalar, x is an m element vector, y is an n element
     ! vector and A is an m by n matrix.
     subroutine stdlib_zgeru(m,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           integer(int32) incx,incy,lda,m,n
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jy,kx
           ! ..
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (m<0) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 5
           else if (incy==0) then
               info = 7
           else if (lda<max(1,m)) then
               info = 9
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zgeru ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or. (alpha==zero)) return
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
     
           if (incy>0) then
               jy = 1
           else
               jy = 1 - (n-1)*incy
           end if
           if (incx==1) then
               loop_20: do j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       loop_10: do i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
                       end do loop_10
                   end if
                   jy = jy + incy
               end do loop_20
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               loop_40: do j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       ix = kx
                       loop_30: do i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
                       end do loop_30
                   end if
                   jy = jy + incy
               end do loop_40
           end if
     
           return
     
           ! end of stdlib_zgeru
     
     end subroutine stdlib_zgeru
     
     
     ! ZHBMV  performs the matrix-vector  operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n hermitian band matrix, with k super-diagonals.
     subroutine stdlib_zhbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) incx,incy,k,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (k<0) then
               info = 3
           else if (lda< (k+1)) then
               info = 6
           else if (incx==0) then
               info = 8
           else if (incy==0) then
               info = 11
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zhbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. ((alpha==zero).and. (beta==one))) return
     
           ! set up the start points in  x  and  y.
     
           if (incx>0) then
               kx = 1
           else
               kx = 1 - (n-1)*incx
           end if
           if (incy>0) then
               ky = 1
           else
               ky = 1 - (n-1)*incy
           end if
     
           ! start the operations. in this version the elements of the array a
           ! are accessed sequentially with one pass through a.
     
           ! first form  y := beta*y.
     
           if (beta/=one) then
               if (incy==1) then
                   if (beta==zero) then
                       loop_10: do i = 1,n
                           y(i) = zero
                       end do loop_10
                   else
                       loop_20: do i = 1,n
                           y(i) = beta*y(i)
                       end do loop_20
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       loop_30: do i = 1,n
                           y(iy) = zero
                           iy = iy + incy
                       end do loop_30
                   else
                       loop_40: do i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do loop_40
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when upper triangle of a is stored.
     
               kplus1 = k + 1
               if ((incx==1) .and. (incy==1)) then
                   loop_60: do j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       l = kplus1 - j
                       loop_50: do i = max(1,j-k),j - 1
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(i)
                       end do loop_50
                       y(j) = y(j) + temp1*dble(a(kplus1,j)) + alpha*temp2
                   end do loop_60
               else
                   jx = kx
                   jy = ky
                   loop_80: do j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       loop_70: do i = max(1,j-k),j - 1
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do loop_70
                       y(jy) = y(jy) + temp1*dble(a(kplus1,j)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       if (j>k) then
                           kx = kx + incx
                           ky = ky + incy
                       end if
                   end do loop_80
               end if
           else
     
              ! form  y  when lower triangle of a is stored.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_100: do j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*dble(a(1,j))
                       l = 1 - j
                       loop_90: do i = j + 1,min(n,j+k)
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(i)
                       end do loop_90
                       y(j) = y(j) + alpha*temp2
                   end do loop_100
               else
                   jx = kx
                   jy = ky
                   loop_120: do j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*dble(a(1,j))
                       l = 1 - j
                       ix = jx
                       iy = jy
                       loop_110: do i = j + 1,min(n,j+k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(ix)
                       end do loop_110
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do loop_120
               end if
           end if
     
           return
     
           ! end of stdlib_zhbmv
     
     end subroutine stdlib_zhbmv
     
     
     ! ZHEMM  performs one of the matrix-matrix operations
     ! C := alpha*A*B + beta*C,
     ! or
     ! C := alpha*B*A + beta*C,
     ! where alpha and beta are scalars, A is an hermitian matrix and  B and
     ! C are m by n matrices.
     subroutine stdlib_zhemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) lda,ldb,ldc,m,n
           character side,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg,max
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,j,k,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! set nrowa as the number of rows of a.
     
           if (stdlib_lsame(side,'l')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = stdlib_lsame(uplo,'u')
     
           ! test the input parameters.
     
           info = 0
           if ((.not.stdlib_lsame(side,'l')) .and. (.not.stdlib_lsame(side,'r'))) then
               info = 1
           else if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 2
           else if (m<0) then
               info = 3
           else if (n<0) then
               info = 4
           else if (lda<max(1,nrowa)) then
               info = 7
           else if (ldb<max(1,m)) then
               info = 9
           else if (ldc<max(1,m)) then
               info = 12
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zhemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   loop_20: do j = 1,n
                       loop_10: do i = 1,m
                           c(i,j) = zero
                       end do loop_10
                   end do loop_20
               else
                   loop_40: do j = 1,n
                       loop_30: do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do loop_30
                   end do loop_40
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   loop_70: do j = 1,n
                       loop_60: do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           loop_50: do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*dconjg(a(k,i))
                           end do loop_50
                           if (beta==zero) then
                               c(i,j) = temp1*dble(a(i,i)) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*dble(a(i,i)) +alpha*temp2
                           end if
                       end do loop_60
                   end do loop_70
               else
                   loop_100: do j = 1,n
                       loop_90: do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           loop_80: do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*dconjg(a(k,i))
                           end do loop_80
                           if (beta==zero) then
                               c(i,j) = temp1*dble(a(i,i)) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*dble(a(i,i)) +alpha*temp2
                           end if
                       end do loop_90
                   end do loop_100
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               loop_170: do j = 1,n
                   temp1 = alpha*dble(a(j,j))
                   if (beta==zero) then
                       loop_110: do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do loop_110
                   else
                       loop_120: do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do loop_120
                   end if
                   loop_140: do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*dconjg(a(j,k))
                       end if
                       loop_130: do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do loop_130
                   end do loop_140
                   loop_160: do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*dconjg(a(j,k))
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       loop_150: do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do loop_150
                   end do loop_160
               end do loop_170
           end if
     
           return
     
           ! end of stdlib_zhemm
     
     end subroutine stdlib_zhemm
     
     
     ! ZHEMV  performs the matrix-vector  operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n hermitian matrix.
     subroutine stdlib_zhemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (lda<max(1,n)) then
               info = 5
           else if (incx==0) then
               info = 7
           else if (incy==0) then
               info = 10
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zhemv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. ((alpha==zero).and. (beta==one))) return
     
           ! set up the start points in  x  and  y.
     
           if (incx>0) then
               kx = 1
           else
               kx = 1 - (n-1)*incx
           end if
           if (incy>0) then
               ky = 1
           else
               ky = 1 - (n-1)*incy
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the triangular part
           ! of a.
     
           ! first form  y := beta*y.
     
           if (beta/=one) then
               if (incy==1) then
                   if (beta==zero) then
                       loop_10: do i = 1,n
                           y(i) = zero
                       end do loop_10
                   else
                       loop_20: do i = 1,n
                           y(i) = beta*y(i)
                       end do loop_20
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       loop_30: do i = 1,n
                           y(iy) = zero
                           iy = iy + incy
                       end do loop_30
                   else
                       loop_40: do i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do loop_40
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when a is stored in upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_60: do j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       loop_50: do i = 1,j - 1
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(i)
                       end do loop_50
                       y(j) = y(j) + temp1*dble(a(j,j)) + alpha*temp2
                   end do loop_60
               else
                   jx = kx
                   jy = ky
                   loop_80: do j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       loop_70: do i = 1,j - 1
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do loop_70
                       y(jy) = y(jy) + temp1*dble(a(j,j)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do loop_80
               end if
           else
     
              ! form  y  when a is stored in lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_100: do j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*dble(a(j,j))
                       loop_90: do i = j + 1,n
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(i)
                       end do loop_90
                       y(j) = y(j) + alpha*temp2
                   end do loop_100
               else
                   jx = kx
                   jy = ky
                   loop_120: do j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*dble(a(j,j))
                       ix = jx
                       iy = jy
                       loop_110: do i = j + 1,n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(ix)
                       end do loop_110
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                   end do loop_120
               end if
           end if
     
           return
     
           ! end of stdlib_zhemv
     
     end subroutine stdlib_zhemv
     
     
     ! ZHER   performs the hermitian rank 1 operation
     ! A := alpha*x*x**H + A,
     ! where alpha is a real scalar, x is an n element vector and A is an
     ! n by n hermitian matrix.
     subroutine stdlib_zher(uplo,n,alpha,x,incx,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,kx
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 5
           else if (lda<max(1,n)) then
               info = 7
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zher  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==dble(zero))) return
     
           ! set the start point in x if the increment is not unity.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the triangular part
           ! of a.
     
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  a  when a is stored in upper triangle.
     
               if (incx==1) then
                   loop_20: do j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           loop_10: do i = 1,j - 1
                               a(i,j) = a(i,j) + x(i)*temp
                           end do loop_10
                           a(j,j) = dble(a(j,j)) + dble(x(j)*temp)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                   end do loop_20
               else
                   jx = kx
                   loop_40: do j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           ix = kx
                           loop_30: do i = 1,j - 1
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
                           end do loop_30
                           a(j,j) = dble(a(j,j)) + dble(x(jx)*temp)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
                   end do loop_40
               end if
           else
     
              ! form  a  when a is stored in lower triangle.
     
               if (incx==1) then
                   loop_60: do j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           a(j,j) = dble(a(j,j)) + dble(temp*x(j))
                           loop_50: do i = j + 1,n
                               a(i,j) = a(i,j) + x(i)*temp
                           end do loop_50
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                   end do loop_60
               else
                   jx = kx
                   loop_80: do j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           a(j,j) = dble(a(j,j)) + dble(temp*x(jx))
                           ix = jx
                           loop_70: do i = j + 1,n
                               ix = ix + incx
                               a(i,j) = a(i,j) + x(ix)*temp
                           end do loop_70
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
                   end do loop_80
               end if
           end if
     
           return
     
           ! end of stdlib_zher
     
     end subroutine stdlib_zher
     
     
     ! ZHER2  performs the hermitian rank 2 operation
     ! A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
     ! where alpha is a scalar, x and y are n element vectors and A is an n
     ! by n hermitian matrix.
     subroutine stdlib_zher2(uplo,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 5
           else if (incy==0) then
               info = 7
           else if (lda<max(1,n)) then
               info = 9
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zher2 ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==zero)) return
     
           ! set up the start points in x and y if the increments are not both
           ! unity.
     
           if ((incx/=1) .or. (incy/=1)) then
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (n-1)*incx
               end if
               if (incy>0) then
                   ky = 1
               else
                   ky = 1 - (n-1)*incy
               end if
               jx = kx
               jy = ky
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through the triangular part
           ! of a.
     
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  a  when a is stored in the upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_20: do j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           loop_10: do i = 1,j - 1
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
                           end do loop_10
                           a(j,j) = dble(a(j,j)) +dble(x(j)*temp1+y(j)*temp2)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                   end do loop_20
               else
                   loop_40: do j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           loop_30: do i = 1,j - 1
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do loop_30
                           a(j,j) = dble(a(j,j)) +dble(x(jx)*temp1+y(jy)*temp2)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                   end do loop_40
               end if
           else
     
              ! form  a  when a is stored in the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_60: do j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           a(j,j) = dble(a(j,j)) +dble(x(j)*temp1+y(j)*temp2)
                           loop_50: do i = j + 1,n
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
                           end do loop_50
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                   end do loop_60
               else
                   loop_80: do j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           a(j,j) = dble(a(j,j)) +dble(x(jx)*temp1+y(jy)*temp2)
                           ix = jx
                           iy = jy
                           loop_70: do i = j + 1,n
                               ix = ix + incx
                               iy = iy + incy
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                           end do loop_70
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                   end do loop_80
               end if
           end if
     
           return
     
           ! end of stdlib_zher2
     
     end subroutine stdlib_zher2
     
     
     ! ZHER2K  performs one of the hermitian rank 2k operations
     ! C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
     ! or
     ! C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
     ! where  alpha and beta  are scalars with  beta  real,  C is an  n by n
     ! hermitian matrix and  A and B  are  n by k matrices in the first case
     ! and  k by n  matrices in the second case.
     subroutine stdlib_zher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           real(dp) beta
           integer(int32) k,lda,ldb,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg,max
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(dp) one
           parameter (one=1.0_dp)
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! test the input parameters.
     
           if (stdlib_lsame(trans,'n')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = stdlib_lsame(uplo,'u')
     
           info = 0
           if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 1
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'c'))) &
                     then
               info = 2
           else if (n<0) then
               info = 3
           else if (k<0) then
               info = 4
           else if (lda<max(1,nrowa)) then
               info = 7
           else if (ldb<max(1,nrowa)) then
               info = 9
           else if (ldc<max(1,n)) then
               info = 12
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zher2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==dble(zero)) then
                       loop_20: do j = 1,n
                           loop_10: do i = 1,j
                               c(i,j) = zero
                           end do loop_10
                       end do loop_20
                   else
                       loop_40: do j = 1,n
                           loop_30: do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do loop_30
                           c(j,j) = beta*dble(c(j,j))
                       end do loop_40
                   end if
               else
                   if (beta==dble(zero)) then
                       loop_60: do j = 1,n
                           loop_50: do i = j,n
                               c(i,j) = zero
                           end do loop_50
                       end do loop_60
                   else
                       loop_80: do j = 1,n
                           c(j,j) = beta*dble(c(j,j))
                           loop_70: do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do loop_70
                       end do loop_80
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
     
               if (upper) then
                   loop_130: do j = 1,n
                       if (beta==dble(zero)) then
                           loop_90: do i = 1,j
                               c(i,j) = zero
                           end do loop_90
                       else if (beta/=one) then
                           loop_100: do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do loop_100
                           c(j,j) = beta*dble(c(j,j))
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       loop_120: do l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*dconjg(b(j,l))
                               temp2 = dconjg(alpha*a(j,l))
                               loop_110: do i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
                               end do loop_110
                               c(j,j) = dble(c(j,j)) +dble(a(j,l)*temp1+b(j,l)*temp2)
                           end if
                       end do loop_120
                   end do loop_130
               else
                   loop_180: do j = 1,n
                       if (beta==dble(zero)) then
                           loop_140: do i = j,n
                               c(i,j) = zero
                           end do loop_140
                       else if (beta/=one) then
                           loop_150: do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do loop_150
                           c(j,j) = beta*dble(c(j,j))
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       loop_170: do l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*dconjg(b(j,l))
                               temp2 = dconjg(alpha*a(j,l))
                               loop_160: do i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
                               end do loop_160
                               c(j,j) = dble(c(j,j)) +dble(a(j,l)*temp1+b(j,l)*temp2)
                           end if
                       end do loop_170
                   end do loop_180
               end if
           else
     
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
     
               if (upper) then
                   loop_210: do j = 1,n
                       loop_200: do i = 1,j
                           temp1 = zero
                           temp2 = zero
                           loop_190: do l = 1,k
                               temp1 = temp1 + dconjg(a(l,i))*b(l,j)
                               temp2 = temp2 + dconjg(b(l,i))*a(l,j)
                           end do loop_190
                           if (i==j) then
                               if (beta==dble(zero)) then
                                   c(j,j) = dble(alpha*temp1+dconjg(alpha)*temp2)
                               else
                                   c(j,j) = beta*dble(c(j,j)) +dble(alpha*temp1+dconjg(alpha)&
                                             *temp2)
                               end if
                           else
                               if (beta==dble(zero)) then
                                   c(i,j) = alpha*temp1 + dconjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 +dconjg(alpha)*temp2
                               end if
                           end if
                       end do loop_200
                   end do loop_210
               else
                   loop_240: do j = 1,n
                       loop_230: do i = j,n
                           temp1 = zero
                           temp2 = zero
                           loop_220: do l = 1,k
                               temp1 = temp1 + dconjg(a(l,i))*b(l,j)
                               temp2 = temp2 + dconjg(b(l,i))*a(l,j)
                           end do loop_220
                           if (i==j) then
                               if (beta==dble(zero)) then
                                   c(j,j) = dble(alpha*temp1+dconjg(alpha)*temp2)
                               else
                                   c(j,j) = beta*dble(c(j,j)) +dble(alpha*temp1+dconjg(alpha)&
                                             *temp2)
                               end if
                           else
                               if (beta==dble(zero)) then
                                   c(i,j) = alpha*temp1 + dconjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 +dconjg(alpha)*temp2
                               end if
                           end if
                       end do loop_230
                   end do loop_240
               end if
           end if
     
           return
     
           ! end of stdlib_zher2k
     
     end subroutine stdlib_zher2k
     
     
     ! ZHERK  performs one of the hermitian rank k operations
     ! C := alpha*A*A**H + beta*C,
     ! or
     ! C := alpha*A**H*A + beta*C,
     ! where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
     ! matrix and  A  is an  n by k  matrix in the  first case and a  k by n
     ! matrix in the second case.
     subroutine stdlib_zherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) k,lda,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dcmplx,dconjg,max
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           real(dp) rtemp
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
     
           ! test the input parameters.
     
           if (stdlib_lsame(trans,'n')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = stdlib_lsame(uplo,'u')
     
           info = 0
           if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 1
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'c'))) &
                     then
               info = 2
           else if (n<0) then
               info = 3
           else if (k<0) then
               info = 4
           else if (lda<max(1,nrowa)) then
               info = 7
           else if (ldc<max(1,n)) then
               info = 10
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zherk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       loop_20: do j = 1,n
                           loop_10: do i = 1,j
                               c(i,j) = zero
                           end do loop_10
                       end do loop_20
                   else
                       loop_40: do j = 1,n
                           loop_30: do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do loop_30
                           c(j,j) = beta*dble(c(j,j))
                       end do loop_40
                   end if
               else
                   if (beta==zero) then
                       loop_60: do j = 1,n
                           loop_50: do i = j,n
                               c(i,j) = zero
                           end do loop_50
                       end do loop_60
                   else
                       loop_80: do j = 1,n
                           c(j,j) = beta*dble(c(j,j))
                           loop_70: do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do loop_70
                       end do loop_80
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**h + beta*c.
     
               if (upper) then
                   loop_130: do j = 1,n
                       if (beta==zero) then
                           loop_90: do i = 1,j
                               c(i,j) = zero
                           end do loop_90
                       else if (beta/=one) then
                           loop_100: do i = 1,j - 1
                               c(i,j) = beta*c(i,j)
                           end do loop_100
                           c(j,j) = beta*dble(c(j,j))
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       loop_120: do l = 1,k
                           if (a(j,l)/=dcmplx(zero)) then
                               temp = alpha*dconjg(a(j,l))
                               loop_110: do i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do loop_110
                               c(j,j) = dble(c(j,j)) + dble(temp*a(i,l))
                           end if
                       end do loop_120
                   end do loop_130
               else
                   loop_180: do j = 1,n
                       if (beta==zero) then
                           loop_140: do i = j,n
                               c(i,j) = zero
                           end do loop_140
                       else if (beta/=one) then
                           c(j,j) = beta*dble(c(j,j))
                           loop_150: do i = j + 1,n
                               c(i,j) = beta*c(i,j)
                           end do loop_150
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       loop_170: do l = 1,k
                           if (a(j,l)/=dcmplx(zero)) then
                               temp = alpha*dconjg(a(j,l))
                               c(j,j) = dble(c(j,j)) + dble(temp*a(j,l))
                               loop_160: do i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do loop_160
                           end if
                       end do loop_170
                   end do loop_180
               end if
           else
     
              ! form  c := alpha*a**h*a + beta*c.
     
               if (upper) then
                   loop_220: do j = 1,n
                       loop_200: do i = 1,j - 1
                           temp = zero
                           loop_190: do l = 1,k
                               temp = temp + dconjg(a(l,i))*a(l,j)
                           end do loop_190
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_200
                       rtemp = zero
                       loop_210: do l = 1,k
                           rtemp = rtemp + dconjg(a(l,j))*a(l,j)
                       end do loop_210
                       if (beta==zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*dble(c(j,j))
                       end if
                   end do loop_220
               else
                   loop_260: do j = 1,n
                       rtemp = zero
                       loop_230: do l = 1,k
                           rtemp = rtemp + dconjg(a(l,j))*a(l,j)
                       end do loop_230
                       if (beta==zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*dble(c(j,j))
                       end if
                       loop_250: do i = j + 1,n
                           temp = zero
                           loop_240: do l = 1,k
                               temp = temp + dconjg(a(l,i))*a(l,j)
                           end do loop_240
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_250
                   end do loop_260
               end if
           end if
     
           return
     
           ! end of stdlib_zherk
     
     end subroutine stdlib_zherk
     
     
     ! ZHPMV  performs the matrix-vector operation
     ! y := alpha*A*x + beta*y,
     ! where alpha and beta are scalars, x and y are n element vectors and
     ! A is an n by n hermitian matrix, supplied in packed form.
     subroutine stdlib_zhpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 6
           else if (incy==0) then
               info = 9
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zhpmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. ((alpha==zero).and. (beta==one))) return
     
           ! set up the start points in  x  and  y.
     
           if (incx>0) then
               kx = 1
           else
               kx = 1 - (n-1)*incx
           end if
           if (incy>0) then
               ky = 1
           else
               ky = 1 - (n-1)*incy
           end if
     
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with one pass through ap.
     
           ! first form  y := beta*y.
     
           if (beta/=one) then
               if (incy==1) then
                   if (beta==zero) then
                       loop_10: do i = 1,n
                           y(i) = zero
                       end do loop_10
                   else
                       loop_20: do i = 1,n
                           y(i) = beta*y(i)
                       end do loop_20
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       loop_30: do i = 1,n
                           y(iy) = zero
                           iy = iy + incy
                       end do loop_30
                   else
                       loop_40: do i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
                       end do loop_40
                   end if
               end if
           end if
           if (alpha==zero) return
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when ap contains the upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_60: do j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       k = kk
                       loop_50: do i = 1,j - 1
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(i)
                           k = k + 1
                       end do loop_50
                       y(j) = y(j) + temp1*dble(ap(kk+j-1)) + alpha*temp2
                       kk = kk + j
                   end do loop_60
               else
                   jx = kx
                   jy = ky
                   loop_80: do j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       loop_70: do k = kk,kk + j - 2
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
                       end do loop_70
                       y(jy) = y(jy) + temp1*dble(ap(kk+j-1)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
                   end do loop_80
               end if
           else
     
              ! form  y  when ap contains the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_100: do j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*dble(ap(kk))
                       k = kk + 1
                       loop_90: do i = j + 1,n
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(i)
                           k = k + 1
                       end do loop_90
                       y(j) = y(j) + alpha*temp2
                       kk = kk + (n-j+1)
                   end do loop_100
               else
                   jx = kx
                   jy = ky
                   loop_120: do j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*dble(ap(kk))
                       ix = jx
                       iy = jy
                       loop_110: do k = kk + 1,kk + n - j
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(ix)
                       end do loop_110
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + (n-j+1)
                   end do loop_120
               end if
           end if
     
           return
     
           ! end of stdlib_zhpmv
     
     end subroutine stdlib_zhpmv
     
     
     ! ZHPR    performs the hermitian rank 1 operation
     ! A := alpha*x*x**H + A,
     ! where alpha is a real scalar, x is an n element vector and A is an
     ! n by n hermitian matrix, supplied in packed form.
     subroutine stdlib_zhpr(uplo,n,alpha,x,incx,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 5
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zhpr  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==dble(zero))) return
     
           ! set the start point in x if the increment is not unity.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with one pass through ap.
     
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  a  when upper triangle is stored in ap.
     
               if (incx==1) then
                   loop_20: do j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           k = kk
                           loop_10: do i = 1,j - 1
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
                           end do loop_10
                           ap(kk+j-1) = dble(ap(kk+j-1)) + dble(x(j)*temp)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       kk = kk + j
                   end do loop_20
               else
                   jx = kx
                   loop_40: do j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           ix = kx
                           loop_30: do k = kk,kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
                           end do loop_30
                           ap(kk+j-1) = dble(ap(kk+j-1)) + dble(x(jx)*temp)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       jx = jx + incx
                       kk = kk + j
                   end do loop_40
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if (incx==1) then
                   loop_60: do j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           ap(kk) = dble(ap(kk)) + dble(temp*x(j))
                           k = kk + 1
                           loop_50: do i = j + 1,n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
                           end do loop_50
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       kk = kk + n - j + 1
                   end do loop_60
               else
                   jx = kx
                   loop_80: do j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           ap(kk) = dble(ap(kk)) + dble(temp*x(jx))
                           ix = jx
                           loop_70: do k = kk + 1,kk + n - j
                               ix = ix + incx
                               ap(k) = ap(k) + x(ix)*temp
                           end do loop_70
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
                   end do loop_80
               end if
           end if
     
           return
     
           ! end of stdlib_zhpr
     
     end subroutine stdlib_zhpr
     
     
     ! ZHPR2  performs the hermitian rank 2 operation
     ! A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
     ! where alpha is a scalar, x and y are n element vectors and A is an
     ! n by n hermitian matrix, supplied in packed form.
     subroutine stdlib_zhpr2(uplo,n,alpha,x,incx,y,incy,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(dp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dble,dconjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (n<0) then
               info = 2
           else if (incx==0) then
               info = 5
           else if (incy==0) then
               info = 7
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zhpr2 ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==zero)) return
     
           ! set up the start points in x and y if the increments are not both
           ! unity.
     
           if ((incx/=1) .or. (incy/=1)) then
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (n-1)*incx
               end if
               if (incy>0) then
                   ky = 1
               else
                   ky = 1 - (n-1)*incy
               end if
               jx = kx
               jy = ky
           end if
     
           ! start the operations. in this version the elements of the array ap
           ! are accessed sequentially with one pass through ap.
     
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  a  when upper triangle is stored in ap.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_20: do j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           k = kk
                           loop_10: do i = 1,j - 1
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
                           end do loop_10
                           ap(kk+j-1) = dble(ap(kk+j-1)) +dble(x(j)*temp1+y(j)*temp2)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       kk = kk + j
                   end do loop_20
               else
                   loop_40: do j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           loop_30: do k = kk,kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
                           end do loop_30
                           ap(kk+j-1) = dble(ap(kk+j-1)) +dble(x(jx)*temp1+y(jy)*temp2)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
                   end do loop_40
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if ((incx==1) .and. (incy==1)) then
                   loop_60: do j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           ap(kk) = dble(ap(kk)) +dble(x(j)*temp1+y(j)*temp2)
                           k = kk + 1
                           loop_50: do i = j + 1,n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
                           end do loop_50
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       kk = kk + n - j + 1
                   end do loop_60
               else
                   loop_80: do j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           ap(kk) = dble(ap(kk)) +dble(x(jx)*temp1+y(jy)*temp2)
                           ix = jx
                           iy = jy
                           loop_70: do k = kk + 1,kk + n - j
                               ix = ix + incx
                               iy = iy + incy
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                           end do loop_70
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
                   end do loop_80
               end if
           end if
     
           return
     
           ! end of stdlib_zhpr2
     
     end subroutine stdlib_zhpr2
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
     subroutine stdlib_zrotg( a, b, c, s )
        integer, parameter :: wp = kind(1.d0)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
        ! .. constants ..
        real(dp), parameter ::dzero = 0.0_dp
        real(dp), parameter :: done  = 1.0_dp
        complex(dp), parameter ::zzero  = 0.0_dp
        ! ..
        ! .. scaling constants ..
     real(dp), parameter :: safmin = real(radix(0._dp),wp)**max(minexponent(0._dp)-1,1-&
               maxexponent(0._dp)   )
     real(dp), parameter :: safmax = real(radix(0._dp),wp)**max(1-minexponent(0._dp),maxexponent(&
               0._dp)-1   )
     real(dp), parameter :: rtmin = sqrt( real(radix(0._dp),wp)**max(minexponent(0._dp)-1,1-&
               maxexponent(0._dp)   ) / epsilon(0._dp) )
     real(dp), parameter :: rtmax = sqrt( real(radix(0._dp),wp)**max(1-minexponent(0._dp),&
               maxexponent(0._dp)-1   ) * epsilon(0._dp) )
        ! ..
        ! .. scalar arguments ..
        real(dp) :: c
        complex(dp) :: a, b, s
        ! ..
        ! .. local scalars ..
        real(dp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
        complex(dp) :: f, fs, g, gs, r, t
        ! ..
        ! .. intrinsic functions ..
        intrinsic :: abs, aimag, conjg, max, min, real, sqrt
        ! ..
        ! .. statement functions ..
        real(dp) :: abssq
        ! ..
        ! .. statement function definitions ..
        abssq( t ) = real( t )**2 + aimag( t )**2
        ! ..
        ! .. executable statements ..
     
        f = a
        g = b
        if( g ==zzero ) then
           c = done
           s =zzero
           r = f
        else if( f ==zzero ) then
           c =dzero
           g1 = max( abs(real(g)), abs(aimag(g)) )
           if( g1 > rtmin .and. g1 < rtmax ) then
     
              ! use unscaled algorithm
     
              g2 = abssq( g )
              d = sqrt( g2 )
              s = conjg( g ) / d
              r = d
           else
     
              ! use scaled algorithm
     
              u = min( safmax, max( safmin, g1 ) )
              uu = done / u
              gs = g*uu
              g2 = abssq( gs )
              d = sqrt( g2 )
              s = conjg( gs ) / d
              r = d*u
           end if
        else
           f1 = max( abs(real(f)), abs(aimag(f)) )
           g1 = max( abs(real(g)), abs(aimag(g)) )
     if( f1 > rtmin .and. f1 < rtmax .and.          g1 > rtmin .and. g1 < rtmax ) then
     
              ! use unscaled algorithm
     
              f2 = abssq( f )
              g2 = abssq( g )
              h2 = f2 + g2
              if( f2 > rtmin .and. h2 < rtmax ) then
                 d = sqrt( f2*h2 )
              else
                 d = sqrt( f2 )*sqrt( h2 )
              end if
              p = 1 / d
              c = f2*p
              s = conjg( g )*( f*p )
              r = f*( h2*p )
           else
     
              ! use scaled algorithm
     
              u = min( safmax, max( safmin, f1, g1 ) )
              uu = done / u
              gs = g*uu
              g2 = abssq( gs )
              if( f1*uu < rtmin ) then
     
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
     
                 v = min( safmax, max( safmin, f1 ) )
                 vv = done / v
                 w = v * uu
                 fs = f*vv
                 f2 = abssq( fs )
                 h2 = f2*w**2 + g2
              else
     
                 ! otherwise use the same scaling for f and g.
     
                 w = done
                 fs = f*uu
                 f2 = abssq( fs )
                 h2 = f2 + g2
              end if
              if( f2 > rtmin .and. h2 < rtmax ) then
                 d = sqrt( f2*h2 )
              else
                 d = sqrt( f2 )*sqrt( h2 )
              end if
              p = 1 / d
              c = ( f2*p )*w
              s = conjg( gs )*( fs*p )
              r = ( fs*( h2*p ) )*u
           end if
        end if
        a = r
        return
     end subroutine stdlib_zrotg
     
     
     ! ZSCAL scales a vector by a constant.
     subroutine stdlib_zscal(n,za,zx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) za
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,nincx
           ! ..
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              do i = 1,n
                 zx(i) = za*zx(i)
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 zx(i) = za*zx(i)
              end do
           end if
           return
     
           ! end of stdlib_zscal
     
     end subroutine stdlib_zscal
     
     
     ! ZSWAP interchanges two vectors.
     subroutine stdlib_zswap(n,zx,incx,zy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*),zy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           complex(dp) ztemp
           integer(int32) i,ix,iy
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
             ! code for both increments equal to 1
              do i = 1,n
                 ztemp = zx(i)
                 zx(i) = zy(i)
                 zy(i) = ztemp
              end do
           else
     
             ! code for unequal increments or equal increments not equal
               ! to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 ztemp = zx(ix)
                 zx(ix) = zy(iy)
                 zy(iy) = ztemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_zswap
     
     end subroutine stdlib_zswap
     
     
     ! ZSYMM  performs one of the matrix-matrix operations
     ! C := alpha*A*B + beta*C,
     ! or
     ! C := alpha*B*A + beta*C,
     ! where  alpha and beta are scalars, A is a symmetric matrix and  B and
     ! C are m by n matrices.
     subroutine stdlib_zsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) lda,ldb,ldc,m,n
           character side,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,j,k,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! set nrowa as the number of rows of a.
     
           if (stdlib_lsame(side,'l')) then
               nrowa = m
           else
               nrowa = n
           end if
           upper = stdlib_lsame(uplo,'u')
     
           ! test the input parameters.
     
           info = 0
           if ((.not.stdlib_lsame(side,'l')) .and. (.not.stdlib_lsame(side,'r'))) then
               info = 1
           else if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 2
           else if (m<0) then
               info = 3
           else if (n<0) then
               info = 4
           else if (lda<max(1,nrowa)) then
               info = 7
           else if (ldb<max(1,m)) then
               info = 9
           else if (ldc<max(1,m)) then
               info = 12
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zsymm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   loop_20: do j = 1,n
                       loop_10: do i = 1,m
                           c(i,j) = zero
                       end do loop_10
                   end do loop_20
               else
                   loop_40: do j = 1,n
                       loop_30: do i = 1,m
                           c(i,j) = beta*c(i,j)
                       end do loop_30
                   end do loop_40
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   loop_70: do j = 1,n
                       loop_60: do i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           loop_50: do k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do loop_50
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
                       end do loop_60
                   end do loop_70
               else
                   loop_100: do j = 1,n
                       loop_90: do i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           loop_80: do k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
                           end do loop_80
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
                       end do loop_90
                   end do loop_100
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               loop_170: do j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta==zero) then
                       loop_110: do i = 1,m
                           c(i,j) = temp1*b(i,j)
                       end do loop_110
                   else
                       loop_120: do i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
                       end do loop_120
                   end if
                   loop_140: do k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       loop_130: do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do loop_130
                   end do loop_140
                   loop_160: do k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       loop_150: do i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
                       end do loop_150
                   end do loop_160
               end do loop_170
           end if
     
           return
     
           ! end of stdlib_zsymm
     
     end subroutine stdlib_zsymm
     
     
     ! ZSYR2K  performs one of the symmetric rank 2k operations
     ! C := alpha*A*B**T + alpha*B*A**T + beta*C,
     ! or
     ! C := alpha*A**T*B + alpha*B**T*A + beta*C,
     ! where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     ! and  A and B  are  n by k  matrices  in the  first  case  and  k by n
     ! matrices in the second case.
     subroutine stdlib_zsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) k,lda,ldb,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           complex(dp) temp1,temp2
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! test the input parameters.
     
           if (stdlib_lsame(trans,'n')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = stdlib_lsame(uplo,'u')
     
           info = 0
           if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 1
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t'))) &
                     then
               info = 2
           else if (n<0) then
               info = 3
           else if (k<0) then
               info = 4
           else if (lda<max(1,nrowa)) then
               info = 7
           else if (ldb<max(1,nrowa)) then
               info = 9
           else if (ldc<max(1,n)) then
               info = 12
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zsyr2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       loop_20: do j = 1,n
                           loop_10: do i = 1,j
                               c(i,j) = zero
                           end do loop_10
                       end do loop_20
                   else
                       loop_40: do j = 1,n
                           loop_30: do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do loop_30
                       end do loop_40
                   end if
               else
                   if (beta==zero) then
                       loop_60: do j = 1,n
                           loop_50: do i = j,n
                               c(i,j) = zero
                           end do loop_50
                       end do loop_60
                   else
                       loop_80: do j = 1,n
                           loop_70: do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do loop_70
                       end do loop_80
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
     
               if (upper) then
                   loop_130: do j = 1,n
                       if (beta==zero) then
                           loop_90: do i = 1,j
                               c(i,j) = zero
                           end do loop_90
                       else if (beta/=one) then
                           loop_100: do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do loop_100
                       end if
                       loop_120: do l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               loop_110: do i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
                               end do loop_110
                           end if
                       end do loop_120
                   end do loop_130
               else
                   loop_180: do j = 1,n
                       if (beta==zero) then
                           loop_140: do i = j,n
                               c(i,j) = zero
                           end do loop_140
                       else if (beta/=one) then
                           loop_150: do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do loop_150
                       end if
                       loop_170: do l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               loop_160: do i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
                               end do loop_160
                           end if
                       end do loop_170
                   end do loop_180
               end if
           else
     
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
     
               if (upper) then
                   loop_210: do j = 1,n
                       loop_200: do i = 1,j
                           temp1 = zero
                           temp2 = zero
                           loop_190: do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do loop_190
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
                       end do loop_200
                   end do loop_210
               else
                   loop_240: do j = 1,n
                       loop_230: do i = j,n
                           temp1 = zero
                           temp2 = zero
                           loop_220: do l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
                           end do loop_220
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
                       end do loop_230
                   end do loop_240
               end if
           end if
     
           return
     
           ! end of stdlib_zsyr2k
     
     end subroutine stdlib_zsyr2k
     
     
     ! ZSYRK  performs one of the symmetric rank k operations
     ! C := alpha*A*A**T + beta*C,
     ! or
     ! C := alpha*A**T*A + beta*C,
     ! where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
     ! and  A  is an  n by k  matrix in the first case and a  k by n  matrix
     ! in the second case.
     subroutine stdlib_zsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha,beta
           integer(int32) k,lda,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! test the input parameters.
     
           if (stdlib_lsame(trans,'n')) then
               nrowa = n
           else
               nrowa = k
           end if
           upper = stdlib_lsame(uplo,'u')
     
           info = 0
           if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 1
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t'))) &
                     then
               info = 2
           else if (n<0) then
               info = 3
           else if (k<0) then
               info = 4
           else if (lda<max(1,nrowa)) then
               info = 7
           else if (ldc<max(1,n)) then
               info = 10
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_zsyrk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       loop_20: do j = 1,n
                           loop_10: do i = 1,j
                               c(i,j) = zero
                           end do loop_10
                       end do loop_20
                   else
                       loop_40: do j = 1,n
                           loop_30: do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do loop_30
                       end do loop_40
                   end if
               else
                   if (beta==zero) then
                       loop_60: do j = 1,n
                           loop_50: do i = j,n
                               c(i,j) = zero
                           end do loop_50
                       end do loop_60
                   else
                       loop_80: do j = 1,n
                           loop_70: do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do loop_70
                       end do loop_80
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**t + beta*c.
     
               if (upper) then
                   loop_130: do j = 1,n
                       if (beta==zero) then
                           loop_90: do i = 1,j
                               c(i,j) = zero
                           end do loop_90
                       else if (beta/=one) then
                           loop_100: do i = 1,j
                               c(i,j) = beta*c(i,j)
                           end do loop_100
                       end if
                       loop_120: do l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               loop_110: do i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do loop_110
                           end if
                       end do loop_120
                   end do loop_130
               else
                   loop_180: do j = 1,n
                       if (beta==zero) then
                           loop_140: do i = j,n
                               c(i,j) = zero
                           end do loop_140
                       else if (beta/=one) then
                           loop_150: do i = j,n
                               c(i,j) = beta*c(i,j)
                           end do loop_150
                       end if
                       loop_170: do l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               loop_160: do i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
                               end do loop_160
                           end if
                       end do loop_170
                   end do loop_180
               end if
           else
     
              ! form  c := alpha*a**t*a + beta*c.
     
               if (upper) then
                   loop_210: do j = 1,n
                       loop_200: do i = 1,j
                           temp = zero
                           loop_190: do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do loop_190
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_200
                   end do loop_210
               else
                   loop_240: do j = 1,n
                       loop_230: do i = j,n
                           temp = zero
                           loop_220: do l = 1,k
                               temp = temp + a(l,i)*a(l,j)
                           end do loop_220
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
                       end do loop_230
                   end do loop_240
               end if
           end if
     
           return
     
           ! end of stdlib_zsyrk
     
     end subroutine stdlib_zsyrk
     
     
     ! ZTBMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular band matrix, with ( k + 1 ) diagonals.
     subroutine stdlib_ztbmv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 2
           else if (.not.stdlib_lsame(diag,'u') .and. .not.stdlib_lsame(diag,'n')) then
               info = 3
           else if (n<0) then
               info = 4
           else if (k<0) then
               info = 5
           else if (lda< (k+1)) then
               info = 7
           else if (incx==0) then
               info = 9
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
           noconj = stdlib_lsame(trans,'t')
           nounit = stdlib_lsame(diag,'n')
     
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx   too small for descending loops.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
     
           if (stdlib_lsame(trans,'n')) then
     
               ! form  x := a*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       loop_20: do j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = kplus1 - j
                               loop_10: do i = max(1,j-k),j - 1
                                   x(i) = x(i) + temp*a(l+i,j)
                               end do loop_10
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
                       end do loop_20
                   else
                       jx = kx
                       loop_40: do j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               loop_30: do i = max(1,j-k),j - 1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix + incx
                               end do loop_30
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j>k) kx = kx + incx
                       end do loop_40
                   end if
               else
                   if (incx==1) then
                       loop_60: do j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = 1 - j
                               loop_50: do i = min(n,j+k),j + 1,-1
                                   x(i) = x(i) + temp*a(l+i,j)
                               end do loop_50
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
                       end do loop_60
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_80: do j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               loop_70: do i = min(n,j+k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix - incx
                               end do loop_70
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
                       end do loop_80
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       loop_110: do j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               loop_90: do i = j - 1,max(1,j-k),-1
                                   temp = temp + a(l+i,j)*x(i)
                               end do loop_90
                           else
                               if (nounit) temp = temp*dconjg(a(kplus1,j))
                               loop_100: do i = j - 1,max(1,j-k),-1
                                   temp = temp + dconjg(a(l+i,j))*x(i)
                               end do loop_100
                           end if
                           x(j) = temp
                       end do loop_110
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_140: do j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               loop_120: do i = j - 1,max(1,j-k),-1
                                   temp = temp + a(l+i,j)*x(ix)
                                   ix = ix - incx
                               end do loop_120
                           else
                               if (nounit) temp = temp*dconjg(a(kplus1,j))
                               loop_130: do i = j - 1,max(1,j-k),-1
                                   temp = temp + dconjg(a(l+i,j))*x(ix)
                                   ix = ix - incx
                               end do loop_130
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do loop_140
                   end if
               else
                   if (incx==1) then
                       loop_170: do j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               loop_150: do i = j + 1,min(n,j+k)
                                   temp = temp + a(l+i,j)*x(i)
                               end do loop_150
                           else
                               if (nounit) temp = temp*dconjg(a(1,j))
                               loop_160: do i = j + 1,min(n,j+k)
                                   temp = temp + dconjg(a(l+i,j))*x(i)
                               end do loop_160
                           end if
                           x(j) = temp
                       end do loop_170
                   else
                       jx = kx
                       loop_200: do j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               loop_180: do i = j + 1,min(n,j+k)
                                   temp = temp + a(l+i,j)*x(ix)
                                   ix = ix + incx
                               end do loop_180
                           else
                               if (nounit) temp = temp*dconjg(a(1,j))
                               loop_190: do i = j + 1,min(n,j+k)
                                   temp = temp + dconjg(a(l+i,j))*x(ix)
                                   ix = ix + incx
                               end do loop_190
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do loop_200
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztbmv
     
     end subroutine stdlib_ztbmv
     
     
     ! ZTBSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular band matrix, with ( k + 1 )
     ! diagonals.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.
     subroutine stdlib_ztbsv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 2
           else if (.not.stdlib_lsame(diag,'u') .and. .not.stdlib_lsame(diag,'n')) then
               info = 3
           else if (n<0) then
               info = 4
           else if (k<0) then
               info = 5
           else if (lda< (k+1)) then
               info = 7
           else if (incx==0) then
               info = 9
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztbsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
           noconj = stdlib_lsame(trans,'t')
           nounit = stdlib_lsame(diag,'n')
     
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed by sequentially with one pass through a.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  x := inv( a )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       loop_20: do j = n,1,-1
                           if (x(j)/=zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               loop_10: do i = j - 1,max(1,j-k),-1
                                   x(i) = x(i) - temp*a(l+i,j)
                               end do loop_10
                           end if
                       end do loop_20
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_40: do j = n,1,-1
                           kx = kx - incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               loop_30: do i = j - 1,max(1,j-k),-1
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix - incx
                               end do loop_30
                           end if
                           jx = jx - incx
                       end do loop_40
                   end if
               else
                   if (incx==1) then
                       loop_60: do j = 1,n
                           if (x(j)/=zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               loop_50: do i = j + 1,min(n,j+k)
                                   x(i) = x(i) - temp*a(l+i,j)
                               end do loop_50
                           end if
                       end do loop_60
                   else
                       jx = kx
                       loop_80: do j = 1,n
                           kx = kx + incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               loop_70: do i = j + 1,min(n,j+k)
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix + incx
                               end do loop_70
                           end if
                           jx = jx + incx
                       end do loop_80
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       loop_110: do j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               loop_90: do i = max(1,j-k),j - 1
                                   temp = temp - a(l+i,j)*x(i)
                               end do loop_90
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               loop_100: do i = max(1,j-k),j - 1
                                   temp = temp - dconjg(a(l+i,j))*x(i)
                               end do loop_100
                               if (nounit) temp = temp/dconjg(a(kplus1,j))
                           end if
                           x(j) = temp
                       end do loop_110
                   else
                       jx = kx
                       loop_140: do j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               loop_120: do i = max(1,j-k),j - 1
                                   temp = temp - a(l+i,j)*x(ix)
                                   ix = ix + incx
                               end do loop_120
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               loop_130: do i = max(1,j-k),j - 1
                                   temp = temp - dconjg(a(l+i,j))*x(ix)
                                   ix = ix + incx
                               end do loop_130
                               if (nounit) temp = temp/dconjg(a(kplus1,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           if (j>k) kx = kx + incx
                       end do loop_140
                   end if
               else
                   if (incx==1) then
                       loop_170: do j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               loop_150: do i = min(n,j+k),j + 1,-1
                                   temp = temp - a(l+i,j)*x(i)
                               end do loop_150
                               if (nounit) temp = temp/a(1,j)
                           else
                               loop_160: do i = min(n,j+k),j + 1,-1
                                   temp = temp - dconjg(a(l+i,j))*x(i)
                               end do loop_160
                               if (nounit) temp = temp/dconjg(a(1,j))
                           end if
                           x(j) = temp
                       end do loop_170
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_200: do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               loop_180: do i = min(n,j+k),j + 1,-1
                                   temp = temp - a(l+i,j)*x(ix)
                                   ix = ix - incx
                               end do loop_180
                               if (nounit) temp = temp/a(1,j)
                           else
                               loop_190: do i = min(n,j+k),j + 1,-1
                                   temp = temp - dconjg(a(l+i,j))*x(ix)
                                   ix = ix - incx
                               end do loop_190
                               if (nounit) temp = temp/dconjg(a(1,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
                       end do loop_200
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztbsv
     
     end subroutine stdlib_ztbsv
     
     
     ! ZTPMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular matrix, supplied in packed form.
     subroutine stdlib_ztpmv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 2
           else if (.not.stdlib_lsame(diag,'u') .and. .not.stdlib_lsame(diag,'n')) then
               info = 3
           else if (n<0) then
               info = 4
           else if (incx==0) then
               info = 7
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztpmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
           noconj = stdlib_lsame(trans,'t')
           nounit = stdlib_lsame(diag,'n')
     
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  x:= a*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = 1
                   if (incx==1) then
                       loop_20: do j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               loop_10: do i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
                               end do loop_10
                               if (nounit) x(j) = x(j)*ap(kk+j-1)
                           end if
                           kk = kk + j
                       end do loop_20
                   else
                       jx = kx
                       loop_40: do j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               loop_30: do k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
                               end do loop_30
                               if (nounit) x(jx) = x(jx)*ap(kk+j-1)
                           end if
                           jx = jx + incx
                           kk = kk + j
                       end do loop_40
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       loop_60: do j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               loop_50: do i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
                               end do loop_50
                               if (nounit) x(j) = x(j)*ap(kk-n+j)
                           end if
                           kk = kk - (n-j+1)
                       end do loop_60
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_80: do j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               loop_70: do k = kk,kk - (n- (j+1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
                               end do loop_70
                               if (nounit) x(jx) = x(jx)*ap(kk-n+j)
                           end if
                           jx = jx - incx
                           kk = kk - (n-j+1)
                       end do loop_80
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       loop_110: do j = n,1,-1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               loop_90: do i = j - 1,1,-1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
                               end do loop_90
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               loop_100: do i = j - 1,1,-1
                                   temp = temp + dconjg(ap(k))*x(i)
                                   k = k - 1
                               end do loop_100
                           end if
                           x(j) = temp
                           kk = kk - j
                       end do loop_110
                   else
                       jx = kx + (n-1)*incx
                       loop_140: do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               loop_120: do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
                               end do loop_120
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               loop_130: do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + dconjg(ap(k))*x(ix)
                               end do loop_130
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
                       end do loop_140
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       loop_170: do j = 1,n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               loop_150: do i = j + 1,n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
                               end do loop_150
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               loop_160: do i = j + 1,n
                                   temp = temp + dconjg(ap(k))*x(i)
                                   k = k + 1
                               end do loop_160
                           end if
                           x(j) = temp
                           kk = kk + (n-j+1)
                       end do loop_170
                   else
                       jx = kx
                       loop_200: do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               loop_180: do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
                               end do loop_180
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               loop_190: do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + dconjg(ap(k))*x(ix)
                               end do loop_190
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n-j+1)
                       end do loop_200
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztpmv
     
     end subroutine stdlib_ztpmv
     
     
     ! ZTPSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular matrix, supplied in packed form.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.
     subroutine stdlib_ztpsv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 2
           else if (.not.stdlib_lsame(diag,'u') .and. .not.stdlib_lsame(diag,'n')) then
               info = 3
           else if (n<0) then
               info = 4
           else if (incx==0) then
               info = 7
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztpsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
           noconj = stdlib_lsame(trans,'t')
           nounit = stdlib_lsame(diag,'n')
     
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of ap are
           ! accessed sequentially with one pass through ap.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  x := inv( a )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       loop_20: do j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               loop_10: do i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
                               end do loop_10
                           end if
                           kk = kk - j
                       end do loop_20
                   else
                       jx = kx + (n-1)*incx
                       loop_40: do j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               loop_30: do k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do loop_30
                           end if
                           jx = jx - incx
                           kk = kk - j
                       end do loop_40
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       loop_60: do j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               loop_50: do i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
                               end do loop_50
                           end if
                           kk = kk + (n-j+1)
                       end do loop_60
                   else
                       jx = kx
                       loop_80: do j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               loop_70: do k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
                               end do loop_70
                           end if
                           jx = jx + incx
                           kk = kk + (n-j+1)
                       end do loop_80
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = 1
                   if (incx==1) then
                       loop_110: do j = 1,n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               loop_90: do i = 1,j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
                               end do loop_90
                               if (nounit) temp = temp/ap(kk+j-1)
                           else
                               loop_100: do i = 1,j - 1
                                   temp = temp - dconjg(ap(k))*x(i)
                                   k = k + 1
                               end do loop_100
                               if (nounit) temp = temp/dconjg(ap(kk+j-1))
                           end if
                           x(j) = temp
                           kk = kk + j
                       end do loop_110
                   else
                       jx = kx
                       loop_140: do j = 1,n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               loop_120: do k = kk,kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
                               end do loop_120
                               if (nounit) temp = temp/ap(kk+j-1)
                           else
                               loop_130: do k = kk,kk + j - 2
                                   temp = temp - dconjg(ap(k))*x(ix)
                                   ix = ix + incx
                               end do loop_130
                               if (nounit) temp = temp/dconjg(ap(kk+j-1))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
                       end do loop_140
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       loop_170: do j = n,1,-1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               loop_150: do i = n,j + 1,-1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
                               end do loop_150
                               if (nounit) temp = temp/ap(kk-n+j)
                           else
                               loop_160: do i = n,j + 1,-1
                                   temp = temp - dconjg(ap(k))*x(i)
                                   k = k - 1
                               end do loop_160
                               if (nounit) temp = temp/dconjg(ap(kk-n+j))
                           end if
                           x(j) = temp
                           kk = kk - (n-j+1)
                       end do loop_170
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_200: do j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               loop_180: do k = kk,kk - (n- (j+1)),-1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
                               end do loop_180
                               if (nounit) temp = temp/ap(kk-n+j)
                           else
                               loop_190: do k = kk,kk - (n- (j+1)),-1
                                   temp = temp - dconjg(ap(k))*x(ix)
                                   ix = ix - incx
                               end do loop_190
                               if (nounit) temp = temp/dconjg(ap(kk-n+j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n-j+1)
                       end do loop_200
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztpsv
     
     end subroutine stdlib_ztpsv
     
     
     ! ZTRMM  performs one of the matrix-matrix operations
     ! B := alpha*op( A )*B,   or   B := alpha*B*op( A )
     ! where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
     ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     ! op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
     subroutine stdlib_ztrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,noconj,nounit,upper
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! test the input parameters.
     
           lside = stdlib_lsame(side,'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
           noconj = stdlib_lsame(transa,'t')
           nounit = stdlib_lsame(diag,'n')
           upper = stdlib_lsame(uplo,'u')
     
           info = 0
           if ((.not.lside) .and. (.not.stdlib_lsame(side,'r'))) then
               info = 1
           else if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 2
           else if ((.not.stdlib_lsame(transa,'n')) .and.(.not.stdlib_lsame(transa,'t')) .and.(&
                     .not.stdlib_lsame(transa,'c'))) then
               info = 3
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n'))) &
                     then
               info = 4
           else if (m<0) then
               info = 5
           else if (n<0) then
               info = 6
           else if (lda<max(1,nrowa)) then
               info = 9
           else if (ldb<max(1,m)) then
               info = 11
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztrmm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               loop_20: do j = 1,n
                   loop_10: do i = 1,m
                       b(i,j) = zero
                   end do loop_10
               end do loop_20
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*a*b.
     
                   if (upper) then
                       loop_50: do j = 1,n
                           loop_40: do k = 1,m
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   loop_30: do i = 1,k - 1
                                       b(i,j) = b(i,j) + temp*a(i,k)
                                   end do loop_30
                                   if (nounit) temp = temp*a(k,k)
                                   b(k,j) = temp
                               end if
                           end do loop_40
                       end do loop_50
                   else
                       loop_80: do j = 1,n
                           loop_70: do k = m,1,-1
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   b(k,j) = temp
                                   if (nounit) b(k,j) = b(k,j)*a(k,k)
                                   loop_60: do i = k + 1,m
                                       b(i,j) = b(i,j) + temp*a(i,k)
                                   end do loop_60
                               end if
                           end do loop_70
                       end do loop_80
                   end if
               else
     
                 ! form  b := alpha*a**t*b   or   b := alpha*a**h*b.
     
                   if (upper) then
                       loop_120: do j = 1,n
                           loop_110: do i = m,1,-1
                               temp = b(i,j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i,i)
                                   loop_90: do k = 1,i - 1
                                       temp = temp + a(k,i)*b(k,j)
                                   end do loop_90
                               else
                                   if (nounit) temp = temp*dconjg(a(i,i))
                                   loop_100: do k = 1,i - 1
                                       temp = temp + dconjg(a(k,i))*b(k,j)
                                   end do loop_100
                               end if
                               b(i,j) = alpha*temp
                           end do loop_110
                       end do loop_120
                   else
                       loop_160: do j = 1,n
                           loop_150: do i = 1,m
                               temp = b(i,j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i,i)
                                   loop_130: do k = i + 1,m
                                       temp = temp + a(k,i)*b(k,j)
                                   end do loop_130
                               else
                                   if (nounit) temp = temp*dconjg(a(i,i))
                                   loop_140: do k = i + 1,m
                                       temp = temp + dconjg(a(k,i))*b(k,j)
                                   end do loop_140
                               end if
                               b(i,j) = alpha*temp
                           end do loop_150
                       end do loop_160
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*a.
     
                   if (upper) then
                       loop_200: do j = n,1,-1
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           loop_170: do i = 1,m
                               b(i,j) = temp*b(i,j)
                           end do loop_170
                           loop_190: do k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   loop_180: do i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
                                   end do loop_180
                               end if
                           end do loop_190
                       end do loop_200
                   else
                       loop_240: do j = 1,n
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           loop_210: do i = 1,m
                               b(i,j) = temp*b(i,j)
                           end do loop_210
                           loop_230: do k = j + 1,n
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   loop_220: do i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
                                   end do loop_220
                               end if
                           end do loop_230
                       end do loop_240
                   end if
               else
     
                 ! form  b := alpha*b*a**t   or   b := alpha*b*a**h.
     
                   if (upper) then
                       loop_280: do k = 1,n
                           loop_260: do j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = alpha*a(j,k)
                                   else
                                       temp = alpha*dconjg(a(j,k))
                                   end if
                                   loop_250: do i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
                                   end do loop_250
                               end if
                           end do loop_260
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k,k)
                               else
                                   temp = temp*dconjg(a(k,k))
                               end if
                           end if
                           if (temp/=one) then
                               loop_270: do i = 1,m
                                   b(i,k) = temp*b(i,k)
                               end do loop_270
                           end if
                       end do loop_280
                   else
                       loop_320: do k = n,1,-1
                           loop_300: do j = k + 1,n
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = alpha*a(j,k)
                                   else
                                       temp = alpha*dconjg(a(j,k))
                                   end if
                                   loop_290: do i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
                                   end do loop_290
                               end if
                           end do loop_300
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k,k)
                               else
                                   temp = temp*dconjg(a(k,k))
                               end if
                           end if
                           if (temp/=one) then
                               loop_310: do i = 1,m
                                   b(i,k) = temp*b(i,k)
                               end do loop_310
                           end if
                       end do loop_320
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrmm
     
     end subroutine stdlib_ztrmm
     
     
     ! ZTRMV  performs one of the matrix-vector operations
     ! x := A*x,   or   x := A**T*x,   or   x := A**H*x,
     ! where x is an n element vector and  A is an n by n unit, or non-unit,
     ! upper or lower triangular matrix.
     subroutine stdlib_ztrmv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 2
           else if (.not.stdlib_lsame(diag,'u') .and. .not.stdlib_lsame(diag,'n')) then
               info = 3
           else if (n<0) then
               info = 4
           else if (lda<max(1,n)) then
               info = 6
           else if (incx==0) then
               info = 8
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztrmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
           noconj = stdlib_lsame(trans,'t')
           nounit = stdlib_lsame(diag,'n')
     
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  x := a*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       loop_20: do j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               loop_10: do i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
                               end do loop_10
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do loop_20
                   else
                       jx = kx
                       loop_40: do j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               loop_30: do i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
                               end do loop_30
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
                       end do loop_40
                   end if
               else
                   if (incx==1) then
                       loop_60: do j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               loop_50: do i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
                               end do loop_50
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
                       end do loop_60
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_80: do j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               loop_70: do i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
                               end do loop_70
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
                       end do loop_80
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       loop_110: do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               loop_90: do i = j - 1,1,-1
                                   temp = temp + a(i,j)*x(i)
                               end do loop_90
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               loop_100: do i = j - 1,1,-1
                                   temp = temp + dconjg(a(i,j))*x(i)
                               end do loop_100
                           end if
                           x(j) = temp
                       end do loop_110
                   else
                       jx = kx + (n-1)*incx
                       loop_140: do j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               loop_120: do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + a(i,j)*x(ix)
                               end do loop_120
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               loop_130: do i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + dconjg(a(i,j))*x(ix)
                               end do loop_130
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do loop_140
                   end if
               else
                   if (incx==1) then
                       loop_170: do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               loop_150: do i = j + 1,n
                                   temp = temp + a(i,j)*x(i)
                               end do loop_150
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               loop_160: do i = j + 1,n
                                   temp = temp + dconjg(a(i,j))*x(i)
                               end do loop_160
                           end if
                           x(j) = temp
                       end do loop_170
                   else
                       jx = kx
                       loop_200: do j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               loop_180: do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + a(i,j)*x(ix)
                               end do loop_180
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               loop_190: do i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + dconjg(a(i,j))*x(ix)
                               end do loop_190
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do loop_200
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrmv
     
     end subroutine stdlib_ztrmv
     
     
     ! ZTRSM  solves one of the matrix equations
     ! op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
     ! where alpha is a scalar, X and B are m by n matrices, A is a unit, or
     ! non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
     ! op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
     ! The matrix X is overwritten on B.
     subroutine stdlib_ztrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,noconj,nounit,upper
           ! ..
           ! .. parameters ..
           complex(dp) one
           parameter (one= (1.0_dp,0.0_dp))
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
     
           ! test the input parameters.
     
           lside = stdlib_lsame(side,'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
           noconj = stdlib_lsame(transa,'t')
           nounit = stdlib_lsame(diag,'n')
           upper = stdlib_lsame(uplo,'u')
     
           info = 0
           if ((.not.lside) .and. (.not.stdlib_lsame(side,'r'))) then
               info = 1
           else if ((.not.upper) .and. (.not.stdlib_lsame(uplo,'l'))) then
               info = 2
           else if ((.not.stdlib_lsame(transa,'n')) .and.(.not.stdlib_lsame(transa,'t')) .and.(&
                     .not.stdlib_lsame(transa,'c'))) then
               info = 3
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n'))) &
                     then
               info = 4
           else if (m<0) then
               info = 5
           else if (n<0) then
               info = 6
           else if (lda<max(1,nrowa)) then
               info = 9
           else if (ldb<max(1,m)) then
               info = 11
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztrsm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               loop_20: do j = 1,n
                   loop_10: do i = 1,m
                       b(i,j) = zero
                   end do loop_10
               end do loop_20
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*inv( a )*b.
     
                   if (upper) then
                       loop_60: do j = 1,n
                           if (alpha/=one) then
                               loop_30: do i = 1,m
                                   b(i,j) = alpha*b(i,j)
                               end do loop_30
                           end if
                           loop_50: do k = m,1,-1
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   loop_40: do i = 1,k - 1
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
                                   end do loop_40
                               end if
                           end do loop_50
                       end do loop_60
                   else
                       loop_100: do j = 1,n
                           if (alpha/=one) then
                               loop_70: do i = 1,m
                                   b(i,j) = alpha*b(i,j)
                               end do loop_70
                           end if
                           loop_90: do k = 1,m
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   loop_80: do i = k + 1,m
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
                                   end do loop_80
                               end if
                           end do loop_90
                       end do loop_100
                   end if
               else
     
                 ! form  b := alpha*inv( a**t )*b
                 ! or    b := alpha*inv( a**h )*b.
     
                   if (upper) then
                       loop_140: do j = 1,n
                           loop_130: do i = 1,m
                               temp = alpha*b(i,j)
                               if (noconj) then
                                   loop_110: do k = 1,i - 1
                                       temp = temp - a(k,i)*b(k,j)
                                   end do loop_110
                                   if (nounit) temp = temp/a(i,i)
                               else
                                   loop_120: do k = 1,i - 1
                                       temp = temp - dconjg(a(k,i))*b(k,j)
                                   end do loop_120
                                   if (nounit) temp = temp/dconjg(a(i,i))
                               end if
                               b(i,j) = temp
                           end do loop_130
                       end do loop_140
                   else
                       loop_180: do j = 1,n
                           loop_170: do i = m,1,-1
                               temp = alpha*b(i,j)
                               if (noconj) then
                                   loop_150: do k = i + 1,m
                                       temp = temp - a(k,i)*b(k,j)
                                   end do loop_150
                                   if (nounit) temp = temp/a(i,i)
                               else
                                   loop_160: do k = i + 1,m
                                       temp = temp - dconjg(a(k,i))*b(k,j)
                                   end do loop_160
                                   if (nounit) temp = temp/dconjg(a(i,i))
                               end if
                               b(i,j) = temp
                           end do loop_170
                       end do loop_180
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*inv( a ).
     
                   if (upper) then
                       loop_230: do j = 1,n
                           if (alpha/=one) then
                               loop_190: do i = 1,m
                                   b(i,j) = alpha*b(i,j)
                               end do loop_190
                           end if
                           loop_210: do k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   loop_200: do i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
                                   end do loop_200
                               end if
                           end do loop_210
                           if (nounit) then
                               temp = one/a(j,j)
                               loop_220: do i = 1,m
                                   b(i,j) = temp*b(i,j)
                               end do loop_220
                           end if
                       end do loop_230
                   else
                       loop_280: do j = n,1,-1
                           if (alpha/=one) then
                               loop_240: do i = 1,m
                                   b(i,j) = alpha*b(i,j)
                               end do loop_240
                           end if
                           loop_260: do k = j + 1,n
                               if (a(k,j)/=zero) then
                                   loop_250: do i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
                                   end do loop_250
                               end if
                           end do loop_260
                           if (nounit) then
                               temp = one/a(j,j)
                               loop_270: do i = 1,m
                                   b(i,j) = temp*b(i,j)
                               end do loop_270
                           end if
                       end do loop_280
                   end if
               else
     
                 ! form  b := alpha*b*inv( a**t )
                 ! or    b := alpha*b*inv( a**h ).
     
                   if (upper) then
                       loop_330: do k = n,1,-1
                           if (nounit) then
                               if (noconj) then
                                   temp = one/a(k,k)
                               else
                                   temp = one/dconjg(a(k,k))
                               end if
                               loop_290: do i = 1,m
                                   b(i,k) = temp*b(i,k)
                               end do loop_290
                           end if
                           loop_310: do j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = a(j,k)
                                   else
                                       temp = dconjg(a(j,k))
                                   end if
                                   loop_300: do i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
                                   end do loop_300
                               end if
                           end do loop_310
                           if (alpha/=one) then
                               loop_320: do i = 1,m
                                   b(i,k) = alpha*b(i,k)
                               end do loop_320
                           end if
                       end do loop_330
                   else
                       loop_380: do k = 1,n
                           if (nounit) then
                               if (noconj) then
                                   temp = one/a(k,k)
                               else
                                   temp = one/dconjg(a(k,k))
                               end if
                               loop_340: do i = 1,m
                                   b(i,k) = temp*b(i,k)
                               end do loop_340
                           end if
                           loop_360: do j = k + 1,n
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = a(j,k)
                                   else
                                       temp = dconjg(a(j,k))
                                   end if
                                   loop_350: do i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
                                   end do loop_350
                               end if
                           end do loop_360
                           if (alpha/=one) then
                               loop_370: do i = 1,m
                                   b(i,k) = alpha*b(i,k)
                               end do loop_370
                           end if
                       end do loop_380
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrsm
     
     end subroutine stdlib_ztrsm
     
     
     ! ZTRSV  solves one of the systems of equations
     ! A*x = b,   or   A**T*x = b,   or   A**H*x = b,
     ! where b and x are n element vectors and A is an n by n unit, or
     ! non-unit, upper or lower triangular matrix.
     ! No test for singularity or near-singularity is included in this
     ! routine. Such tests must be performed before calling this routine.
     subroutine stdlib_ztrsv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(dp) zero
           parameter (zero= (0.0_dp,0.0_dp))
           ! ..
           ! .. local scalars ..
           complex(dp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic dconjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t') &
                     .and..not.stdlib_lsame(trans,'c')) then
               info = 2
           else if (.not.stdlib_lsame(diag,'u') .and. .not.stdlib_lsame(diag,'n')) then
               info = 3
           else if (n<0) then
               info = 4
           else if (lda<max(1,n)) then
               info = 6
           else if (incx==0) then
               info = 8
           end if
           if (info/=0) then
               call stdlib_xerbla('stdlib_ztrsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
           noconj = stdlib_lsame(trans,'t')
           nounit = stdlib_lsame(diag,'n')
     
           ! set up the start point in x if the increment is not unity. this
           ! will be  ( n - 1 )*incx  too small for descending loops.
     
           if (incx<=0) then
               kx = 1 - (n-1)*incx
           else if (incx/=1) then
               kx = 1
           end if
     
           ! start the operations. in this version the elements of a are
           ! accessed sequentially with one pass through a.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  x := inv( a )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       loop_20: do j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               loop_10: do i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
                               end do loop_10
                           end if
                       end do loop_20
                   else
                       jx = kx + (n-1)*incx
                       loop_40: do j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               loop_30: do i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do loop_30
                           end if
                           jx = jx - incx
                       end do loop_40
                   end if
               else
                   if (incx==1) then
                       loop_60: do j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               loop_50: do i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
                               end do loop_50
                           end if
                       end do loop_60
                   else
                       jx = kx
                       loop_80: do j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               loop_70: do i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
                               end do loop_70
                           end if
                           jx = jx + incx
                       end do loop_80
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       loop_110: do j = 1,n
                           temp = x(j)
                           if (noconj) then
                               loop_90: do i = 1,j - 1
                                   temp = temp - a(i,j)*x(i)
                               end do loop_90
                               if (nounit) temp = temp/a(j,j)
                           else
                               loop_100: do i = 1,j - 1
                                   temp = temp - dconjg(a(i,j))*x(i)
                               end do loop_100
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(j) = temp
                       end do loop_110
                   else
                       jx = kx
                       loop_140: do j = 1,n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               loop_120: do i = 1,j - 1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix + incx
                               end do loop_120
                               if (nounit) temp = temp/a(j,j)
                           else
                               loop_130: do i = 1,j - 1
                                   temp = temp - dconjg(a(i,j))*x(ix)
                                   ix = ix + incx
                               end do loop_130
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                       end do loop_140
                   end if
               else
                   if (incx==1) then
                       loop_170: do j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               loop_150: do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(i)
                               end do loop_150
                               if (nounit) temp = temp/a(j,j)
                           else
                               loop_160: do i = n,j + 1,-1
                                   temp = temp - dconjg(a(i,j))*x(i)
                               end do loop_160
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(j) = temp
                       end do loop_170
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       loop_200: do j = n,1,-1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               loop_180: do i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix - incx
                               end do loop_180
                               if (nounit) temp = temp/a(j,j)
                           else
                               loop_190: do i = n,j + 1,-1
                                   temp = temp - dconjg(a(i,j))*x(ix)
                                   ix = ix - incx
                               end do loop_190
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                       end do loop_200
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrsv
     
     end subroutine stdlib_ztrsv



end module stdlib_linalg_blas_z
