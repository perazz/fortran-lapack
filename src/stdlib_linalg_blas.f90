module stdlib_linalg_blas
     use stdlib_linalg_constants
     implicit none(type,external)
     private






     public :: sp,dp,lk,int32,int64
     public :: stdlib_caxpy
     public :: stdlib_ccopy
     public :: stdlib_cdotc
     public :: stdlib_cdotu
     public :: stdlib_cgbmv
     public :: stdlib_cgemm
     public :: stdlib_cgemv
     public :: stdlib_cgerc
     public :: stdlib_cgeru
     public :: stdlib_chbmv
     public :: stdlib_chemm
     public :: stdlib_chemv
     public :: stdlib_cher
     public :: stdlib_cher2
     public :: stdlib_cher2k
     public :: stdlib_cherk
     public :: stdlib_chpmv
     public :: stdlib_chpr
     public :: stdlib_chpr2
     public :: stdlib_crotg
     public :: stdlib_cscal
     public :: stdlib_csrot
     public :: stdlib_csscal
     public :: stdlib_cswap
     public :: stdlib_csymm
     public :: stdlib_csyr2k
     public :: stdlib_csyrk
     public :: stdlib_ctbmv
     public :: stdlib_ctbsv
     public :: stdlib_ctpmv
     public :: stdlib_ctpsv
     public :: stdlib_ctrmm
     public :: stdlib_ctrmv
     public :: stdlib_ctrsm
     public :: stdlib_ctrsv
     public :: stdlib_dasum
     public :: stdlib_daxpy
     public :: stdlib_dcabs1
     public :: stdlib_dcopy
     public :: stdlib_ddot
     public :: stdlib_dgbmv
     public :: stdlib_dgemm
     public :: stdlib_dgemv
     public :: stdlib_dger
     public :: stdlib_dnrm2
     public :: stdlib_drot
     public :: stdlib_drotg
     public :: stdlib_drotm
     public :: stdlib_drotmg
     public :: stdlib_dsbmv
     public :: stdlib_dscal
     public :: stdlib_dsdot
     public :: stdlib_dspmv
     public :: stdlib_dspr
     public :: stdlib_dspr2
     public :: stdlib_dswap
     public :: stdlib_dsymm
     public :: stdlib_dsymv
     public :: stdlib_dsyr
     public :: stdlib_dsyr2
     public :: stdlib_dsyr2k
     public :: stdlib_dsyrk
     public :: stdlib_dtbmv
     public :: stdlib_dtbsv
     public :: stdlib_dtpmv
     public :: stdlib_dtpsv
     public :: stdlib_dtrmm
     public :: stdlib_dtrmv
     public :: stdlib_dtrsm
     public :: stdlib_dtrsv
     public :: stdlib_dzasum
     public :: stdlib_dznrm2
     public :: stdlib_icamax
     public :: stdlib_idamax
     public :: stdlib_isamax
     public :: stdlib_izamax
     public :: stdlib_lsame
     public :: stdlib_sasum
     public :: stdlib_saxpy
     public :: stdlib_scabs1
     public :: stdlib_scasum
     public :: stdlib_scnrm2
     public :: stdlib_scopy
     public :: stdlib_sdot
     public :: stdlib_sdsdot
     public :: stdlib_sgbmv
     public :: stdlib_sgemm
     public :: stdlib_sgemv
     public :: stdlib_sger
     public :: stdlib_snrm2
     public :: stdlib_srot
     public :: stdlib_srotg
     public :: stdlib_srotm
     public :: stdlib_srotmg
     public :: stdlib_ssbmv
     public :: stdlib_sscal
     public :: stdlib_sspmv
     public :: stdlib_sspr
     public :: stdlib_sspr2
     public :: stdlib_sswap
     public :: stdlib_ssymm
     public :: stdlib_ssymv
     public :: stdlib_ssyr
     public :: stdlib_ssyr2
     public :: stdlib_ssyr2k
     public :: stdlib_ssyrk
     public :: stdlib_stbmv
     public :: stdlib_stbsv
     public :: stdlib_stpmv
     public :: stdlib_stpsv
     public :: stdlib_strmm
     public :: stdlib_strmv
     public :: stdlib_strsm
     public :: stdlib_strsv
     public :: stdlib_xerbla
     public :: stdlib_xerbla_array
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
     
     
     subroutine stdlib_ccopy(n,cx,incx,cy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*),cy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 cy(i) = cx(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 cy(iy) = cx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_ccopy
     
     end subroutine stdlib_ccopy

     
     
     complex(sp) function stdlib_cdotc(n,cx,incx,cy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*),cy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           complex(sp) ctemp
           integer(int32) i,ix,iy
           ! ..
           ! .. intrinsic functions ..
           intrinsic conjg
           ! ..
           ctemp = (0.0,0.0)
           stdlib_cdotc = (0.0,0.0)
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 ctemp = ctemp + conjg(cx(i))*cy(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 ctemp = ctemp + conjg(cx(ix))*cy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_cdotc = ctemp
           return
     
           ! end of stdlib_cdotc
     
     end function stdlib_cdotc

     
     
     complex(sp) function stdlib_cdotu(n,cx,incx,cy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*),cy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           complex(sp) ctemp
           integer(int32) i,ix,iy
           ! ..
           ctemp = (0.0,0.0)
           stdlib_cdotu = (0.0,0.0)
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 ctemp = ctemp + cx(i)*cy(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 ctemp = ctemp + cx(ix)*cy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_cdotu = ctemp
           return
     
           ! end of stdlib_cdotu
     
     end function stdlib_cdotu

     
     
     subroutine stdlib_crotg( a, b, c, s )
        integer, parameter :: wp = kind(1.e0)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
        ! .. constants ..
        real(sp), parameter ::szero = 0.0_sp
        real(sp), parameter :: sone  = 1.0_sp
        complex(sp), parameter :: czero  = 0.0_sp
        ! ..
        ! .. scaling constants ..
     real(sp), parameter :: safmin = real(radix(0._sp),wp)**max(minexponent(0._sp)-1,1-&
          maxexponent(0._sp)   )
     real(sp), parameter :: safmax = real(radix(0._sp),wp)**max(1-minexponent(0._sp),maxexponent(&
          0._sp)-1   )
     real(sp), parameter :: rtmin = sqrt( real(radix(0._sp),wp)**max(minexponent(0._sp)-1,1-&
          maxexponent(0._sp)   ) / epsilon(0._sp) )
     real(sp), parameter :: rtmax = sqrt( real(radix(0._sp),wp)**max(1-minexponent(0._sp),&
          maxexponent(0._sp)-1   ) * epsilon(0._sp) )
        ! ..
        ! .. scalar arguments ..
        real(sp) :: c
        complex(sp) :: a, b, s
        ! ..
        ! .. local scalars ..
        real(sp) :: d, f1, f2, g1, g2, h2, p, u, uu, v, vv, w
        complex(sp) :: f, fs, g, gs, r, t
        ! ..
        ! .. intrinsic functions ..
        intrinsic :: abs, aimag, conjg, max, min, real, sqrt
        ! ..
        ! .. statement functions ..
        real(sp) :: abssq
        ! ..
        ! .. statement function definitions ..
        abssq( t ) = real( t )**2 + aimag( t )**2
        ! ..
        ! .. executable statements ..
     
        f = a
        g = b
        if( g == czero ) then
           c = sone
           s = czero
           r = f
        else if( f == czero ) then
           c =szero
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
              uu = sone / u
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
              uu = sone / u
              gs = g*uu
              g2 = abssq( gs )
              if( f1*uu < rtmin ) then
     
                 ! f is not well-scaled when scaled by g1.
                 ! use a different scaling for f.
     
                 v = min( safmax, max( safmin, f1 ) )
                 vv = sone / v
                 w = v * uu
                 fs = f*vv
                 f2 = abssq( fs )
                 h2 = f2*w**2 + g2
              else
     
                 ! otherwise use the same scaling for f and g.
     
                 w = sone
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
     end subroutine stdlib_crotg

     
     
     subroutine stdlib_cscal(n,ca,cx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) ca
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,nincx
           ! ..
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              do i = 1,n
                 cx(i) = ca*cx(i)
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 cx(i) = ca*cx(i)
              end do
           end if
           return
     
           ! end of stdlib_cscal
     
     end subroutine stdlib_cscal

     
     
     subroutine stdlib_csrot( n, cx, incx, cy, incy, c, s )
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32)           incx, incy, n
           real(sp)              c, s
           ! ..
           ! .. array arguments ..
           complex(sp)           cx( * ), cy( * )
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32)           i, ix, iy
           complex(sp)           ctemp
           ! ..
           ! .. executable statements ..
     
           if( n<=0 )return
           if( incx==1 .and. incy==1 ) then
     
              ! code for both increments equal to 1
     
              do i = 1, n
                 ctemp = c*cx( i ) + s*cy( i )
                 cy( i ) = c*cy( i ) - s*cx( i )
                 cx( i ) = ctemp
              end do
           else
     
              ! code for unequal increments or equal increments not equal
                ! to 1
     
              ix = 1
              iy = 1
              if( incx<0 )ix = ( -n+1 )*incx + 1
              if( incy<0 )iy = ( -n+1 )*incy + 1
              do i = 1, n
                 ctemp = c*cx( ix ) + s*cy( iy )
                 cy( iy ) = c*cy( iy ) - s*cx( ix )
                 cx( ix ) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_csrot
     
     end subroutine stdlib_csrot

     
     
     subroutine stdlib_csscal(n,sa,cx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) sa
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic aimag,cmplx,real
           ! ..
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              do i = 1,n
                 cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
              end do
           end if
           return
     
           ! end of stdlib_csscal
     
     end subroutine stdlib_csscal

     
     
     subroutine stdlib_cswap(n,cx,incx,cy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*),cy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           complex(sp) ctemp
           integer(int32) i,ix,iy
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
             ! code for both increments equal to 1
              do i = 1,n
                 ctemp = cx(i)
                 cx(i) = cy(i)
                 cy(i) = ctemp
              end do
           else
     
             ! code for unequal increments or equal increments not equal
               ! to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 ctemp = cx(ix)
                 cx(ix) = cy(iy)
                 cy(iy) = ctemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_cswap
     
     end subroutine stdlib_cswap

     
     
     real(dp) function stdlib_dasum(n,dx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dtemp
           integer(int32) i,m,mp1,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic dabs,mod
           ! ..
           stdlib_dasum = 0.0d0
           dtemp = 0.0d0
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
              ! code for increment equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,6)
              if (m/=0) then
                 do i = 1,m
                    dtemp = dtemp + dabs(dx(i))
                 end do
                 if (n<6) then
                    stdlib_dasum = dtemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,6
                 dtemp = dtemp + dabs(dx(i)) + dabs(dx(i+1)) +dabs(dx(i+2)) + dabs(dx(i+3)) +dabs(&
          dx(i+4)) + dabs(dx(i+5))
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 dtemp = dtemp + dabs(dx(i))
              end do
           end if
           stdlib_dasum = dtemp
           return
     
           ! end of stdlib_dasum
     
     end function stdlib_dasum

     
     
     subroutine stdlib_daxpy(n,da,dx,incx,dy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) da
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*),dy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0) return
           if (da==0.0d0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,4)
              if (m/=0) then
                 do i = 1,m
                    dy(i) = dy(i) + da*dx(i)
                 end do
              end if
              if (n<4) return
              mp1 = m + 1
              do i = mp1,n,4
                 dy(i) = dy(i) + da*dx(i)
                 dy(i+1) = dy(i+1) + da*dx(i+1)
                 dy(i+2) = dy(i+2) + da*dx(i+2)
                 dy(i+3) = dy(i+3) + da*dx(i+3)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
               dy(iy) = dy(iy) + da*dx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_daxpy
     
     end subroutine stdlib_daxpy

     
     
     real(dp) function stdlib_dcabs1(z)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(dp) z
           ! ..
           ! ..
        ! =====================================================================
     
           ! .. intrinsic functions ..
           intrinsic abs,dble,dimag
     
           stdlib_dcabs1 = abs(dble(z)) + abs(dimag(z))
           return
     
           ! end of stdlib_dcabs1
     
     end function stdlib_dcabs1

     
     
     subroutine stdlib_dcopy(n,dx,incx,dy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*),dy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,7)
              if (m/=0) then
                 do i = 1,m
                    dy(i) = dx(i)
                 end do
                 if (n<7) return
              end if
              mp1 = m + 1
              do i = mp1,n,7
                 dy(i) = dx(i)
                 dy(i+1) = dx(i+1)
                 dy(i+2) = dx(i+2)
                 dy(i+3) = dx(i+3)
                 dy(i+4) = dx(i+4)
                 dy(i+5) = dx(i+5)
                 dy(i+6) = dx(i+6)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 dy(iy) = dx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_dcopy
     
     end subroutine stdlib_dcopy

     
     
     real(dp) function stdlib_ddot(n,dx,incx,dy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*),dy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dtemp
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           stdlib_ddot = 0.0d0
           dtemp = 0.0d0
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,5)
              if (m/=0) then
                 do i = 1,m
                    dtemp = dtemp + dx(i)*dy(i)
                 end do
                 if (n<5) then
                    stdlib_ddot=dtemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,5
               dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) +dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) +&
           dx(i+4)*dy(i+4)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 dtemp = dtemp + dx(ix)*dy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_ddot = dtemp
           return
     
           ! end of stdlib_ddot
     
     end function stdlib_ddot

     
     
     function stdlib_dnrm2( n, x, incx )
        integer, parameter :: wp = kind(1.d0)
        real(dp) :: stdlib_dnrm2
     
        ! -- reference blas level1 routine (version 3.9.1) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
     
        ! .. constants ..
        real(dp), parameter ::dzero = 0.0_dp
        real(dp), parameter :: done  = 1.0_dp
        real(dp), parameter :: maxn = huge(0.0_dp)
        ! ..
        ! .. blue's scaling constants ..
     real(dp), parameter :: tsml = real(radix(0._dp), wp)**ceiling(       (minexponent(0._dp) - 1)&
           * 0.5_dp)
     real(dp), parameter :: tbig = real(radix(0._dp), wp)**floor(       (maxexponent(0._dp) -&
           digits(0._dp) + 1) * 0.5_dp)
     real(dp), parameter :: ssml = real(radix(0._dp), wp)**( - floor(       (minexponent(0._dp) -&
           digits(0._dp)) * 0.5_dp))
     real(dp), parameter :: sbig = real(radix(0._dp), wp)**( - ceiling(       (maxexponent(0._dp)&
           + digits(0._dp) - 1) * 0.5_dp))
        ! ..
        ! .. scalar arguments ..
        integer :: incx, n
        ! ..
        ! .. array arguments ..
        real(dp) :: x(*)
        ! ..
        ! .. local scalars ..
        integer :: i, ix
        logical :: notbig
        real(dp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
     
        ! quick return if possible
     
        stdlib_dnrm2 =dzero
        if( n <= 0 ) return
     
        scl = done
        sumsq =dzero
     
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
     
        notbig = .true.
        asml =dzero
        amed =dzero
        abig =dzero
        ix = 1
        if( incx < 0 ) ix = 1 - (n-1)*incx
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
     
        ! combine abig and amed or amed and asml if more than done
        ! accumulator was used.
     
        if (abig >dzero) then
     
           ! combine abig and amed if abig > 0.
     
           if ( (amed >dzero) .or. (amed > maxn) .or. (amed /= amed) ) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = done / sbig
           sumsq = abig
        else if (asml >dzero) then
     
           ! combine amed and asml if asml > 0.
     
           if ( (amed >dzero) .or. (amed > maxn) .or. (amed /= amed) ) then
              amed = sqrt(amed)
              asml = sqrt(asml) / ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = done
              sumsq = ymax**2*( done + (ymin/ymax)**2 )
           else
              scl = done / ssml
              sumsq = asml
           end if
        else
     
           ! otherwise all values are mid-range
     
           scl = done
           sumsq = amed
        end if
        stdlib_dnrm2 = scl*sqrt( sumsq )
        return
     end function stdlib_dnrm2

     
     
     subroutine stdlib_drot(n,dx,incx,dy,incy,c,s)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) c,s
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*),dy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dtemp
           integer(int32) i,ix,iy
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
             ! code for both increments equal to 1
     
              do i = 1,n
                 dtemp = c*dx(i) + s*dy(i)
                 dy(i) = c*dy(i) - s*dx(i)
                 dx(i) = dtemp
              end do
           else
     
             ! code for unequal increments or equal increments not equal
               ! to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 dtemp = c*dx(ix) + s*dy(iy)
                 dy(iy) = c*dy(iy) - s*dx(ix)
                 dx(ix) = dtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_drot
     
     end subroutine stdlib_drot

     
     
     subroutine stdlib_drotg( a, b, c, s )
        integer, parameter :: wp = kind(1.d0)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
        ! .. constants ..
        real(dp), parameter ::dzero = 0.0_dp
        real(dp), parameter :: done  = 1.0_dp
        ! ..
        ! .. scaling constants ..
     real(dp), parameter :: safmin = real(radix(0._dp),wp)**max(minexponent(0._dp)-1,1-&
          maxexponent(0._dp)   )
     real(dp), parameter :: safmax = real(radix(0._dp),wp)**max(1-minexponent(0._dp),maxexponent(&
          0._dp)-1   )
        ! ..
        ! .. scalar arguments ..
        real(dp) :: a, b, c, s
        ! ..
        ! .. local scalars ..
        real(dp) :: anorm, bnorm, scl, sigma, r, z
        ! ..
        anorm = abs(a)
        bnorm = abs(b)
        if( bnorm ==dzero ) then
           c = done
           s =dzero
           b =dzero
        else if( anorm ==dzero ) then
           c =dzero
           s = done
           a = b
           b = done
        else
           scl = min( safmax, max( safmin, anorm, bnorm ) )
           if( anorm > bnorm ) then
              sigma = sign(done,a)
           else
              sigma = sign(done,b)
           end if
           r = sigma*( scl*sqrt((a/scl)**2 + (b/scl)**2) )
           c = a/r
           s = b/r
           if( anorm > bnorm ) then
              z = s
           else if( c /=dzero ) then
              z = done/c
           else
              z = done
           end if
           a = r
           b = z
        end if
        return
     end subroutine stdlib_drotg

     
     
     subroutine stdlib_drotm(n,dx,incx,dy,incy,dparam)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(dp) dparam(5),dx(*),dy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dflag,dh11,dh12,dh21,dh22,two,w,z,zero
           integer(int32) i,kx,ky,nsteps
           ! ..
           ! .. data statements ..
           data zero,two/0.d0,2.d0/
           ! ..
     
           dflag = dparam(1)
           if (n<=0 .or. (dflag+two==zero)) return
           if (incx==incy.and.incx>0) then
     
              nsteps = n*incx
              if (dflag<zero) then
                 dh11 = dparam(2)
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 dh22 = dparam(5)
                 do i = 1,nsteps,incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w*dh11 + z*dh12
                    dy(i) = w*dh21 + z*dh22
                 end do
              else if (dflag==zero) then
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 do i = 1,nsteps,incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w + z*dh12
                    dy(i) = w*dh21 + z
                 end do
              else
                 dh11 = dparam(2)
                 dh22 = dparam(5)
                 do i = 1,nsteps,incx
                    w = dx(i)
                    z = dy(i)
                    dx(i) = w*dh11 + z
                    dy(i) = -w + dh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx<0) kx = 1 + (1-n)*incx
              if (incy<0) ky = 1 + (1-n)*incy
     
              if (dflag<zero) then
                 dh11 = dparam(2)
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 dh22 = dparam(5)
                 do i = 1,n
                    w = dx(kx)
                    z = dy(ky)
                    dx(kx) = w*dh11 + z*dh12
                    dy(ky) = w*dh21 + z*dh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (dflag==zero) then
                 dh12 = dparam(4)
                 dh21 = dparam(3)
                 do i = 1,n
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
                  do i = 1,n
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
     
           ! end of stdlib_drotm
     
     end subroutine stdlib_drotm

     
     
     subroutine stdlib_drotmg(dd1,dd2,dx1,dy1,dparam)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) dd1,dd2,dx1,dy1
           ! ..
           ! .. array arguments ..
           real(dp) dparam(5)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dflag,dh11,dh12,dh21,dh22,dp1,dp2,dq1,dq2,dtemp,du,gam,gamsq,one,rgamsq,two,&
          zero
           ! ..
           ! .. intrinsic functions ..
           intrinsic dabs
           ! ..
           ! .. data statements ..
     
           data zero,one,two/0.d0,1.d0,2.d0/
           data gam,gamsq,rgamsq/4096.d0,16777216.d0,5.9604645d-8/
           ! ..
           if (dd1<zero) then
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
              if (dp2==zero) then
                 dflag = -two
                 dparam(1) = dflag
                 return
              end if
              ! regular-case..
              dp1 = dd1*dx1
              dq2 = dp2*dy1
              dq1 = dp1*dx1
     
              if (dabs(dq1)>dabs(dq2)) then
                 dh21 = -dy1/dx1
                 dh12 = dp2/dp1
     
                 du = one - dh12*dh21
     
                if (du>zero) then
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
                 if (dq2<zero) then
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
              if (dd1/=zero) then
                 do while ((dd1<=rgamsq) .or. (dd1>=gamsq))
                    if (dflag==zero) then
                       dh11 = one
                       dh22 = one
                       dflag = -one
                    else
                       dh21 = -one
                       dh12 = one
                       dflag = -one
                    end if
                    if (dd1<=rgamsq) then
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
                 enddo
              end if
              if (dd2/=zero) then
                 do while ( (dabs(dd2)<=rgamsq) .or. (dabs(dd2)>=gamsq) )
                    if (dflag==zero) then
                       dh11 = one
                       dh22 = one
                       dflag = -one
                    else
                       dh21 = -one
                       dh12 = one
                       dflag = -one
                    end if
                    if (dabs(dd2)<=rgamsq) then
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
           if (dflag<zero) then
              dparam(2) = dh11
              dparam(3) = dh21
              dparam(4) = dh12
              dparam(5) = dh22
           else if (dflag==zero) then
              dparam(3) = dh21
              dparam(4) = dh12
           else
              dparam(2) = dh11
              dparam(5) = dh22
           end if
           dparam(1) = dflag
           return
     
           ! end of stdlib_drotmg
     
     end subroutine stdlib_drotmg

     
     
     subroutine stdlib_dscal(n,da,dx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) da
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,m,mp1,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,5)
              if (m/=0) then
                 do i = 1,m
                    dx(i) = da*dx(i)
                 end do
                 if (n<5) return
              end if
              mp1 = m + 1
              do i = mp1,n,5
                 dx(i) = da*dx(i)
                 dx(i+1) = da*dx(i+1)
                 dx(i+2) = da*dx(i+2)
                 dx(i+3) = da*dx(i+3)
                 dx(i+4) = da*dx(i+4)
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 dx(i) = da*dx(i)
              end do
           end if
           return
     
           ! end of stdlib_dscal
     
     end subroutine stdlib_dscal

     
     
     real(dp) function stdlib_dsdot(n,sx,incx,sy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! ..
     
        ! authors:
        ! ========
        ! lawson, c. l., (jpl), hanson, r. j., (snla),
        ! kincaid, d. r., (u. of texas), krogh, f. t., (jpl)
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,kx,ky,ns
           ! ..
           ! .. intrinsic functions ..
           intrinsic dble
           ! ..
           stdlib_dsdot = 0.0d0
           if (n<=0) return
           if (incx==incy .and. incx>0) then
     
           ! code for equal, positive, non-unit increments.
     
              ns = n*incx
              do i = 1,ns,incx
                 stdlib_dsdot = stdlib_dsdot + dble(sx(i))*dble(sy(i))
              end do
           else
     
           ! code for unequal or nonpositive increments.
     
              kx = 1
              ky = 1
              if (incx<0) kx = 1 + (1-n)*incx
              if (incy<0) ky = 1 + (1-n)*incy
              do i = 1,n
                 stdlib_dsdot = stdlib_dsdot + dble(sx(kx))*dble(sy(ky))
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           return
     
           ! end of stdlib_dsdot
     
     end function stdlib_dsdot

     
     
     subroutine stdlib_dswap(n,dx,incx,dy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*),dy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dtemp
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
             ! code for both increments equal to 1
     
     
             ! clean-up loop
     
              m = mod(n,3)
              if (m/=0) then
                 do i = 1,m
                    dtemp = dx(i)
                    dx(i) = dy(i)
                    dy(i) = dtemp
                 end do
                 if (n<3) return
              end if
              mp1 = m + 1
              do i = mp1,n,3
                 dtemp = dx(i)
                 dx(i) = dy(i)
                 dy(i) = dtemp
                 dtemp = dx(i+1)
                 dx(i+1) = dy(i+1)
                 dy(i+1) = dtemp
                 dtemp = dx(i+2)
                 dx(i+2) = dy(i+2)
                 dy(i+2) = dtemp
              end do
           else
     
             ! code for unequal increments or equal increments not equal
               ! to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 dtemp = dx(ix)
                 dx(ix) = dy(iy)
                 dy(iy) = dtemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_dswap
     
     end subroutine stdlib_dswap

     
     
     real(dp) function stdlib_dzasum(n,zx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) stemp
           integer(int32) i,nincx
           ! ..
     
     
           stdlib_dzasum = 0.0d0
           stemp = 0.0d0
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              do i = 1,n
                 stemp = stemp + stdlib_dcabs1(zx(i))
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + stdlib_dcabs1(zx(i))
              end do
           end if
           stdlib_dzasum = stemp
           return
     
           ! end of stdlib_dzasum
     
     end function stdlib_dzasum

     
     
     function stdlib_dznrm2( n, x, incx )
        integer, parameter :: wp = kind(1.d0)
        real(dp) :: stdlib_dznrm2
     
        ! -- reference blas level1 routine (version 3.9.1) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
     
        ! .. constants ..
        real(dp), parameter ::dzero = 0.0_dp
        real(dp), parameter :: done  = 1.0_dp
        real(dp), parameter :: maxn = huge(0.0_dp)
        ! ..
        ! .. blue's scaling constants ..
     real(dp), parameter :: tsml = real(radix(0._dp), wp)**ceiling(       (minexponent(0._dp) - 1)&
           * 0.5_dp)
     real(dp), parameter :: tbig = real(radix(0._dp), wp)**floor(       (maxexponent(0._dp) -&
           digits(0._dp) + 1) * 0.5_dp)
     real(dp), parameter :: ssml = real(radix(0._dp), wp)**( - floor(       (minexponent(0._dp) -&
           digits(0._dp)) * 0.5_dp))
     real(dp), parameter :: sbig = real(radix(0._dp), wp)**( - ceiling(       (maxexponent(0._dp)&
           + digits(0._dp) - 1) * 0.5_dp))
        ! ..
        ! .. scalar arguments ..
        integer :: incx, n
        ! ..
        ! .. array arguments ..
        complex(dp) :: x(*)
        ! ..
        ! .. local scalars ..
        integer :: i, ix
        logical :: notbig
        real(dp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
     
        ! quick return if possible
     
        stdlib_dznrm2 =dzero
        if( n <= 0 ) return
     
        scl = done
        sumsq =dzero
     
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
     
        notbig = .true.
        asml =dzero
        amed =dzero
        abig =dzero
        ix = 1
        if( incx < 0 ) ix = 1 - (n-1)*incx
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
     
        ! combine abig and amed or amed and asml if more than done
        ! accumulator was used.
     
        if (abig >dzero) then
     
           ! combine abig and amed if abig > 0.
     
           if ( (amed >dzero) .or. (amed > maxn) .or. (amed /= amed) ) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = done / sbig
           sumsq = abig
        else if (asml >dzero) then
     
           ! combine amed and asml if asml > 0.
     
           if ( (amed >dzero) .or. (amed > maxn) .or. (amed /= amed) ) then
              amed = sqrt(amed)
              asml = sqrt(asml) / ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = done
              sumsq = ymax**2*( done + (ymin/ymax)**2 )
           else
              scl = done / ssml
              sumsq = asml
           end if
        else
     
           ! otherwise all values are mid-range
     
           scl = done
           sumsq = amed
        end if
        stdlib_dznrm2 = scl*sqrt( sumsq )
        return
     end function stdlib_dznrm2

     
     
     integer(int32) function stdlib_idamax(n,dx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           real(dp) dx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dmax
           integer(int32) i,ix
           ! ..
           ! .. intrinsic functions ..
           intrinsic dabs
           ! ..
           stdlib_idamax = 0
           if (n<1 .or. incx<=0) return
           stdlib_idamax = 1
           if (n==1) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              dmax = dabs(dx(1))
              do i = 2,n
                 if (dabs(dx(i))>dmax) then
                    stdlib_idamax = i
                    dmax = dabs(dx(i))
                 end if
              end do
           else
     
              ! code for increment not equal to 1
     
              ix = 1
              dmax = dabs(dx(1))
              ix = ix + incx
              do i = 2,n
                 if (dabs(dx(ix))>dmax) then
                    stdlib_idamax = i
                    dmax = dabs(dx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     
           ! end of stdlib_idamax
     
     end function stdlib_idamax

     
     
     integer(int32) function stdlib_isamax(n,sx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) smax
           integer(int32) i,ix
           ! ..
           ! .. intrinsic functions ..
           intrinsic abs
           ! ..
           stdlib_isamax = 0
           if (n<1 .or. incx<=0) return
           stdlib_isamax = 1
           if (n==1) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              smax = abs(sx(1))
              do i = 2,n
                 if (abs(sx(i))>smax) then
                    stdlib_isamax = i
                    smax = abs(sx(i))
                 end if
              end do
           else
     
              ! code for increment not equal to 1
     
              ix = 1
              smax = abs(sx(1))
              ix = ix + incx
              do i = 2,n
                 if (abs(sx(ix))>smax) then
                    stdlib_isamax = i
                    smax = abs(sx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     
           ! end of stdlib_isamax
     
     end function stdlib_isamax

     
     
     integer(int32) function stdlib_izamax(n,zx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(dp) zx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(dp) dmax
           integer(int32) i,ix
           ! ..
     
     
           stdlib_izamax = 0
           if (n<1 .or. incx<=0) return
           stdlib_izamax = 1
           if (n==1) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              dmax = stdlib_dcabs1(zx(1))
              do i = 2,n
                 if (stdlib_dcabs1(zx(i))>dmax) then
                    stdlib_izamax = i
                    dmax = stdlib_dcabs1(zx(i))
                 end if
              end do
           else
     
              ! code for increment not equal to 1
     
              ix = 1
              dmax = stdlib_dcabs1(zx(1))
              ix = ix + incx
              do i = 2,n
                 if (stdlib_dcabs1(zx(ix))>dmax) then
                    stdlib_izamax = i
                    dmax = stdlib_dcabs1(zx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     
           ! end of stdlib_izamax
     
     end function stdlib_izamax

     
     
     logical(lk) function stdlib_lsame(ca,cb)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           character ca,cb
           ! ..
     
       ! =====================================================================
     
           ! .. intrinsic functions ..
           intrinsic ichar
           ! ..
           ! .. local scalars ..
           integer(int32) inta,intb,zcode
           ! ..
     
           ! test if the characters are equal
     
           stdlib_lsame = ca == cb
           if (stdlib_lsame) return
     
           ! now test for equivalence if both characters are alphabetic.
     
           zcode = ichar('z')
     
           ! use 'z' rather than 'a' so that ascii can be detected on prime
           ! machines, on which ichar returns a value with bit 8 set.
           ! ichar('a') on prime machines returns 193 which is the same as
           ! ichar('a') on an ebcdic machine.
     
           inta = ichar(ca)
           intb = ichar(cb)
     
           if (zcode==90 .or. zcode==122) then
     
              ! ascii is assumed - zcode is the ascii code of either lower or
              ! upper case 'z'.
     
               if (inta>=97 .and. inta<=122) inta = inta - 32
               if (intb>=97 .and. intb<=122) intb = intb - 32
     
           else if (zcode==233 .or. zcode==169) then
     
              ! ebcdic is assumed - zcode is the ebcdic code of either lower or
              ! upper case 'z'.
     
               if (&
          inta>=129 .and. inta<=137 .or.inta>=145 .and. inta<=153 .or.inta>=162 .and. inta<=169)&
           inta = inta + 64
               if (&
          intb>=129 .and. intb<=137 .or.intb>=145 .and. intb<=153 .or.intb>=162 .and. intb<=169)&
           intb = intb + 64
     
           else if (zcode==218 .or. zcode==250) then
     
              ! ascii is assumed, on prime machines - zcode is the ascii code
              ! plus 128 of either lower or upper case 'z'.
     
               if (inta>=225 .and. inta<=250) inta = inta - 32
               if (intb>=225 .and. intb<=250) intb = intb - 32
           end if
           stdlib_lsame = inta == intb
     
           ! return
     
           ! end of stdlib_lsame
     
     end function stdlib_lsame

     
     
     real(sp) function stdlib_sasum(n,sx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) stemp
           integer(int32) i,m,mp1,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic abs,mod
           ! ..
           stdlib_sasum = 0.0e0
           stemp = 0.0e0
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
              ! code for increment equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,6)
              if (m/=0) then
                 do i = 1,m
                    stemp = stemp + abs(sx(i))
                 end do
                 if (n<6) then
                    stdlib_sasum = stemp
                    return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,6
                 stemp = stemp + abs(sx(i)) + abs(sx(i+1)) +abs(sx(i+2)) + abs(sx(i+3)) +abs(sx(i+&
          4)) + abs(sx(i+5))
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + abs(sx(i))
              end do
           end if
           stdlib_sasum = stemp
           return
     
           ! end of stdlib_sasum
     
     end function stdlib_sasum

     
     
     subroutine stdlib_saxpy(n,sa,sx,incx,sy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) sa
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0) return
           if (sa==0.0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,4)
              if (m/=0) then
                 do i = 1,m
                    sy(i) = sy(i) + sa*sx(i)
                 end do
              end if
              if (n<4) return
              mp1 = m + 1
              do i = mp1,n,4
                 sy(i) = sy(i) + sa*sx(i)
                 sy(i+1) = sy(i+1) + sa*sx(i+1)
                 sy(i+2) = sy(i+2) + sa*sx(i+2)
                 sy(i+3) = sy(i+3) + sa*sx(i+3)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
               sy(iy) = sy(iy) + sa*sx(ix)
               ix = ix + incx
               iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_saxpy
     
     end subroutine stdlib_saxpy

     
     
     real(sp) function stdlib_scabs1(z)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) z
           ! ..
     
        ! =====================================================================
     
           ! .. intrinsic functions ..
           intrinsic abs,aimag,real
           ! ..
           stdlib_scabs1 = abs(real(z)) + abs(aimag(z))
           return
     
           ! end of stdlib_scabs1
     
     end function stdlib_scabs1

     
     
     real(sp) function stdlib_scasum(n,cx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) stemp
           integer(int32) i,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic abs,aimag,real
           ! ..
           stdlib_scasum = 0.0e0
           stemp = 0.0e0
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              do i = 1,n
                 stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
              end do
           end if
           stdlib_scasum = stemp
           return
     
           ! end of stdlib_scasum
     
     end function stdlib_scasum

     
     
     function stdlib_scnrm2( n, x, incx )
        integer, parameter :: wp = kind(1.e0)
        real(sp) :: stdlib_scnrm2
     
        ! -- reference blas level1 routine (version 3.9.1) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
     
        ! .. constants ..
        real(sp), parameter ::szero = 0.0_sp
        real(sp), parameter :: sone  = 1.0_sp
        real(sp), parameter :: maxn = huge(0.0_sp)
        ! ..
        ! .. blue's scaling constants ..
     real(sp), parameter :: tsml = real(radix(0._sp), wp)**ceiling(       (minexponent(0._sp) - 1)&
           * 0.5_sp)
     real(sp), parameter :: tbig = real(radix(0._sp), wp)**floor(       (maxexponent(0._sp) -&
           digits(0._sp) + 1) * 0.5_sp)
     real(sp), parameter :: ssml = real(radix(0._sp), wp)**( - floor(       (minexponent(0._sp) -&
           digits(0._sp)) * 0.5_sp))
     real(sp), parameter :: sbig = real(radix(0._sp), wp)**( - ceiling(       (maxexponent(0._sp)&
           + digits(0._sp) - 1) * 0.5_sp))
        ! ..
        ! .. scalar arguments ..
        integer :: incx, n
        ! ..
        ! .. array arguments ..
        complex(sp) :: x(*)
        ! ..
        ! .. local scalars ..
        integer :: i, ix
        logical :: notbig
        real(sp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
     
        ! quick return if possible
     
        stdlib_scnrm2 =szero
        if( n <= 0 ) return
     
        scl = sone
        sumsq =szero
     
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
     
        notbig = .true.
        asml =szero
        amed =szero
        abig =szero
        ix = 1
        if( incx < 0 ) ix = 1 - (n-1)*incx
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
     
        ! combine abig and amed or amed and asml if more than sone
        ! accumulator was used.
     
        if (abig >szero) then
     
           ! combine abig and amed if abig > 0.
     
           if ( (amed >szero) .or. (amed > maxn) .or. (amed /= amed) ) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = sone / sbig
           sumsq = abig
        else if (asml >szero) then
     
           ! combine amed and asml if asml > 0.
     
           if ( (amed >szero) .or. (amed > maxn) .or. (amed /= amed) ) then
              amed = sqrt(amed)
              asml = sqrt(asml) / ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = sone
              sumsq = ymax**2*( sone + (ymin/ymax)**2 )
           else
              scl = sone / ssml
              sumsq = asml
           end if
        else
     
           ! otherwise all values are mid-range
     
           scl = sone
           sumsq = amed
        end if
        stdlib_scnrm2 = scl*sqrt( sumsq )
        return
     end function stdlib_scnrm2

     
     
     subroutine stdlib_scopy(n,sx,incx,sy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,7)
              if (m/=0) then
                 do i = 1,m
                    sy(i) = sx(i)
                 end do
                 if (n<7) return
              end if
              mp1 = m + 1
              do i = mp1,n,7
                 sy(i) = sx(i)
                 sy(i+1) = sx(i+1)
                 sy(i+2) = sx(i+2)
                 sy(i+3) = sx(i+3)
                 sy(i+4) = sx(i+4)
                 sy(i+5) = sx(i+5)
                 sy(i+6) = sx(i+6)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 sy(iy) = sx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_scopy
     
     end subroutine stdlib_scopy

     
     
     real(sp) function stdlib_sdot(n,sx,incx,sy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) stemp
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           stemp = 0.0e0
           stdlib_sdot = 0.0e0
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,5)
              if (m/=0) then
                 do i = 1,m
                    stemp = stemp + sx(i)*sy(i)
                 end do
                 if (n<5) then
                    stdlib_sdot=stemp
                 return
                 end if
              end if
              mp1 = m + 1
              do i = mp1,n,5
               stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) +&
           sx(i+4)*sy(i+4)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 stemp = stemp + sx(ix)*sy(iy)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           stdlib_sdot = stemp
           return
     
           ! end of stdlib_sdot
     
     end function stdlib_sdot

     
     
     real(sp) function stdlib_sdsdot(n,sb,sx,incx,sy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) sb
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! .. local scalars ..
           real(dp) stdlib_dsdot
           integer(int32) i,kx,ky,ns
           ! ..
           ! .. intrinsic functions ..
           intrinsic dble
           ! ..
           stdlib_dsdot = sb
           if (n<=0) then
              stdlib_sdsdot = stdlib_dsdot
              return
           end if
           if (incx==incy .and. incx>0) then
     
           ! code for equal and positive increments.
     
              ns = n*incx
              do i = 1,ns,incx
                 stdlib_dsdot = stdlib_dsdot + dble(sx(i))*dble(sy(i))
              end do
           else
     
           ! code for unequal or nonpositive increments.
     
              kx = 1
              ky = 1
              if (incx<0) kx = 1 + (1-n)*incx
              if (incy<0) ky = 1 + (1-n)*incy
              do i = 1,n
                 stdlib_dsdot = stdlib_dsdot + dble(sx(kx))*dble(sy(ky))
                 kx = kx + incx
                 ky = ky + incy
              end do
           end if
           stdlib_sdsdot = stdlib_dsdot
           return
     
           ! end of stdlib_sdsdot
     
     end function stdlib_sdsdot

     
     
     function stdlib_snrm2( n, x, incx )
        integer, parameter :: wp = kind(1.e0)
        real(sp) :: stdlib_snrm2
     
        ! -- reference blas level1 routine (version 3.9.1) --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! march 2021
     
        ! .. constants ..
        real(sp), parameter ::szero = 0.0_sp
        real(sp), parameter :: sone  = 1.0_sp
        real(sp), parameter :: maxn = huge(0.0_sp)
        ! ..
        ! .. blue's scaling constants ..
     real(sp), parameter :: tsml = real(radix(0._sp), wp)**ceiling(       (minexponent(0._sp) - 1)&
           * 0.5_sp)
     real(sp), parameter :: tbig = real(radix(0._sp), wp)**floor(       (maxexponent(0._sp) -&
           digits(0._sp) + 1) * 0.5_sp)
     real(sp), parameter :: ssml = real(radix(0._sp), wp)**( - floor(       (minexponent(0._sp) -&
           digits(0._sp)) * 0.5_sp))
     real(sp), parameter :: sbig = real(radix(0._sp), wp)**( - ceiling(       (maxexponent(0._sp)&
           + digits(0._sp) - 1) * 0.5_sp))
        ! ..
        ! .. scalar arguments ..
        integer :: incx, n
        ! ..
        ! .. array arguments ..
        real(sp) :: x(*)
        ! ..
        ! .. local scalars ..
        integer :: i, ix
        logical :: notbig
        real(sp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
     
        ! quick return if possible
     
        stdlib_snrm2 =szero
        if( n <= 0 ) return
     
        scl = sone
        sumsq =szero
     
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
     
        notbig = .true.
        asml =szero
        amed =szero
        abig =szero
        ix = 1
        if( incx < 0 ) ix = 1 - (n-1)*incx
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
     
        ! combine abig and amed or amed and asml if more than sone
        ! accumulator was used.
     
        if (abig >szero) then
     
           ! combine abig and amed if abig > 0.
     
           if ( (amed >szero) .or. (amed > maxn) .or. (amed /= amed) ) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = sone / sbig
           sumsq = abig
        else if (asml >szero) then
     
           ! combine amed and asml if asml > 0.
     
           if ( (amed >szero) .or. (amed > maxn) .or. (amed /= amed) ) then
              amed = sqrt(amed)
              asml = sqrt(asml) / ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = sone
              sumsq = ymax**2*( sone + (ymin/ymax)**2 )
           else
              scl = sone / ssml
              sumsq = asml
           end if
        else
     
           ! otherwise all values are mid-range
     
           scl = sone
           sumsq = amed
        end if
        stdlib_snrm2 = scl*sqrt( sumsq )
        return
     end function stdlib_snrm2

     
     
     subroutine stdlib_srot(n,sx,incx,sy,incy,c,s)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) c,s
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) stemp
           integer(int32) i,ix,iy
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
             ! code for both increments equal to 1
     
              do i = 1,n
                 stemp = c*sx(i) + s*sy(i)
                 sy(i) = c*sy(i) - s*sx(i)
                 sx(i) = stemp
              end do
           else
     
             ! code for unequal increments or equal increments not equal
               ! to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 stemp = c*sx(ix) + s*sy(iy)
                 sy(iy) = c*sy(iy) - s*sx(ix)
                 sx(ix) = stemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_srot
     
     end subroutine stdlib_srot

     
     
     subroutine stdlib_srotg( a, b, c, s )
        integer, parameter :: wp = kind(1.e0)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
        ! .. constants ..
        real(sp), parameter ::szero = 0.0_sp
        real(sp), parameter :: sone  = 1.0_sp
        ! ..
        ! .. scaling constants ..
     real(sp), parameter :: safmin = real(radix(0._sp),wp)**max(minexponent(0._sp)-1,1-&
          maxexponent(0._sp)   )
     real(sp), parameter :: safmax = real(radix(0._sp),wp)**max(1-minexponent(0._sp),maxexponent(&
          0._sp)-1   )
        ! ..
        ! .. scalar arguments ..
        real(sp) :: a, b, c, s
        ! ..
        ! .. local scalars ..
        real(sp) :: anorm, bnorm, scl, sigma, r, z
        ! ..
        anorm = abs(a)
        bnorm = abs(b)
        if( bnorm ==szero ) then
           c = sone
           s =szero
           b =szero
        else if( anorm ==szero ) then
           c =szero
           s = sone
           a = b
           b = sone
        else
           scl = min( safmax, max( safmin, anorm, bnorm ) )
           if( anorm > bnorm ) then
              sigma = sign(sone,a)
           else
              sigma = sign(sone,b)
           end if
           r = sigma*( scl*sqrt((a/scl)**2 + (b/scl)**2) )
           c = a/r
           s = b/r
           if( anorm > bnorm ) then
              z = s
           else if( c /=szero ) then
              z = sone/c
           else
              z = sone
           end if
           a = r
           b = z
        end if
        return
     end subroutine stdlib_srotg

     
     
     subroutine stdlib_srotm(n,sx,incx,sy,incy,sparam)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sparam(5),sx(*),sy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) sflag,sh11,sh12,sh21,sh22,two,w,z,zero
           integer(int32) i,kx,ky,nsteps
           ! ..
           ! .. data statements ..
           data zero,two/0.e0,2.e0/
           ! ..
     
           sflag = sparam(1)
           if (n<=0 .or. (sflag+two==zero)) return
           if (incx==incy.and.incx>0) then
     
              nsteps = n*incx
              if (sflag<zero) then
                 sh11 = sparam(2)
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 sh22 = sparam(5)
                 do i = 1,nsteps,incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z*sh12
                    sy(i) = w*sh21 + z*sh22
                 end do
              else if (sflag==zero) then
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 do i = 1,nsteps,incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w + z*sh12
                    sy(i) = w*sh21 + z
                 end do
              else
                 sh11 = sparam(2)
                 sh22 = sparam(5)
                 do i = 1,nsteps,incx
                    w = sx(i)
                    z = sy(i)
                    sx(i) = w*sh11 + z
                    sy(i) = -w + sh22*z
                 end do
              end if
           else
              kx = 1
              ky = 1
              if (incx<0) kx = 1 + (1-n)*incx
              if (incy<0) ky = 1 + (1-n)*incy
     
              if (sflag<zero) then
                 sh11 = sparam(2)
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 sh22 = sparam(5)
                 do i = 1,n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w*sh11 + z*sh12
                    sy(ky) = w*sh21 + z*sh22
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else if (sflag==zero) then
                 sh12 = sparam(4)
                 sh21 = sparam(3)
                 do i = 1,n
                    w = sx(kx)
                    z = sy(ky)
                    sx(kx) = w + z*sh12
                    sy(ky) = w*sh21 + z
                    kx = kx + incx
                    ky = ky + incy
                 end do
              else
                  sh11 = sparam(2)
                  sh22 = sparam(5)
                  do i = 1,n
                     w = sx(kx)
                     z = sy(ky)
                     sx(kx) = w*sh11 + z
                     sy(ky) = -w + sh22*z
                     kx = kx + incx
                     ky = ky + incy
                 end do
              end if
           end if
           return
     
           ! end of stdlib_srotm
     
     end subroutine stdlib_srotm

     
     
     subroutine stdlib_srotmg(sd1,sd2,sx1,sy1,sparam)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) sd1,sd2,sx1,sy1
           ! ..
           ! .. array arguments ..
           real(sp) sparam(5)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) gam,gamsq,one,rgamsq,sflag,sh11,sh12,sh21,sh22,sp1,sp2,sq1,sq2,stemp,su,two,&
          zero
           ! ..
           ! .. intrinsic functions ..
           intrinsic abs
           ! ..
           ! .. data statements ..
     
           data zero,one,two/0.e0,1.e0,2.e0/
           data gam,gamsq,rgamsq/4096.e0,1.67772e7,5.96046e-8/
           ! ..
           if (sd1<zero) then
              ! go zero-h-d-and-sx1..
              sflag = -one
              sh11 = zero
              sh12 = zero
              sh21 = zero
              sh22 = zero
     
              sd1 = zero
              sd2 = zero
              sx1 = zero
           else
              ! case-sd1-nonnegative
              sp2 = sd2*sy1
              if (sp2==zero) then
                 sflag = -two
                 sparam(1) = sflag
                 return
              end if
              ! regular-case..
              sp1 = sd1*sx1
              sq2 = sp2*sy1
              sq1 = sp1*sx1
     
              if (abs(sq1)>abs(sq2)) then
                 sh21 = -sy1/sx1
                 sh12 = sp2/sp1
     
                 su = one - sh12*sh21
     
                if (su>zero) then
                  sflag = zero
                  sd1 = sd1/su
                  sd2 = sd2/su
                  sx1 = sx1*su
                else
                  ! this code path if here for safety. we do not expect this
                  ! condition to ever hold except in edge cases with rounding
                  ! errors. see doi: 10.1145/355841.355847
                  sflag = -one
                  sh11 = zero
                  sh12 = zero
                  sh21 = zero
                  sh22 = zero
     
                  sd1 = zero
                  sd2 = zero
                  sx1 = zero
                end if
              else
                 if (sq2<zero) then
                    ! go zero-h-d-and-sx1..
                    sflag = -one
                    sh11 = zero
                    sh12 = zero
                    sh21 = zero
                    sh22 = zero
     
                    sd1 = zero
                    sd2 = zero
                    sx1 = zero
                 else
                    sflag = one
                    sh11 = sp1/sp2
                    sh22 = sx1/sy1
                    su = one + sh11*sh22
                    stemp = sd2/su
                    sd2 = sd1/su
                    sd1 = stemp
                    sx1 = sy1*su
                 end if
              end if
           ! procedure..scale-check
              if (sd1/=zero) then
                 do while ((sd1<=rgamsq) .or. (sd1>=gamsq))
                    if (sflag==zero) then
                       sh11 = one
                       sh22 = one
                       sflag = -one
                    else
                       sh21 = -one
                       sh12 = one
                       sflag = -one
                    end if
                    if (sd1<=rgamsq) then
                       sd1 = sd1*gam**2
                       sx1 = sx1/gam
                       sh11 = sh11/gam
                       sh12 = sh12/gam
                    else
                       sd1 = sd1/gam**2
                       sx1 = sx1*gam
                       sh11 = sh11*gam
                       sh12 = sh12*gam
                    end if
                 enddo
              end if
              if (sd2/=zero) then
                 do while ( (abs(sd2)<=rgamsq) .or. (abs(sd2)>=gamsq) )
                    if (sflag==zero) then
                       sh11 = one
                       sh22 = one
                       sflag = -one
                    else
                       sh21 = -one
                       sh12 = one
                       sflag = -one
                    end if
                    if (abs(sd2)<=rgamsq) then
                       sd2 = sd2*gam**2
                       sh21 = sh21/gam
                       sh22 = sh22/gam
                    else
                       sd2 = sd2/gam**2
                       sh21 = sh21*gam
                       sh22 = sh22*gam
                    end if
                 end do
              end if
           end if
           if (sflag<zero) then
              sparam(2) = sh11
              sparam(3) = sh21
              sparam(4) = sh12
              sparam(5) = sh22
           else if (sflag==zero) then
              sparam(3) = sh21
              sparam(4) = sh12
           else
              sparam(2) = sh11
              sparam(5) = sh22
           end if
           sparam(1) = sflag
           return
     
           ! end of stdlib_srotmg
     
     end subroutine stdlib_srotmg

     
     
     subroutine stdlib_sscal(n,sa,sx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) sa
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,m,mp1,nincx
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0 .or. incx<=0) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
     
              ! clean-up loop
     
              m = mod(n,5)
              if (m/=0) then
                 do i = 1,m
                    sx(i) = sa*sx(i)
                 end do
                 if (n<5) return
              end if
              mp1 = m + 1
              do i = mp1,n,5
                 sx(i) = sa*sx(i)
                 sx(i+1) = sa*sx(i+1)
                 sx(i+2) = sa*sx(i+2)
                 sx(i+3) = sa*sx(i+3)
                 sx(i+4) = sa*sx(i+4)
              end do
           else
     
              ! code for increment not equal to 1
     
              nincx = n*incx
              do i = 1,nincx,incx
                 sx(i) = sa*sx(i)
              end do
           end if
           return
     
           ! end of stdlib_sscal
     
     end subroutine stdlib_sscal

     
     
     subroutine stdlib_sswap(n,sx,incx,sy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           real(sp) sx(*),sy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) stemp
           integer(int32) i,ix,iy,m,mp1
           ! ..
           ! .. intrinsic functions ..
           intrinsic mod
           ! ..
           if (n<=0) return
           if (incx==1 .and. incy==1) then
     
             ! code for both increments equal to 1
     
     
             ! clean-up loop
     
              m = mod(n,3)
              if (m/=0) then
                 do i = 1,m
                    stemp = sx(i)
                    sx(i) = sy(i)
                    sy(i) = stemp
                 end do
                 if (n<3) return
              end if
              mp1 = m + 1
              do i = mp1,n,3
                 stemp = sx(i)
                 sx(i) = sy(i)
                 sy(i) = stemp
                 stemp = sx(i+1)
                 sx(i+1) = sy(i+1)
                 sy(i+1) = stemp
                 stemp = sx(i+2)
                 sx(i+2) = sy(i+2)
                 sy(i+2) = stemp
              end do
           else
     
             ! code for unequal increments or equal increments not equal
               ! to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 stemp = sx(ix)
                 sx(ix) = sy(iy)
                 sy(iy) = stemp
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
           return
     
           ! end of stdlib_sswap
     
     end subroutine stdlib_sswap

     
     
     subroutine stdlib_xerbla( srname, info )
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           character*(*)      srname
           integer(int32)            info
           ! ..
     
       ! =====================================================================
     
           ! .. intrinsic functions ..
           intrinsic          len_trim
           ! ..
           ! .. executable statements ..
     
           write( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
     
           stop
     
      9999 format( ' ** on entry to ', a, ' parameter number ', i2, ' had ','an illegal value' )
     
           ! end of stdlib_xerbla
     
     end subroutine stdlib_xerbla

     
     
     subroutine stdlib_xerbla_array(srname_array, srname_len, info)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) srname_len, info
           ! ..
           ! .. array arguments ..
           character(1) srname_array(srname_len)
           ! ..
     
       ! =====================================================================
     
           ! ..
           ! .. local scalars ..
           integer(int32) i
           ! ..
           ! .. local arrays ..
           character*32 srname
           ! ..
           ! .. intrinsic functions ..
           intrinsic min, len
           ! ..
     
           ! .. executable statements ..
           srname = ''
           do i = 1, min( srname_len, len( srname ) )
              srname( i:i ) = srname_array( i )
           end do
           call stdlib_xerbla( srname, info )
           return
     
           ! end of stdlib_xerbla_array
     
     end subroutine stdlib_xerbla_array

     
     
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
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('zgbmv ',info)
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kup1 = ku + 1
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       k = kup1 - j
                       do 50 i = max(1,j-ku),min(m,j+kl)
                           y(i) = y(i) + temp*a(k+i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       k = kup1 - j
                       do 70 i = max(1,j-ku),min(m,j+kl)
                           y(iy) = y(iy) + temp*a(k+i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
                       if (j>ku) ky = ky + incy
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
     
               jy = ky
               if (incx==1) then
                   do 110 j = 1,n
                       temp = zero
                       k = kup1 - j
                       if (noconj) then
                           do 90 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + a(k+i,j)*x(i)
        90                 continue
                       else
                           do 100 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + dconjg(a(k+i,j))*x(i)
       100                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       110         continue
               else
                   do 140 j = 1,n
                       temp = zero
                       ix = kx
                       k = kup1 - j
                       if (noconj) then
                           do 120 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + a(k+i,j)*x(ix)
                               ix = ix + incx
       120                 continue
                       else
                           do 130 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + dconjg(a(k+i,j))*x(ix)
                               ix = ix + incx
       130                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j>ku) kx = kx + incx
       140         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zgbmv
     
     end subroutine stdlib_zgbmv

     
     
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
           else if ((.not.notb) .and. (.not.conjb) .and.(.not.stdlib_lsame(transb,'t')))&
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
               call stdlib_xerbla('zgemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.(((alpha==zero).or. (k==0)).and. (beta==one)))&
           return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (notb) then
               if (nota) then
     
                 ! form  c := alpha*a*b + beta*c.
     
                   do 90 j = 1,n
                       if (beta==zero) then
                           do 50 i = 1,m
                               c(i,j) = zero
        50                 continue
                       else if (beta/=one) then
                           do 60 i = 1,m
                               c(i,j) = beta*c(i,j)
        60                 continue
                       end if
                       do 80 l = 1,k
                           temp = alpha*b(l,j)
                           do 70 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
        70                 continue
        80             continue
        90         continue
               else if (conja) then
     
                 ! form  c := alpha*a**h*b + beta*c.
     
                   do 120 j = 1,n
                       do 110 i = 1,m
                           temp = zero
                           do 100 l = 1,k
                               temp = temp + dconjg(a(l,i))*b(l,j)
       100                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       110             continue
       120         continue
               else
     
                 ! form  c := alpha*a**t*b + beta*c
     
                   do 150 j = 1,n
                       do 140 i = 1,m
                           temp = zero
                           do 130 l = 1,k
                               temp = temp + a(l,i)*b(l,j)
       130                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       140             continue
       150         continue
               end if
           else if (nota) then
               if (conjb) then
     
                 ! form  c := alpha*a*b**h + beta*c.
     
                   do 200 j = 1,n
                       if (beta==zero) then
                           do 160 i = 1,m
                               c(i,j) = zero
       160                 continue
                       else if (beta/=one) then
                           do 170 i = 1,m
                               c(i,j) = beta*c(i,j)
       170                 continue
                       end if
                       do 190 l = 1,k
                           temp = alpha*dconjg(b(j,l))
                           do 180 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
       180                 continue
       190             continue
       200         continue
               else
     
                 ! form  c := alpha*a*b**t + beta*c
     
                   do 250 j = 1,n
                       if (beta==zero) then
                           do 210 i = 1,m
                               c(i,j) = zero
       210                 continue
                       else if (beta/=one) then
                           do 220 i = 1,m
                               c(i,j) = beta*c(i,j)
       220                 continue
                       end if
                       do 240 l = 1,k
                           temp = alpha*b(j,l)
                           do 230 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
       230                 continue
       240             continue
       250         continue
               end if
           else if (conja) then
               if (conjb) then
     
                 ! form  c := alpha*a**h*b**h + beta*c.
     
                   do 280 j = 1,n
                       do 270 i = 1,m
                           temp = zero
                           do 260 l = 1,k
                               temp = temp + dconjg(a(l,i))*dconjg(b(j,l))
       260                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       270             continue
       280         continue
               else
     
                 ! form  c := alpha*a**h*b**t + beta*c
     
                   do 310 j = 1,n
                       do 300 i = 1,m
                           temp = zero
                           do 290 l = 1,k
                               temp = temp + dconjg(a(l,i))*b(j,l)
       290                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       300             continue
       310         continue
               end if
           else
               if (conjb) then
     
                 ! form  c := alpha*a**t*b**h + beta*c
     
                   do 340 j = 1,n
                       do 330 i = 1,m
                           temp = zero
                           do 320 l = 1,k
                               temp = temp + a(l,i)*dconjg(b(j,l))
       320                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       330             continue
       340         continue
               else
     
                 ! form  c := alpha*a**t*b**t + beta*c
     
                   do 370 j = 1,n
                       do 360 i = 1,m
                           temp = zero
                           do 350 l = 1,k
                               temp = temp + a(l,i)*b(j,l)
       350                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       360             continue
       370         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zgemm
     
     end subroutine stdlib_zgemm

     
     
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
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('zgemv ',info)
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       do 50 i = 1,m
                           y(i) = y(i) + temp*a(i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       do 70 i = 1,m
                           y(iy) = y(iy) + temp*a(i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
     
               jy = ky
               if (incx==1) then
                   do 110 j = 1,n
                       temp = zero
                       if (noconj) then
                           do 90 i = 1,m
                               temp = temp + a(i,j)*x(i)
        90                 continue
                       else
                           do 100 i = 1,m
                               temp = temp + dconjg(a(i,j))*x(i)
       100                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       110         continue
               else
                   do 140 j = 1,n
                       temp = zero
                       ix = kx
                       if (noconj) then
                           do 120 i = 1,m
                               temp = temp + a(i,j)*x(ix)
                               ix = ix + incx
       120                 continue
                       else
                           do 130 i = 1,m
                               temp = temp + dconjg(a(i,j))*x(ix)
                               ix = ix + incx
       130                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       140         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zgemv
     
     end subroutine stdlib_zgemv

     
     
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
               call stdlib_xerbla('zgerc ',info)
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
               do 20 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*dconjg(y(jy))
                       do 10 i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
        10             continue
                   end if
                   jy = jy + incy
        20     continue
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               do 40 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*dconjg(y(jy))
                       ix = kx
                       do 30 i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
        30             continue
                   end if
                   jy = jy + incy
        40     continue
           end if
     
           return
     
           ! end of stdlib_zgerc
     
     end subroutine stdlib_zgerc

     
     
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
               call stdlib_xerbla('zgeru ',info)
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
               do 20 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       do 10 i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
        10             continue
                   end if
                   jy = jy + incy
        20     continue
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               do 40 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       ix = kx
                       do 30 i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
        30             continue
                   end if
                   jy = jy + incy
        40     continue
           end if
     
           return
     
           ! end of stdlib_zgeru
     
     end subroutine stdlib_zgeru

     
     
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
               call stdlib_xerbla('zhbmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when upper triangle of a is stored.
     
               kplus1 = k + 1
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       l = kplus1 - j
                       do 50 i = max(1,j-k),j - 1
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(i)
        50             continue
                       y(j) = y(j) + temp1*dble(a(kplus1,j)) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       do 70 i = max(1,j-k),j - 1
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*dble(a(kplus1,j)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       if (j>k) then
                           kx = kx + incx
                           ky = ky + incy
                       end if
        80         continue
               end if
           else
     
              ! form  y  when lower triangle of a is stored.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*dble(a(1,j))
                       l = 1 - j
                       do 90 i = j + 1,min(n,j+k)
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*dble(a(1,j))
                       l = 1 - j
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,min(n,j+k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + dconjg(a(l+i,j))*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zhbmv
     
     end subroutine stdlib_zhbmv

     
     
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
               call stdlib_xerbla('zhemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   do 70 j = 1,n
                       do 60 i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 50 k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*dconjg(a(k,i))
        50                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*dble(a(i,i)) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*dble(a(i,i)) +alpha*temp2
                           end if
        60             continue
        70         continue
               else
                   do 100 j = 1,n
                       do 90 i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 80 k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*dconjg(a(k,i))
        80                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*dble(a(i,i)) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*dble(a(i,i)) +alpha*temp2
                           end if
        90             continue
       100         continue
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               do 170 j = 1,n
                   temp1 = alpha*dble(a(j,j))
                   if (beta==zero) then
                       do 110 i = 1,m
                           c(i,j) = temp1*b(i,j)
       110             continue
                   else
                       do 120 i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
       120             continue
                   end if
                   do 140 k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*dconjg(a(j,k))
                       end if
                       do 130 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       130             continue
       140         continue
                   do 160 k = j + 1,n
                       if (upper) then
                           temp1 = alpha*dconjg(a(j,k))
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do 150 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       150             continue
       160         continue
       170     continue
           end if
     
           return
     
           ! end of stdlib_zhemm
     
     end subroutine stdlib_zhemm

     
     
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
               call stdlib_xerbla('zhemv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when a is stored in upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(i)
        50             continue
                       y(j) = y(j) + temp1*dble(a(j,j)) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 i = 1,j - 1
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*dble(a(j,j)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           else
     
              ! form  y  when a is stored in lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*dble(a(j,j))
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*dble(a(j,j))
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + dconjg(a(i,j))*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zhemv
     
     end subroutine stdlib_zhemv

     
     
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
               call stdlib_xerbla('zher  ',info)
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           do 10 i = 1,j - 1
                               a(i,j) = a(i,j) + x(i)*temp
        10                 continue
                           a(j,j) = dble(a(j,j)) + dble(x(j)*temp)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           ix = kx
                           do 30 i = 1,j - 1
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                           a(j,j) = dble(a(j,j)) + dble(x(jx)*temp)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
        40         continue
               end if
           else
     
              ! form  a  when a is stored in lower triangle.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           a(j,j) = dble(a(j,j)) + dble(temp*x(j))
                           do 50 i = j + 1,n
                               a(i,j) = a(i,j) + x(i)*temp
        50                 continue
                       else
                           a(j,j) = dble(a(j,j))
                       end if
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           a(j,j) = dble(a(j,j)) + dble(temp*x(jx))
                           ix = jx
                           do 70 i = j + 1,n
                               ix = ix + incx
                               a(i,j) = a(i,j) + x(ix)*temp
        70                 continue
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zher
     
     end subroutine stdlib_zher

     
     
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
               call stdlib_xerbla('zher2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           do 10 i = 1,j - 1
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        10                 continue
                           a(j,j) = dble(a(j,j)) +dble(x(j)*temp1+y(j)*temp2)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           do 30 i = 1,j - 1
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                           a(j,j) = dble(a(j,j)) +dble(x(jx)*temp1+y(jy)*temp2)
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
                       jy = jy + incy
        40         continue
               end if
           else
     
              ! form  a  when a is stored in the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           a(j,j) = dble(a(j,j)) +dble(x(j)*temp1+y(j)*temp2)
                           do 50 i = j + 1,n
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        50                 continue
                       else
                           a(j,j) = dble(a(j,j))
                       end if
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           a(j,j) = dble(a(j,j)) +dble(x(jx)*temp1+y(jy)*temp2)
                           ix = jx
                           iy = jy
                           do 70 i = j + 1,n
                               ix = ix + incx
                               iy = iy + incy
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
        70                 continue
                       else
                           a(j,j) = dble(a(j,j))
                       end if
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zher2
     
     end subroutine stdlib_zher2

     
     
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'c')))&
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
               call stdlib_xerbla('zher2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==dble(zero)) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
        30                 continue
                           c(j,j) = beta*dble(c(j,j))
        40             continue
                   end if
               else
                   if (beta==dble(zero)) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           c(j,j) = beta*dble(c(j,j))
                           do 70 i = j + 1,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==dble(zero)) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
       100                 continue
                           c(j,j) = beta*dble(c(j,j))
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       do 120 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*dconjg(b(j,l))
                               temp2 = dconjg(alpha*a(j,l))
                               do 110 i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       110                     continue
                               c(j,j) = dble(c(j,j)) +dble(a(j,l)*temp1+b(j,l)*temp2)
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==dble(zero)) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j + 1,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                           c(j,j) = beta*dble(c(j,j))
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       do 170 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*dconjg(b(j,l))
                               temp2 = dconjg(alpha*a(j,l))
                               do 160 i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       160                     continue
                               c(j,j) = dble(c(j,j)) +dble(a(j,l)*temp1+b(j,l)*temp2)
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do 190 l = 1,k
                               temp1 = temp1 + dconjg(a(l,i))*b(l,j)
                               temp2 = temp2 + dconjg(b(l,i))*a(l,j)
       190                 continue
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
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp1 = zero
                           temp2 = zero
                           do 220 l = 1,k
                               temp1 = temp1 + dconjg(a(l,i))*b(l,j)
                               temp2 = temp2 + dconjg(b(l,i))*a(l,j)
       220                 continue
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
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zher2k
     
     end subroutine stdlib_zher2k

     
     
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'c')))&
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
               call stdlib_xerbla('zherk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
        30                 continue
                           c(j,j) = beta*dble(c(j,j))
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           c(j,j) = beta*dble(c(j,j))
                           do 70 i = j + 1,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**h + beta*c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
       100                 continue
                           c(j,j) = beta*dble(c(j,j))
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       do 120 l = 1,k
                           if (a(j,l)/=dcmplx(zero)) then
                               temp = alpha*dconjg(a(j,l))
                               do 110 i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
       110                     continue
                               c(j,j) = dble(c(j,j)) + dble(temp*a(i,l))
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           c(j,j) = beta*dble(c(j,j))
                           do 150 i = j + 1,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       else
                           c(j,j) = dble(c(j,j))
                       end if
                       do 170 l = 1,k
                           if (a(j,l)/=dcmplx(zero)) then
                               temp = alpha*dconjg(a(j,l))
                               c(j,j) = dble(c(j,j)) + dble(temp*a(j,l))
                               do 160 i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**h*a + beta*c.
     
               if (upper) then
                   do 220 j = 1,n
                       do 200 i = 1,j - 1
                           temp = zero
                           do 190 l = 1,k
                               temp = temp + dconjg(a(l,i))*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       200             continue
                       rtemp = zero
                       do 210 l = 1,k
                           rtemp = rtemp + dconjg(a(l,j))*a(l,j)
       210             continue
                       if (beta==zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*dble(c(j,j))
                       end if
       220         continue
               else
                   do 260 j = 1,n
                       rtemp = zero
                       do 230 l = 1,k
                           rtemp = rtemp + dconjg(a(l,j))*a(l,j)
       230             continue
                       if (beta==zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*dble(c(j,j))
                       end if
                       do 250 i = j + 1,n
                           temp = zero
                           do 240 l = 1,k
                               temp = temp + dconjg(a(l,i))*a(l,j)
       240                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       250             continue
       260         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zherk
     
     end subroutine stdlib_zherk

     
     
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
               call stdlib_xerbla('zhpmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when ap contains the upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       k = kk
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(i)
                           k = k + 1
        50             continue
                       y(j) = y(j) + temp1*dble(ap(kk+j-1)) + alpha*temp2
                       kk = kk + j
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 k = kk,kk + j - 2
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*dble(ap(kk+j-1)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        80         continue
               end if
           else
     
              ! form  y  when ap contains the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*dble(ap(kk))
                       k = kk + 1
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(i)
                           k = k + 1
        90             continue
                       y(j) = y(j) + alpha*temp2
                       kk = kk + (n-j+1)
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*dble(ap(kk))
                       ix = jx
                       iy = jy
                       do 110 k = kk + 1,kk + n - j
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + dconjg(ap(k))*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + (n-j+1)
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zhpmv
     
     end subroutine stdlib_zhpmv

     
     
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
               call stdlib_xerbla('zhpr  ',info)
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           k = kk
                           do 10 i = 1,j - 1
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        10                 continue
                           ap(kk+j-1) = dble(ap(kk+j-1)) + dble(x(j)*temp)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       kk = kk + j
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           ix = kx
                           do 30 k = kk,kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                           ap(kk+j-1) = dble(ap(kk+j-1)) + dble(x(jx)*temp)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       jx = jx + incx
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*dconjg(x(j))
                           ap(kk) = dble(ap(kk)) + dble(temp*x(j))
                           k = kk + 1
                           do 50 i = j + 1,n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        50                 continue
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*dconjg(x(jx))
                           ap(kk) = dble(ap(kk)) + dble(temp*x(jx))
                           ix = jx
                           do 70 k = kk + 1,kk + n - j
                               ix = ix + incx
                               ap(k) = ap(k) + x(ix)*temp
        70                 continue
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zhpr
     
     end subroutine stdlib_zhpr

     
     
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
               call stdlib_xerbla('zhpr2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           k = kk
                           do 10 i = 1,j - 1
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        10                 continue
                           ap(kk+j-1) = dble(ap(kk+j-1)) +dble(x(j)*temp1+y(j)*temp2)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       kk = kk + j
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           do 30 k = kk,kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                           ap(kk+j-1) = dble(ap(kk+j-1)) +dble(x(jx)*temp1+y(jy)*temp2)
                       else
                           ap(kk+j-1) = dble(ap(kk+j-1))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*dconjg(y(j))
                           temp2 = dconjg(alpha*x(j))
                           ap(kk) = dble(ap(kk)) +dble(x(j)*temp1+y(j)*temp2)
                           k = kk + 1
                           do 50 i = j + 1,n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        50                 continue
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*dconjg(y(jy))
                           temp2 = dconjg(alpha*x(jx))
                           ap(kk) = dble(ap(kk)) +dble(x(jx)*temp1+y(jy)*temp2)
                           ix = jx
                           iy = jy
                           do 70 k = kk + 1,kk + n - j
                               ix = ix + incx
                               iy = iy + incy
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
        70                 continue
                       else
                           ap(kk) = dble(ap(kk))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zhpr2
     
     end subroutine stdlib_zhpr2

     
     
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
               call stdlib_xerbla('zsymm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   do 70 j = 1,n
                       do 60 i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 50 k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        50                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        60             continue
        70         continue
               else
                   do 100 j = 1,n
                       do 90 i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 80 k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        80                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        90             continue
       100         continue
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               do 170 j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta==zero) then
                       do 110 i = 1,m
                           c(i,j) = temp1*b(i,j)
       110             continue
                   else
                       do 120 i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
       120             continue
                   end if
                   do 140 k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do 130 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       130             continue
       140         continue
                   do 160 k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do 150 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       150             continue
       160         continue
       170     continue
           end if
     
           return
     
           ! end of stdlib_zsymm
     
     end subroutine stdlib_zsymm

     
     
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')))&
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
               call stdlib_xerbla('zsyr2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do 190 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp1 = zero
                           temp2 = zero
                           do 220 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zsyr2k
     
     end subroutine stdlib_zsyr2k

     
     
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')))&
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
               call stdlib_xerbla('zsyrk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**t + beta*c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*a + beta*c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp = zero
                           do 190 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp = zero
                           do 220 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_zsyrk
     
     end subroutine stdlib_zsyrk

     
     
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
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ztbmv ',info)
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do 10 i = max(1,j-k),j - 1
                                   x(i) = x(i) + temp*a(l+i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do 30 i = max(1,j-k),j - 1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j>k) kx = kx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = 1 - j
                               do 50 i = min(n,j+k),j + 1,-1
                                   x(i) = x(i) + temp*a(l+i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do 70 i = min(n,j+k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 110 j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do 90 i = j - 1,max(1,j-k),-1
                                   temp = temp + a(l+i,j)*x(i)
        90                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(kplus1,j))
                               do 100 i = j - 1,max(1,j-k),-1
                                   temp = temp + dconjg(a(l+i,j))*x(i)
       100                     continue
                           end if
                           x(j) = temp
       110             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 140 j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do 120 i = j - 1,max(1,j-k),-1
                                   temp = temp + a(l+i,j)*x(ix)
                                   ix = ix - incx
       120                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(kplus1,j))
                               do 130 i = j - 1,max(1,j-k),-1
                                   temp = temp + dconjg(a(l+i,j))*x(ix)
                                   ix = ix - incx
       130                     continue
                           end if
                           x(jx) = temp
                           jx = jx - incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do 150 i = j + 1,min(n,j+k)
                                   temp = temp + a(l+i,j)*x(i)
       150                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(1,j))
                               do 160 i = j + 1,min(n,j+k)
                                   temp = temp + dconjg(a(l+i,j))*x(i)
       160                     continue
                           end if
                           x(j) = temp
       170             continue
                   else
                       jx = kx
                       do 200 j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do 180 i = j + 1,min(n,j+k)
                                   temp = temp + a(l+i,j)*x(ix)
                                   ix = ix + incx
       180                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(1,j))
                               do 190 i = j + 1,min(n,j+k)
                                   temp = temp + dconjg(a(l+i,j))*x(ix)
                                   ix = ix + incx
       190                     continue
                           end if
                           x(jx) = temp
                           jx = jx + incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztbmv
     
     end subroutine stdlib_ztbmv

     
     
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
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ztbsv ',info)
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do 10 i = j - 1,max(1,j-k),-1
                                   x(i) = x(i) - temp*a(l+i,j)
        10                     continue
                           end if
        20             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 40 j = n,1,-1
                           kx = kx - incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do 30 i = j - 1,max(1,j-k),-1
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix - incx
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do 50 i = j + 1,min(n,j+k)
                                   x(i) = x(i) - temp*a(l+i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           kx = kx + incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do 70 i = j + 1,min(n,j+k)
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix + incx
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 110 j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               do 90 i = max(1,j-k),j - 1
                                   temp = temp - a(l+i,j)*x(i)
        90                     continue
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do 100 i = max(1,j-k),j - 1
                                   temp = temp - dconjg(a(l+i,j))*x(i)
       100                     continue
                               if (nounit) temp = temp/dconjg(a(kplus1,j))
                           end if
                           x(j) = temp
       110             continue
                   else
                       jx = kx
                       do 140 j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               do 120 i = max(1,j-k),j - 1
                                   temp = temp - a(l+i,j)*x(ix)
                                   ix = ix + incx
       120                     continue
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do 130 i = max(1,j-k),j - 1
                                   temp = temp - dconjg(a(l+i,j))*x(ix)
                                   ix = ix + incx
       130                     continue
                               if (nounit) temp = temp/dconjg(a(kplus1,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           if (j>k) kx = kx + incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               do 150 i = min(n,j+k),j + 1,-1
                                   temp = temp - a(l+i,j)*x(i)
       150                     continue
                               if (nounit) temp = temp/a(1,j)
                           else
                               do 160 i = min(n,j+k),j + 1,-1
                                   temp = temp - dconjg(a(l+i,j))*x(i)
       160                     continue
                               if (nounit) temp = temp/dconjg(a(1,j))
                           end if
                           x(j) = temp
       170             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 200 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               do 180 i = min(n,j+k),j + 1,-1
                                   temp = temp - a(l+i,j)*x(ix)
                                   ix = ix - incx
       180                     continue
                               if (nounit) temp = temp/a(1,j)
                           else
                               do 190 i = min(n,j+k),j + 1,-1
                                   temp = temp - dconjg(a(l+i,j))*x(ix)
                                   ix = ix - incx
       190                     continue
                               if (nounit) temp = temp/dconjg(a(1,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztbsv
     
     end subroutine stdlib_ztbsv

     
     
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
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ztpmv ',info)
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
        10                     continue
                               if (nounit) x(j) = x(j)*ap(kk+j-1)
                           end if
                           kk = kk + j
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk+j-1)
                           end if
                           jx = jx + incx
                           kk = kk + j
        40             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
        50                     continue
                               if (nounit) x(j) = x(j)*ap(kk-n+j)
                           end if
                           kk = kk - (n-j+1)
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 k = kk,kk - (n- (j+1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk-n+j)
                           end if
                           jx = jx - incx
                           kk = kk - (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 110 j = n,1,-1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 90 i = j - 1,1,-1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
        90                     continue
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               do 100 i = j - 1,1,-1
                                   temp = temp + dconjg(ap(k))*x(i)
                                   k = k - 1
       100                     continue
                           end if
                           x(j) = temp
                           kk = kk - j
       110             continue
                   else
                       jx = kx + (n-1)*incx
                       do 140 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 120 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
       120                     continue
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               do 130 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + dconjg(ap(k))*x(ix)
       130                     continue
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
       140             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 170 j = 1,n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 150 i = j + 1,n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
       150                     continue
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               do 160 i = j + 1,n
                                   temp = temp + dconjg(ap(k))*x(i)
                                   k = k + 1
       160                     continue
                           end if
                           x(j) = temp
                           kk = kk + (n-j+1)
       170             continue
                   else
                       jx = kx
                       do 200 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 180 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
       180                     continue
                           else
                               if (nounit) temp = temp*dconjg(ap(kk))
                               do 190 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + dconjg(ap(k))*x(ix)
       190                     continue
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n-j+1)
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztpmv
     
     end subroutine stdlib_ztpmv

     
     
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
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ztpsv ',info)
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
        10                     continue
                           end if
                           kk = kk - j
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 30 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
        30                     continue
                           end if
                           jx = jx - incx
                           kk = kk - j
        40             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
        50                     continue
                           end if
                           kk = kk + (n-j+1)
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 70 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
        70                     continue
                           end if
                           jx = jx + incx
                           kk = kk + (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = 1
                   if (incx==1) then
                       do 110 j = 1,n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do 90 i = 1,j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
        90                     continue
                               if (nounit) temp = temp/ap(kk+j-1)
                           else
                               do 100 i = 1,j - 1
                                   temp = temp - dconjg(ap(k))*x(i)
                                   k = k + 1
       100                     continue
                               if (nounit) temp = temp/dconjg(ap(kk+j-1))
                           end if
                           x(j) = temp
                           kk = kk + j
       110             continue
                   else
                       jx = kx
                       do 140 j = 1,n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do 120 k = kk,kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
       120                     continue
                               if (nounit) temp = temp/ap(kk+j-1)
                           else
                               do 130 k = kk,kk + j - 2
                                   temp = temp - dconjg(ap(k))*x(ix)
                                   ix = ix + incx
       130                     continue
                               if (nounit) temp = temp/dconjg(ap(kk+j-1))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
       140             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 170 j = n,1,-1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do 150 i = n,j + 1,-1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
       150                     continue
                               if (nounit) temp = temp/ap(kk-n+j)
                           else
                               do 160 i = n,j + 1,-1
                                   temp = temp - dconjg(ap(k))*x(i)
                                   k = k - 1
       160                     continue
                               if (nounit) temp = temp/dconjg(ap(kk-n+j))
                           end if
                           x(j) = temp
                           kk = kk - (n-j+1)
       170             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 200 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do 180 k = kk,kk - (n- (j+1)),-1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
       180                     continue
                               if (nounit) temp = temp/ap(kk-n+j)
                           else
                               do 190 k = kk,kk - (n- (j+1)),-1
                                   temp = temp - dconjg(ap(k))*x(ix)
                                   ix = ix - incx
       190                     continue
                               if (nounit) temp = temp/dconjg(ap(kk-n+j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n-j+1)
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztpsv
     
     end subroutine stdlib_ztpsv

     
     
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('ztrmm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*a*b.
     
                   if (upper) then
                       do 50 j = 1,n
                           do 40 k = 1,m
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   do 30 i = 1,k - 1
                                       b(i,j) = b(i,j) + temp*a(i,k)
        30                         continue
                                   if (nounit) temp = temp*a(k,k)
                                   b(k,j) = temp
                               end if
        40                 continue
        50             continue
                   else
                       do 80 j = 1,n
                           do 70 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   b(k,j) = temp
                                   if (nounit) b(k,j) = b(k,j)*a(k,k)
                                   do 60 i = k + 1,m
                                       b(i,j) = b(i,j) + temp*a(i,k)
        60                         continue
                               end if
        70                 continue
        80             continue
                   end if
               else
     
                 ! form  b := alpha*a**t*b   or   b := alpha*a**h*b.
     
                   if (upper) then
                       do 120 j = 1,n
                           do 110 i = m,1,-1
                               temp = b(i,j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i,i)
                                   do 90 k = 1,i - 1
                                       temp = temp + a(k,i)*b(k,j)
        90                         continue
                               else
                                   if (nounit) temp = temp*dconjg(a(i,i))
                                   do 100 k = 1,i - 1
                                       temp = temp + dconjg(a(k,i))*b(k,j)
       100                         continue
                               end if
                               b(i,j) = alpha*temp
       110                 continue
       120             continue
                   else
                       do 160 j = 1,n
                           do 150 i = 1,m
                               temp = b(i,j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i,i)
                                   do 130 k = i + 1,m
                                       temp = temp + a(k,i)*b(k,j)
       130                         continue
                               else
                                   if (nounit) temp = temp*dconjg(a(i,i))
                                   do 140 k = i + 1,m
                                       temp = temp + dconjg(a(k,i))*b(k,j)
       140                         continue
                               end if
                               b(i,j) = alpha*temp
       150                 continue
       160             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*a.
     
                   if (upper) then
                       do 200 j = n,1,-1
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 170 i = 1,m
                               b(i,j) = temp*b(i,j)
       170                 continue
                           do 190 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 180 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       180                         continue
                               end if
       190                 continue
       200             continue
                   else
                       do 240 j = 1,n
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 210 i = 1,m
                               b(i,j) = temp*b(i,j)
       210                 continue
                           do 230 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 220 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       220                         continue
                               end if
       230                 continue
       240             continue
                   end if
               else
     
                 ! form  b := alpha*b*a**t   or   b := alpha*b*a**h.
     
                   if (upper) then
                       do 280 k = 1,n
                           do 260 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = alpha*a(j,k)
                                   else
                                       temp = alpha*dconjg(a(j,k))
                                   end if
                                   do 250 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       250                         continue
                               end if
       260                 continue
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k,k)
                               else
                                   temp = temp*dconjg(a(k,k))
                               end if
                           end if
                           if (temp/=one) then
                               do 270 i = 1,m
                                   b(i,k) = temp*b(i,k)
       270                     continue
                           end if
       280             continue
                   else
                       do 320 k = n,1,-1
                           do 300 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = alpha*a(j,k)
                                   else
                                       temp = alpha*dconjg(a(j,k))
                                   end if
                                   do 290 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       290                         continue
                               end if
       300                 continue
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k,k)
                               else
                                   temp = temp*dconjg(a(k,k))
                               end if
                           end if
                           if (temp/=one) then
                               do 310 i = 1,m
                                   b(i,k) = temp*b(i,k)
       310                     continue
                           end if
       320             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrmm
     
     end subroutine stdlib_ztrmm

     
     
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
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ztrmv ',info)
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 110 j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 90 i = j - 1,1,-1
                                   temp = temp + a(i,j)*x(i)
        90                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               do 100 i = j - 1,1,-1
                                   temp = temp + dconjg(a(i,j))*x(i)
       100                     continue
                           end if
                           x(j) = temp
       110             continue
                   else
                       jx = kx + (n-1)*incx
                       do 140 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 120 i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + a(i,j)*x(ix)
       120                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               do 130 i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + dconjg(a(i,j))*x(ix)
       130                     continue
                           end if
                           x(jx) = temp
                           jx = jx - incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = 1,n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 150 i = j + 1,n
                                   temp = temp + a(i,j)*x(i)
       150                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               do 160 i = j + 1,n
                                   temp = temp + dconjg(a(i,j))*x(i)
       160                     continue
                           end if
                           x(j) = temp
       170             continue
                   else
                       jx = kx
                       do 200 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 180 i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + a(i,j)*x(ix)
       180                     continue
                           else
                               if (nounit) temp = temp*dconjg(a(j,j))
                               do 190 i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + dconjg(a(i,j))*x(ix)
       190                     continue
                           end if
                           x(jx) = temp
                           jx = jx + incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrmv
     
     end subroutine stdlib_ztrmv

     
     
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('ztrsm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*inv( a )*b.
     
                   if (upper) then
                       do 60 j = 1,n
                           if (alpha/=one) then
                               do 30 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        30                     continue
                           end if
                           do 50 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 40 i = 1,k - 1
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        40                         continue
                               end if
        50                 continue
        60             continue
                   else
                       do 100 j = 1,n
                           if (alpha/=one) then
                               do 70 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        70                     continue
                           end if
                           do 90 k = 1,m
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 80 i = k + 1,m
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        80                         continue
                               end if
        90                 continue
       100             continue
                   end if
               else
     
                 ! form  b := alpha*inv( a**t )*b
                 ! or    b := alpha*inv( a**h )*b.
     
                   if (upper) then
                       do 140 j = 1,n
                           do 130 i = 1,m
                               temp = alpha*b(i,j)
                               if (noconj) then
                                   do 110 k = 1,i - 1
                                       temp = temp - a(k,i)*b(k,j)
       110                         continue
                                   if (nounit) temp = temp/a(i,i)
                               else
                                   do 120 k = 1,i - 1
                                       temp = temp - dconjg(a(k,i))*b(k,j)
       120                         continue
                                   if (nounit) temp = temp/dconjg(a(i,i))
                               end if
                               b(i,j) = temp
       130                 continue
       140             continue
                   else
                       do 180 j = 1,n
                           do 170 i = m,1,-1
                               temp = alpha*b(i,j)
                               if (noconj) then
                                   do 150 k = i + 1,m
                                       temp = temp - a(k,i)*b(k,j)
       150                         continue
                                   if (nounit) temp = temp/a(i,i)
                               else
                                   do 160 k = i + 1,m
                                       temp = temp - dconjg(a(k,i))*b(k,j)
       160                         continue
                                   if (nounit) temp = temp/dconjg(a(i,i))
                               end if
                               b(i,j) = temp
       170                 continue
       180             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*inv( a ).
     
                   if (upper) then
                       do 230 j = 1,n
                           if (alpha/=one) then
                               do 190 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       190                     continue
                           end if
                           do 210 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   do 200 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       200                         continue
                               end if
       210                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 220 i = 1,m
                                   b(i,j) = temp*b(i,j)
       220                     continue
                           end if
       230             continue
                   else
                       do 280 j = n,1,-1
                           if (alpha/=one) then
                               do 240 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       240                     continue
                           end if
                           do 260 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   do 250 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       250                         continue
                               end if
       260                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 270 i = 1,m
                                   b(i,j) = temp*b(i,j)
       270                     continue
                           end if
       280             continue
                   end if
               else
     
                 ! form  b := alpha*b*inv( a**t )
                 ! or    b := alpha*b*inv( a**h ).
     
                   if (upper) then
                       do 330 k = n,1,-1
                           if (nounit) then
                               if (noconj) then
                                   temp = one/a(k,k)
                               else
                                   temp = one/dconjg(a(k,k))
                               end if
                               do 290 i = 1,m
                                   b(i,k) = temp*b(i,k)
       290                     continue
                           end if
                           do 310 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = a(j,k)
                                   else
                                       temp = dconjg(a(j,k))
                                   end if
                                   do 300 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       300                         continue
                               end if
       310                 continue
                           if (alpha/=one) then
                               do 320 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       320                     continue
                           end if
       330             continue
                   else
                       do 380 k = 1,n
                           if (nounit) then
                               if (noconj) then
                                   temp = one/a(k,k)
                               else
                                   temp = one/dconjg(a(k,k))
                               end if
                               do 340 i = 1,m
                                   b(i,k) = temp*b(i,k)
       340                     continue
                           end if
                           do 360 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = a(j,k)
                                   else
                                       temp = dconjg(a(j,k))
                                   end if
                                   do 350 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       350                         continue
                               end if
       360                 continue
                           if (alpha/=one) then
                               do 370 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       370                     continue
                           end if
       380             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrsm
     
     end subroutine stdlib_ztrsm

     
     
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
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ztrsv ',info)
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
        10                     continue
                           end if
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 30 i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 70 i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 110 j = 1,n
                           temp = x(j)
                           if (noconj) then
                               do 90 i = 1,j - 1
                                   temp = temp - a(i,j)*x(i)
        90                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 100 i = 1,j - 1
                                   temp = temp - dconjg(a(i,j))*x(i)
       100                     continue
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(j) = temp
       110             continue
                   else
                       jx = kx
                       do 140 j = 1,n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do 120 i = 1,j - 1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix + incx
       120                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 130 i = 1,j - 1
                                   temp = temp - dconjg(a(i,j))*x(ix)
                                   ix = ix + incx
       130                     continue
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               do 150 i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(i)
       150                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 160 i = n,j + 1,-1
                                   temp = temp - dconjg(a(i,j))*x(i)
       160                     continue
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(j) = temp
       170             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 200 j = n,1,-1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do 180 i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix - incx
       180                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 190 i = n,j + 1,-1
                                   temp = temp - dconjg(a(i,j))*x(ix)
                                   ix = ix - incx
       190                     continue
                               if (nounit) temp = temp/dconjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ztrsv
     
     end subroutine stdlib_ztrsv

     
     
     subroutine stdlib_caxpy(n,ca,cx,incx,cy,incy)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) ca
           integer(int32) incx,incy,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*),cy(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           integer(int32) i,ix,iy
           ! ..
     
     
           if (n<=0) return
           if (stdlib_scabs1(ca)==0.0_sp) return
           if (incx==1 .and. incy==1) then
     
              ! code for both increments equal to 1
     
              do i = 1,n
                 cy(i) = cy(i) + ca*cx(i)
              end do
           else
     
              ! code for unequal increments or equal increments
                ! not equal to 1
     
              ix = 1
              iy = 1
              if (incx<0) ix = (-n+1)*incx + 1
              if (incy<0) iy = (-n+1)*incy + 1
              do i = 1,n
                 cy(iy) = cy(iy) + ca*cx(ix)
                 ix = ix + incx
                 iy = iy + incy
              end do
           end if
     
           return
     
           ! end of stdlib_caxpy
     
     end subroutine stdlib_caxpy

     
     
     subroutine stdlib_cgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) incx,incy,kl,ku,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,k,kup1,kx,ky,lenx,leny
           logical(lk) noconj
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('cgbmv ',info)
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kup1 = ku + 1
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       k = kup1 - j
                       do 50 i = max(1,j-ku),min(m,j+kl)
                           y(i) = y(i) + temp*a(k+i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       k = kup1 - j
                       do 70 i = max(1,j-ku),min(m,j+kl)
                           y(iy) = y(iy) + temp*a(k+i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
                       if (j>ku) ky = ky + incy
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
     
               jy = ky
               if (incx==1) then
                   do 110 j = 1,n
                       temp = zero
                       k = kup1 - j
                       if (noconj) then
                           do 90 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + a(k+i,j)*x(i)
        90                 continue
                       else
                           do 100 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + conjg(a(k+i,j))*x(i)
       100                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       110         continue
               else
                   do 140 j = 1,n
                       temp = zero
                       ix = kx
                       k = kup1 - j
                       if (noconj) then
                           do 120 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + a(k+i,j)*x(ix)
                               ix = ix + incx
       120                 continue
                       else
                           do 130 i = max(1,j-ku),min(m,j+kl)
                               temp = temp + conjg(a(k+i,j))*x(ix)
                               ix = ix + incx
       130                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j>ku) kx = kx + incx
       140         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cgbmv
     
     end subroutine stdlib_cgbmv

     
     
     subroutine stdlib_cgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) k,lda,ldb,ldc,m,n
           character transa,transb
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,j,l,nrowa,nrowb
           logical(lk) conja,conjb,nota,notb
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
     
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! conjugated or transposed, set  conja and conjb  as true if  a  and
           ! b  respectively are to be  transposed but  not conjugated  and set
           ! nrowa and  nrowb  as the number of rows of  a  and  b  respectively.
     
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
           else if ((.not.notb) .and. (.not.conjb) .and.(.not.stdlib_lsame(transb,'t')))&
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
               call stdlib_xerbla('cgemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.(((alpha==zero).or. (k==0)).and. (beta==one)))&
           return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (notb) then
               if (nota) then
     
                 ! form  c := alpha*a*b + beta*c.
     
                   do 90 j = 1,n
                       if (beta==zero) then
                           do 50 i = 1,m
                               c(i,j) = zero
        50                 continue
                       else if (beta/=one) then
                           do 60 i = 1,m
                               c(i,j) = beta*c(i,j)
        60                 continue
                       end if
                       do 80 l = 1,k
                           temp = alpha*b(l,j)
                           do 70 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
        70                 continue
        80             continue
        90         continue
               else if (conja) then
     
                 ! form  c := alpha*a**h*b + beta*c.
     
                   do 120 j = 1,n
                       do 110 i = 1,m
                           temp = zero
                           do 100 l = 1,k
                               temp = temp + conjg(a(l,i))*b(l,j)
       100                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       110             continue
       120         continue
               else
     
                 ! form  c := alpha*a**t*b + beta*c
     
                   do 150 j = 1,n
                       do 140 i = 1,m
                           temp = zero
                           do 130 l = 1,k
                               temp = temp + a(l,i)*b(l,j)
       130                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       140             continue
       150         continue
               end if
           else if (nota) then
               if (conjb) then
     
                 ! form  c := alpha*a*b**h + beta*c.
     
                   do 200 j = 1,n
                       if (beta==zero) then
                           do 160 i = 1,m
                               c(i,j) = zero
       160                 continue
                       else if (beta/=one) then
                           do 170 i = 1,m
                               c(i,j) = beta*c(i,j)
       170                 continue
                       end if
                       do 190 l = 1,k
                           temp = alpha*conjg(b(j,l))
                           do 180 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
       180                 continue
       190             continue
       200         continue
               else
     
                 ! form  c := alpha*a*b**t + beta*c
     
                   do 250 j = 1,n
                       if (beta==zero) then
                           do 210 i = 1,m
                               c(i,j) = zero
       210                 continue
                       else if (beta/=one) then
                           do 220 i = 1,m
                               c(i,j) = beta*c(i,j)
       220                 continue
                       end if
                       do 240 l = 1,k
                           temp = alpha*b(j,l)
                           do 230 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
       230                 continue
       240             continue
       250         continue
               end if
           else if (conja) then
               if (conjb) then
     
                 ! form  c := alpha*a**h*b**h + beta*c.
     
                   do 280 j = 1,n
                       do 270 i = 1,m
                           temp = zero
                           do 260 l = 1,k
                               temp = temp + conjg(a(l,i))*conjg(b(j,l))
       260                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       270             continue
       280         continue
               else
     
                 ! form  c := alpha*a**h*b**t + beta*c
     
                   do 310 j = 1,n
                       do 300 i = 1,m
                           temp = zero
                           do 290 l = 1,k
                               temp = temp + conjg(a(l,i))*b(j,l)
       290                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       300             continue
       310         continue
               end if
           else
               if (conjb) then
     
                 ! form  c := alpha*a**t*b**h + beta*c
     
                   do 340 j = 1,n
                       do 330 i = 1,m
                           temp = zero
                           do 320 l = 1,k
                               temp = temp + a(l,i)*conjg(b(j,l))
       320                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       330             continue
       340         continue
               else
     
                 ! form  c := alpha*a**t*b**t + beta*c
     
                   do 370 j = 1,n
                       do 360 i = 1,m
                           temp = zero
                           do 350 l = 1,k
                               temp = temp + a(l,i)*b(j,l)
       350                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       360             continue
       370         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cgemm
     
     end subroutine stdlib_cgemm

     
     
     subroutine stdlib_cgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) incx,incy,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
           logical(lk) noconj
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('cgemv ',info)
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       do 50 i = 1,m
                           y(i) = y(i) + temp*a(i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       do 70 i = 1,m
                           y(iy) = y(iy) + temp*a(i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y  or  y := alpha*a**h*x + y.
     
               jy = ky
               if (incx==1) then
                   do 110 j = 1,n
                       temp = zero
                       if (noconj) then
                           do 90 i = 1,m
                               temp = temp + a(i,j)*x(i)
        90                 continue
                       else
                           do 100 i = 1,m
                               temp = temp + conjg(a(i,j))*x(i)
       100                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       110         continue
               else
                   do 140 j = 1,n
                       temp = zero
                       ix = kx
                       if (noconj) then
                           do 120 i = 1,m
                               temp = temp + a(i,j)*x(ix)
                               ix = ix + incx
       120                 continue
                       else
                           do 130 i = 1,m
                               temp = temp + conjg(a(i,j))*x(ix)
                               ix = ix + incx
       130                 continue
                       end if
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       140         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cgemv
     
     end subroutine stdlib_cgemv

     
     
     subroutine stdlib_cgerc(m,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           integer(int32) incx,incy,lda,m,n
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jy,kx
           ! ..
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
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
               call stdlib_xerbla('cgerc ',info)
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
               do 20 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*conjg(y(jy))
                       do 10 i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
        10             continue
                   end if
                   jy = jy + incy
        20     continue
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               do 40 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*conjg(y(jy))
                       ix = kx
                       do 30 i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
        30             continue
                   end if
                   jy = jy + incy
        40     continue
           end if
     
           return
     
           ! end of stdlib_cgerc
     
     end subroutine stdlib_cgerc

     
     
     subroutine stdlib_cgeru(m,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           integer(int32) incx,incy,lda,m,n
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
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
               call stdlib_xerbla('cgeru ',info)
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
               do 20 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       do 10 i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
        10             continue
                   end if
                   jy = jy + incy
        20     continue
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               do 40 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       ix = kx
                       do 30 i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
        30             continue
                   end if
                   jy = jy + incy
        40     continue
           end if
     
           return
     
           ! end of stdlib_cgeru
     
     end subroutine stdlib_cgeru

     
     
     subroutine stdlib_chbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) incx,incy,k,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,min,real
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
               call stdlib_xerbla('chbmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when upper triangle of a is stored.
     
               kplus1 = k + 1
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       l = kplus1 - j
                       do 50 i = max(1,j-k),j - 1
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + conjg(a(l+i,j))*x(i)
        50             continue
                       y(j) = y(j) + temp1*real(a(kplus1,j)) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       do 70 i = max(1,j-k),j - 1
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + conjg(a(l+i,j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*real(a(kplus1,j)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       if (j>k) then
                           kx = kx + incx
                           ky = ky + incy
                       end if
        80         continue
               end if
           else
     
              ! form  y  when lower triangle of a is stored.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*real(a(1,j))
                       l = 1 - j
                       do 90 i = j + 1,min(n,j+k)
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + conjg(a(l+i,j))*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*real(a(1,j))
                       l = 1 - j
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,min(n,j+k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + conjg(a(l+i,j))*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_chbmv
     
     end subroutine stdlib_chbmv

     
     
     subroutine stdlib_chemm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) lda,ldb,ldc,m,n
           character side,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,real
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,j,k,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
               call stdlib_xerbla('chemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   do 70 j = 1,n
                       do 60 i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 50 k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
        50                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*real(a(i,i)) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i)) +alpha*temp2
                           end if
        60             continue
        70         continue
               else
                   do 100 j = 1,n
                       do 90 i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 80 k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*conjg(a(k,i))
        80                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*real(a(i,i)) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*real(a(i,i)) +alpha*temp2
                           end if
        90             continue
       100         continue
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               do 170 j = 1,n
                   temp1 = alpha*real(a(j,j))
                   if (beta==zero) then
                       do 110 i = 1,m
                           c(i,j) = temp1*b(i,j)
       110             continue
                   else
                       do 120 i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
       120             continue
                   end if
                   do 140 k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*conjg(a(j,k))
                       end if
                       do 130 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       130             continue
       140         continue
                   do 160 k = j + 1,n
                       if (upper) then
                           temp1 = alpha*conjg(a(j,k))
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do 150 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       150             continue
       160         continue
       170     continue
           end if
     
           return
     
           ! end of stdlib_chemm
     
     end subroutine stdlib_chemm

     
     
     subroutine stdlib_chemv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,real
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
               call stdlib_xerbla('chemv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when a is stored in upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + conjg(a(i,j))*x(i)
        50             continue
                       y(j) = y(j) + temp1*real(a(j,j)) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 i = 1,j - 1
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + conjg(a(i,j))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*real(a(j,j)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           else
     
              ! form  y  when a is stored in lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*real(a(j,j))
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + conjg(a(i,j))*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*real(a(j,j))
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + conjg(a(i,j))*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_chemv
     
     end subroutine stdlib_chemv

     
     
     subroutine stdlib_cher(uplo,n,alpha,x,incx,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,kx
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,real
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
               call stdlib_xerbla('cher  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==real(zero))) return
     
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*conjg(x(j))
                           do 10 i = 1,j - 1
                               a(i,j) = a(i,j) + x(i)*temp
        10                 continue
                           a(j,j) = real(a(j,j)) + real(x(j)*temp)
                       else
                           a(j,j) = real(a(j,j))
                       end if
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*conjg(x(jx))
                           ix = kx
                           do 30 i = 1,j - 1
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                           a(j,j) = real(a(j,j)) + real(x(jx)*temp)
                       else
                           a(j,j) = real(a(j,j))
                       end if
                       jx = jx + incx
        40         continue
               end if
           else
     
              ! form  a  when a is stored in lower triangle.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*conjg(x(j))
                           a(j,j) = real(a(j,j)) + real(temp*x(j))
                           do 50 i = j + 1,n
                               a(i,j) = a(i,j) + x(i)*temp
        50                 continue
                       else
                           a(j,j) = real(a(j,j))
                       end if
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*conjg(x(jx))
                           a(j,j) = real(a(j,j)) + real(temp*x(jx))
                           ix = jx
                           do 70 i = j + 1,n
                               ix = ix + incx
                               a(i,j) = a(i,j) + x(ix)*temp
        70                 continue
                       else
                           a(j,j) = real(a(j,j))
                       end if
                       jx = jx + incx
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cher
     
     end subroutine stdlib_cher

     
     
     subroutine stdlib_cher2(uplo,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,real
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
               call stdlib_xerbla('cher2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           do 10 i = 1,j - 1
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        10                 continue
                           a(j,j) = real(a(j,j)) +real(x(j)*temp1+y(j)*temp2)
                       else
                           a(j,j) = real(a(j,j))
                       end if
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           do 30 i = 1,j - 1
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                           a(j,j) = real(a(j,j)) +real(x(jx)*temp1+y(jy)*temp2)
                       else
                           a(j,j) = real(a(j,j))
                       end if
                       jx = jx + incx
                       jy = jy + incy
        40         continue
               end if
           else
     
              ! form  a  when a is stored in the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           a(j,j) = real(a(j,j)) +real(x(j)*temp1+y(j)*temp2)
                           do 50 i = j + 1,n
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        50                 continue
                       else
                           a(j,j) = real(a(j,j))
                       end if
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           a(j,j) = real(a(j,j)) +real(x(jx)*temp1+y(jy)*temp2)
                           ix = jx
                           iy = jy
                           do 70 i = j + 1,n
                               ix = ix + incx
                               iy = iy + incy
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
        70                 continue
                       else
                           a(j,j) = real(a(j,j))
                       end if
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cher2
     
     end subroutine stdlib_cher2

     
     
     subroutine stdlib_cher2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           real(sp) beta
           integer(int32) k,lda,ldb,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,real
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(sp) one
           parameter (one=1.0_sp)
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'c')))&
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
               call stdlib_xerbla('cher2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==real(zero)) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
        30                 continue
                           c(j,j) = beta*real(c(j,j))
        40             continue
                   end if
               else
                   if (beta==real(zero)) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           c(j,j) = beta*real(c(j,j))
                           do 70 i = j + 1,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**h + conjg( alpha )*b*a**h +
                         ! c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==real(zero)) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
       100                 continue
                           c(j,j) = beta*real(c(j,j))
                       else
                           c(j,j) = real(c(j,j))
                       end if
                       do 120 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do 110 i = 1,j - 1
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       110                     continue
                               c(j,j) = real(c(j,j)) +real(a(j,l)*temp1+b(j,l)*temp2)
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==real(zero)) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j + 1,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                           c(j,j) = beta*real(c(j,j))
                       else
                           c(j,j) = real(c(j,j))
                       end if
                       do 170 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*conjg(b(j,l))
                               temp2 = conjg(alpha*a(j,l))
                               do 160 i = j + 1,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       160                     continue
                               c(j,j) = real(c(j,j)) +real(a(j,l)*temp1+b(j,l)*temp2)
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**h*b + conjg( alpha )*b**h*a +
                         ! c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do 190 l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
       190                 continue
                           if (i==j) then
                               if (beta==real(zero)) then
                                   c(j,j) = real(alpha*temp1+conjg(alpha)*temp2)
                               else
                                   c(j,j) = beta*real(c(j,j)) +real(alpha*temp1+conjg(alpha)&
          *temp2)
                               end if
                           else
                               if (beta==real(zero)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 +conjg(alpha)*temp2
                               end if
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp1 = zero
                           temp2 = zero
                           do 220 l = 1,k
                               temp1 = temp1 + conjg(a(l,i))*b(l,j)
                               temp2 = temp2 + conjg(b(l,i))*a(l,j)
       220                 continue
                           if (i==j) then
                               if (beta==real(zero)) then
                                   c(j,j) = real(alpha*temp1+conjg(alpha)*temp2)
                               else
                                   c(j,j) = beta*real(c(j,j)) +real(alpha*temp1+conjg(alpha)&
          *temp2)
                               end if
                           else
                               if (beta==real(zero)) then
                                   c(i,j) = alpha*temp1 + conjg(alpha)*temp2
                               else
                                   c(i,j) = beta*c(i,j) + alpha*temp1 +conjg(alpha)*temp2
                               end if
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cher2k
     
     end subroutine stdlib_cher2k

     
     
     subroutine stdlib_cherk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) k,lda,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic cmplx,conjg,max,real
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           real(sp) rtemp
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'c')))&
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
               call stdlib_xerbla('cherk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
        30                 continue
                           c(j,j) = beta*real(c(j,j))
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           c(j,j) = beta*real(c(j,j))
                           do 70 i = j + 1,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**h + beta*c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j - 1
                               c(i,j) = beta*c(i,j)
       100                 continue
                           c(j,j) = beta*real(c(j,j))
                       else
                           c(j,j) = real(c(j,j))
                       end if
                       do 120 l = 1,k
                           if (a(j,l)/=cmplx(zero)) then
                               temp = alpha*conjg(a(j,l))
                               do 110 i = 1,j - 1
                                   c(i,j) = c(i,j) + temp*a(i,l)
       110                     continue
                               c(j,j) = real(c(j,j)) + real(temp*a(i,l))
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           c(j,j) = beta*real(c(j,j))
                           do 150 i = j + 1,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       else
                           c(j,j) = real(c(j,j))
                       end if
                       do 170 l = 1,k
                           if (a(j,l)/=cmplx(zero)) then
                               temp = alpha*conjg(a(j,l))
                               c(j,j) = real(c(j,j)) + real(temp*a(j,l))
                               do 160 i = j + 1,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**h*a + beta*c.
     
               if (upper) then
                   do 220 j = 1,n
                       do 200 i = 1,j - 1
                           temp = zero
                           do 190 l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       200             continue
                       rtemp = zero
                       do 210 l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
       210             continue
                       if (beta==zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j))
                       end if
       220         continue
               else
                   do 260 j = 1,n
                       rtemp = zero
                       do 230 l = 1,k
                           rtemp = rtemp + conjg(a(l,j))*a(l,j)
       230             continue
                       if (beta==zero) then
                           c(j,j) = alpha*rtemp
                       else
                           c(j,j) = alpha*rtemp + beta*real(c(j,j))
                       end if
                       do 250 i = j + 1,n
                           temp = zero
                           do 240 l = 1,k
                               temp = temp + conjg(a(l,i))*a(l,j)
       240                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       250             continue
       260         continue
               end if
           end if
     
           return
     
           ! end of stdlib_cherk
     
     end subroutine stdlib_cherk

     
     
     subroutine stdlib_chpmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,real
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
               call stdlib_xerbla('chpmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when ap contains the upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       k = kk
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(i)
                           k = k + 1
        50             continue
                       y(j) = y(j) + temp1*real(ap(kk+j-1)) + alpha*temp2
                       kk = kk + j
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 k = kk,kk + j - 2
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*real(ap(kk+j-1)) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        80         continue
               end if
           else
     
              ! form  y  when ap contains the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*real(ap(kk))
                       k = kk + 1
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(i)
                           k = k + 1
        90             continue
                       y(j) = y(j) + alpha*temp2
                       kk = kk + (n-j+1)
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*real(ap(kk))
                       ix = jx
                       iy = jy
                       do 110 k = kk + 1,kk + n - j
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + conjg(ap(k))*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + (n-j+1)
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_chpmv
     
     end subroutine stdlib_chpmv

     
     
     subroutine stdlib_chpr(uplo,n,alpha,x,incx,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,real
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
               call stdlib_xerbla('chpr  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==real(zero))) return
     
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*conjg(x(j))
                           k = kk
                           do 10 i = 1,j - 1
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        10                 continue
                           ap(kk+j-1) = real(ap(kk+j-1)) + real(x(j)*temp)
                       else
                           ap(kk+j-1) = real(ap(kk+j-1))
                       end if
                       kk = kk + j
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*conjg(x(jx))
                           ix = kx
                           do 30 k = kk,kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                           ap(kk+j-1) = real(ap(kk+j-1)) + real(x(jx)*temp)
                       else
                           ap(kk+j-1) = real(ap(kk+j-1))
                       end if
                       jx = jx + incx
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*conjg(x(j))
                           ap(kk) = real(ap(kk)) + real(temp*x(j))
                           k = kk + 1
                           do 50 i = j + 1,n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        50                 continue
                       else
                           ap(kk) = real(ap(kk))
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*conjg(x(jx))
                           ap(kk) = real(ap(kk)) + real(temp*x(jx))
                           ix = jx
                           do 70 k = kk + 1,kk + n - j
                               ix = ix + incx
                               ap(k) = ap(k) + x(ix)*temp
        70                 continue
                       else
                           ap(kk) = real(ap(kk))
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_chpr
     
     end subroutine stdlib_chpr

     
     
     subroutine stdlib_chpr2(uplo,n,alpha,x,incx,y,incy,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           complex(sp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,real
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
               call stdlib_xerbla('chpr2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           k = kk
                           do 10 i = 1,j - 1
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        10                 continue
                           ap(kk+j-1) = real(ap(kk+j-1)) +real(x(j)*temp1+y(j)*temp2)
                       else
                           ap(kk+j-1) = real(ap(kk+j-1))
                       end if
                       kk = kk + j
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           ix = kx
                           iy = ky
                           do 30 k = kk,kk + j - 2
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                           ap(kk+j-1) = real(ap(kk+j-1)) +real(x(jx)*temp1+y(jy)*temp2)
                       else
                           ap(kk+j-1) = real(ap(kk+j-1))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*conjg(y(j))
                           temp2 = conjg(alpha*x(j))
                           ap(kk) = real(ap(kk)) +real(x(j)*temp1+y(j)*temp2)
                           k = kk + 1
                           do 50 i = j + 1,n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        50                 continue
                       else
                           ap(kk) = real(ap(kk))
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*conjg(y(jy))
                           temp2 = conjg(alpha*x(jx))
                           ap(kk) = real(ap(kk)) +real(x(jx)*temp1+y(jy)*temp2)
                           ix = jx
                           iy = jy
                           do 70 k = kk + 1,kk + n - j
                               ix = ix + incx
                               iy = iy + incy
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
        70                 continue
                       else
                           ap(kk) = real(ap(kk))
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_chpr2
     
     end subroutine stdlib_chpr2

     
     
     subroutine stdlib_csymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) lda,ldb,ldc,m,n
           character side,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,j,k,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
               call stdlib_xerbla('csymm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   do 70 j = 1,n
                       do 60 i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 50 k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        50                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        60             continue
        70         continue
               else
                   do 100 j = 1,n
                       do 90 i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 80 k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        80                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        90             continue
       100         continue
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               do 170 j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta==zero) then
                       do 110 i = 1,m
                           c(i,j) = temp1*b(i,j)
       110             continue
                   else
                       do 120 i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
       120             continue
                   end if
                   do 140 k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do 130 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       130             continue
       140         continue
                   do 160 k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do 150 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       150             continue
       160         continue
       170     continue
           end if
     
           return
     
           ! end of stdlib_csymm
     
     end subroutine stdlib_csymm

     
     
     subroutine stdlib_csyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) k,lda,ldb,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           complex(sp) temp1,temp2
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')))&
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
               call stdlib_xerbla('csyr2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do 190 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp1 = zero
                           temp2 = zero
                           do 220 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_csyr2k
     
     end subroutine stdlib_csyr2k

     
     
     subroutine stdlib_csyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha,beta
           integer(int32) k,lda,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')))&
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
               call stdlib_xerbla('csyrk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**t + beta*c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*a + beta*c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp = zero
                           do 190 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp = zero
                           do 220 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_csyrk
     
     end subroutine stdlib_csyrk

     
     
     subroutine stdlib_ctbmv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ctbmv ',info)
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do 10 i = max(1,j-k),j - 1
                                   x(i) = x(i) + temp*a(l+i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do 30 i = max(1,j-k),j - 1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j>k) kx = kx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = 1 - j
                               do 50 i = min(n,j+k),j + 1,-1
                                   x(i) = x(i) + temp*a(l+i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do 70 i = min(n,j+k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 110 j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do 90 i = j - 1,max(1,j-k),-1
                                   temp = temp + a(l+i,j)*x(i)
        90                     continue
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do 100 i = j - 1,max(1,j-k),-1
                                   temp = temp + conjg(a(l+i,j))*x(i)
       100                     continue
                           end if
                           x(j) = temp
       110             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 140 j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(kplus1,j)
                               do 120 i = j - 1,max(1,j-k),-1
                                   temp = temp + a(l+i,j)*x(ix)
                                   ix = ix - incx
       120                     continue
                           else
                               if (nounit) temp = temp*conjg(a(kplus1,j))
                               do 130 i = j - 1,max(1,j-k),-1
                                   temp = temp + conjg(a(l+i,j))*x(ix)
                                   ix = ix - incx
       130                     continue
                           end if
                           x(jx) = temp
                           jx = jx - incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do 150 i = j + 1,min(n,j+k)
                                   temp = temp + a(l+i,j)*x(i)
       150                     continue
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do 160 i = j + 1,min(n,j+k)
                                   temp = temp + conjg(a(l+i,j))*x(i)
       160                     continue
                           end if
                           x(j) = temp
       170             continue
                   else
                       jx = kx
                       do 200 j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               if (nounit) temp = temp*a(1,j)
                               do 180 i = j + 1,min(n,j+k)
                                   temp = temp + a(l+i,j)*x(ix)
                                   ix = ix + incx
       180                     continue
                           else
                               if (nounit) temp = temp*conjg(a(1,j))
                               do 190 i = j + 1,min(n,j+k)
                                   temp = temp + conjg(a(l+i,j))*x(ix)
                                   ix = ix + incx
       190                     continue
                           end if
                           x(jx) = temp
                           jx = jx + incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctbmv
     
     end subroutine stdlib_ctbmv

     
     
     subroutine stdlib_ctbsv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ctbsv ',info)
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do 10 i = j - 1,max(1,j-k),-1
                                   x(i) = x(i) - temp*a(l+i,j)
        10                     continue
                           end if
        20             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 40 j = n,1,-1
                           kx = kx - incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do 30 i = j - 1,max(1,j-k),-1
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix - incx
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do 50 i = j + 1,min(n,j+k)
                                   x(i) = x(i) - temp*a(l+i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           kx = kx + incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do 70 i = j + 1,min(n,j+k)
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix + incx
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 110 j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           if (noconj) then
                               do 90 i = max(1,j-k),j - 1
                                   temp = temp - a(l+i,j)*x(i)
        90                     continue
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do 100 i = max(1,j-k),j - 1
                                   temp = temp - conjg(a(l+i,j))*x(i)
       100                     continue
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(j) = temp
       110             continue
                   else
                       jx = kx
                       do 140 j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           if (noconj) then
                               do 120 i = max(1,j-k),j - 1
                                   temp = temp - a(l+i,j)*x(ix)
                                   ix = ix + incx
       120                     continue
                               if (nounit) temp = temp/a(kplus1,j)
                           else
                               do 130 i = max(1,j-k),j - 1
                                   temp = temp - conjg(a(l+i,j))*x(ix)
                                   ix = ix + incx
       130                     continue
                               if (nounit) temp = temp/conjg(a(kplus1,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           if (j>k) kx = kx + incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           if (noconj) then
                               do 150 i = min(n,j+k),j + 1,-1
                                   temp = temp - a(l+i,j)*x(i)
       150                     continue
                               if (nounit) temp = temp/a(1,j)
                           else
                               do 160 i = min(n,j+k),j + 1,-1
                                   temp = temp - conjg(a(l+i,j))*x(i)
       160                     continue
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(j) = temp
       170             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 200 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           if (noconj) then
                               do 180 i = min(n,j+k),j + 1,-1
                                   temp = temp - a(l+i,j)*x(ix)
                                   ix = ix - incx
       180                     continue
                               if (nounit) temp = temp/a(1,j)
                           else
                               do 190 i = min(n,j+k),j + 1,-1
                                   temp = temp - conjg(a(l+i,j))*x(ix)
                                   ix = ix - incx
       190                     continue
                               if (nounit) temp = temp/conjg(a(1,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctbsv
     
     end subroutine stdlib_ctbsv

     
     
     subroutine stdlib_ctpmv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ctpmv ',info)
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
        10                     continue
                               if (nounit) x(j) = x(j)*ap(kk+j-1)
                           end if
                           kk = kk + j
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk+j-1)
                           end if
                           jx = jx + incx
                           kk = kk + j
        40             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
        50                     continue
                               if (nounit) x(j) = x(j)*ap(kk-n+j)
                           end if
                           kk = kk - (n-j+1)
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 k = kk,kk - (n- (j+1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk-n+j)
                           end if
                           jx = jx - incx
                           kk = kk - (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 110 j = n,1,-1
                           temp = x(j)
                           k = kk - 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 90 i = j - 1,1,-1
                                   temp = temp + ap(k)*x(i)
                                   k = k - 1
        90                     continue
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do 100 i = j - 1,1,-1
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k - 1
       100                     continue
                           end if
                           x(j) = temp
                           kk = kk - j
       110             continue
                   else
                       jx = kx + (n-1)*incx
                       do 140 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 120 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + ap(k)*x(ix)
       120                     continue
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do 130 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(ap(k))*x(ix)
       130                     continue
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
       140             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 170 j = 1,n
                           temp = x(j)
                           k = kk + 1
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 150 i = j + 1,n
                                   temp = temp + ap(k)*x(i)
                                   k = k + 1
       150                     continue
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do 160 i = j + 1,n
                                   temp = temp + conjg(ap(k))*x(i)
                                   k = k + 1
       160                     continue
                           end if
                           x(j) = temp
                           kk = kk + (n-j+1)
       170             continue
                   else
                       jx = kx
                       do 200 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*ap(kk)
                               do 180 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + ap(k)*x(ix)
       180                     continue
                           else
                               if (nounit) temp = temp*conjg(ap(kk))
                               do 190 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   temp = temp + conjg(ap(k))*x(ix)
       190                     continue
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n-j+1)
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctpmv
     
     end subroutine stdlib_ctpmv

     
     
     subroutine stdlib_ctpsv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ctpsv ',info)
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
        10                     continue
                           end if
                           kk = kk - j
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 30 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
        30                     continue
                           end if
                           jx = jx - incx
                           kk = kk - j
        40             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
        50                     continue
                           end if
                           kk = kk + (n-j+1)
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 70 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
        70                     continue
                           end if
                           jx = jx + incx
                           kk = kk + (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = 1
                   if (incx==1) then
                       do 110 j = 1,n
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do 90 i = 1,j - 1
                                   temp = temp - ap(k)*x(i)
                                   k = k + 1
        90                     continue
                               if (nounit) temp = temp/ap(kk+j-1)
                           else
                               do 100 i = 1,j - 1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k + 1
       100                     continue
                               if (nounit) temp = temp/conjg(ap(kk+j-1))
                           end if
                           x(j) = temp
                           kk = kk + j
       110             continue
                   else
                       jx = kx
                       do 140 j = 1,n
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do 120 k = kk,kk + j - 2
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix + incx
       120                     continue
                               if (nounit) temp = temp/ap(kk+j-1)
                           else
                               do 130 k = kk,kk + j - 2
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix + incx
       130                     continue
                               if (nounit) temp = temp/conjg(ap(kk+j-1))
                           end if
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
       140             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 170 j = n,1,-1
                           temp = x(j)
                           k = kk
                           if (noconj) then
                               do 150 i = n,j + 1,-1
                                   temp = temp - ap(k)*x(i)
                                   k = k - 1
       150                     continue
                               if (nounit) temp = temp/ap(kk-n+j)
                           else
                               do 160 i = n,j + 1,-1
                                   temp = temp - conjg(ap(k))*x(i)
                                   k = k - 1
       160                     continue
                               if (nounit) temp = temp/conjg(ap(kk-n+j))
                           end if
                           x(j) = temp
                           kk = kk - (n-j+1)
       170             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 200 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           if (noconj) then
                               do 180 k = kk,kk - (n- (j+1)),-1
                                   temp = temp - ap(k)*x(ix)
                                   ix = ix - incx
       180                     continue
                               if (nounit) temp = temp/ap(kk-n+j)
                           else
                               do 190 k = kk,kk - (n- (j+1)),-1
                                   temp = temp - conjg(ap(k))*x(ix)
                                   ix = ix - incx
       190                     continue
                               if (nounit) temp = temp/conjg(ap(kk-n+j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n-j+1)
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctpsv
     
     end subroutine stdlib_ctpsv

     
     
     subroutine stdlib_ctrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,noconj,nounit,upper
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('ctrmm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*a*b.
     
                   if (upper) then
                       do 50 j = 1,n
                           do 40 k = 1,m
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   do 30 i = 1,k - 1
                                       b(i,j) = b(i,j) + temp*a(i,k)
        30                         continue
                                   if (nounit) temp = temp*a(k,k)
                                   b(k,j) = temp
                               end if
        40                 continue
        50             continue
                   else
                       do 80 j = 1,n
                           do 70 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   b(k,j) = temp
                                   if (nounit) b(k,j) = b(k,j)*a(k,k)
                                   do 60 i = k + 1,m
                                       b(i,j) = b(i,j) + temp*a(i,k)
        60                         continue
                               end if
        70                 continue
        80             continue
                   end if
               else
     
                 ! form  b := alpha*a**t*b   or   b := alpha*a**h*b.
     
                   if (upper) then
                       do 120 j = 1,n
                           do 110 i = m,1,-1
                               temp = b(i,j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i,i)
                                   do 90 k = 1,i - 1
                                       temp = temp + a(k,i)*b(k,j)
        90                         continue
                               else
                                   if (nounit) temp = temp*conjg(a(i,i))
                                   do 100 k = 1,i - 1
                                       temp = temp + conjg(a(k,i))*b(k,j)
       100                         continue
                               end if
                               b(i,j) = alpha*temp
       110                 continue
       120             continue
                   else
                       do 160 j = 1,n
                           do 150 i = 1,m
                               temp = b(i,j)
                               if (noconj) then
                                   if (nounit) temp = temp*a(i,i)
                                   do 130 k = i + 1,m
                                       temp = temp + a(k,i)*b(k,j)
       130                         continue
                               else
                                   if (nounit) temp = temp*conjg(a(i,i))
                                   do 140 k = i + 1,m
                                       temp = temp + conjg(a(k,i))*b(k,j)
       140                         continue
                               end if
                               b(i,j) = alpha*temp
       150                 continue
       160             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*a.
     
                   if (upper) then
                       do 200 j = n,1,-1
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 170 i = 1,m
                               b(i,j) = temp*b(i,j)
       170                 continue
                           do 190 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 180 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       180                         continue
                               end if
       190                 continue
       200             continue
                   else
                       do 240 j = 1,n
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 210 i = 1,m
                               b(i,j) = temp*b(i,j)
       210                 continue
                           do 230 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 220 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       220                         continue
                               end if
       230                 continue
       240             continue
                   end if
               else
     
                 ! form  b := alpha*b*a**t   or   b := alpha*b*a**h.
     
                   if (upper) then
                       do 280 k = 1,n
                           do 260 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = alpha*a(j,k)
                                   else
                                       temp = alpha*conjg(a(j,k))
                                   end if
                                   do 250 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       250                         continue
                               end if
       260                 continue
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k,k)
                               else
                                   temp = temp*conjg(a(k,k))
                               end if
                           end if
                           if (temp/=one) then
                               do 270 i = 1,m
                                   b(i,k) = temp*b(i,k)
       270                     continue
                           end if
       280             continue
                   else
                       do 320 k = n,1,-1
                           do 300 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = alpha*a(j,k)
                                   else
                                       temp = alpha*conjg(a(j,k))
                                   end if
                                   do 290 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       290                         continue
                               end if
       300                 continue
                           temp = alpha
                           if (nounit) then
                               if (noconj) then
                                   temp = temp*a(k,k)
                               else
                                   temp = temp*conjg(a(k,k))
                               end if
                           end if
                           if (temp/=one) then
                               do 310 i = 1,m
                                   b(i,k) = temp*b(i,k)
       310                     continue
                           end if
       320             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctrmm
     
     end subroutine stdlib_ctrmm

     
     
     subroutine stdlib_ctrmv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ctrmv ',info)
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x  or  x := a**h*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 110 j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 90 i = j - 1,1,-1
                                   temp = temp + a(i,j)*x(i)
        90                     continue
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do 100 i = j - 1,1,-1
                                   temp = temp + conjg(a(i,j))*x(i)
       100                     continue
                           end if
                           x(j) = temp
       110             continue
                   else
                       jx = kx + (n-1)*incx
                       do 140 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 120 i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + a(i,j)*x(ix)
       120                     continue
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do 130 i = j - 1,1,-1
                                   ix = ix - incx
                                   temp = temp + conjg(a(i,j))*x(ix)
       130                     continue
                           end if
                           x(jx) = temp
                           jx = jx - incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = 1,n
                           temp = x(j)
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 150 i = j + 1,n
                                   temp = temp + a(i,j)*x(i)
       150                     continue
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do 160 i = j + 1,n
                                   temp = temp + conjg(a(i,j))*x(i)
       160                     continue
                           end if
                           x(j) = temp
       170             continue
                   else
                       jx = kx
                       do 200 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (noconj) then
                               if (nounit) temp = temp*a(j,j)
                               do 180 i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + a(i,j)*x(ix)
       180                     continue
                           else
                               if (nounit) temp = temp*conjg(a(j,j))
                               do 190 i = j + 1,n
                                   ix = ix + incx
                                   temp = temp + conjg(a(i,j))*x(ix)
       190                     continue
                           end if
                           x(jx) = temp
                           jx = jx + incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctrmv
     
     end subroutine stdlib_ctrmv

     
     
     subroutine stdlib_ctrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           complex(sp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,noconj,nounit,upper
           ! ..
           ! .. parameters ..
           complex(sp) one
           parameter (one= (1.0_sp,0.0_sp))
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('ctrsm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*inv( a )*b.
     
                   if (upper) then
                       do 60 j = 1,n
                           if (alpha/=one) then
                               do 30 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        30                     continue
                           end if
                           do 50 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 40 i = 1,k - 1
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        40                         continue
                               end if
        50                 continue
        60             continue
                   else
                       do 100 j = 1,n
                           if (alpha/=one) then
                               do 70 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        70                     continue
                           end if
                           do 90 k = 1,m
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 80 i = k + 1,m
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        80                         continue
                               end if
        90                 continue
       100             continue
                   end if
               else
     
                 ! form  b := alpha*inv( a**t )*b
                 ! or    b := alpha*inv( a**h )*b.
     
                   if (upper) then
                       do 140 j = 1,n
                           do 130 i = 1,m
                               temp = alpha*b(i,j)
                               if (noconj) then
                                   do 110 k = 1,i - 1
                                       temp = temp - a(k,i)*b(k,j)
       110                         continue
                                   if (nounit) temp = temp/a(i,i)
                               else
                                   do 120 k = 1,i - 1
                                       temp = temp - conjg(a(k,i))*b(k,j)
       120                         continue
                                   if (nounit) temp = temp/conjg(a(i,i))
                               end if
                               b(i,j) = temp
       130                 continue
       140             continue
                   else
                       do 180 j = 1,n
                           do 170 i = m,1,-1
                               temp = alpha*b(i,j)
                               if (noconj) then
                                   do 150 k = i + 1,m
                                       temp = temp - a(k,i)*b(k,j)
       150                         continue
                                   if (nounit) temp = temp/a(i,i)
                               else
                                   do 160 k = i + 1,m
                                       temp = temp - conjg(a(k,i))*b(k,j)
       160                         continue
                                   if (nounit) temp = temp/conjg(a(i,i))
                               end if
                               b(i,j) = temp
       170                 continue
       180             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*inv( a ).
     
                   if (upper) then
                       do 230 j = 1,n
                           if (alpha/=one) then
                               do 190 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       190                     continue
                           end if
                           do 210 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   do 200 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       200                         continue
                               end if
       210                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 220 i = 1,m
                                   b(i,j) = temp*b(i,j)
       220                     continue
                           end if
       230             continue
                   else
                       do 280 j = n,1,-1
                           if (alpha/=one) then
                               do 240 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       240                     continue
                           end if
                           do 260 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   do 250 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       250                         continue
                               end if
       260                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 270 i = 1,m
                                   b(i,j) = temp*b(i,j)
       270                     continue
                           end if
       280             continue
                   end if
               else
     
                 ! form  b := alpha*b*inv( a**t )
                 ! or    b := alpha*b*inv( a**h ).
     
                   if (upper) then
                       do 330 k = n,1,-1
                           if (nounit) then
                               if (noconj) then
                                   temp = one/a(k,k)
                               else
                                   temp = one/conjg(a(k,k))
                               end if
                               do 290 i = 1,m
                                   b(i,k) = temp*b(i,k)
       290                     continue
                           end if
                           do 310 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = a(j,k)
                                   else
                                       temp = conjg(a(j,k))
                                   end if
                                   do 300 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       300                         continue
                               end if
       310                 continue
                           if (alpha/=one) then
                               do 320 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       320                     continue
                           end if
       330             continue
                   else
                       do 380 k = 1,n
                           if (nounit) then
                               if (noconj) then
                                   temp = one/a(k,k)
                               else
                                   temp = one/conjg(a(k,k))
                               end if
                               do 340 i = 1,m
                                   b(i,k) = temp*b(i,k)
       340                     continue
                           end if
                           do 360 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   if (noconj) then
                                       temp = a(j,k)
                                   else
                                       temp = conjg(a(j,k))
                                   end if
                                   do 350 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       350                         continue
                               end if
       360                 continue
                           if (alpha/=one) then
                               do 370 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       370                     continue
                           end if
       380             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctrsm
     
     end subroutine stdlib_ctrsm

     
     
     subroutine stdlib_ctrsv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           complex(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           complex(sp) zero
           parameter (zero= (0.0_sp,0.0_sp))
           ! ..
           ! .. local scalars ..
           complex(sp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) noconj,nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic conjg,max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('ctrsv ',info)
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
        10                     continue
                           end if
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 30 i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 70 i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x  or  x := inv( a**h )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 110 j = 1,n
                           temp = x(j)
                           if (noconj) then
                               do 90 i = 1,j - 1
                                   temp = temp - a(i,j)*x(i)
        90                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 100 i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(i)
       100                     continue
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
       110             continue
                   else
                       jx = kx
                       do 140 j = 1,n
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do 120 i = 1,j - 1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix + incx
       120                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 130 i = 1,j - 1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix + incx
       130                     continue
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx + incx
       140             continue
                   end if
               else
                   if (incx==1) then
                       do 170 j = n,1,-1
                           temp = x(j)
                           if (noconj) then
                               do 150 i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(i)
       150                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 160 i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(i)
       160                     continue
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(j) = temp
       170             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 200 j = n,1,-1
                           ix = kx
                           temp = x(jx)
                           if (noconj) then
                               do 180 i = n,j + 1,-1
                                   temp = temp - a(i,j)*x(ix)
                                   ix = ix - incx
       180                     continue
                               if (nounit) temp = temp/a(j,j)
                           else
                               do 190 i = n,j + 1,-1
                                   temp = temp - conjg(a(i,j))*x(ix)
                                   ix = ix - incx
       190                     continue
                               if (nounit) temp = temp/conjg(a(j,j))
                           end if
                           x(jx) = temp
                           jx = jx - incx
       200             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_ctrsv
     
     end subroutine stdlib_ctrsv

     
     
     subroutine stdlib_dgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) incx,incy,kl,ku,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,k,kup1,kx,ky,lenx,leny
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dgbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kup1 = ku + 1
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       k = kup1 - j
                       do 50 i = max(1,j-ku),min(m,j+kl)
                           y(i) = y(i) + temp*a(k+i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       k = kup1 - j
                       do 70 i = max(1,j-ku),min(m,j+kl)
                           y(iy) = y(iy) + temp*a(k+i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
                       if (j>ku) ky = ky + incy
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y.
     
               jy = ky
               if (incx==1) then
                   do 100 j = 1,n
                       temp = zero
                       k = kup1 - j
                       do 90 i = max(1,j-ku),min(m,j+kl)
                           temp = temp + a(k+i,j)*x(i)
        90             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       100         continue
               else
                   do 120 j = 1,n
                       temp = zero
                       ix = kx
                       k = kup1 - j
                       do 110 i = max(1,j-ku),min(m,j+kl)
                           temp = temp + a(k+i,j)*x(ix)
                           ix = ix + incx
       110             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j>ku) kx = kx + incx
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dgbmv
     
     end subroutine stdlib_dgbmv

     
     
     subroutine stdlib_dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) k,lda,ldb,ldc,m,n
           character transa,transb
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,j,l,nrowa,nrowb
           logical(lk) nota,notb
           ! ..
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
     
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
     
           nota = stdlib_lsame(transa,'n')
           notb = stdlib_lsame(transb,'n')
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
           if ((.not.nota) .and. (.not.stdlib_lsame(transa,'c')) .and.(.not.stdlib_lsame(transa,&
          't'))) then
               info = 1
           else if ((.not.notb) .and. (.not.stdlib_lsame(transb,'c')) .and.(.not.stdlib_lsame(&
          transb,'t'))) then
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
               call stdlib_xerbla('dgemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.(((alpha==zero).or. (k==0)).and. (beta==one)))&
           return
     
           ! and if  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (notb) then
               if (nota) then
     
                 ! form  c := alpha*a*b + beta*c.
     
                   do 90 j = 1,n
                       if (beta==zero) then
                           do 50 i = 1,m
                               c(i,j) = zero
        50                 continue
                       else if (beta/=one) then
                           do 60 i = 1,m
                               c(i,j) = beta*c(i,j)
        60                 continue
                       end if
                       do 80 l = 1,k
                           temp = alpha*b(l,j)
                           do 70 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
        70                 continue
        80             continue
        90         continue
               else
     
                 ! form  c := alpha*a**t*b + beta*c
     
                   do 120 j = 1,n
                       do 110 i = 1,m
                           temp = zero
                           do 100 l = 1,k
                               temp = temp + a(l,i)*b(l,j)
       100                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       110             continue
       120         continue
               end if
           else
               if (nota) then
     
                 ! form  c := alpha*a*b**t + beta*c
     
                   do 170 j = 1,n
                       if (beta==zero) then
                           do 130 i = 1,m
                               c(i,j) = zero
       130                 continue
                       else if (beta/=one) then
                           do 140 i = 1,m
                               c(i,j) = beta*c(i,j)
       140                 continue
                       end if
                       do 160 l = 1,k
                           temp = alpha*b(j,l)
                           do 150 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
       150                 continue
       160             continue
       170         continue
               else
     
                 ! form  c := alpha*a**t*b**t + beta*c
     
                   do 200 j = 1,n
                       do 190 i = 1,m
                           temp = zero
                           do 180 l = 1,k
                               temp = temp + a(l,i)*b(j,l)
       180                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       190             continue
       200         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dgemm
     
     end subroutine stdlib_dgemm

     
     
     subroutine stdlib_dgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) incx,incy,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dgemv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       do 50 i = 1,m
                           y(i) = y(i) + temp*a(i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       do 70 i = 1,m
                           y(iy) = y(iy) + temp*a(i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y.
     
               jy = ky
               if (incx==1) then
                   do 100 j = 1,n
                       temp = zero
                       do 90 i = 1,m
                           temp = temp + a(i,j)*x(i)
        90             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       100         continue
               else
                   do 120 j = 1,n
                       temp = zero
                       ix = kx
                       do 110 i = 1,m
                           temp = temp + a(i,j)*x(ix)
                           ix = ix + incx
       110             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dgemv
     
     end subroutine stdlib_dgemv

     
     
     subroutine stdlib_dger(m,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,incy,lda,m,n
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
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
               call stdlib_xerbla('dger  ',info)
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
               do 20 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       do 10 i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
        10             continue
                   end if
                   jy = jy + incy
        20     continue
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               do 40 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       ix = kx
                       do 30 i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
        30             continue
                   end if
                   jy = jy + incy
        40     continue
           end if
     
           return
     
           ! end of stdlib_dger
     
     end subroutine stdlib_dger

     
     
     subroutine stdlib_dsbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) incx,incy,k,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
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
               call stdlib_xerbla('dsbmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when upper triangle of a is stored.
     
               kplus1 = k + 1
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       l = kplus1 - j
                       do 50 i = max(1,j-k),j - 1
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(i)
        50             continue
                       y(j) = y(j) + temp1*a(kplus1,j) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       do 70 i = max(1,j-k),j - 1
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*a(kplus1,j) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       if (j>k) then
                           kx = kx + incx
                           ky = ky + incy
                       end if
        80         continue
               end if
           else
     
              ! form  y  when lower triangle of a is stored.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*a(1,j)
                       l = 1 - j
                       do 90 i = j + 1,min(n,j+k)
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*a(1,j)
                       l = 1 - j
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,min(n,j+k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dsbmv
     
     end subroutine stdlib_dsbmv

     
     
     subroutine stdlib_dspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
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
               call stdlib_xerbla('dspmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when ap contains the upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       k = kk
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(i)
                           k = k + 1
        50             continue
                       y(j) = y(j) + temp1*ap(kk+j-1) + alpha*temp2
                       kk = kk + j
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 k = kk,kk + j - 2
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*ap(kk+j-1) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        80         continue
               end if
           else
     
              ! form  y  when ap contains the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*ap(kk)
                       k = kk + 1
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(i)
                           k = k + 1
        90             continue
                       y(j) = y(j) + alpha*temp2
                       kk = kk + (n-j+1)
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*ap(kk)
                       ix = jx
                       iy = jy
                       do 110 k = kk + 1,kk + n - j
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + (n-j+1)
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dspmv
     
     end subroutine stdlib_dspmv

     
     
     subroutine stdlib_dspr(uplo,n,alpha,x,incx,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
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
               call stdlib_xerbla('dspr  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==zero)) return
     
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           k = kk
                           do 10 i = 1,j
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        10                 continue
                       end if
                       kk = kk + j
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = kx
                           do 30 k = kk,kk + j - 1
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                       end if
                       jx = jx + incx
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           k = kk
                           do 50 i = j,n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        50                 continue
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = jx
                           do 70 k = kk,kk + n - j
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
        70                 continue
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dspr
     
     end subroutine stdlib_dspr

     
     
     subroutine stdlib_dspr2(uplo,n,alpha,x,incx,y,incy,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
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
               call stdlib_xerbla('dspr2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           k = kk
                           do 10 i = 1,j
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        10                 continue
                       end if
                       kk = kk + j
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = kx
                           iy = ky
                           do 30 k = kk,kk + j - 1
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           k = kk
                           do 50 i = j,n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        50                 continue
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = jx
                           iy = jy
                           do 70 k = kk,kk + n - j
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        70                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dspr2
     
     end subroutine stdlib_dspr2

     
     
     subroutine stdlib_dsymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) lda,ldb,ldc,m,n
           character side,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
           integer(int32) i,info,j,k,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
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
               call stdlib_xerbla('dsymm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   do 70 j = 1,n
                       do 60 i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 50 k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        50                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        60             continue
        70         continue
               else
                   do 100 j = 1,n
                       do 90 i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 80 k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        80                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        90             continue
       100         continue
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               do 170 j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta==zero) then
                       do 110 i = 1,m
                           c(i,j) = temp1*b(i,j)
       110             continue
                   else
                       do 120 i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
       120             continue
                   end if
                   do 140 k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do 130 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       130             continue
       140         continue
                   do 160 k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do 150 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       150             continue
       160         continue
       170     continue
           end if
     
           return
     
           ! end of stdlib_dsymm
     
     end subroutine stdlib_dsymm

     
     
     subroutine stdlib_dsymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
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
               call stdlib_xerbla('dsymv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when a is stored in upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(i)
        50             continue
                       y(j) = y(j) + temp1*a(j,j) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 i = 1,j - 1
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           else
     
              ! form  y  when a is stored in lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*a(j,j)
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*a(j,j)
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dsymv
     
     end subroutine stdlib_dsymv

     
     
     subroutine stdlib_dsyr(uplo,n,alpha,x,incx,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,kx
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
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
               call stdlib_xerbla('dsyr  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==zero)) return
     
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           do 10 i = 1,j
                               a(i,j) = a(i,j) + x(i)*temp
        10                 continue
                       end if
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = kx
                           do 30 i = 1,j
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                       end if
                       jx = jx + incx
        40         continue
               end if
           else
     
              ! form  a  when a is stored in lower triangle.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           do 50 i = j,n
                               a(i,j) = a(i,j) + x(i)*temp
        50                 continue
                       end if
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = jx
                           do 70 i = j,n
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
        70                 continue
                       end if
                       jx = jx + incx
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dsyr
     
     end subroutine stdlib_dsyr

     
     
     subroutine stdlib_dsyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
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
               call stdlib_xerbla('dsyr2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           do 10 i = 1,j
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        10                 continue
                       end if
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = kx
                           iy = ky
                           do 30 i = 1,j
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
        40         continue
               end if
           else
     
              ! form  a  when a is stored in the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           do 50 i = j,n
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        50                 continue
                       end if
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = jx
                           iy = jy
                           do 70 i = j,n
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        70                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dsyr2
     
     end subroutine stdlib_dsyr2

     
     
     subroutine stdlib_dsyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) k,lda,ldb,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(dp) temp1,temp2
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')) .and.(&
          .not.stdlib_lsame(trans,'c'))) then
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
               call stdlib_xerbla('dsyr2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do 190 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp1 = zero
                           temp2 = zero
                           do 220 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dsyr2k
     
     end subroutine stdlib_dsyr2k

     
     
     subroutine stdlib_dsyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha,beta
           integer(int32) k,lda,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(dp) temp
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')) .and.(&
          .not.stdlib_lsame(trans,'c'))) then
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
               call stdlib_xerbla('dsyrk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**t + beta*c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*a + beta*c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp = zero
                           do 190 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp = zero
                           do 220 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_dsyrk
     
     end subroutine stdlib_dsyrk

     
     
     subroutine stdlib_dtbmv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dtbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do 10 i = max(1,j-k),j - 1
                                   x(i) = x(i) + temp*a(l+i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do 30 i = max(1,j-k),j - 1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j>k) kx = kx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = 1 - j
                               do 50 i = min(n,j+k),j + 1,-1
                                   x(i) = x(i) + temp*a(l+i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do 70 i = min(n,j+k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 100 j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do 90 i = j - 1,max(1,j-k),-1
                               temp = temp + a(l+i,j)*x(i)
        90                 continue
                           x(j) = temp
       100             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 120 j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do 110 i = j - 1,max(1,j-k),-1
                               temp = temp + a(l+i,j)*x(ix)
                               ix = ix - incx
       110                 continue
                           x(jx) = temp
                           jx = jx - incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do 130 i = j + 1,min(n,j+k)
                               temp = temp + a(l+i,j)*x(i)
       130                 continue
                           x(j) = temp
       140             continue
                   else
                       jx = kx
                       do 160 j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do 150 i = j + 1,min(n,j+k)
                               temp = temp + a(l+i,j)*x(ix)
                               ix = ix + incx
       150                 continue
                           x(jx) = temp
                           jx = jx + incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtbmv
     
     end subroutine stdlib_dtbmv

     
     
     subroutine stdlib_dtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dtbsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do 10 i = j - 1,max(1,j-k),-1
                                   x(i) = x(i) - temp*a(l+i,j)
        10                     continue
                           end if
        20             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 40 j = n,1,-1
                           kx = kx - incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do 30 i = j - 1,max(1,j-k),-1
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix - incx
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do 50 i = j + 1,min(n,j+k)
                                   x(i) = x(i) - temp*a(l+i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           kx = kx + incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do 70 i = j + 1,min(n,j+k)
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix + incx
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t)*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 100 j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           do 90 i = max(1,j-k),j - 1
                               temp = temp - a(l+i,j)*x(i)
        90                 continue
                           if (nounit) temp = temp/a(kplus1,j)
                           x(j) = temp
       100             continue
                   else
                       jx = kx
                       do 120 j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           do 110 i = max(1,j-k),j - 1
                               temp = temp - a(l+i,j)*x(ix)
                               ix = ix + incx
       110                 continue
                           if (nounit) temp = temp/a(kplus1,j)
                           x(jx) = temp
                           jx = jx + incx
                           if (j>k) kx = kx + incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           do 130 i = min(n,j+k),j + 1,-1
                               temp = temp - a(l+i,j)*x(i)
       130                 continue
                           if (nounit) temp = temp/a(1,j)
                           x(j) = temp
       140             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 160 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           do 150 i = min(n,j+k),j + 1,-1
                               temp = temp - a(l+i,j)*x(ix)
                               ix = ix - incx
       150                 continue
                           if (nounit) temp = temp/a(1,j)
                           x(jx) = temp
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtbsv
     
     end subroutine stdlib_dtbsv

     
     
     subroutine stdlib_dtpmv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) nounit
           ! ..
     
     
     
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dtpmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
        10                     continue
                               if (nounit) x(j) = x(j)*ap(kk+j-1)
                           end if
                           kk = kk + j
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk+j-1)
                           end if
                           jx = jx + incx
                           kk = kk + j
        40             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
        50                     continue
                               if (nounit) x(j) = x(j)*ap(kk-n+j)
                           end if
                           kk = kk - (n-j+1)
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 k = kk,kk - (n- (j+1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk-n+j)
                           end if
                           jx = jx - incx
                           kk = kk - (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 100 j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk - 1
                           do 90 i = j - 1,1,-1
                               temp = temp + ap(k)*x(i)
                               k = k - 1
        90                 continue
                           x(j) = temp
                           kk = kk - j
       100             continue
                   else
                       jx = kx + (n-1)*incx
                       do 120 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do 110 k = kk - 1,kk - j + 1,-1
                               ix = ix - incx
                               temp = temp + ap(k)*x(ix)
       110                 continue
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
       120             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 140 j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk + 1
                           do 130 i = j + 1,n
                               temp = temp + ap(k)*x(i)
                               k = k + 1
       130                 continue
                           x(j) = temp
                           kk = kk + (n-j+1)
       140             continue
                   else
                       jx = kx
                       do 160 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do 150 k = kk + 1,kk + n - j
                               ix = ix + incx
                               temp = temp + ap(k)*x(ix)
       150                 continue
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n-j+1)
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtpmv
     
     end subroutine stdlib_dtpmv

     
     
     subroutine stdlib_dtpsv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) nounit
           ! ..
     
     
     
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dtpsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
        10                     continue
                           end if
                           kk = kk - j
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 30 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
        30                     continue
                           end if
                           jx = jx - incx
                           kk = kk - j
        40             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
        50                     continue
                           end if
                           kk = kk + (n-j+1)
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 70 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
        70                     continue
                           end if
                           jx = jx + incx
                           kk = kk + (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = 1
                   if (incx==1) then
                       do 100 j = 1,n
                           temp = x(j)
                           k = kk
                           do 90 i = 1,j - 1
                               temp = temp - ap(k)*x(i)
                               k = k + 1
        90                 continue
                           if (nounit) temp = temp/ap(kk+j-1)
                           x(j) = temp
                           kk = kk + j
       100             continue
                   else
                       jx = kx
                       do 120 j = 1,n
                           temp = x(jx)
                           ix = kx
                           do 110 k = kk,kk + j - 2
                               temp = temp - ap(k)*x(ix)
                               ix = ix + incx
       110                 continue
                           if (nounit) temp = temp/ap(kk+j-1)
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
       120             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 140 j = n,1,-1
                           temp = x(j)
                           k = kk
                           do 130 i = n,j + 1,-1
                               temp = temp - ap(k)*x(i)
                               k = k - 1
       130                 continue
                           if (nounit) temp = temp/ap(kk-n+j)
                           x(j) = temp
                           kk = kk - (n-j+1)
       140             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 160 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do 150 k = kk,kk - (n- (j+1)),-1
                               temp = temp - ap(k)*x(ix)
                               ix = ix - incx
       150                 continue
                           if (nounit) temp = temp/ap(kk-n+j)
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n-j+1)
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtpsv
     
     end subroutine stdlib_dtpsv

     
     
     subroutine stdlib_dtrmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,nounit,upper
           ! ..
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
     
           ! test the input parameters.
     
           lside = stdlib_lsame(side,'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('dtrmm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*a*b.
     
                   if (upper) then
                       do 50 j = 1,n
                           do 40 k = 1,m
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   do 30 i = 1,k - 1
                                       b(i,j) = b(i,j) + temp*a(i,k)
        30                         continue
                                   if (nounit) temp = temp*a(k,k)
                                   b(k,j) = temp
                               end if
        40                 continue
        50             continue
                   else
                       do 80 j = 1,n
                           do 70 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   b(k,j) = temp
                                   if (nounit) b(k,j) = b(k,j)*a(k,k)
                                   do 60 i = k + 1,m
                                       b(i,j) = b(i,j) + temp*a(i,k)
        60                         continue
                               end if
        70                 continue
        80             continue
                   end if
               else
     
                 ! form  b := alpha*a**t*b.
     
                   if (upper) then
                       do 110 j = 1,n
                           do 100 i = m,1,-1
                               temp = b(i,j)
                               if (nounit) temp = temp*a(i,i)
                               do 90 k = 1,i - 1
                                   temp = temp + a(k,i)*b(k,j)
        90                     continue
                               b(i,j) = alpha*temp
       100                 continue
       110             continue
                   else
                       do 140 j = 1,n
                           do 130 i = 1,m
                               temp = b(i,j)
                               if (nounit) temp = temp*a(i,i)
                               do 120 k = i + 1,m
                                   temp = temp + a(k,i)*b(k,j)
       120                     continue
                               b(i,j) = alpha*temp
       130                 continue
       140             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*a.
     
                   if (upper) then
                       do 180 j = n,1,-1
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 150 i = 1,m
                               b(i,j) = temp*b(i,j)
       150                 continue
                           do 170 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 160 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       160                         continue
                               end if
       170                 continue
       180             continue
                   else
                       do 220 j = 1,n
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 190 i = 1,m
                               b(i,j) = temp*b(i,j)
       190                 continue
                           do 210 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 200 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       200                         continue
                               end if
       210                 continue
       220             continue
                   end if
               else
     
                 ! form  b := alpha*b*a**t.
     
                   if (upper) then
                       do 260 k = 1,n
                           do 240 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   temp = alpha*a(j,k)
                                   do 230 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       230                         continue
                               end if
       240                 continue
                           temp = alpha
                           if (nounit) temp = temp*a(k,k)
                           if (temp/=one) then
                               do 250 i = 1,m
                                   b(i,k) = temp*b(i,k)
       250                     continue
                           end if
       260             continue
                   else
                       do 300 k = n,1,-1
                           do 280 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   temp = alpha*a(j,k)
                                   do 270 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       270                         continue
                               end if
       280                 continue
                           temp = alpha
                           if (nounit) temp = temp*a(k,k)
                           if (temp/=one) then
                               do 290 i = 1,m
                                   b(i,k) = temp*b(i,k)
       290                     continue
                           end if
       300             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtrmm
     
     end subroutine stdlib_dtrmm

     
     
     subroutine stdlib_dtrmv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dtrmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 100 j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do 90 i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
        90                 continue
                           x(j) = temp
       100             continue
                   else
                       jx = kx + (n-1)*incx
                       do 120 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do 110 i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
       110                 continue
                           x(jx) = temp
                           jx = jx - incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do 130 i = j + 1,n
                               temp = temp + a(i,j)*x(i)
       130                 continue
                           x(j) = temp
       140             continue
                   else
                       jx = kx
                       do 160 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do 150 i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
       150                 continue
                           x(jx) = temp
                           jx = jx + incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtrmv
     
     end subroutine stdlib_dtrmv

     
     
     subroutine stdlib_dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(dp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,nounit,upper
           ! ..
           ! .. parameters ..
           real(dp) one,zero
           parameter (one=1.0_dp,zero=0.0_dp)
           ! ..
     
           ! test the input parameters.
     
           lside = stdlib_lsame(side,'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('dtrsm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*inv( a )*b.
     
                   if (upper) then
                       do 60 j = 1,n
                           if (alpha/=one) then
                               do 30 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        30                     continue
                           end if
                           do 50 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 40 i = 1,k - 1
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        40                         continue
                               end if
        50                 continue
        60             continue
                   else
                       do 100 j = 1,n
                           if (alpha/=one) then
                               do 70 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        70                     continue
                           end if
                           do 90 k = 1,m
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 80 i = k + 1,m
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        80                         continue
                               end if
        90                 continue
       100             continue
                   end if
               else
     
                 ! form  b := alpha*inv( a**t )*b.
     
                   if (upper) then
                       do 130 j = 1,n
                           do 120 i = 1,m
                               temp = alpha*b(i,j)
                               do 110 k = 1,i - 1
                                   temp = temp - a(k,i)*b(k,j)
       110                     continue
                               if (nounit) temp = temp/a(i,i)
                               b(i,j) = temp
       120                 continue
       130             continue
                   else
                       do 160 j = 1,n
                           do 150 i = m,1,-1
                               temp = alpha*b(i,j)
                               do 140 k = i + 1,m
                                   temp = temp - a(k,i)*b(k,j)
       140                     continue
                               if (nounit) temp = temp/a(i,i)
                               b(i,j) = temp
       150                 continue
       160             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*inv( a ).
     
                   if (upper) then
                       do 210 j = 1,n
                           if (alpha/=one) then
                               do 170 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       170                     continue
                           end if
                           do 190 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   do 180 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       180                         continue
                               end if
       190                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 200 i = 1,m
                                   b(i,j) = temp*b(i,j)
       200                     continue
                           end if
       210             continue
                   else
                       do 260 j = n,1,-1
                           if (alpha/=one) then
                               do 220 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       220                     continue
                           end if
                           do 240 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   do 230 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       230                         continue
                               end if
       240                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 250 i = 1,m
                                   b(i,j) = temp*b(i,j)
       250                     continue
                           end if
       260             continue
                   end if
               else
     
                 ! form  b := alpha*b*inv( a**t ).
     
                   if (upper) then
                       do 310 k = n,1,-1
                           if (nounit) then
                               temp = one/a(k,k)
                               do 270 i = 1,m
                                   b(i,k) = temp*b(i,k)
       270                     continue
                           end if
                           do 290 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   temp = a(j,k)
                                   do 280 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       280                         continue
                               end if
       290                 continue
                           if (alpha/=one) then
                               do 300 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       300                     continue
                           end if
       310             continue
                   else
                       do 360 k = 1,n
                           if (nounit) then
                               temp = one/a(k,k)
                               do 320 i = 1,m
                                   b(i,k) = temp*b(i,k)
       320                     continue
                           end if
                           do 340 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   temp = a(j,k)
                                   do 330 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       330                         continue
                               end if
       340                 continue
                           if (alpha/=one) then
                               do 350 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       350                     continue
                           end if
       360             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtrsm
     
     end subroutine stdlib_dtrsm

     
     
     subroutine stdlib_dtrsv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(dp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(dp) zero
           parameter (zero=0.0_dp)
           ! ..
           ! .. local scalars ..
           real(dp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('dtrsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
        10                     continue
                           end if
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 30 i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 70 i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 100 j = 1,n
                           temp = x(j)
                           do 90 i = 1,j - 1
                               temp = temp - a(i,j)*x(i)
        90                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
       100             continue
                   else
                       jx = kx
                       do 120 j = 1,n
                           temp = x(jx)
                           ix = kx
                           do 110 i = 1,j - 1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix + incx
       110                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx + incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = n,1,-1
                           temp = x(j)
                           do 130 i = n,j + 1,-1
                               temp = temp - a(i,j)*x(i)
       130                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
       140             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 160 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do 150 i = n,j + 1,-1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix - incx
       150                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx - incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_dtrsv
     
     end subroutine stdlib_dtrsv

     
     
     integer(int32) function stdlib_icamax(n,cx,incx)
     
        ! -- reference blas level1 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           ! ..
           ! .. array arguments ..
           complex(sp) cx(*)
           ! ..
     
        ! =====================================================================
     
           ! .. local scalars ..
           real(sp) smax
           integer(int32) i,ix
           ! ..
     
     
           stdlib_icamax = 0
           if (n<1 .or. incx<=0) return
           stdlib_icamax = 1
           if (n==1) return
           if (incx==1) then
     
              ! code for increment equal to 1
     
              smax = stdlib_scabs1(cx(1))
              do i = 2,n
                 if (stdlib_scabs1(cx(i))>smax) then
                    stdlib_icamax = i
                    smax = stdlib_scabs1(cx(i))
                 end if
              end do
           else
     
              ! code for increment not equal to 1
     
              ix = 1
              smax = stdlib_scabs1(cx(1))
              ix = ix + incx
              do i = 2,n
                 if (stdlib_scabs1(cx(ix))>smax) then
                    stdlib_icamax = i
                    smax = stdlib_scabs1(cx(ix))
                 end if
                 ix = ix + incx
              end do
           end if
           return
     
           ! end of stdlib_icamax
     
     end function stdlib_icamax

     
     
     subroutine stdlib_sgbmv(trans,m,n,kl,ku,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) incx,incy,kl,ku,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,k,kup1,kx,ky,lenx,leny
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('sgbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kup1 = ku + 1
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       k = kup1 - j
                       do 50 i = max(1,j-ku),min(m,j+kl)
                           y(i) = y(i) + temp*a(k+i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       k = kup1 - j
                       do 70 i = max(1,j-ku),min(m,j+kl)
                           y(iy) = y(iy) + temp*a(k+i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
                       if (j>ku) ky = ky + incy
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y.
     
               jy = ky
               if (incx==1) then
                   do 100 j = 1,n
                       temp = zero
                       k = kup1 - j
                       do 90 i = max(1,j-ku),min(m,j+kl)
                           temp = temp + a(k+i,j)*x(i)
        90             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       100         continue
               else
                   do 120 j = 1,n
                       temp = zero
                       ix = kx
                       k = kup1 - j
                       do 110 i = max(1,j-ku),min(m,j+kl)
                           temp = temp + a(k+i,j)*x(ix)
                           ix = ix + incx
       110             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
                       if (j>ku) kx = kx + incx
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_sgbmv
     
     end subroutine stdlib_sgbmv

     
     
     subroutine stdlib_sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) k,lda,ldb,ldc,m,n
           character transa,transb
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,j,l,nrowa,nrowb
           logical(lk) nota,notb
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
     
           ! set  nota  and  notb  as  true if  a  and  b  respectively are not
           ! transposed and set  nrowa and nrowb  as the number of rows of  a
           ! and  b  respectively.
     
           nota = stdlib_lsame(transa,'n')
           notb = stdlib_lsame(transb,'n')
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
           if ((.not.nota) .and. (.not.stdlib_lsame(transa,'c')) .and.(.not.stdlib_lsame(transa,&
          't'))) then
               info = 1
           else if ((.not.notb) .and. (.not.stdlib_lsame(transb,'c')) .and.(.not.stdlib_lsame(&
          transb,'t'))) then
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
               call stdlib_xerbla('sgemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.(((alpha==zero).or. (k==0)).and. (beta==one)))&
           return
     
           ! and if  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (notb) then
               if (nota) then
     
                 ! form  c := alpha*a*b + beta*c.
     
                   do 90 j = 1,n
                       if (beta==zero) then
                           do 50 i = 1,m
                               c(i,j) = zero
        50                 continue
                       else if (beta/=one) then
                           do 60 i = 1,m
                               c(i,j) = beta*c(i,j)
        60                 continue
                       end if
                       do 80 l = 1,k
                           temp = alpha*b(l,j)
                           do 70 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
        70                 continue
        80             continue
        90         continue
               else
     
                 ! form  c := alpha*a**t*b + beta*c
     
                   do 120 j = 1,n
                       do 110 i = 1,m
                           temp = zero
                           do 100 l = 1,k
                               temp = temp + a(l,i)*b(l,j)
       100                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       110             continue
       120         continue
               end if
           else
               if (nota) then
     
                 ! form  c := alpha*a*b**t + beta*c
     
                   do 170 j = 1,n
                       if (beta==zero) then
                           do 130 i = 1,m
                               c(i,j) = zero
       130                 continue
                       else if (beta/=one) then
                           do 140 i = 1,m
                               c(i,j) = beta*c(i,j)
       140                 continue
                       end if
                       do 160 l = 1,k
                           temp = alpha*b(j,l)
                           do 150 i = 1,m
                               c(i,j) = c(i,j) + temp*a(i,l)
       150                 continue
       160             continue
       170         continue
               else
     
                 ! form  c := alpha*a**t*b**t + beta*c
     
                   do 200 j = 1,n
                       do 190 i = 1,m
                           temp = zero
                           do 180 l = 1,k
                               temp = temp + a(l,i)*b(j,l)
       180                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       190             continue
       200         continue
               end if
           end if
     
           return
     
           ! end of stdlib_sgemm
     
     end subroutine stdlib_sgemm

     
     
     subroutine stdlib_sgemv(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) incx,incy,lda,m,n
           character trans
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky,lenx,leny
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('sgemv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
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
                       do 10 i = 1,leny
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,leny
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,leny
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,leny
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(trans,'n')) then
     
              ! form  y := alpha*a*x + y.
     
               jx = kx
               if (incy==1) then
                   do 60 j = 1,n
                       temp = alpha*x(jx)
                       do 50 i = 1,m
                           y(i) = y(i) + temp*a(i,j)
        50             continue
                       jx = jx + incx
        60         continue
               else
                   do 80 j = 1,n
                       temp = alpha*x(jx)
                       iy = ky
                       do 70 i = 1,m
                           y(iy) = y(iy) + temp*a(i,j)
                           iy = iy + incy
        70             continue
                       jx = jx + incx
        80         continue
               end if
           else
     
              ! form  y := alpha*a**t*x + y.
     
               jy = ky
               if (incx==1) then
                   do 100 j = 1,n
                       temp = zero
                       do 90 i = 1,m
                           temp = temp + a(i,j)*x(i)
        90             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       100         continue
               else
                   do 120 j = 1,n
                       temp = zero
                       ix = kx
                       do 110 i = 1,m
                           temp = temp + a(i,j)*x(ix)
                           ix = ix + incx
       110             continue
                       y(jy) = y(jy) + alpha*temp
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_sgemv
     
     end subroutine stdlib_sgemv

     
     
     subroutine stdlib_sger(m,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,incy,lda,m,n
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
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
               call stdlib_xerbla('sger  ',info)
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
               do 20 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       do 10 i = 1,m
                           a(i,j) = a(i,j) + x(i)*temp
        10             continue
                   end if
                   jy = jy + incy
        20     continue
           else
               if (incx>0) then
                   kx = 1
               else
                   kx = 1 - (m-1)*incx
               end if
               do 40 j = 1,n
                   if (y(jy)/=zero) then
                       temp = alpha*y(jy)
                       ix = kx
                       do 30 i = 1,m
                           a(i,j) = a(i,j) + x(ix)*temp
                           ix = ix + incx
        30             continue
                   end if
                   jy = jy + incy
        40     continue
           end if
     
           return
     
           ! end of stdlib_sger
     
     end subroutine stdlib_sger

     
     
     subroutine stdlib_ssbmv(uplo,n,k,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) incx,incy,k,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kplus1,kx,ky,l
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
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
               call stdlib_xerbla('ssbmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when upper triangle of a is stored.
     
               kplus1 = k + 1
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       l = kplus1 - j
                       do 50 i = max(1,j-k),j - 1
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(i)
        50             continue
                       y(j) = y(j) + temp1*a(kplus1,j) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       l = kplus1 - j
                       do 70 i = max(1,j-k),j - 1
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*a(kplus1,j) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       if (j>k) then
                           kx = kx + incx
                           ky = ky + incy
                       end if
        80         continue
               end if
           else
     
              ! form  y  when lower triangle of a is stored.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*a(1,j)
                       l = 1 - j
                       do 90 i = j + 1,min(n,j+k)
                           y(i) = y(i) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*a(1,j)
                       l = 1 - j
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,min(n,j+k)
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(l+i,j)
                           temp2 = temp2 + a(l+i,j)*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_ssbmv
     
     end subroutine stdlib_ssbmv

     
     
     subroutine stdlib_sspmv(uplo,n,alpha,ap,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
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
               call stdlib_xerbla('sspmv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           kk = 1
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when ap contains the upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       k = kk
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(i)
                           k = k + 1
        50             continue
                       y(j) = y(j) + temp1*ap(kk+j-1) + alpha*temp2
                       kk = kk + j
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 k = kk,kk + j - 2
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*ap(kk+j-1) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        80         continue
               end if
           else
     
              ! form  y  when ap contains the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*ap(kk)
                       k = kk + 1
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(i)
                           k = k + 1
        90             continue
                       y(j) = y(j) + alpha*temp2
                       kk = kk + (n-j+1)
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*ap(kk)
                       ix = jx
                       iy = jy
                       do 110 k = kk + 1,kk + n - j
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*ap(k)
                           temp2 = temp2 + ap(k)*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + (n-j+1)
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_sspmv
     
     end subroutine stdlib_sspmv

     
     
     subroutine stdlib_sspr(uplo,n,alpha,x,incx,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
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
               call stdlib_xerbla('sspr  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==zero)) return
     
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           k = kk
                           do 10 i = 1,j
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        10                 continue
                       end if
                       kk = kk + j
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = kx
                           do 30 k = kk,kk + j - 1
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                       end if
                       jx = jx + incx
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           k = kk
                           do 50 i = j,n
                               ap(k) = ap(k) + x(i)*temp
                               k = k + 1
        50                 continue
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = jx
                           do 70 k = kk,kk + n - j
                               ap(k) = ap(k) + x(ix)*temp
                               ix = ix + incx
        70                 continue
                       end if
                       jx = jx + incx
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_sspr
     
     end subroutine stdlib_sspr

     
     
     subroutine stdlib_sspr2(uplo,n,alpha,x,incx,y,incy,ap)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,incy,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) ap(*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,k,kk,kx,ky
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
               call stdlib_xerbla('sspr2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           k = kk
                           do 10 i = 1,j
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        10                 continue
                       end if
                       kk = kk + j
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = kx
                           iy = ky
                           do 30 k = kk,kk + j - 1
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + j
        40         continue
               end if
           else
     
              ! form  a  when lower triangle is stored in ap.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           k = kk
                           do 50 i = j,n
                               ap(k) = ap(k) + x(i)*temp1 + y(i)*temp2
                               k = k + 1
        50                 continue
                       end if
                       kk = kk + n - j + 1
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = jx
                           iy = jy
                           do 70 k = kk,kk + n - j
                               ap(k) = ap(k) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        70                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
                       kk = kk + n - j + 1
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_sspr2
     
     end subroutine stdlib_sspr2

     
     
     subroutine stdlib_ssymm(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) lda,ldb,ldc,m,n
           character side,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,j,k,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
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
               call stdlib_xerbla('ssymm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.((alpha==zero).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (beta==zero) then
                   do 20 j = 1,n
                       do 10 i = 1,m
                           c(i,j) = zero
        10             continue
        20         continue
               else
                   do 40 j = 1,n
                       do 30 i = 1,m
                           c(i,j) = beta*c(i,j)
        30             continue
        40         continue
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(side,'l')) then
     
              ! form  c := alpha*a*b + beta*c.
     
               if (upper) then
                   do 70 j = 1,n
                       do 60 i = 1,m
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 50 k = 1,i - 1
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        50                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        60             continue
        70         continue
               else
                   do 100 j = 1,n
                       do 90 i = m,1,-1
                           temp1 = alpha*b(i,j)
                           temp2 = zero
                           do 80 k = i + 1,m
                               c(k,j) = c(k,j) + temp1*a(k,i)
                               temp2 = temp2 + b(k,j)*a(k,i)
        80                 continue
                           if (beta==zero) then
                               c(i,j) = temp1*a(i,i) + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + temp1*a(i,i) +alpha*temp2
                           end if
        90             continue
       100         continue
               end if
           else
     
              ! form  c := alpha*b*a + beta*c.
     
               do 170 j = 1,n
                   temp1 = alpha*a(j,j)
                   if (beta==zero) then
                       do 110 i = 1,m
                           c(i,j) = temp1*b(i,j)
       110             continue
                   else
                       do 120 i = 1,m
                           c(i,j) = beta*c(i,j) + temp1*b(i,j)
       120             continue
                   end if
                   do 140 k = 1,j - 1
                       if (upper) then
                           temp1 = alpha*a(k,j)
                       else
                           temp1 = alpha*a(j,k)
                       end if
                       do 130 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       130             continue
       140         continue
                   do 160 k = j + 1,n
                       if (upper) then
                           temp1 = alpha*a(j,k)
                       else
                           temp1 = alpha*a(k,j)
                       end if
                       do 150 i = 1,m
                           c(i,j) = c(i,j) + temp1*b(i,k)
       150             continue
       160         continue
       170     continue
           end if
     
           return
     
           ! end of stdlib_ssymm
     
     end subroutine stdlib_ssymm

     
     
     subroutine stdlib_ssymv(uplo,n,alpha,a,lda,x,incx,beta,y,incy)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
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
               call stdlib_xerbla('ssymv ',info)
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
                       do 10 i = 1,n
                           y(i) = zero
        10             continue
                   else
                       do 20 i = 1,n
                           y(i) = beta*y(i)
        20             continue
                   end if
               else
                   iy = ky
                   if (beta==zero) then
                       do 30 i = 1,n
                           y(iy) = zero
                           iy = iy + incy
        30             continue
                   else
                       do 40 i = 1,n
                           y(iy) = beta*y(iy)
                           iy = iy + incy
        40             continue
                   end if
               end if
           end if
           if (alpha==zero) return
           if (stdlib_lsame(uplo,'u')) then
     
              ! form  y  when a is stored in upper triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       do 50 i = 1,j - 1
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(i)
        50             continue
                       y(j) = y(j) + temp1*a(j,j) + alpha*temp2
        60         continue
               else
                   jx = kx
                   jy = ky
                   do 80 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       ix = kx
                       iy = ky
                       do 70 i = 1,j - 1
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(ix)
                           ix = ix + incx
                           iy = iy + incy
        70             continue
                       y(jy) = y(jy) + temp1*a(j,j) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           else
     
              ! form  y  when a is stored in lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 100 j = 1,n
                       temp1 = alpha*x(j)
                       temp2 = zero
                       y(j) = y(j) + temp1*a(j,j)
                       do 90 i = j + 1,n
                           y(i) = y(i) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(i)
        90             continue
                       y(j) = y(j) + alpha*temp2
       100         continue
               else
                   jx = kx
                   jy = ky
                   do 120 j = 1,n
                       temp1 = alpha*x(jx)
                       temp2 = zero
                       y(jy) = y(jy) + temp1*a(j,j)
                       ix = jx
                       iy = jy
                       do 110 i = j + 1,n
                           ix = ix + incx
                           iy = iy + incy
                           y(iy) = y(iy) + temp1*a(i,j)
                           temp2 = temp2 + a(i,j)*x(ix)
       110             continue
                       y(jy) = y(jy) + alpha*temp2
                       jx = jx + incx
                       jy = jy + incy
       120         continue
               end if
           end if
     
           return
     
           ! end of stdlib_ssymv
     
     end subroutine stdlib_ssymv

     
     
     subroutine stdlib_ssyr(uplo,n,alpha,x,incx,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,kx
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
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
               call stdlib_xerbla('ssyr  ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (alpha==zero)) return
     
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
                   do 20 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           do 10 i = 1,j
                               a(i,j) = a(i,j) + x(i)*temp
        10                 continue
                       end if
        20         continue
               else
                   jx = kx
                   do 40 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = kx
                           do 30 i = 1,j
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
        30                 continue
                       end if
                       jx = jx + incx
        40         continue
               end if
           else
     
              ! form  a  when a is stored in lower triangle.
     
               if (incx==1) then
                   do 60 j = 1,n
                       if (x(j)/=zero) then
                           temp = alpha*x(j)
                           do 50 i = j,n
                               a(i,j) = a(i,j) + x(i)*temp
        50                 continue
                       end if
        60         continue
               else
                   jx = kx
                   do 80 j = 1,n
                       if (x(jx)/=zero) then
                           temp = alpha*x(jx)
                           ix = jx
                           do 70 i = j,n
                               a(i,j) = a(i,j) + x(ix)*temp
                               ix = ix + incx
        70                 continue
                       end if
                       jx = jx + incx
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_ssyr
     
     end subroutine stdlib_ssyr

     
     
     subroutine stdlib_ssyr2(uplo,n,alpha,x,incx,y,incy,a,lda)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) incx,incy,lda,n
           character uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*),y(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,ix,iy,j,jx,jy,kx,ky
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
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
               call stdlib_xerbla('ssyr2 ',info)
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
                   do 20 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           do 10 i = 1,j
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        10                 continue
                       end if
        20         continue
               else
                   do 40 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = kx
                           iy = ky
                           do 30 i = 1,j
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        30                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
        40         continue
               end if
           else
     
              ! form  a  when a is stored in the lower triangle.
     
               if ((incx==1) .and. (incy==1)) then
                   do 60 j = 1,n
                       if ((x(j)/=zero) .or. (y(j)/=zero)) then
                           temp1 = alpha*y(j)
                           temp2 = alpha*x(j)
                           do 50 i = j,n
                               a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
        50                 continue
                       end if
        60         continue
               else
                   do 80 j = 1,n
                       if ((x(jx)/=zero) .or. (y(jy)/=zero)) then
                           temp1 = alpha*y(jy)
                           temp2 = alpha*x(jx)
                           ix = jx
                           iy = jy
                           do 70 i = j,n
                               a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                               ix = ix + incx
                               iy = iy + incy
        70                 continue
                       end if
                       jx = jx + incx
                       jy = jy + incy
        80         continue
               end if
           end if
     
           return
     
           ! end of stdlib_ssyr2
     
     end subroutine stdlib_ssyr2

     
     
     subroutine stdlib_ssyr2k(uplo,trans,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) k,lda,ldb,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),b(ldb,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(sp) temp1,temp2
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')) .and.(&
          .not.stdlib_lsame(trans,'c'))) then
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
               call stdlib_xerbla('ssyr2k',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*b**t + alpha*b*a**t + c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if ((a(j,l)/=zero) .or. (b(j,l)/=zero)) then
                               temp1 = alpha*b(j,l)
                               temp2 = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + a(i,l)*temp1 +b(i,l)*temp2
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*b + alpha*b**t*a + c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp1 = zero
                           temp2 = zero
                           do 190 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp1 = zero
                           temp2 = zero
                           do 220 l = 1,k
                               temp1 = temp1 + a(l,i)*b(l,j)
                               temp2 = temp2 + b(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp1 + alpha*temp2
                           else
                               c(i,j) = beta*c(i,j) + alpha*temp1 +alpha*temp2
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_ssyr2k
     
     end subroutine stdlib_ssyr2k

     
     
     subroutine stdlib_ssyrk(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha,beta
           integer(int32) k,lda,ldc,n
           character trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),c(ldc,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,j,l,nrowa
           logical(lk) upper
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
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
           else if ((.not.stdlib_lsame(trans,'n')) .and.(.not.stdlib_lsame(trans,'t')) .and.(&
          .not.stdlib_lsame(trans,'c'))) then
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
               call stdlib_xerbla('ssyrk ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((n==0) .or. (((alpha==zero).or.(k==0)).and. (beta==one))) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               if (upper) then
                   if (beta==zero) then
                       do 20 j = 1,n
                           do 10 i = 1,j
                               c(i,j) = zero
        10                 continue
        20             continue
                   else
                       do 40 j = 1,n
                           do 30 i = 1,j
                               c(i,j) = beta*c(i,j)
        30                 continue
        40             continue
                   end if
               else
                   if (beta==zero) then
                       do 60 j = 1,n
                           do 50 i = j,n
                               c(i,j) = zero
        50                 continue
        60             continue
                   else
                       do 80 j = 1,n
                           do 70 i = j,n
                               c(i,j) = beta*c(i,j)
        70                 continue
        80             continue
                   end if
               end if
               return
           end if
     
           ! start the operations.
     
           if (stdlib_lsame(trans,'n')) then
     
              ! form  c := alpha*a*a**t + beta*c.
     
               if (upper) then
                   do 130 j = 1,n
                       if (beta==zero) then
                           do 90 i = 1,j
                               c(i,j) = zero
        90                 continue
                       else if (beta/=one) then
                           do 100 i = 1,j
                               c(i,j) = beta*c(i,j)
       100                 continue
                       end if
                       do 120 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 110 i = 1,j
                                   c(i,j) = c(i,j) + temp*a(i,l)
       110                     continue
                           end if
       120             continue
       130         continue
               else
                   do 180 j = 1,n
                       if (beta==zero) then
                           do 140 i = j,n
                               c(i,j) = zero
       140                 continue
                       else if (beta/=one) then
                           do 150 i = j,n
                               c(i,j) = beta*c(i,j)
       150                 continue
                       end if
                       do 170 l = 1,k
                           if (a(j,l)/=zero) then
                               temp = alpha*a(j,l)
                               do 160 i = j,n
                                   c(i,j) = c(i,j) + temp*a(i,l)
       160                     continue
                           end if
       170             continue
       180         continue
               end if
           else
     
              ! form  c := alpha*a**t*a + beta*c.
     
               if (upper) then
                   do 210 j = 1,n
                       do 200 i = 1,j
                           temp = zero
                           do 190 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       190                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       200             continue
       210         continue
               else
                   do 240 j = 1,n
                       do 230 i = j,n
                           temp = zero
                           do 220 l = 1,k
                               temp = temp + a(l,i)*a(l,j)
       220                 continue
                           if (beta==zero) then
                               c(i,j) = alpha*temp
                           else
                               c(i,j) = alpha*temp + beta*c(i,j)
                           end if
       230             continue
       240         continue
               end if
           end if
     
           return
     
           ! end of stdlib_ssyrk
     
     end subroutine stdlib_ssyrk

     
     
     subroutine stdlib_stbmv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('stbmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = kplus1 - j
                               do 10 i = max(1,j-k),j - 1
                                   x(i) = x(i) + temp*a(l+i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(kplus1,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = kplus1 - j
                               do 30 i = max(1,j-k),j - 1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(kplus1,j)
                           end if
                           jx = jx + incx
                           if (j>k) kx = kx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               l = 1 - j
                               do 50 i = min(n,j+k),j + 1,-1
                                   x(i) = x(i) + temp*a(l+i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(1,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               l = 1 - j
                               do 70 i = min(n,j+k),j + 1,-1
                                   x(ix) = x(ix) + temp*a(l+i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(1,j)
                           end if
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 100 j = n,1,-1
                           temp = x(j)
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do 90 i = j - 1,max(1,j-k),-1
                               temp = temp + a(l+i,j)*x(i)
        90                 continue
                           x(j) = temp
       100             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 120 j = n,1,-1
                           temp = x(jx)
                           kx = kx - incx
                           ix = kx
                           l = kplus1 - j
                           if (nounit) temp = temp*a(kplus1,j)
                           do 110 i = j - 1,max(1,j-k),-1
                               temp = temp + a(l+i,j)*x(ix)
                               ix = ix - incx
       110                 continue
                           x(jx) = temp
                           jx = jx - incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = 1,n
                           temp = x(j)
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do 130 i = j + 1,min(n,j+k)
                               temp = temp + a(l+i,j)*x(i)
       130                 continue
                           x(j) = temp
       140             continue
                   else
                       jx = kx
                       do 160 j = 1,n
                           temp = x(jx)
                           kx = kx + incx
                           ix = kx
                           l = 1 - j
                           if (nounit) temp = temp*a(1,j)
                           do 150 i = j + 1,min(n,j+k)
                               temp = temp + a(l+i,j)*x(ix)
                               ix = ix + incx
       150                 continue
                           x(jx) = temp
                           jx = jx + incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_stbmv
     
     end subroutine stdlib_stbmv

     
     
     subroutine stdlib_stbsv(uplo,trans,diag,n,k,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,k,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,kplus1,kx,l
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max,min
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('stbsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               l = kplus1 - j
                               if (nounit) x(j) = x(j)/a(kplus1,j)
                               temp = x(j)
                               do 10 i = j - 1,max(1,j-k),-1
                                   x(i) = x(i) - temp*a(l+i,j)
        10                     continue
                           end if
        20             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 40 j = n,1,-1
                           kx = kx - incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = kplus1 - j
                               if (nounit) x(jx) = x(jx)/a(kplus1,j)
                               temp = x(jx)
                               do 30 i = j - 1,max(1,j-k),-1
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix - incx
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               l = 1 - j
                               if (nounit) x(j) = x(j)/a(1,j)
                               temp = x(j)
                               do 50 i = j + 1,min(n,j+k)
                                   x(i) = x(i) - temp*a(l+i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           kx = kx + incx
                           if (x(jx)/=zero) then
                               ix = kx
                               l = 1 - j
                               if (nounit) x(jx) = x(jx)/a(1,j)
                               temp = x(jx)
                               do 70 i = j + 1,min(n,j+k)
                                   x(ix) = x(ix) - temp*a(l+i,j)
                                   ix = ix + incx
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t)*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kplus1 = k + 1
                   if (incx==1) then
                       do 100 j = 1,n
                           temp = x(j)
                           l = kplus1 - j
                           do 90 i = max(1,j-k),j - 1
                               temp = temp - a(l+i,j)*x(i)
        90                 continue
                           if (nounit) temp = temp/a(kplus1,j)
                           x(j) = temp
       100             continue
                   else
                       jx = kx
                       do 120 j = 1,n
                           temp = x(jx)
                           ix = kx
                           l = kplus1 - j
                           do 110 i = max(1,j-k),j - 1
                               temp = temp - a(l+i,j)*x(ix)
                               ix = ix + incx
       110                 continue
                           if (nounit) temp = temp/a(kplus1,j)
                           x(jx) = temp
                           jx = jx + incx
                           if (j>k) kx = kx + incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = n,1,-1
                           temp = x(j)
                           l = 1 - j
                           do 130 i = min(n,j+k),j + 1,-1
                               temp = temp - a(l+i,j)*x(i)
       130                 continue
                           if (nounit) temp = temp/a(1,j)
                           x(j) = temp
       140             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 160 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           l = 1 - j
                           do 150 i = min(n,j+k),j + 1,-1
                               temp = temp - a(l+i,j)*x(ix)
                               ix = ix - incx
       150                 continue
                           if (nounit) temp = temp/a(1,j)
                           x(jx) = temp
                           jx = jx - incx
                           if ((n-j)>=k) kx = kx - incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_stbsv
     
     end subroutine stdlib_stbsv

     
     
     subroutine stdlib_stpmv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) nounit
           ! ..
     
     
     
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('stpmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k + 1
        10                     continue
                               if (nounit) x(j) = x(j)*ap(kk+j-1)
                           end if
                           kk = kk + j
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 k = kk,kk + j - 2
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk+j-1)
                           end if
                           jx = jx + incx
                           kk = kk + j
        40             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               k = kk
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*ap(k)
                                   k = k - 1
        50                     continue
                               if (nounit) x(j) = x(j)*ap(kk-n+j)
                           end if
                           kk = kk - (n-j+1)
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 k = kk,kk - (n- (j+1)),-1
                                   x(ix) = x(ix) + temp*ap(k)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*ap(kk-n+j)
                           end if
                           jx = jx - incx
                           kk = kk - (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 100 j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk - 1
                           do 90 i = j - 1,1,-1
                               temp = temp + ap(k)*x(i)
                               k = k - 1
        90                 continue
                           x(j) = temp
                           kk = kk - j
       100             continue
                   else
                       jx = kx + (n-1)*incx
                       do 120 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do 110 k = kk - 1,kk - j + 1,-1
                               ix = ix - incx
                               temp = temp + ap(k)*x(ix)
       110                 continue
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - j
       120             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 140 j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*ap(kk)
                           k = kk + 1
                           do 130 i = j + 1,n
                               temp = temp + ap(k)*x(i)
                               k = k + 1
       130                 continue
                           x(j) = temp
                           kk = kk + (n-j+1)
       140             continue
                   else
                       jx = kx
                       do 160 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*ap(kk)
                           do 150 k = kk + 1,kk + n - j
                               ix = ix + incx
                               temp = temp + ap(k)*x(ix)
       150                 continue
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + (n-j+1)
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_stpmv
     
     end subroutine stdlib_stpmv

     
     
     subroutine stdlib_stpsv(uplo,trans,diag,n,ap,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) ap(*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,k,kk,kx
           logical(lk) nounit
           ! ..
     
     
     
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('stpsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk - 1
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*ap(k)
                                   k = k - 1
        10                     continue
                           end if
                           kk = kk - j
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 30 k = kk - 1,kk - j + 1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*ap(k)
        30                     continue
                           end if
                           jx = jx - incx
                           kk = kk - j
        40             continue
                   end if
               else
                   kk = 1
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/ap(kk)
                               temp = x(j)
                               k = kk + 1
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*ap(k)
                                   k = k + 1
        50                     continue
                           end if
                           kk = kk + (n-j+1)
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/ap(kk)
                               temp = x(jx)
                               ix = jx
                               do 70 k = kk + 1,kk + n - j
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*ap(k)
        70                     continue
                           end if
                           jx = jx + incx
                           kk = kk + (n-j+1)
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   kk = 1
                   if (incx==1) then
                       do 100 j = 1,n
                           temp = x(j)
                           k = kk
                           do 90 i = 1,j - 1
                               temp = temp - ap(k)*x(i)
                               k = k + 1
        90                 continue
                           if (nounit) temp = temp/ap(kk+j-1)
                           x(j) = temp
                           kk = kk + j
       100             continue
                   else
                       jx = kx
                       do 120 j = 1,n
                           temp = x(jx)
                           ix = kx
                           do 110 k = kk,kk + j - 2
                               temp = temp - ap(k)*x(ix)
                               ix = ix + incx
       110                 continue
                           if (nounit) temp = temp/ap(kk+j-1)
                           x(jx) = temp
                           jx = jx + incx
                           kk = kk + j
       120             continue
                   end if
               else
                   kk = (n* (n+1))/2
                   if (incx==1) then
                       do 140 j = n,1,-1
                           temp = x(j)
                           k = kk
                           do 130 i = n,j + 1,-1
                               temp = temp - ap(k)*x(i)
                               k = k - 1
       130                 continue
                           if (nounit) temp = temp/ap(kk-n+j)
                           x(j) = temp
                           kk = kk - (n-j+1)
       140             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 160 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do 150 k = kk,kk - (n- (j+1)),-1
                               temp = temp - ap(k)*x(ix)
                               ix = ix - incx
       150                 continue
                           if (nounit) temp = temp/ap(kk-n+j)
                           x(jx) = temp
                           jx = jx - incx
                           kk = kk - (n-j+1)
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_stpsv
     
     end subroutine stdlib_stpsv

     
     
     subroutine stdlib_strmm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,nounit,upper
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
     
           ! test the input parameters.
     
           lside = stdlib_lsame(side,'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('strmm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*a*b.
     
                   if (upper) then
                       do 50 j = 1,n
                           do 40 k = 1,m
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   do 30 i = 1,k - 1
                                       b(i,j) = b(i,j) + temp*a(i,k)
        30                         continue
                                   if (nounit) temp = temp*a(k,k)
                                   b(k,j) = temp
                               end if
        40                 continue
        50             continue
                   else
                       do 80 j = 1,n
                           do 70 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   temp = alpha*b(k,j)
                                   b(k,j) = temp
                                   if (nounit) b(k,j) = b(k,j)*a(k,k)
                                   do 60 i = k + 1,m
                                       b(i,j) = b(i,j) + temp*a(i,k)
        60                         continue
                               end if
        70                 continue
        80             continue
                   end if
               else
     
                 ! form  b := alpha*a**t*b.
     
                   if (upper) then
                       do 110 j = 1,n
                           do 100 i = m,1,-1
                               temp = b(i,j)
                               if (nounit) temp = temp*a(i,i)
                               do 90 k = 1,i - 1
                                   temp = temp + a(k,i)*b(k,j)
        90                     continue
                               b(i,j) = alpha*temp
       100                 continue
       110             continue
                   else
                       do 140 j = 1,n
                           do 130 i = 1,m
                               temp = b(i,j)
                               if (nounit) temp = temp*a(i,i)
                               do 120 k = i + 1,m
                                   temp = temp + a(k,i)*b(k,j)
       120                     continue
                               b(i,j) = alpha*temp
       130                 continue
       140             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*a.
     
                   if (upper) then
                       do 180 j = n,1,-1
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 150 i = 1,m
                               b(i,j) = temp*b(i,j)
       150                 continue
                           do 170 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 160 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       160                         continue
                               end if
       170                 continue
       180             continue
                   else
                       do 220 j = 1,n
                           temp = alpha
                           if (nounit) temp = temp*a(j,j)
                           do 190 i = 1,m
                               b(i,j) = temp*b(i,j)
       190                 continue
                           do 210 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   temp = alpha*a(k,j)
                                   do 200 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       200                         continue
                               end if
       210                 continue
       220             continue
                   end if
               else
     
                 ! form  b := alpha*b*a**t.
     
                   if (upper) then
                       do 260 k = 1,n
                           do 240 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   temp = alpha*a(j,k)
                                   do 230 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       230                         continue
                               end if
       240                 continue
                           temp = alpha
                           if (nounit) temp = temp*a(k,k)
                           if (temp/=one) then
                               do 250 i = 1,m
                                   b(i,k) = temp*b(i,k)
       250                     continue
                           end if
       260             continue
                   else
                       do 300 k = n,1,-1
                           do 280 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   temp = alpha*a(j,k)
                                   do 270 i = 1,m
                                       b(i,j) = b(i,j) + temp*b(i,k)
       270                         continue
                               end if
       280                 continue
                           temp = alpha
                           if (nounit) temp = temp*a(k,k)
                           if (temp/=one) then
                               do 290 i = 1,m
                                   b(i,k) = temp*b(i,k)
       290                     continue
                           end if
       300             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_strmm
     
     end subroutine stdlib_strmm

     
     
     subroutine stdlib_strmv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('strmv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = 1,n
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 10 i = 1,j - 1
                                   x(i) = x(i) + temp*a(i,j)
        10                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        20             continue
                   else
                       jx = kx
                       do 40 j = 1,n
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 30 i = 1,j - 1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix + incx
        30                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx + incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = n,1,-1
                           if (x(j)/=zero) then
                               temp = x(j)
                               do 50 i = n,j + 1,-1
                                   x(i) = x(i) + temp*a(i,j)
        50                     continue
                               if (nounit) x(j) = x(j)*a(j,j)
                           end if
        60             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 80 j = n,1,-1
                           if (x(jx)/=zero) then
                               temp = x(jx)
                               ix = kx
                               do 70 i = n,j + 1,-1
                                   x(ix) = x(ix) + temp*a(i,j)
                                   ix = ix - incx
        70                     continue
                               if (nounit) x(jx) = x(jx)*a(j,j)
                           end if
                           jx = jx - incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := a**t*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 100 j = n,1,-1
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do 90 i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
        90                 continue
                           x(j) = temp
       100             continue
                   else
                       jx = kx + (n-1)*incx
                       do 120 j = n,1,-1
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do 110 i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
       110                 continue
                           x(jx) = temp
                           jx = jx - incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = 1,n
                           temp = x(j)
                           if (nounit) temp = temp*a(j,j)
                           do 130 i = j + 1,n
                               temp = temp + a(i,j)*x(i)
       130                 continue
                           x(j) = temp
       140             continue
                   else
                       jx = kx
                       do 160 j = 1,n
                           temp = x(jx)
                           ix = jx
                           if (nounit) temp = temp*a(j,j)
                           do 150 i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
       150                 continue
                           x(jx) = temp
                           jx = jx + incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_strmv
     
     end subroutine stdlib_strmv

     
     
     subroutine stdlib_strsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
     
        ! -- reference blas level3 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           real(sp) alpha
           integer(int32) lda,ldb,m,n
           character diag,side,transa,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),b(ldb,*)
           ! ..
     
        ! =====================================================================
     
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,j,k,nrowa
           logical(lk) lside,nounit,upper
           ! ..
           ! .. parameters ..
           real(sp) one,zero
           parameter (one=1.0_sp,zero=0.0_sp)
           ! ..
     
           ! test the input parameters.
     
           lside = stdlib_lsame(side,'l')
           if (lside) then
               nrowa = m
           else
               nrowa = n
           end if
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
           else if ((.not.stdlib_lsame(diag,'u')) .and. (.not.stdlib_lsame(diag,'n')))&
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
               call stdlib_xerbla('strsm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (m==0 .or. n==0) return
     
           ! and when  alpha.eq.zero.
     
           if (alpha==zero) then
               do 20 j = 1,n
                   do 10 i = 1,m
                       b(i,j) = zero
        10         continue
        20     continue
               return
           end if
     
           ! start the operations.
     
           if (lside) then
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*inv( a )*b.
     
                   if (upper) then
                       do 60 j = 1,n
                           if (alpha/=one) then
                               do 30 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        30                     continue
                           end if
                           do 50 k = m,1,-1
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 40 i = 1,k - 1
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        40                         continue
                               end if
        50                 continue
        60             continue
                   else
                       do 100 j = 1,n
                           if (alpha/=one) then
                               do 70 i = 1,m
                                   b(i,j) = alpha*b(i,j)
        70                     continue
                           end if
                           do 90 k = 1,m
                               if (b(k,j)/=zero) then
                                   if (nounit) b(k,j) = b(k,j)/a(k,k)
                                   do 80 i = k + 1,m
                                       b(i,j) = b(i,j) - b(k,j)*a(i,k)
        80                         continue
                               end if
        90                 continue
       100             continue
                   end if
               else
     
                 ! form  b := alpha*inv( a**t )*b.
     
                   if (upper) then
                       do 130 j = 1,n
                           do 120 i = 1,m
                               temp = alpha*b(i,j)
                               do 110 k = 1,i - 1
                                   temp = temp - a(k,i)*b(k,j)
       110                     continue
                               if (nounit) temp = temp/a(i,i)
                               b(i,j) = temp
       120                 continue
       130             continue
                   else
                       do 160 j = 1,n
                           do 150 i = m,1,-1
                               temp = alpha*b(i,j)
                               do 140 k = i + 1,m
                                   temp = temp - a(k,i)*b(k,j)
       140                     continue
                               if (nounit) temp = temp/a(i,i)
                               b(i,j) = temp
       150                 continue
       160             continue
                   end if
               end if
           else
               if (stdlib_lsame(transa,'n')) then
     
                 ! form  b := alpha*b*inv( a ).
     
                   if (upper) then
                       do 210 j = 1,n
                           if (alpha/=one) then
                               do 170 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       170                     continue
                           end if
                           do 190 k = 1,j - 1
                               if (a(k,j)/=zero) then
                                   do 180 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       180                         continue
                               end if
       190                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 200 i = 1,m
                                   b(i,j) = temp*b(i,j)
       200                     continue
                           end if
       210             continue
                   else
                       do 260 j = n,1,-1
                           if (alpha/=one) then
                               do 220 i = 1,m
                                   b(i,j) = alpha*b(i,j)
       220                     continue
                           end if
                           do 240 k = j + 1,n
                               if (a(k,j)/=zero) then
                                   do 230 i = 1,m
                                       b(i,j) = b(i,j) - a(k,j)*b(i,k)
       230                         continue
                               end if
       240                 continue
                           if (nounit) then
                               temp = one/a(j,j)
                               do 250 i = 1,m
                                   b(i,j) = temp*b(i,j)
       250                     continue
                           end if
       260             continue
                   end if
               else
     
                 ! form  b := alpha*b*inv( a**t ).
     
                   if (upper) then
                       do 310 k = n,1,-1
                           if (nounit) then
                               temp = one/a(k,k)
                               do 270 i = 1,m
                                   b(i,k) = temp*b(i,k)
       270                     continue
                           end if
                           do 290 j = 1,k - 1
                               if (a(j,k)/=zero) then
                                   temp = a(j,k)
                                   do 280 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       280                         continue
                               end if
       290                 continue
                           if (alpha/=one) then
                               do 300 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       300                     continue
                           end if
       310             continue
                   else
                       do 360 k = 1,n
                           if (nounit) then
                               temp = one/a(k,k)
                               do 320 i = 1,m
                                   b(i,k) = temp*b(i,k)
       320                     continue
                           end if
                           do 340 j = k + 1,n
                               if (a(j,k)/=zero) then
                                   temp = a(j,k)
                                   do 330 i = 1,m
                                       b(i,j) = b(i,j) - temp*b(i,k)
       330                         continue
                               end if
       340                 continue
                           if (alpha/=one) then
                               do 350 i = 1,m
                                   b(i,k) = alpha*b(i,k)
       350                     continue
                           end if
       360             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_strsm
     
     end subroutine stdlib_strsm

     
     
     subroutine stdlib_strsv(uplo,trans,diag,n,a,lda,x,incx)
     
        ! -- reference blas level2 routine --
        ! -- reference blas is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
     
           ! .. scalar arguments ..
           integer(int32) incx,lda,n
           character diag,trans,uplo
           ! ..
           ! .. array arguments ..
           real(sp) a(lda,*),x(*)
           ! ..
     
        ! =====================================================================
     
           ! .. parameters ..
           real(sp) zero
           parameter (zero=0.0_sp)
           ! ..
           ! .. local scalars ..
           real(sp) temp
           integer(int32) i,info,ix,j,jx,kx
           logical(lk) nounit
           ! ..
     
     
     
           ! .. intrinsic functions ..
           intrinsic max
           ! ..
     
           ! test the input parameters.
     
           info = 0
           if (.not.stdlib_lsame(uplo,'u') .and. .not.stdlib_lsame(uplo,'l')) then
               info = 1
           else if (.not.stdlib_lsame(trans,'n') .and. .not.stdlib_lsame(trans,'t')&
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
               call stdlib_xerbla('strsv ',info)
               return
           end if
     
           ! quick return if possible.
     
           if (n==0) return
     
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
                       do 20 j = n,1,-1
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 10 i = j - 1,1,-1
                                   x(i) = x(i) - temp*a(i,j)
        10                     continue
                           end if
        20             continue
                   else
                       jx = kx + (n-1)*incx
                       do 40 j = n,1,-1
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 30 i = j - 1,1,-1
                                   ix = ix - incx
                                   x(ix) = x(ix) - temp*a(i,j)
        30                     continue
                           end if
                           jx = jx - incx
        40             continue
                   end if
               else
                   if (incx==1) then
                       do 60 j = 1,n
                           if (x(j)/=zero) then
                               if (nounit) x(j) = x(j)/a(j,j)
                               temp = x(j)
                               do 50 i = j + 1,n
                                   x(i) = x(i) - temp*a(i,j)
        50                     continue
                           end if
        60             continue
                   else
                       jx = kx
                       do 80 j = 1,n
                           if (x(jx)/=zero) then
                               if (nounit) x(jx) = x(jx)/a(j,j)
                               temp = x(jx)
                               ix = jx
                               do 70 i = j + 1,n
                                   ix = ix + incx
                                   x(ix) = x(ix) - temp*a(i,j)
        70                     continue
                           end if
                           jx = jx + incx
        80             continue
                   end if
               end if
           else
     
              ! form  x := inv( a**t )*x.
     
               if (stdlib_lsame(uplo,'u')) then
                   if (incx==1) then
                       do 100 j = 1,n
                           temp = x(j)
                           do 90 i = 1,j - 1
                               temp = temp - a(i,j)*x(i)
        90                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
       100             continue
                   else
                       jx = kx
                       do 120 j = 1,n
                           temp = x(jx)
                           ix = kx
                           do 110 i = 1,j - 1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix + incx
       110                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx + incx
       120             continue
                   end if
               else
                   if (incx==1) then
                       do 140 j = n,1,-1
                           temp = x(j)
                           do 130 i = n,j + 1,-1
                               temp = temp - a(i,j)*x(i)
       130                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(j) = temp
       140             continue
                   else
                       kx = kx + (n-1)*incx
                       jx = kx
                       do 160 j = n,1,-1
                           temp = x(jx)
                           ix = kx
                           do 150 i = n,j + 1,-1
                               temp = temp - a(i,j)*x(ix)
                               ix = ix - incx
       150                 continue
                           if (nounit) temp = temp/a(j,j)
                           x(jx) = temp
                           jx = jx - incx
       160             continue
                   end if
               end if
           end if
     
           return
     
           ! end of stdlib_strsv
     
     end subroutine stdlib_strsv




end module stdlib_linalg_blas
