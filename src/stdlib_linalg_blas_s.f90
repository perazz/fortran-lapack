module stdlib_linalg_blas_s
     use stdlib_linalg_constants
     use stdlib_linalg_blas_aux
     implicit none(type,external)
     private






     public :: sp,dp,lk,int32,int64
     public :: stdlib_sasum
     public :: stdlib_saxpy
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


     contains
     
     
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
     real(sp), parameter :: tbig = real(radix(0._sp), wp)**floor(       (maxexponent(0._sp) - &
               digits(0._sp) + 1) * 0.5_sp)
     real(sp), parameter :: ssml = real(radix(0._sp), wp)**( - floor(       (minexponent(0._sp) - &
               digits(0._sp)) * 0.5_sp))
     real(sp), parameter :: sbig = real(radix(0._sp), wp)**( - ceiling(       (maxexponent(0._sp) &
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
               stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + &
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
               call stdlib_xerbla('stdlib_sgbmv ',info)
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
               call stdlib_xerbla('stdlib_sgemm ',info)
               return
           end if
     
           ! quick return if possible.
     
           if ((m==0) .or. (n==0) .or.(((alpha==zero).or. (k==0)).and. (beta==one))) &
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
               call stdlib_xerbla('stdlib_sgemv ',info)
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
               call stdlib_xerbla('stdlib_sger  ',info)
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
     real(sp), parameter :: tbig = real(radix(0._sp), wp)**floor(       (maxexponent(0._sp) - &
               digits(0._sp) + 1) * 0.5_sp)
     real(sp), parameter :: ssml = real(radix(0._sp), wp)**( - floor(       (minexponent(0._sp) - &
               digits(0._sp)) * 0.5_sp))
     real(sp), parameter :: sbig = real(radix(0._sp), wp)**( - ceiling(       (maxexponent(0._sp) &
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
               call stdlib_xerbla('stdlib_ssbmv ',info)
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
               call stdlib_xerbla('stdlib_sspmv ',info)
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
               call stdlib_xerbla('stdlib_sspr  ',info)
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
               call stdlib_xerbla('stdlib_sspr2 ',info)
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
               call stdlib_xerbla('stdlib_ssymm ',info)
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
               call stdlib_xerbla('stdlib_ssymv ',info)
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
               call stdlib_xerbla('stdlib_ssyr  ',info)
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
               call stdlib_xerbla('stdlib_ssyr2 ',info)
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
               call stdlib_xerbla('stdlib_ssyr2k',info)
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
               call stdlib_xerbla('stdlib_ssyrk ',info)
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
               call stdlib_xerbla('stdlib_stbmv ',info)
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
               call stdlib_xerbla('stdlib_stbsv ',info)
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
               call stdlib_xerbla('stdlib_stpmv ',info)
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
               call stdlib_xerbla('stdlib_stpsv ',info)
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
               call stdlib_xerbla('stdlib_strmm ',info)
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
               call stdlib_xerbla('stdlib_strmv ',info)
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
               call stdlib_xerbla('stdlib_strsm ',info)
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
               call stdlib_xerbla('stdlib_strsv ',info)
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



end module stdlib_linalg_blas_s
