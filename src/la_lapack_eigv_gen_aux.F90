!> Nonsymmetric eigenproblem helpers: 2-by-2 standardization, Sylvester solves, diagonal block swaps
module la_lapack_eigv_gen_aux
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_mnorm
     use la_lapack_blas_like_scalar
     use la_lapack_givens_jacobi_rot
     use la_lapack_householder_reflectors
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slasy2
     public :: la_slaln2
     public :: la_slanv2
     public :: la_slaexc
     public :: la_strexc
     public :: la_dlasy2
     public :: la_dlaln2
     public :: la_dlanv2
     public :: la_dlaexc
     public :: la_dtrexc
#ifdef LA_WITH_XDP
     public :: la_xlasy2
     public :: la_xlaln2
     public :: la_xlanv2
     public :: la_xlaexc
     public :: la_xtrexc
#endif
#ifdef LA_WITH_QP
     public :: la_qlasy2
     public :: la_qlaln2
     public :: la_qlanv2
     public :: la_qlaexc
     public :: la_qtrexc
#endif
     public :: la_ctrexc
     public :: la_ztrexc
#ifdef LA_WITH_XDP
     public :: la_ytrexc
#endif
#ifdef LA_WITH_QP
     public :: la_wtrexc
#endif

     contains

     !> SLASY2: solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
     !> op(TL)*X + ISGN*X*op(TR) = SCALE*B,
     !> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
     !> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.

     pure subroutine la_slasy2(ltranl,ltranr,isgn,n1,n2,tl,ldtl,tr,ldtr,b,ldb, &
               scale,x,ldx,xnorm,info)
        use la_constants_sp,only:zero,half,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltranl,ltranr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,ldb,ldtl,ldtr,ldx,n1,n2
           real(sp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(sp),intent(in) :: b(ldb,*),tl(ldtl,*),tr(ldtr,*)
           real(sp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: bswap,xswap
           integer(ilp) :: i,ip,ipiv,ipsv,j,jp,jpsv,k
           real(sp) :: bet,eps,gam,l21,sgn,smin,smlnum,tau1,temp,u11,u12,u22, &
                     xmax
           ! Local Arrays
           logical(lk) :: bswpiv(4),xswpiv(4)
           integer(ilp) :: jpiv(4),locl21(4),locu12(4),locu22(4)
           real(sp) :: btmp(4),t16(4,4),tmp(4),x2(2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Data Statements
           locu12 = [3,4,1,2]
           locl21 = [2,1,4,3]
           locu22 = [4,3,2,1]
           xswpiv = [.false.,.false.,.true.,.true.]
           bswpiv = [.false.,.true.,.false.,.true.]
           ! Executable Statements
           ! do not check the input parameters for errors
           info = 0
           ! quick return if possible
           if (n1 == 0 .or. n2 == 0) return
           ! set constants to control overflow
           eps = la_slamch('P')
           smlnum = la_slamch('S')/eps
           sgn = isgn
           k = n1 + n1 + n2 - 2
           go to(10,20,30,50) k
           ! 1 by 1: tl11*x + sgn*x*tr11 = b11
           10 continue
           tau1 = tl(1,1) + sgn*tr(1,1)
           bet = abs(tau1)
           if (bet <= smlnum) then
              tau1 = smlnum
              bet = smlnum
              info = 1
           end if
           scale = one
           gam = abs(b(1,1))
           if (smlnum*gam > bet) scale = one/gam
           x(1,1) = (b(1,1)*scale)/tau1
           xnorm = abs(x(1,1))
           return
           ! 1 by 2:
           ! tl11*[x11 x12] + isgn*[x11 x12]*op[tr11 tr12]  = [b11 b12]
                                             ! [tr21 tr22]
                                             20 continue
           smin = max(eps*max(abs(tl(1,1)),abs(tr(1,1)),abs(tr(1,2)),abs(tr( &
                     2,1)),abs(tr(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(1,1) + sgn*tr(2,2)
           if (ltranr) then
              tmp(2) = sgn*tr(2,1)
              tmp(3) = sgn*tr(1,2)
           else
              tmp(2) = sgn*tr(1,2)
              tmp(3) = sgn*tr(2,1)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(1,2)
           go to 40
           ! 2 by 1:
                ! op[tl11 tl12]*[x11] + isgn* [x11]*tr11  = [b11]
                  ! [tl21 tl22] [x21]         [x21]         [b21]
                  30 continue
           smin = max(eps*max(abs(tr(1,1)),abs(tl(1,1)),abs(tl(1,2)),abs(tl( &
                     2,1)),abs(tl(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(2,2) + sgn*tr(1,1)
           if (ltranl) then
              tmp(2) = tl(1,2)
              tmp(3) = tl(2,1)
           else
              tmp(2) = tl(2,1)
              tmp(3) = tl(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           40 continue
           ! solve 2 by 2 system using complete pivoting.
           ! set pivots less than smin to smin.
           ipiv = la_isamax(4,tmp,1)
           u11 = tmp(ipiv)
           if (abs(u11) <= smin) then
              info = 1
              u11 = smin
           end if
           u12 = tmp(locu12(ipiv))
           l21 = tmp(locl21(ipiv))/u11
           u22 = tmp(locu22(ipiv)) - u12*l21
           xswap = xswpiv(ipiv)
           bswap = bswpiv(ipiv)
           if (abs(u22) <= smin) then
              info = 1
              u22 = smin
           end if
           if (bswap) then
              temp = btmp(2)
              btmp(2) = btmp(1) - l21*temp
              btmp(1) = temp
           else
              btmp(2) = btmp(2) - l21*btmp(1)
           end if
           scale = one
           if ((two*smlnum)*abs(btmp(2)) > abs(u22) .or. (two*smlnum)*abs(btmp(1)) > abs( &
                      u11)) then
              scale = half/max(abs(btmp(1)),abs(btmp(2)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
           end if
           x2(2) = btmp(2)/u22
           x2(1) = btmp(1)/u11 - (u12/u11)*x2(2)
           if (xswap) then
              temp = x2(2)
              x2(2) = x2(1)
              x2(1) = temp
           end if
           x(1,1) = x2(1)
           if (n1 == 1) then
              x(1,2) = x2(2)
              xnorm = abs(x(1,1)) + abs(x(1,2))
           else
              x(2,1) = x2(2)
              xnorm = max(abs(x(1,1)),abs(x(2,1)))
           end if
           return
           ! 2 by 2:
           ! op[tl11 tl12]*[x11 x12] +isgn* [x11 x12]*op[tr11 tr12] = [b11 b12]
             ! [tl21 tl22] [x21 x22]        [x21 x22]   [tr21 tr22]   [b21 b22]
           ! solve equivalent 4 by 4 system using complete pivoting.
           ! set pivots less than smin to smin.
           50 continue
           smin = max(abs(tr(1,1)),abs(tr(1,2)),abs(tr(2,1)),abs(tr(2,2)))

           smin = max(smin,abs(tl(1,1)),abs(tl(1,2)),abs(tl(2,1)),abs(tl(2, &
                     2)))
           smin = max(eps*smin,smlnum)
           btmp(1) = zero
           call la_scopy(16,btmp,0,t16,1)
           t16(1,1) = tl(1,1) + sgn*tr(1,1)
           t16(2,2) = tl(2,2) + sgn*tr(1,1)
           t16(3,3) = tl(1,1) + sgn*tr(2,2)
           t16(4,4) = tl(2,2) + sgn*tr(2,2)
           if (ltranl) then
              t16(1,2) = tl(2,1)
              t16(2,1) = tl(1,2)
              t16(3,4) = tl(2,1)
              t16(4,3) = tl(1,2)
           else
              t16(1,2) = tl(1,2)
              t16(2,1) = tl(2,1)
              t16(3,4) = tl(1,2)
              t16(4,3) = tl(2,1)
           end if
           if (ltranr) then
              t16(1,3) = sgn*tr(1,2)
              t16(2,4) = sgn*tr(1,2)
              t16(3,1) = sgn*tr(2,1)
              t16(4,2) = sgn*tr(2,1)
           else
              t16(1,3) = sgn*tr(2,1)
              t16(2,4) = sgn*tr(2,1)
              t16(3,1) = sgn*tr(1,2)
              t16(4,2) = sgn*tr(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           btmp(3) = b(1,2)
           btmp(4) = b(2,2)
           ! perform elimination
           loop_100: do i = 1,3
              xmax = zero
              do ip = i,4
                 do jp = i,4
                    if (abs(t16(ip,jp)) >= xmax) then
                       xmax = abs(t16(ip,jp))
                       ipsv = ip
                       jpsv = jp
                    end if
                 end do
              end do
              if (ipsv /= i) then
                 call la_sswap(4,t16(ipsv,1),4,t16(i,1),4)
                 temp = btmp(i)
                 btmp(i) = btmp(ipsv)
                 btmp(ipsv) = temp
              end if
              if (jpsv /= i) call la_sswap(4,t16(1,jpsv),1,t16(1,i),1)
              jpiv(i) = jpsv
              if (abs(t16(i,i)) < smin) then
                 info = 1
                 t16(i,i) = smin
              end if
              do j = i + 1,4
                 t16(j,i) = t16(j,i)/t16(i,i)
                 btmp(j) = btmp(j) - t16(j,i)*btmp(i)
                 do k = i + 1,4
                    t16(j,k) = t16(j,k) - t16(j,i)*t16(i,k)
                 end do
              end do
           end do loop_100
           if (abs(t16(4,4)) < smin) then
              info = 1
              t16(4,4) = smin
           end if
           scale = one
           if ((eight*smlnum)*abs(btmp(1)) > abs(t16(1,1)) .or. (eight*smlnum)*abs( &
           btmp(2)) > abs(t16(2,2)) .or. (eight*smlnum)*abs(btmp(3)) > abs(t16(3,3)) &
                      .or. (eight*smlnum)*abs(btmp(4)) > abs(t16(4,4))) then
              scale = (one/eight)/max(abs(btmp(1)),abs(btmp(2)),abs(btmp(3)), &
                        abs(btmp(4)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
              btmp(3) = btmp(3)*scale
              btmp(4) = btmp(4)*scale
           end if
           do i = 1,4
              k = 5 - i
              temp = one/t16(k,k)
              tmp(k) = btmp(k)*temp
              do j = k + 1,4
                 tmp(k) = tmp(k) - (temp*t16(k,j))*tmp(j)
              end do
           end do
           do i = 1,3
              if (jpiv(4 - i) /= 4 - i) then
                 temp = tmp(4 - i)
                 tmp(4 - i) = tmp(jpiv(4 - i))
                 tmp(jpiv(4 - i)) = temp
              end if
           end do
           x(1,1) = tmp(1)
           x(2,1) = tmp(2)
           x(1,2) = tmp(3)
           x(2,2) = tmp(4)
           xnorm = max(abs(tmp(1)) + abs(tmp(3)),abs(tmp(2)) + abs(tmp(4)))
           return
     end subroutine la_slasy2
     !> DLASY2: solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
     !> op(TL)*X + ISGN*X*op(TR) = SCALE*B,
     !> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
     !> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.

     pure subroutine la_dlasy2(ltranl,ltranr,isgn,n1,n2,tl,ldtl,tr,ldtr,b,ldb, &
               scale,x,ldx,xnorm,info)
        use la_constants_dp,only:zero,half,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltranl,ltranr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,ldb,ldtl,ldtr,ldx,n1,n2
           real(dp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(dp),intent(in) :: b(ldb,*),tl(ldtl,*),tr(ldtr,*)
           real(dp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: bswap,xswap
           integer(ilp) :: i,ip,ipiv,ipsv,j,jp,jpsv,k
           real(dp) :: bet,eps,gam,l21,sgn,smin,smlnum,tau1,temp,u11,u12,u22, &
                     xmax
           ! Local Arrays
           logical(lk) :: bswpiv(4),xswpiv(4)
           integer(ilp) :: jpiv(4),locl21(4),locu12(4),locu22(4)
           real(dp) :: btmp(4),t16(4,4),tmp(4),x2(2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Data Statements
           locu12 = [3,4,1,2]
           locl21 = [2,1,4,3]
           locu22 = [4,3,2,1]
           xswpiv = [.false.,.false.,.true.,.true.]
           bswpiv = [.false.,.true.,.false.,.true.]
           ! Executable Statements
           ! do not check the input parameters for errors
           info = 0
           ! quick return if possible
           if (n1 == 0 .or. n2 == 0) return
           ! set constants to control overflow
           eps = la_dlamch('P')
           smlnum = la_dlamch('S')/eps
           sgn = isgn
           k = n1 + n1 + n2 - 2
           go to(10,20,30,50) k
           ! 1 by 1: tl11*x + sgn*x*tr11 = b11
           10 continue
           tau1 = tl(1,1) + sgn*tr(1,1)
           bet = abs(tau1)
           if (bet <= smlnum) then
              tau1 = smlnum
              bet = smlnum
              info = 1
           end if
           scale = one
           gam = abs(b(1,1))
           if (smlnum*gam > bet) scale = one/gam
           x(1,1) = (b(1,1)*scale)/tau1
           xnorm = abs(x(1,1))
           return
           ! 1 by 2:
           ! tl11*[x11 x12] + isgn*[x11 x12]*op[tr11 tr12]  = [b11 b12]
                                             ! [tr21 tr22]
                                             20 continue
           smin = max(eps*max(abs(tl(1,1)),abs(tr(1,1)),abs(tr(1,2)),abs(tr( &
                     2,1)),abs(tr(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(1,1) + sgn*tr(2,2)
           if (ltranr) then
              tmp(2) = sgn*tr(2,1)
              tmp(3) = sgn*tr(1,2)
           else
              tmp(2) = sgn*tr(1,2)
              tmp(3) = sgn*tr(2,1)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(1,2)
           go to 40
           ! 2 by 1:
                ! op[tl11 tl12]*[x11] + isgn* [x11]*tr11  = [b11]
                  ! [tl21 tl22] [x21]         [x21]         [b21]
                  30 continue
           smin = max(eps*max(abs(tr(1,1)),abs(tl(1,1)),abs(tl(1,2)),abs(tl( &
                     2,1)),abs(tl(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(2,2) + sgn*tr(1,1)
           if (ltranl) then
              tmp(2) = tl(1,2)
              tmp(3) = tl(2,1)
           else
              tmp(2) = tl(2,1)
              tmp(3) = tl(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           40 continue
           ! solve 2 by 2 system using complete pivoting.
           ! set pivots less than smin to smin.
           ipiv = la_idamax(4,tmp,1)
           u11 = tmp(ipiv)
           if (abs(u11) <= smin) then
              info = 1
              u11 = smin
           end if
           u12 = tmp(locu12(ipiv))
           l21 = tmp(locl21(ipiv))/u11
           u22 = tmp(locu22(ipiv)) - u12*l21
           xswap = xswpiv(ipiv)
           bswap = bswpiv(ipiv)
           if (abs(u22) <= smin) then
              info = 1
              u22 = smin
           end if
           if (bswap) then
              temp = btmp(2)
              btmp(2) = btmp(1) - l21*temp
              btmp(1) = temp
           else
              btmp(2) = btmp(2) - l21*btmp(1)
           end if
           scale = one
           if ((two*smlnum)*abs(btmp(2)) > abs(u22) .or. (two*smlnum)*abs(btmp(1)) > abs( &
                      u11)) then
              scale = half/max(abs(btmp(1)),abs(btmp(2)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
           end if
           x2(2) = btmp(2)/u22
           x2(1) = btmp(1)/u11 - (u12/u11)*x2(2)
           if (xswap) then
              temp = x2(2)
              x2(2) = x2(1)
              x2(1) = temp
           end if
           x(1,1) = x2(1)
           if (n1 == 1) then
              x(1,2) = x2(2)
              xnorm = abs(x(1,1)) + abs(x(1,2))
           else
              x(2,1) = x2(2)
              xnorm = max(abs(x(1,1)),abs(x(2,1)))
           end if
           return
           ! 2 by 2:
           ! op[tl11 tl12]*[x11 x12] +isgn* [x11 x12]*op[tr11 tr12] = [b11 b12]
             ! [tl21 tl22] [x21 x22]        [x21 x22]   [tr21 tr22]   [b21 b22]
           ! solve equivalent 4 by 4 system using complete pivoting.
           ! set pivots less than smin to smin.
           50 continue
           smin = max(abs(tr(1,1)),abs(tr(1,2)),abs(tr(2,1)),abs(tr(2,2)))

           smin = max(smin,abs(tl(1,1)),abs(tl(1,2)),abs(tl(2,1)),abs(tl(2, &
                     2)))
           smin = max(eps*smin,smlnum)
           btmp(1) = zero
           call la_dcopy(16,btmp,0,t16,1)
           t16(1,1) = tl(1,1) + sgn*tr(1,1)
           t16(2,2) = tl(2,2) + sgn*tr(1,1)
           t16(3,3) = tl(1,1) + sgn*tr(2,2)
           t16(4,4) = tl(2,2) + sgn*tr(2,2)
           if (ltranl) then
              t16(1,2) = tl(2,1)
              t16(2,1) = tl(1,2)
              t16(3,4) = tl(2,1)
              t16(4,3) = tl(1,2)
           else
              t16(1,2) = tl(1,2)
              t16(2,1) = tl(2,1)
              t16(3,4) = tl(1,2)
              t16(4,3) = tl(2,1)
           end if
           if (ltranr) then
              t16(1,3) = sgn*tr(1,2)
              t16(2,4) = sgn*tr(1,2)
              t16(3,1) = sgn*tr(2,1)
              t16(4,2) = sgn*tr(2,1)
           else
              t16(1,3) = sgn*tr(2,1)
              t16(2,4) = sgn*tr(2,1)
              t16(3,1) = sgn*tr(1,2)
              t16(4,2) = sgn*tr(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           btmp(3) = b(1,2)
           btmp(4) = b(2,2)
           ! perform elimination
           loop_100: do i = 1,3
              xmax = zero
              do ip = i,4
                 do jp = i,4
                    if (abs(t16(ip,jp)) >= xmax) then
                       xmax = abs(t16(ip,jp))
                       ipsv = ip
                       jpsv = jp
                    end if
                 end do
              end do
              if (ipsv /= i) then
                 call la_dswap(4,t16(ipsv,1),4,t16(i,1),4)
                 temp = btmp(i)
                 btmp(i) = btmp(ipsv)
                 btmp(ipsv) = temp
              end if
              if (jpsv /= i) call la_dswap(4,t16(1,jpsv),1,t16(1,i),1)
              jpiv(i) = jpsv
              if (abs(t16(i,i)) < smin) then
                 info = 1
                 t16(i,i) = smin
              end if
              do j = i + 1,4
                 t16(j,i) = t16(j,i)/t16(i,i)
                 btmp(j) = btmp(j) - t16(j,i)*btmp(i)
                 do k = i + 1,4
                    t16(j,k) = t16(j,k) - t16(j,i)*t16(i,k)
                 end do
              end do
           end do loop_100
           if (abs(t16(4,4)) < smin) then
              info = 1
              t16(4,4) = smin
           end if
           scale = one
           if ((eight*smlnum)*abs(btmp(1)) > abs(t16(1,1)) .or. (eight*smlnum)*abs( &
           btmp(2)) > abs(t16(2,2)) .or. (eight*smlnum)*abs(btmp(3)) > abs(t16(3,3)) &
                      .or. (eight*smlnum)*abs(btmp(4)) > abs(t16(4,4))) then
              scale = (one/eight)/max(abs(btmp(1)),abs(btmp(2)),abs(btmp(3)), &
                        abs(btmp(4)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
              btmp(3) = btmp(3)*scale
              btmp(4) = btmp(4)*scale
           end if
           do i = 1,4
              k = 5 - i
              temp = one/t16(k,k)
              tmp(k) = btmp(k)*temp
              do j = k + 1,4
                 tmp(k) = tmp(k) - (temp*t16(k,j))*tmp(j)
              end do
           end do
           do i = 1,3
              if (jpiv(4 - i) /= 4 - i) then
                 temp = tmp(4 - i)
                 tmp(4 - i) = tmp(jpiv(4 - i))
                 tmp(jpiv(4 - i)) = temp
              end if
           end do
           x(1,1) = tmp(1)
           x(2,1) = tmp(2)
           x(1,2) = tmp(3)
           x(2,2) = tmp(4)
           xnorm = max(abs(tmp(1)) + abs(tmp(3)),abs(tmp(2)) + abs(tmp(4)))
           return
     end subroutine la_dlasy2
#ifdef LA_WITH_XDP
     !> XLASY2: solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
     !> op(TL)*X + ISGN*X*op(TR) = SCALE*B,
     !> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
     !> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.

     pure subroutine la_xlasy2(ltranl,ltranr,isgn,n1,n2,tl,ldtl,tr,ldtr,b,ldb, &
               scale,x,ldx,xnorm,info)
        use la_constants_xdp,only:zero,half,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltranl,ltranr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,ldb,ldtl,ldtr,ldx,n1,n2
           real(xdp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(xdp),intent(in) :: b(ldb,*),tl(ldtl,*),tr(ldtr,*)
           real(xdp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: bswap,xswap
           integer(ilp) :: i,ip,ipiv,ipsv,j,jp,jpsv,k
           real(xdp) :: bet,eps,gam,l21,sgn,smin,smlnum,tau1,temp,u11,u12,u22, &
                     xmax
           ! Local Arrays
           logical(lk) :: bswpiv(4),xswpiv(4)
           integer(ilp) :: jpiv(4),locl21(4),locu12(4),locu22(4)
           real(xdp) :: btmp(4),t16(4,4),tmp(4),x2(2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Data Statements
           locu12 = [3,4,1,2]
           locl21 = [2,1,4,3]
           locu22 = [4,3,2,1]
           xswpiv = [.false.,.false.,.true.,.true.]
           bswpiv = [.false.,.true.,.false.,.true.]
           ! Executable Statements
           ! do not check the input parameters for errors
           info = 0
           ! quick return if possible
           if (n1 == 0 .or. n2 == 0) return
           ! set constants to control overflow
           eps = la_xlamch('P')
           smlnum = la_xlamch('S')/eps
           sgn = isgn
           k = n1 + n1 + n2 - 2
           go to(10,20,30,50) k
           ! 1 by 1: tl11*x + sgn*x*tr11 = b11
           10 continue
           tau1 = tl(1,1) + sgn*tr(1,1)
           bet = abs(tau1)
           if (bet <= smlnum) then
              tau1 = smlnum
              bet = smlnum
              info = 1
           end if
           scale = one
           gam = abs(b(1,1))
           if (smlnum*gam > bet) scale = one/gam
           x(1,1) = (b(1,1)*scale)/tau1
           xnorm = abs(x(1,1))
           return
           ! 1 by 2:
           ! tl11*[x11 x12] + isgn*[x11 x12]*op[tr11 tr12]  = [b11 b12]
                                             ! [tr21 tr22]
                                             20 continue
           smin = max(eps*max(abs(tl(1,1)),abs(tr(1,1)),abs(tr(1,2)),abs(tr( &
                     2,1)),abs(tr(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(1,1) + sgn*tr(2,2)
           if (ltranr) then
              tmp(2) = sgn*tr(2,1)
              tmp(3) = sgn*tr(1,2)
           else
              tmp(2) = sgn*tr(1,2)
              tmp(3) = sgn*tr(2,1)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(1,2)
           go to 40
           ! 2 by 1:
                ! op[tl11 tl12]*[x11] + isgn* [x11]*tr11  = [b11]
                  ! [tl21 tl22] [x21]         [x21]         [b21]
                  30 continue
           smin = max(eps*max(abs(tr(1,1)),abs(tl(1,1)),abs(tl(1,2)),abs(tl( &
                     2,1)),abs(tl(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(2,2) + sgn*tr(1,1)
           if (ltranl) then
              tmp(2) = tl(1,2)
              tmp(3) = tl(2,1)
           else
              tmp(2) = tl(2,1)
              tmp(3) = tl(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           40 continue
           ! solve 2 by 2 system using complete pivoting.
           ! set pivots less than smin to smin.
           ipiv = la_ixamax(4,tmp,1)
           u11 = tmp(ipiv)
           if (abs(u11) <= smin) then
              info = 1
              u11 = smin
           end if
           u12 = tmp(locu12(ipiv))
           l21 = tmp(locl21(ipiv))/u11
           u22 = tmp(locu22(ipiv)) - u12*l21
           xswap = xswpiv(ipiv)
           bswap = bswpiv(ipiv)
           if (abs(u22) <= smin) then
              info = 1
              u22 = smin
           end if
           if (bswap) then
              temp = btmp(2)
              btmp(2) = btmp(1) - l21*temp
              btmp(1) = temp
           else
              btmp(2) = btmp(2) - l21*btmp(1)
           end if
           scale = one
           if ((two*smlnum)*abs(btmp(2)) > abs(u22) .or. (two*smlnum)*abs(btmp(1)) > abs( &
                      u11)) then
              scale = half/max(abs(btmp(1)),abs(btmp(2)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
           end if
           x2(2) = btmp(2)/u22
           x2(1) = btmp(1)/u11 - (u12/u11)*x2(2)
           if (xswap) then
              temp = x2(2)
              x2(2) = x2(1)
              x2(1) = temp
           end if
           x(1,1) = x2(1)
           if (n1 == 1) then
              x(1,2) = x2(2)
              xnorm = abs(x(1,1)) + abs(x(1,2))
           else
              x(2,1) = x2(2)
              xnorm = max(abs(x(1,1)),abs(x(2,1)))
           end if
           return
           ! 2 by 2:
           ! op[tl11 tl12]*[x11 x12] +isgn* [x11 x12]*op[tr11 tr12] = [b11 b12]
             ! [tl21 tl22] [x21 x22]        [x21 x22]   [tr21 tr22]   [b21 b22]
           ! solve equivalent 4 by 4 system using complete pivoting.
           ! set pivots less than smin to smin.
           50 continue
           smin = max(abs(tr(1,1)),abs(tr(1,2)),abs(tr(2,1)),abs(tr(2,2)))

           smin = max(smin,abs(tl(1,1)),abs(tl(1,2)),abs(tl(2,1)),abs(tl(2, &
                     2)))
           smin = max(eps*smin,smlnum)
           btmp(1) = zero
           call la_xcopy(16,btmp,0,t16,1)
           t16(1,1) = tl(1,1) + sgn*tr(1,1)
           t16(2,2) = tl(2,2) + sgn*tr(1,1)
           t16(3,3) = tl(1,1) + sgn*tr(2,2)
           t16(4,4) = tl(2,2) + sgn*tr(2,2)
           if (ltranl) then
              t16(1,2) = tl(2,1)
              t16(2,1) = tl(1,2)
              t16(3,4) = tl(2,1)
              t16(4,3) = tl(1,2)
           else
              t16(1,2) = tl(1,2)
              t16(2,1) = tl(2,1)
              t16(3,4) = tl(1,2)
              t16(4,3) = tl(2,1)
           end if
           if (ltranr) then
              t16(1,3) = sgn*tr(1,2)
              t16(2,4) = sgn*tr(1,2)
              t16(3,1) = sgn*tr(2,1)
              t16(4,2) = sgn*tr(2,1)
           else
              t16(1,3) = sgn*tr(2,1)
              t16(2,4) = sgn*tr(2,1)
              t16(3,1) = sgn*tr(1,2)
              t16(4,2) = sgn*tr(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           btmp(3) = b(1,2)
           btmp(4) = b(2,2)
           ! perform elimination
           loop_100: do i = 1,3
              xmax = zero
              do ip = i,4
                 do jp = i,4
                    if (abs(t16(ip,jp)) >= xmax) then
                       xmax = abs(t16(ip,jp))
                       ipsv = ip
                       jpsv = jp
                    end if
                 end do
              end do
              if (ipsv /= i) then
                 call la_xswap(4,t16(ipsv,1),4,t16(i,1),4)
                 temp = btmp(i)
                 btmp(i) = btmp(ipsv)
                 btmp(ipsv) = temp
              end if
              if (jpsv /= i) call la_xswap(4,t16(1,jpsv),1,t16(1,i),1)
              jpiv(i) = jpsv
              if (abs(t16(i,i)) < smin) then
                 info = 1
                 t16(i,i) = smin
              end if
              do j = i + 1,4
                 t16(j,i) = t16(j,i)/t16(i,i)
                 btmp(j) = btmp(j) - t16(j,i)*btmp(i)
                 do k = i + 1,4
                    t16(j,k) = t16(j,k) - t16(j,i)*t16(i,k)
                 end do
              end do
           end do loop_100
           if (abs(t16(4,4)) < smin) then
              info = 1
              t16(4,4) = smin
           end if
           scale = one
           if ((eight*smlnum)*abs(btmp(1)) > abs(t16(1,1)) .or. (eight*smlnum)*abs( &
           btmp(2)) > abs(t16(2,2)) .or. (eight*smlnum)*abs(btmp(3)) > abs(t16(3,3)) &
                      .or. (eight*smlnum)*abs(btmp(4)) > abs(t16(4,4))) then
              scale = (one/eight)/max(abs(btmp(1)),abs(btmp(2)),abs(btmp(3)), &
                        abs(btmp(4)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
              btmp(3) = btmp(3)*scale
              btmp(4) = btmp(4)*scale
           end if
           do i = 1,4
              k = 5 - i
              temp = one/t16(k,k)
              tmp(k) = btmp(k)*temp
              do j = k + 1,4
                 tmp(k) = tmp(k) - (temp*t16(k,j))*tmp(j)
              end do
           end do
           do i = 1,3
              if (jpiv(4 - i) /= 4 - i) then
                 temp = tmp(4 - i)
                 tmp(4 - i) = tmp(jpiv(4 - i))
                 tmp(jpiv(4 - i)) = temp
              end if
           end do
           x(1,1) = tmp(1)
           x(2,1) = tmp(2)
           x(1,2) = tmp(3)
           x(2,2) = tmp(4)
           xnorm = max(abs(tmp(1)) + abs(tmp(3)),abs(tmp(2)) + abs(tmp(4)))
           return
     end subroutine la_xlasy2
#endif
#ifdef LA_WITH_QP
     !> QLASY2: solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
     !> op(TL)*X + ISGN*X*op(TR) = SCALE*B,
     !> where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
     !> -1.  op(T) = T or T**T, where T**T denotes the transpose of T.

     pure subroutine la_qlasy2(ltranl,ltranr,isgn,n1,n2,tl,ldtl,tr,ldtr,b,ldb, &
               scale,x,ldx,xnorm,info)
        use la_constants_qp,only:zero,half,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltranl,ltranr
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: isgn,ldb,ldtl,ldtr,ldx,n1,n2
           real(qp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(qp),intent(in) :: b(ldb,*),tl(ldtl,*),tr(ldtr,*)
           real(qp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: bswap,xswap
           integer(ilp) :: i,ip,ipiv,ipsv,j,jp,jpsv,k
           real(qp) :: bet,eps,gam,l21,sgn,smin,smlnum,tau1,temp,u11,u12,u22, &
                     xmax
           ! Local Arrays
           logical(lk) :: bswpiv(4),xswpiv(4)
           integer(ilp) :: jpiv(4),locl21(4),locu12(4),locu22(4)
           real(qp) :: btmp(4),t16(4,4),tmp(4),x2(2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Data Statements
           locu12 = [3,4,1,2]
           locl21 = [2,1,4,3]
           locu22 = [4,3,2,1]
           xswpiv = [.false.,.false.,.true.,.true.]
           bswpiv = [.false.,.true.,.false.,.true.]
           ! Executable Statements
           ! do not check the input parameters for errors
           info = 0
           ! quick return if possible
           if (n1 == 0 .or. n2 == 0) return
           ! set constants to control overflow
           eps = la_qlamch('P')
           smlnum = la_qlamch('S')/eps
           sgn = isgn
           k = n1 + n1 + n2 - 2
           go to(10,20,30,50) k
           ! 1 by 1: tl11*x + sgn*x*tr11 = b11
           10 continue
           tau1 = tl(1,1) + sgn*tr(1,1)
           bet = abs(tau1)
           if (bet <= smlnum) then
              tau1 = smlnum
              bet = smlnum
              info = 1
           end if
           scale = one
           gam = abs(b(1,1))
           if (smlnum*gam > bet) scale = one/gam
           x(1,1) = (b(1,1)*scale)/tau1
           xnorm = abs(x(1,1))
           return
           ! 1 by 2:
           ! tl11*[x11 x12] + isgn*[x11 x12]*op[tr11 tr12]  = [b11 b12]
                                             ! [tr21 tr22]
                                             20 continue
           smin = max(eps*max(abs(tl(1,1)),abs(tr(1,1)),abs(tr(1,2)),abs(tr( &
                     2,1)),abs(tr(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(1,1) + sgn*tr(2,2)
           if (ltranr) then
              tmp(2) = sgn*tr(2,1)
              tmp(3) = sgn*tr(1,2)
           else
              tmp(2) = sgn*tr(1,2)
              tmp(3) = sgn*tr(2,1)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(1,2)
           go to 40
           ! 2 by 1:
                ! op[tl11 tl12]*[x11] + isgn* [x11]*tr11  = [b11]
                  ! [tl21 tl22] [x21]         [x21]         [b21]
                  30 continue
           smin = max(eps*max(abs(tr(1,1)),abs(tl(1,1)),abs(tl(1,2)),abs(tl( &
                     2,1)),abs(tl(2,2))),smlnum)
           tmp(1) = tl(1,1) + sgn*tr(1,1)
           tmp(4) = tl(2,2) + sgn*tr(1,1)
           if (ltranl) then
              tmp(2) = tl(1,2)
              tmp(3) = tl(2,1)
           else
              tmp(2) = tl(2,1)
              tmp(3) = tl(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           40 continue
           ! solve 2 by 2 system using complete pivoting.
           ! set pivots less than smin to smin.
           ipiv = la_iqamax(4,tmp,1)
           u11 = tmp(ipiv)
           if (abs(u11) <= smin) then
              info = 1
              u11 = smin
           end if
           u12 = tmp(locu12(ipiv))
           l21 = tmp(locl21(ipiv))/u11
           u22 = tmp(locu22(ipiv)) - u12*l21
           xswap = xswpiv(ipiv)
           bswap = bswpiv(ipiv)
           if (abs(u22) <= smin) then
              info = 1
              u22 = smin
           end if
           if (bswap) then
              temp = btmp(2)
              btmp(2) = btmp(1) - l21*temp
              btmp(1) = temp
           else
              btmp(2) = btmp(2) - l21*btmp(1)
           end if
           scale = one
           if ((two*smlnum)*abs(btmp(2)) > abs(u22) .or. (two*smlnum)*abs(btmp(1)) > abs( &
                      u11)) then
              scale = half/max(abs(btmp(1)),abs(btmp(2)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
           end if
           x2(2) = btmp(2)/u22
           x2(1) = btmp(1)/u11 - (u12/u11)*x2(2)
           if (xswap) then
              temp = x2(2)
              x2(2) = x2(1)
              x2(1) = temp
           end if
           x(1,1) = x2(1)
           if (n1 == 1) then
              x(1,2) = x2(2)
              xnorm = abs(x(1,1)) + abs(x(1,2))
           else
              x(2,1) = x2(2)
              xnorm = max(abs(x(1,1)),abs(x(2,1)))
           end if
           return
           ! 2 by 2:
           ! op[tl11 tl12]*[x11 x12] +isgn* [x11 x12]*op[tr11 tr12] = [b11 b12]
             ! [tl21 tl22] [x21 x22]        [x21 x22]   [tr21 tr22]   [b21 b22]
           ! solve equivalent 4 by 4 system using complete pivoting.
           ! set pivots less than smin to smin.
           50 continue
           smin = max(abs(tr(1,1)),abs(tr(1,2)),abs(tr(2,1)),abs(tr(2,2)))

           smin = max(smin,abs(tl(1,1)),abs(tl(1,2)),abs(tl(2,1)),abs(tl(2, &
                     2)))
           smin = max(eps*smin,smlnum)
           btmp(1) = zero
           call la_qcopy(16,btmp,0,t16,1)
           t16(1,1) = tl(1,1) + sgn*tr(1,1)
           t16(2,2) = tl(2,2) + sgn*tr(1,1)
           t16(3,3) = tl(1,1) + sgn*tr(2,2)
           t16(4,4) = tl(2,2) + sgn*tr(2,2)
           if (ltranl) then
              t16(1,2) = tl(2,1)
              t16(2,1) = tl(1,2)
              t16(3,4) = tl(2,1)
              t16(4,3) = tl(1,2)
           else
              t16(1,2) = tl(1,2)
              t16(2,1) = tl(2,1)
              t16(3,4) = tl(1,2)
              t16(4,3) = tl(2,1)
           end if
           if (ltranr) then
              t16(1,3) = sgn*tr(1,2)
              t16(2,4) = sgn*tr(1,2)
              t16(3,1) = sgn*tr(2,1)
              t16(4,2) = sgn*tr(2,1)
           else
              t16(1,3) = sgn*tr(2,1)
              t16(2,4) = sgn*tr(2,1)
              t16(3,1) = sgn*tr(1,2)
              t16(4,2) = sgn*tr(1,2)
           end if
           btmp(1) = b(1,1)
           btmp(2) = b(2,1)
           btmp(3) = b(1,2)
           btmp(4) = b(2,2)
           ! perform elimination
           loop_100: do i = 1,3
              xmax = zero
              do ip = i,4
                 do jp = i,4
                    if (abs(t16(ip,jp)) >= xmax) then
                       xmax = abs(t16(ip,jp))
                       ipsv = ip
                       jpsv = jp
                    end if
                 end do
              end do
              if (ipsv /= i) then
                 call la_qswap(4,t16(ipsv,1),4,t16(i,1),4)
                 temp = btmp(i)
                 btmp(i) = btmp(ipsv)
                 btmp(ipsv) = temp
              end if
              if (jpsv /= i) call la_qswap(4,t16(1,jpsv),1,t16(1,i),1)
              jpiv(i) = jpsv
              if (abs(t16(i,i)) < smin) then
                 info = 1
                 t16(i,i) = smin
              end if
              do j = i + 1,4
                 t16(j,i) = t16(j,i)/t16(i,i)
                 btmp(j) = btmp(j) - t16(j,i)*btmp(i)
                 do k = i + 1,4
                    t16(j,k) = t16(j,k) - t16(j,i)*t16(i,k)
                 end do
              end do
           end do loop_100
           if (abs(t16(4,4)) < smin) then
              info = 1
              t16(4,4) = smin
           end if
           scale = one
           if ((eight*smlnum)*abs(btmp(1)) > abs(t16(1,1)) .or. (eight*smlnum)*abs( &
           btmp(2)) > abs(t16(2,2)) .or. (eight*smlnum)*abs(btmp(3)) > abs(t16(3,3)) &
                      .or. (eight*smlnum)*abs(btmp(4)) > abs(t16(4,4))) then
              scale = (one/eight)/max(abs(btmp(1)),abs(btmp(2)),abs(btmp(3)), &
                        abs(btmp(4)))
              btmp(1) = btmp(1)*scale
              btmp(2) = btmp(2)*scale
              btmp(3) = btmp(3)*scale
              btmp(4) = btmp(4)*scale
           end if
           do i = 1,4
              k = 5 - i
              temp = one/t16(k,k)
              tmp(k) = btmp(k)*temp
              do j = k + 1,4
                 tmp(k) = tmp(k) - (temp*t16(k,j))*tmp(j)
              end do
           end do
           do i = 1,3
              if (jpiv(4 - i) /= 4 - i) then
                 temp = tmp(4 - i)
                 tmp(4 - i) = tmp(jpiv(4 - i))
                 tmp(jpiv(4 - i)) = temp
              end if
           end do
           x(1,1) = tmp(1)
           x(2,1) = tmp(2)
           x(1,2) = tmp(3)
           x(2,2) = tmp(4)
           xnorm = max(abs(tmp(1)) + abs(tmp(3)),abs(tmp(2)) + abs(tmp(4)))
           return
     end subroutine la_qlasy2
#endif

     !> SLALN2: solves a system of the form  (ca A - w D ) X = s B
     !> or (ca A**T - w D) X = s B   with possible scaling ("s") and
     !> perturbation of A.  (A**T means A-transpose.)
     !> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
     !> real diagonal matrix, w is a real or complex value, and X and B are
     !> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
     !> may be 1 or 2.
     !> If w is complex, X and B are represented as NA x 2 matrices,
     !> the first column of each being the real part and the second
     !> being the imaginary part.
     !> "s" is a scaling factor (<= 1), computed by SLALN2, which is
     !> so chosen that X can be computed without overflow.  X is further
     !> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
     !> than overflow.
     !> If both singular values of (ca A - w D) are less than SMIN,
     !> SMIN*identity will be used instead of (ca A - w D).  If only one
     !> singular value is less than SMIN, one element of (ca A - w D) will be
     !> perturbed enough to make the smallest singular value roughly SMIN.
     !> If both singular values are at least SMIN, (ca A - w D) will not be
     !> perturbed.  In any case, the perturbation will be at most some small
     !> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
     !> are computed by infinity-norm approximations, and thus will only be
     !> correct to a factor of 2 or so.
     !> Note: all input quantities are assumed to be smaller than overflow
     !> by a reasonable factor.  (See BIGNUM.)

     pure subroutine la_slaln2(ltrans,na,nw,smin,ca,a,lda,d1,d2,b,ldb,wr,wi,x, &
               ldx,scale,xnorm,info)
        use la_constants_sp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltrans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,na,nw
           real(sp),intent(in) :: ca,d1,d2,smin,wi,wr
           real(sp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: icmax,j
           real(sp) :: bbnd,bi1,bi2,bignum,bnorm,br1,br2,ci21,ci22,cmax,cnorm,cr21, &
           cr22,csi,csr,li21,lr21,smini,smlnum,temp,u22abs,ui11,ui11r,ui12,ui12s, &
                     ui22,ur11,ur11r,ur12,ur12s,ur22,xi1,xi2,xr1,xr2
           ! Local Arrays
           logical(lk) :: cswap(4),rswap(4)
           integer(ilp) :: ipivot(4,4)
           real(sp) :: ci(2,2),civ(4),cr(2,2),crv(4)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Equivalences
           equivalence(ci(1,1),civ(1)), (cr(1,1),crv(1))
           ! Data Statements
           cswap = [.false.,.false.,.true.,.true.]
           rswap = [.false.,.true.,.false.,.true.]
           ipivot = reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1], [4,4])
           ! Executable Statements
           ! compute bignum
           smlnum = two*la_slamch('SAFE MINIMUM')
           bignum = one/smlnum
           smini = max(smin,smlnum)
           ! don't check for input errors
           info = 0
           ! standard initializations
           scale = one
           if (na == 1) then
              ! 1 x 1  (i.e., scalar) system   c x = b
              if (nw == 1) then
                 ! real 1x1 system.
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 cnorm = abs(csr)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 x(1,1) = (b(1,1)*scale)/csr
                 xnorm = abs(x(1,1))
              else
                 ! complex 1x1 system (w is complex)
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 csi = -wi*d1
                 cnorm = abs(csr) + abs(csi)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    csi = zero
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1)) + abs(b(1,2))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 call la_sladiv(scale*b(1,1),scale*b(1,2),csr,csi,x(1,1),x(1, &
                           2))
                 xnorm = abs(x(1,1)) + abs(x(1,2))
              end if
           else
              ! 2x2 system
              ! compute the realpart of  c = ca a - w d  (or  ca a**t - w d,KIND=sp)
              cr(1,1) = ca*a(1,1) - wr*d1
              cr(2,2) = ca*a(2,2) - wr*d2
              if (ltrans) then
                 cr(1,2) = ca*a(2,1)
                 cr(2,1) = ca*a(1,2)
              else
                 cr(2,1) = ca*a(2,1)
                 cr(1,2) = ca*a(1,2)
              end if
              if (nw == 1) then
                 ! real2x2 system  (w is real,KIND=sp)
                 ! find the largest element in c
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) > cmax) then
                       cmax = abs(crv(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)),abs(b(2,1)))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ur11r = one/ur11
                 lr21 = ur11r*cr21
                 ur22 = cr22 - ur12*lr21
                 ! if smaller pivot < smini, use smini
                 if (abs(ur22) < smini) then
                    ur22 = smini
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br1 = b(2,1)
                    br2 = b(1,1)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                 end if
                 br2 = br2 - lr21*br1
                 bbnd = max(abs(br1*(ur22*ur11r)),abs(br2))
                 if (bbnd > one .and. abs(ur22) < one) then
                    if (bbnd >= bignum*abs(ur22)) scale = one/bbnd
                 end if
                 xr2 = (br2*scale)/ur22
                 xr1 = (scale*br1)*ur11r - xr2*(ur11r*ur12)
                 if (cswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                 end if
                 xnorm = max(abs(xr1),abs(xr2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              else
                 ! complex 2x2 system  (w is complex)
                 ! find the largest element in c
                 ci(1,1) = -wi*d1
                 ci(2,1) = zero
                 ci(1,2) = zero
                 ci(2,2) = -wi*d2
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) + abs(civ(j)) > cmax) then
                       cmax = abs(crv(j)) + abs(civ(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2) &
                               ))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    x(1,2) = temp*b(1,2)
                    x(2,2) = temp*b(2,2)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 ui11 = civ(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ci21 = civ(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 ui12 = civ(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ci22 = civ(ipivot(4,icmax))
                 if (icmax == 1 .or. icmax == 4) then
                    ! code when off-diagonals of pivoted c are real
                    if (abs(ur11) > abs(ui11)) then
                       temp = ui11/ur11
                       ur11r = one/(ur11*(one + temp**2))
                       ui11r = -temp*ur11r
                    else
                       temp = ur11/ui11
                       ui11r = -one/(ui11*(one + temp**2))
                       ur11r = -temp*ui11r
                    end if
                    lr21 = cr21*ur11r
                    li21 = cr21*ui11r
                    ur12s = ur12*ur11r
                    ui12s = ur12*ui11r
                    ur22 = cr22 - ur12*lr21
                    ui22 = ci22 - ur12*li21
                 else
                    ! code when diagonals of pivoted c are real
                    ur11r = one/ur11
                    ui11r = zero
                    lr21 = cr21*ur11r
                    li21 = ci21*ur11r
                    ur12s = ur12*ur11r
                    ui12s = ui12*ur11r
                    ur22 = cr22 - ur12*lr21 + ui12*li21
                    ui22 = -ur12*li21 - ui12*lr21
                 end if
                 u22abs = abs(ur22) + abs(ui22)
                 ! if smaller pivot < smini, use smini
                 if (u22abs < smini) then
                    ur22 = smini
                    ui22 = zero
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br2 = b(1,1)
                    br1 = b(2,1)
                    bi2 = b(1,2)
                    bi1 = b(2,2)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                    bi1 = b(1,2)
                    bi2 = b(2,2)
                 end if
                 br2 = br2 - lr21*br1 + li21*bi1
                 bi2 = bi2 - li21*br1 - lr21*bi1
                 bbnd = max((abs(br1) + abs(bi1))*(u22abs*(abs(ur11r) + abs(ui11r))), &
                           abs(br2) + abs(bi2))
                 if (bbnd > one .and. u22abs < one) then
                    if (bbnd >= bignum*u22abs) then
                       scale = one/bbnd
                       br1 = scale*br1
                       bi1 = scale*bi1
                       br2 = scale*br2
                       bi2 = scale*bi2
                    end if
                 end if
                 call la_sladiv(br2,bi2,ur22,ui22,xr2,xi2)
                 xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
                 xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
                 if (cswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                    x(1,2) = xi2
                    x(2,2) = xi1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                    x(1,2) = xi1
                    x(2,2) = xi2
                 end if
                 xnorm = max(abs(xr1) + abs(xi1),abs(xr2) + abs(xi2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       x(1,2) = temp*x(1,2)
                       x(2,2) = temp*x(2,2)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_slaln2
     !> DLALN2: solves a system of the form  (ca A - w D ) X = s B
     !> or (ca A**T - w D) X = s B   with possible scaling ("s") and
     !> perturbation of A.  (A**T means A-transpose.)
     !> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
     !> real diagonal matrix, w is a real or complex value, and X and B are
     !> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
     !> may be 1 or 2.
     !> If w is complex, X and B are represented as NA x 2 matrices,
     !> the first column of each being the real part and the second
     !> being the imaginary part.
     !> "s" is a scaling factor (<= 1), computed by DLALN2, which is
     !> so chosen that X can be computed without overflow.  X is further
     !> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
     !> than overflow.
     !> If both singular values of (ca A - w D) are less than SMIN,
     !> SMIN*identity will be used instead of (ca A - w D).  If only one
     !> singular value is less than SMIN, one element of (ca A - w D) will be
     !> perturbed enough to make the smallest singular value roughly SMIN.
     !> If both singular values are at least SMIN, (ca A - w D) will not be
     !> perturbed.  In any case, the perturbation will be at most some small
     !> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
     !> are computed by infinity-norm approximations, and thus will only be
     !> correct to a factor of 2 or so.
     !> Note: all input quantities are assumed to be smaller than overflow
     !> by a reasonable factor.  (See BIGNUM.)

     pure subroutine la_dlaln2(ltrans,na,nw,smin,ca,a,lda,d1,d2,b,ldb,wr,wi,x, &
               ldx,scale,xnorm,info)
        use la_constants_dp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltrans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,na,nw
           real(dp),intent(in) :: ca,d1,d2,smin,wi,wr
           real(dp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: icmax,j
           real(dp) :: bbnd,bi1,bi2,bignum,bnorm,br1,br2,ci21,ci22,cmax,cnorm,cr21, &
           cr22,csi,csr,li21,lr21,smini,smlnum,temp,u22abs,ui11,ui11r,ui12,ui12s, &
                     ui22,ur11,ur11r,ur12,ur12s,ur22,xi1,xi2,xr1,xr2
           ! Local Arrays
           logical(lk) :: rswap(4),zswap(4)
           integer(ilp) :: ipivot(4,4)
           real(dp) :: ci(2,2),civ(4),cr(2,2),crv(4)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Equivalences
           equivalence(ci(1,1),civ(1)), (cr(1,1),crv(1))
           ! Data Statements
           zswap = [.false.,.false.,.true.,.true.]
           rswap = [.false.,.true.,.false.,.true.]
           ipivot = reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1], [4,4])
           ! Executable Statements
           ! compute bignum
           smlnum = two*la_dlamch('SAFE MINIMUM')
           bignum = one/smlnum
           smini = max(smin,smlnum)
           ! don't check for input errors
           info = 0
           ! standard initializations
           scale = one
           if (na == 1) then
              ! 1 x 1  (i.e., scalar) system   c x = b
              if (nw == 1) then
                 ! real 1x1 system.
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 cnorm = abs(csr)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 x(1,1) = (b(1,1)*scale)/csr
                 xnorm = abs(x(1,1))
              else
                 ! complex 1x1 system (w is complex)
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 csi = -wi*d1
                 cnorm = abs(csr) + abs(csi)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    csi = zero
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1)) + abs(b(1,2))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 call la_dladiv(scale*b(1,1),scale*b(1,2),csr,csi,x(1,1),x(1, &
                           2))
                 xnorm = abs(x(1,1)) + abs(x(1,2))
              end if
           else
              ! 2x2 system
              ! compute the realpart of  c = ca a - w d  (or  ca a**t - w d,KIND=dp)
              cr(1,1) = ca*a(1,1) - wr*d1
              cr(2,2) = ca*a(2,2) - wr*d2
              if (ltrans) then
                 cr(1,2) = ca*a(2,1)
                 cr(2,1) = ca*a(1,2)
              else
                 cr(2,1) = ca*a(2,1)
                 cr(1,2) = ca*a(1,2)
              end if
              if (nw == 1) then
                 ! real2x2 system  (w is real,KIND=dp)
                 ! find the largest element in c
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) > cmax) then
                       cmax = abs(crv(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)),abs(b(2,1)))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ur11r = one/ur11
                 lr21 = ur11r*cr21
                 ur22 = cr22 - ur12*lr21
                 ! if smaller pivot < smini, use smini
                 if (abs(ur22) < smini) then
                    ur22 = smini
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br1 = b(2,1)
                    br2 = b(1,1)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                 end if
                 br2 = br2 - lr21*br1
                 bbnd = max(abs(br1*(ur22*ur11r)),abs(br2))
                 if (bbnd > one .and. abs(ur22) < one) then
                    if (bbnd >= bignum*abs(ur22)) scale = one/bbnd
                 end if
                 xr2 = (br2*scale)/ur22
                 xr1 = (scale*br1)*ur11r - xr2*(ur11r*ur12)
                 if (zswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                 end if
                 xnorm = max(abs(xr1),abs(xr2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              else
                 ! complex 2x2 system  (w is complex)
                 ! find the largest element in c
                 ci(1,1) = -wi*d1
                 ci(2,1) = zero
                 ci(1,2) = zero
                 ci(2,2) = -wi*d2
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) + abs(civ(j)) > cmax) then
                       cmax = abs(crv(j)) + abs(civ(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2) &
                               ))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    x(1,2) = temp*b(1,2)
                    x(2,2) = temp*b(2,2)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 ui11 = civ(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ci21 = civ(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 ui12 = civ(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ci22 = civ(ipivot(4,icmax))
                 if (icmax == 1 .or. icmax == 4) then
                    ! code when off-diagonals of pivoted c are real
                    if (abs(ur11) > abs(ui11)) then
                       temp = ui11/ur11
                       ur11r = one/(ur11*(one + temp**2))
                       ui11r = -temp*ur11r
                    else
                       temp = ur11/ui11
                       ui11r = -one/(ui11*(one + temp**2))
                       ur11r = -temp*ui11r
                    end if
                    lr21 = cr21*ur11r
                    li21 = cr21*ui11r
                    ur12s = ur12*ur11r
                    ui12s = ur12*ui11r
                    ur22 = cr22 - ur12*lr21
                    ui22 = ci22 - ur12*li21
                 else
                    ! code when diagonals of pivoted c are real
                    ur11r = one/ur11
                    ui11r = zero
                    lr21 = cr21*ur11r
                    li21 = ci21*ur11r
                    ur12s = ur12*ur11r
                    ui12s = ui12*ur11r
                    ur22 = cr22 - ur12*lr21 + ui12*li21
                    ui22 = -ur12*li21 - ui12*lr21
                 end if
                 u22abs = abs(ur22) + abs(ui22)
                 ! if smaller pivot < smini, use smini
                 if (u22abs < smini) then
                    ur22 = smini
                    ui22 = zero
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br2 = b(1,1)
                    br1 = b(2,1)
                    bi2 = b(1,2)
                    bi1 = b(2,2)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                    bi1 = b(1,2)
                    bi2 = b(2,2)
                 end if
                 br2 = br2 - lr21*br1 + li21*bi1
                 bi2 = bi2 - li21*br1 - lr21*bi1
                 bbnd = max((abs(br1) + abs(bi1))*(u22abs*(abs(ur11r) + abs(ui11r))), &
                           abs(br2) + abs(bi2))
                 if (bbnd > one .and. u22abs < one) then
                    if (bbnd >= bignum*u22abs) then
                       scale = one/bbnd
                       br1 = scale*br1
                       bi1 = scale*bi1
                       br2 = scale*br2
                       bi2 = scale*bi2
                    end if
                 end if
                 call la_dladiv(br2,bi2,ur22,ui22,xr2,xi2)
                 xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
                 xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
                 if (zswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                    x(1,2) = xi2
                    x(2,2) = xi1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                    x(1,2) = xi1
                    x(2,2) = xi2
                 end if
                 xnorm = max(abs(xr1) + abs(xi1),abs(xr2) + abs(xi2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       x(1,2) = temp*x(1,2)
                       x(2,2) = temp*x(2,2)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_dlaln2
#ifdef LA_WITH_XDP
     !> XLALN2: solves a system of the form  (ca A - w D ) X = s B
     !> or (ca A**T - w D) X = s B   with possible scaling ("s") and
     !> perturbation of A.  (A**T means A-transpose.)
     !> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
     !> real diagonal matrix, w is a real or complex value, and X and B are
     !> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
     !> may be 1 or 2.
     !> If w is complex, X and B are represented as NA x 2 matrices,
     !> the first column of each being the real part and the second
     !> being the imaginary part.
     !> "s" is a scaling factor (<= 1), computed by XLALN2, which is
     !> so chosen that X can be computed without overflow.  X is further
     !> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
     !> than overflow.
     !> If both singular values of (ca A - w D) are less than SMIN,
     !> SMIN*identity will be used instead of (ca A - w D).  If only one
     !> singular value is less than SMIN, one element of (ca A - w D) will be
     !> perturbed enough to make the smallest singular value roughly SMIN.
     !> If both singular values are at least SMIN, (ca A - w D) will not be
     !> perturbed.  In any case, the perturbation will be at most some small
     !> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
     !> are computed by infinity-norm approximations, and thus will only be
     !> correct to a factor of 2 or so.
     !> Note: all input quantities are assumed to be smaller than overflow
     !> by a reasonable factor.  (See BIGNUM.)

     pure subroutine la_xlaln2(ltrans,na,nw,smin,ca,a,lda,d1,d2,b,ldb,wr,wi,x, &
               ldx,scale,xnorm,info)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltrans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,na,nw
           real(xdp),intent(in) :: ca,d1,d2,smin,wi,wr
           real(xdp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: icmax,j
           real(xdp) :: bbnd,bi1,bi2,bignum,bnorm,br1,br2,ci21,ci22,cmax,cnorm,cr21, &
           cr22,csi,csr,li21,lr21,smini,smlnum,temp,u22abs,ui11,ui11r,ui12,ui12s, &
                     ui22,ur11,ur11r,ur12,ur12s,ur22,xi1,xi2,xr1,xr2
           ! Local Arrays
           logical(lk) :: rswap(4),yswap(4)
           integer(ilp) :: ipivot(4,4)
           real(xdp) :: ci(2,2),civ(4),cr(2,2),crv(4)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Equivalences
           equivalence(ci(1,1),civ(1)), (cr(1,1),crv(1))
           ! Data Statements
           yswap = [.false.,.false.,.true.,.true.]
           rswap = [.false.,.true.,.false.,.true.]
           ipivot = reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1], [4,4])
           ! Executable Statements
           ! compute bignum
           smlnum = two*la_xlamch('SAFE MINIMUM')
           bignum = one/smlnum
           smini = max(smin,smlnum)
           ! don't check for input errors
           info = 0
           ! standard initializations
           scale = one
           if (na == 1) then
              ! 1 x 1  (i.e., scalar) system   c x = b
              if (nw == 1) then
                 ! real 1x1 system.
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 cnorm = abs(csr)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 x(1,1) = (b(1,1)*scale)/csr
                 xnorm = abs(x(1,1))
              else
                 ! complex 1x1 system (w is complex)
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 csi = -wi*d1
                 cnorm = abs(csr) + abs(csi)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    csi = zero
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1)) + abs(b(1,2))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 call la_xladiv(scale*b(1,1),scale*b(1,2),csr,csi,x(1,1),x(1, &
                           2))
                 xnorm = abs(x(1,1)) + abs(x(1,2))
              end if
           else
              ! 2x2 system
              ! compute the realpart of  c = ca a - w d  (or  ca a**t - w d,KIND=xdp)
              cr(1,1) = ca*a(1,1) - wr*d1
              cr(2,2) = ca*a(2,2) - wr*d2
              if (ltrans) then
                 cr(1,2) = ca*a(2,1)
                 cr(2,1) = ca*a(1,2)
              else
                 cr(2,1) = ca*a(2,1)
                 cr(1,2) = ca*a(1,2)
              end if
              if (nw == 1) then
                 ! real2x2 system  (w is real,KIND=xdp)
                 ! find the largest element in c
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) > cmax) then
                       cmax = abs(crv(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)),abs(b(2,1)))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ur11r = one/ur11
                 lr21 = ur11r*cr21
                 ur22 = cr22 - ur12*lr21
                 ! if smaller pivot < smini, use smini
                 if (abs(ur22) < smini) then
                    ur22 = smini
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br1 = b(2,1)
                    br2 = b(1,1)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                 end if
                 br2 = br2 - lr21*br1
                 bbnd = max(abs(br1*(ur22*ur11r)),abs(br2))
                 if (bbnd > one .and. abs(ur22) < one) then
                    if (bbnd >= bignum*abs(ur22)) scale = one/bbnd
                 end if
                 xr2 = (br2*scale)/ur22
                 xr1 = (scale*br1)*ur11r - xr2*(ur11r*ur12)
                 if (yswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                 end if
                 xnorm = max(abs(xr1),abs(xr2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              else
                 ! complex 2x2 system  (w is complex)
                 ! find the largest element in c
                 ci(1,1) = -wi*d1
                 ci(2,1) = zero
                 ci(1,2) = zero
                 ci(2,2) = -wi*d2
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) + abs(civ(j)) > cmax) then
                       cmax = abs(crv(j)) + abs(civ(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2) &
                               ))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    x(1,2) = temp*b(1,2)
                    x(2,2) = temp*b(2,2)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 ui11 = civ(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ci21 = civ(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 ui12 = civ(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ci22 = civ(ipivot(4,icmax))
                 if (icmax == 1 .or. icmax == 4) then
                    ! code when off-diagonals of pivoted c are real
                    if (abs(ur11) > abs(ui11)) then
                       temp = ui11/ur11
                       ur11r = one/(ur11*(one + temp**2))
                       ui11r = -temp*ur11r
                    else
                       temp = ur11/ui11
                       ui11r = -one/(ui11*(one + temp**2))
                       ur11r = -temp*ui11r
                    end if
                    lr21 = cr21*ur11r
                    li21 = cr21*ui11r
                    ur12s = ur12*ur11r
                    ui12s = ur12*ui11r
                    ur22 = cr22 - ur12*lr21
                    ui22 = ci22 - ur12*li21
                 else
                    ! code when diagonals of pivoted c are real
                    ur11r = one/ur11
                    ui11r = zero
                    lr21 = cr21*ur11r
                    li21 = ci21*ur11r
                    ur12s = ur12*ur11r
                    ui12s = ui12*ur11r
                    ur22 = cr22 - ur12*lr21 + ui12*li21
                    ui22 = -ur12*li21 - ui12*lr21
                 end if
                 u22abs = abs(ur22) + abs(ui22)
                 ! if smaller pivot < smini, use smini
                 if (u22abs < smini) then
                    ur22 = smini
                    ui22 = zero
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br2 = b(1,1)
                    br1 = b(2,1)
                    bi2 = b(1,2)
                    bi1 = b(2,2)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                    bi1 = b(1,2)
                    bi2 = b(2,2)
                 end if
                 br2 = br2 - lr21*br1 + li21*bi1
                 bi2 = bi2 - li21*br1 - lr21*bi1
                 bbnd = max((abs(br1) + abs(bi1))*(u22abs*(abs(ur11r) + abs(ui11r))), &
                           abs(br2) + abs(bi2))
                 if (bbnd > one .and. u22abs < one) then
                    if (bbnd >= bignum*u22abs) then
                       scale = one/bbnd
                       br1 = scale*br1
                       bi1 = scale*bi1
                       br2 = scale*br2
                       bi2 = scale*bi2
                    end if
                 end if
                 call la_xladiv(br2,bi2,ur22,ui22,xr2,xi2)
                 xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
                 xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
                 if (yswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                    x(1,2) = xi2
                    x(2,2) = xi1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                    x(1,2) = xi1
                    x(2,2) = xi2
                 end if
                 xnorm = max(abs(xr1) + abs(xi1),abs(xr2) + abs(xi2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       x(1,2) = temp*x(1,2)
                       x(2,2) = temp*x(2,2)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_xlaln2
#endif
#ifdef LA_WITH_QP
     !> QLALN2: solves a system of the form  (ca A - w D ) X = s B
     !> or (ca A**T - w D) X = s B   with possible scaling ("s") and
     !> perturbation of A.  (A**T means A-transpose.)
     !> A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
     !> real diagonal matrix, w is a real or complex value, and X and B are
     !> NA x 1 matrices -- real if w is real, complex if w is complex.  NA
     !> may be 1 or 2.
     !> If w is complex, X and B are represented as NA x 2 matrices,
     !> the first column of each being the real part and the second
     !> being the imaginary part.
     !> "s" is a scaling factor (<= 1), computed by QLALN2, which is
     !> so chosen that X can be computed without overflow.  X is further
     !> scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
     !> than overflow.
     !> If both singular values of (ca A - w D) are less than SMIN,
     !> SMIN*identity will be used instead of (ca A - w D).  If only one
     !> singular value is less than SMIN, one element of (ca A - w D) will be
     !> perturbed enough to make the smallest singular value roughly SMIN.
     !> If both singular values are at least SMIN, (ca A - w D) will not be
     !> perturbed.  In any case, the perturbation will be at most some small
     !> multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
     !> are computed by infinity-norm approximations, and thus will only be
     !> correct to a factor of 2 or so.
     !> Note: all input quantities are assumed to be smaller than overflow
     !> by a reasonable factor.  (See BIGNUM.)

     pure subroutine la_qlaln2(ltrans,na,nw,smin,ca,a,lda,d1,d2,b,ldb,wr,wi,x, &
               ldx,scale,xnorm,info)
        use la_constants_qp,only:zero,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ltrans
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,na,nw
           real(qp),intent(in) :: ca,d1,d2,smin,wi,wr
           real(qp),intent(out) :: scale,xnorm
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: x(ldx,*)
       ! =====================================================================

           ! Local Scalars
           integer(ilp) :: icmax,j
           real(qp) :: bbnd,bi1,bi2,bignum,bnorm,br1,br2,ci21,ci22,cmax,cnorm,cr21, &
           cr22,csi,csr,li21,lr21,smini,smlnum,temp,u22abs,ui11,ui11r,ui12,ui12s, &
                     ui22,ur11,ur11r,ur12,ur12s,ur22,xi1,xi2,xr1,xr2
           ! Local Arrays
           logical(lk) :: rswap(4),wswap(4)
           integer(ilp) :: ipivot(4,4)
           real(qp) :: ci(2,2),civ(4),cr(2,2),crv(4)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Equivalences
           equivalence(ci(1,1),civ(1)), (cr(1,1),crv(1))
           ! Data Statements
           wswap = [.false.,.false.,.true.,.true.]
           rswap = [.false.,.true.,.false.,.true.]
           ipivot = reshape([1,2,3,4,2,1,4,3,3,4,1,2,4,3,2,1], [4,4])
           ! Executable Statements
           ! compute bignum
           smlnum = two*la_qlamch('SAFE MINIMUM')
           bignum = one/smlnum
           smini = max(smin,smlnum)
           ! don't check for input errors
           info = 0
           ! standard initializations
           scale = one
           if (na == 1) then
              ! 1 x 1  (i.e., scalar) system   c x = b
              if (nw == 1) then
                 ! real 1x1 system.
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 cnorm = abs(csr)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 x(1,1) = (b(1,1)*scale)/csr
                 xnorm = abs(x(1,1))
              else
                 ! complex 1x1 system (w is complex)
                 ! c = ca a - w d
                 csr = ca*a(1,1) - wr*d1
                 csi = -wi*d1
                 cnorm = abs(csr) + abs(csi)
                 ! if | c | < smini, use c = smini
                 if (cnorm < smini) then
                    csr = smini
                    csi = zero
                    cnorm = smini
                    info = 1
                 end if
                 ! check scaling for  x = b / c
                 bnorm = abs(b(1,1)) + abs(b(1,2))
                 if (cnorm < one .and. bnorm > one) then
                    if (bnorm > bignum*cnorm) scale = one/bnorm
                 end if
                 ! compute x
                 call la_qladiv(scale*b(1,1),scale*b(1,2),csr,csi,x(1,1),x(1, &
                           2))
                 xnorm = abs(x(1,1)) + abs(x(1,2))
              end if
           else
              ! 2x2 system
              ! compute the realpart of  c = ca a - w d  (or  ca a**t - w d,KIND=qp)
              cr(1,1) = ca*a(1,1) - wr*d1
              cr(2,2) = ca*a(2,2) - wr*d2
              if (ltrans) then
                 cr(1,2) = ca*a(2,1)
                 cr(2,1) = ca*a(1,2)
              else
                 cr(2,1) = ca*a(2,1)
                 cr(1,2) = ca*a(1,2)
              end if
              if (nw == 1) then
                 ! real2x2 system  (w is real,KIND=qp)
                 ! find the largest element in c
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) > cmax) then
                       cmax = abs(crv(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)),abs(b(2,1)))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ur11r = one/ur11
                 lr21 = ur11r*cr21
                 ur22 = cr22 - ur12*lr21
                 ! if smaller pivot < smini, use smini
                 if (abs(ur22) < smini) then
                    ur22 = smini
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br1 = b(2,1)
                    br2 = b(1,1)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                 end if
                 br2 = br2 - lr21*br1
                 bbnd = max(abs(br1*(ur22*ur11r)),abs(br2))
                 if (bbnd > one .and. abs(ur22) < one) then
                    if (bbnd >= bignum*abs(ur22)) scale = one/bbnd
                 end if
                 xr2 = (br2*scale)/ur22
                 xr1 = (scale*br1)*ur11r - xr2*(ur11r*ur12)
                 if (wswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                 end if
                 xnorm = max(abs(xr1),abs(xr2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              else
                 ! complex 2x2 system  (w is complex)
                 ! find the largest element in c
                 ci(1,1) = -wi*d1
                 ci(2,1) = zero
                 ci(1,2) = zero
                 ci(2,2) = -wi*d2
                 cmax = zero
                 icmax = 0
                 do j = 1,4
                    if (abs(crv(j)) + abs(civ(j)) > cmax) then
                       cmax = abs(crv(j)) + abs(civ(j))
                       icmax = j
                    end if
                 end do
                 ! if norm(c) < smini, use smini*identity.
                 if (cmax < smini) then
                    bnorm = max(abs(b(1,1)) + abs(b(1,2)),abs(b(2,1)) + abs(b(2,2) &
                               ))
                    if (smini < one .and. bnorm > one) then
                       if (bnorm > bignum*smini) scale = one/bnorm
                    end if
                    temp = scale/smini
                    x(1,1) = temp*b(1,1)
                    x(2,1) = temp*b(2,1)
                    x(1,2) = temp*b(1,2)
                    x(2,2) = temp*b(2,2)
                    xnorm = temp*bnorm
                    info = 1
                    return
                 end if
                 ! gaussian elimination with complete pivoting.
                 ur11 = crv(icmax)
                 ui11 = civ(icmax)
                 cr21 = crv(ipivot(2,icmax))
                 ci21 = civ(ipivot(2,icmax))
                 ur12 = crv(ipivot(3,icmax))
                 ui12 = civ(ipivot(3,icmax))
                 cr22 = crv(ipivot(4,icmax))
                 ci22 = civ(ipivot(4,icmax))
                 if (icmax == 1 .or. icmax == 4) then
                    ! code when off-diagonals of pivoted c are real
                    if (abs(ur11) > abs(ui11)) then
                       temp = ui11/ur11
                       ur11r = one/(ur11*(one + temp**2))
                       ui11r = -temp*ur11r
                    else
                       temp = ur11/ui11
                       ui11r = -one/(ui11*(one + temp**2))
                       ur11r = -temp*ui11r
                    end if
                    lr21 = cr21*ur11r
                    li21 = cr21*ui11r
                    ur12s = ur12*ur11r
                    ui12s = ur12*ui11r
                    ur22 = cr22 - ur12*lr21
                    ui22 = ci22 - ur12*li21
                 else
                    ! code when diagonals of pivoted c are real
                    ur11r = one/ur11
                    ui11r = zero
                    lr21 = cr21*ur11r
                    li21 = ci21*ur11r
                    ur12s = ur12*ur11r
                    ui12s = ui12*ur11r
                    ur22 = cr22 - ur12*lr21 + ui12*li21
                    ui22 = -ur12*li21 - ui12*lr21
                 end if
                 u22abs = abs(ur22) + abs(ui22)
                 ! if smaller pivot < smini, use smini
                 if (u22abs < smini) then
                    ur22 = smini
                    ui22 = zero
                    info = 1
                 end if
                 if (rswap(icmax)) then
                    br2 = b(1,1)
                    br1 = b(2,1)
                    bi2 = b(1,2)
                    bi1 = b(2,2)
                 else
                    br1 = b(1,1)
                    br2 = b(2,1)
                    bi1 = b(1,2)
                    bi2 = b(2,2)
                 end if
                 br2 = br2 - lr21*br1 + li21*bi1
                 bi2 = bi2 - li21*br1 - lr21*bi1
                 bbnd = max((abs(br1) + abs(bi1))*(u22abs*(abs(ur11r) + abs(ui11r))), &
                           abs(br2) + abs(bi2))
                 if (bbnd > one .and. u22abs < one) then
                    if (bbnd >= bignum*u22abs) then
                       scale = one/bbnd
                       br1 = scale*br1
                       bi1 = scale*bi1
                       br2 = scale*br2
                       bi2 = scale*bi2
                    end if
                 end if
                 call la_qladiv(br2,bi2,ur22,ui22,xr2,xi2)
                 xr1 = ur11r*br1 - ui11r*bi1 - ur12s*xr2 + ui12s*xi2
                 xi1 = ui11r*br1 + ur11r*bi1 - ui12s*xr2 - ur12s*xi2
                 if (wswap(icmax)) then
                    x(1,1) = xr2
                    x(2,1) = xr1
                    x(1,2) = xi2
                    x(2,2) = xi1
                 else
                    x(1,1) = xr1
                    x(2,1) = xr2
                    x(1,2) = xi1
                    x(2,2) = xi2
                 end if
                 xnorm = max(abs(xr1) + abs(xi1),abs(xr2) + abs(xi2))
                 ! further scaling if  norm(a) norm(x) > overflow
                 if (xnorm > one .and. cmax > one) then
                    if (xnorm > bignum/cmax) then
                       temp = cmax/bignum
                       x(1,1) = temp*x(1,1)
                       x(2,1) = temp*x(2,1)
                       x(1,2) = temp*x(1,2)
                       x(2,2) = temp*x(2,2)
                       xnorm = temp*xnorm
                       scale = temp*scale
                    end if
                 end if
              end if
           end if
           return
     end subroutine la_qlaln2
#endif

     !> SLANV2: computes the Schur factorization of a real 2-by-2 nonsymmetric
     !> matrix in standard form:
     !> [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
     !> [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
     !> where either
     !> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
     !> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
     !> conjugate eigenvalues.

     pure subroutine la_slanv2(a,b,c,d,rt1r,rt1i,rt2r,rt2i,cs,sn)
        use la_constants_sp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(inout) :: a,b,c,d
           real(sp),intent(out) :: cs,rt1i,rt1r,rt2i,rt2r,sn
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: multpl = 4.0e+0_sp

           ! Local Scalars
           real(sp) :: aa,bb,bcmax,bcmis,cc,cs1,dd,eps,p,sab,sac,scale,sigma,sn1, &
                     tau,temp,z,safmin,safmn2,safmx2
           integer(ilp) :: count
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sign,sqrt
           ! Executable Statements
           safmin = la_slamch('S')
           eps = la_slamch('P')
           safmn2 = la_slamch('B')**int(log(safmin/eps)/log(la_slamch('B'))/ &
                     two,KIND=ilp)
           safmx2 = one/safmn2
           if (c == zero) then
              cs = one
              sn = zero
           else if (b == zero) then
              ! swap rows and columns
              cs = zero
              sn = one
              temp = d
              d = a
              a = temp
              b = -c
              c = zero
           else if ((a - d) == zero .and. sign(one,b) /= sign(one,c)) then
              cs = one
              sn = zero
           else
              temp = a - d
              p = half*temp
              bcmax = max(abs(b),abs(c))
              bcmis = min(abs(b),abs(c))*sign(one,b)*sign(one,c)
              scale = max(abs(p),bcmax)
              z = (p/scale)*p + (bcmax/scale)*bcmis
              ! if z is of the order of the machine accuracy, postpone the
              ! decision on the nature of eigenvalues
              if (z >= multpl*eps) then
                 ! real eigenvalues. compute a and d.
                 z = p + sign(sqrt(scale)*sqrt(z),p)
                 a = d + z
                 d = d - (bcmax/z)*bcmis
                 ! compute b and the rotation matrix
                 tau = la_slapy2(c,z)
                 cs = z/tau
                 sn = c/tau
                 b = b - c
                 c = zero
              else
                 ! complex eigenvalues, or real(almost,KIND=sp) equal eigenvalues.
                 ! make diagonal elements equal.
                 count = 0
                 sigma = b + c
                 10 continue
                 count = count + 1
                 scale = max(abs(temp),abs(sigma))
                 if (scale >= safmx2) then
                    sigma = sigma*safmn2
                    temp = temp*safmn2
                    if (count <= 20) goto 10
                 end if
                 if (scale <= safmn2) then
                    sigma = sigma*safmx2
                    temp = temp*safmx2
                    if (count <= 20) goto 10
                 end if
                 p = half*temp
                 tau = la_slapy2(sigma,temp)
                 cs = sqrt(half*(one + abs(sigma)/tau))
                 sn = -(p/(tau*cs))*sign(one,sigma)
                 ! compute [ aa  bb ] = [ a  b ] [ cs -sn ]
                         ! [ cc  dd ]   [ c  d ] [ sn  cs ]
                 aa = a*cs + b*sn
                 bb = -a*sn + b*cs
                 cc = c*cs + d*sn
                 dd = -c*sn + d*cs
                 ! compute [ a  b ] = [ cs  sn ] [ aa  bb ]
                         ! [ c  d ]   [-sn  cs ] [ cc  dd ]
                 a = aa*cs + cc*sn
                 b = bb*cs + dd*sn
                 c = -aa*sn + cc*cs
                 d = -bb*sn + dd*cs
                 temp = half*(a + d)
                 a = temp
                 d = temp
                 if (c /= zero) then
                    if (b /= zero) then
                       if (sign(one,b) == sign(one,c)) then
                          ! real eigenvalues: reduce to upper triangular form
                          sab = sqrt(abs(b))
                          sac = sqrt(abs(c))
                          p = sign(sab*sac,c)
                          tau = one/sqrt(abs(b + c))
                          a = temp + p
                          d = temp - p
                          b = b - c
                          c = zero
                          cs1 = sab*tau
                          sn1 = sac*tau
                          temp = cs*cs1 - sn*sn1
                          sn = cs*sn1 + sn*cs1
                          cs = temp
                       end if
                    else
                       b = -c
                       c = zero
                       temp = cs
                       cs = -sn
                       sn = temp
                    end if
                 end if
              end if
           end if
           ! store eigenvalues in (rt1r,rt1i) and (rt2r,rt2i).
           rt1r = a
           rt2r = d
           if (c == zero) then
              rt1i = zero
              rt2i = zero
           else
              rt1i = sqrt(abs(b))*sqrt(abs(c))
              rt2i = -rt1i
           end if
           return
     end subroutine la_slanv2
     !> DLANV2: computes the Schur factorization of a real 2-by-2 nonsymmetric
     !> matrix in standard form:
     !> [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
     !> [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
     !> where either
     !> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
     !> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
     !> conjugate eigenvalues.

     pure subroutine la_dlanv2(a,b,c,d,rt1r,rt1i,rt2r,rt2i,cs,sn)
        use la_constants_dp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(inout) :: a,b,c,d
           real(dp),intent(out) :: cs,rt1i,rt1r,rt2i,rt2r,sn
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: multpl = 4.0e+0_dp

           ! Local Scalars
           real(dp) :: aa,bb,bcmax,bcmis,cc,cs1,dd,eps,p,sab,sac,scale,sigma,sn1, &
                     tau,temp,z,safmin,safmn2,safmx2
           integer(ilp) :: count
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sign,sqrt
           ! Executable Statements
           safmin = la_dlamch('S')
           eps = la_dlamch('P')
           safmn2 = la_dlamch('B')**int(log(safmin/eps)/log(la_dlamch('B'))/ &
                     two,KIND=ilp)
           safmx2 = one/safmn2
           if (c == zero) then
              cs = one
              sn = zero
           else if (b == zero) then
              ! swap rows and columns
              cs = zero
              sn = one
              temp = d
              d = a
              a = temp
              b = -c
              c = zero
           else if ((a - d) == zero .and. sign(one,b) /= sign(one,c)) then
              cs = one
              sn = zero
           else
              temp = a - d
              p = half*temp
              bcmax = max(abs(b),abs(c))
              bcmis = min(abs(b),abs(c))*sign(one,b)*sign(one,c)
              scale = max(abs(p),bcmax)
              z = (p/scale)*p + (bcmax/scale)*bcmis
              ! if z is of the order of the machine accuracy, postpone the
              ! decision on the nature of eigenvalues
              if (z >= multpl*eps) then
                 ! real eigenvalues. compute a and d.
                 z = p + sign(sqrt(scale)*sqrt(z),p)
                 a = d + z
                 d = d - (bcmax/z)*bcmis
                 ! compute b and the rotation matrix
                 tau = la_dlapy2(c,z)
                 cs = z/tau
                 sn = c/tau
                 b = b - c
                 c = zero
              else
                 ! complex eigenvalues, or real(almost,KIND=dp) equal eigenvalues.
                 ! make diagonal elements equal.
                 count = 0
                 sigma = b + c
                 10 continue
                 count = count + 1
                 scale = max(abs(temp),abs(sigma))
                 if (scale >= safmx2) then
                    sigma = sigma*safmn2
                    temp = temp*safmn2
                    if (count <= 20) goto 10
                 end if
                 if (scale <= safmn2) then
                    sigma = sigma*safmx2
                    temp = temp*safmx2
                    if (count <= 20) goto 10
                 end if
                 p = half*temp
                 tau = la_dlapy2(sigma,temp)
                 cs = sqrt(half*(one + abs(sigma)/tau))
                 sn = -(p/(tau*cs))*sign(one,sigma)
                 ! compute [ aa  bb ] = [ a  b ] [ cs -sn ]
                         ! [ cc  dd ]   [ c  d ] [ sn  cs ]
                 aa = a*cs + b*sn
                 bb = -a*sn + b*cs
                 cc = c*cs + d*sn
                 dd = -c*sn + d*cs
                 ! compute [ a  b ] = [ cs  sn ] [ aa  bb ]
                         ! [ c  d ]   [-sn  cs ] [ cc  dd ]
                 a = aa*cs + cc*sn
                 b = bb*cs + dd*sn
                 c = -aa*sn + cc*cs
                 d = -bb*sn + dd*cs
                 temp = half*(a + d)
                 a = temp
                 d = temp
                 if (c /= zero) then
                    if (b /= zero) then
                       if (sign(one,b) == sign(one,c)) then
                          ! real eigenvalues: reduce to upper triangular form
                          sab = sqrt(abs(b))
                          sac = sqrt(abs(c))
                          p = sign(sab*sac,c)
                          tau = one/sqrt(abs(b + c))
                          a = temp + p
                          d = temp - p
                          b = b - c
                          c = zero
                          cs1 = sab*tau
                          sn1 = sac*tau
                          temp = cs*cs1 - sn*sn1
                          sn = cs*sn1 + sn*cs1
                          cs = temp
                       end if
                    else
                       b = -c
                       c = zero
                       temp = cs
                       cs = -sn
                       sn = temp
                    end if
                 end if
              end if
           end if
           ! store eigenvalues in (rt1r,rt1i) and (rt2r,rt2i).
           rt1r = a
           rt2r = d
           if (c == zero) then
              rt1i = zero
              rt2i = zero
           else
              rt1i = sqrt(abs(b))*sqrt(abs(c))
              rt2i = -rt1i
           end if
           return
     end subroutine la_dlanv2
#ifdef LA_WITH_XDP
     !> XLANV2: computes the Schur factorization of a real 2-by-2 nonsymmetric
     !> matrix in standard form:
     !> [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
     !> [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
     !> where either
     !> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
     !> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
     !> conjugate eigenvalues.

     pure subroutine la_xlanv2(a,b,c,d,rt1r,rt1i,rt2r,rt2i,cs,sn)
        use la_constants_xdp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(inout) :: a,b,c,d
           real(xdp),intent(out) :: cs,rt1i,rt1r,rt2i,rt2r,sn
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: multpl = 4.0e+0_xdp

           ! Local Scalars
           real(xdp) :: aa,bb,bcmax,bcmis,cc,cs1,dd,eps,p,sab,sac,scale,sigma,sn1, &
                     tau,temp,z,safmin,safmn2,safmx2
           integer(ilp) :: count
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sign,sqrt
           ! Executable Statements
           safmin = la_xlamch('S')
           eps = la_xlamch('P')
           safmn2 = la_xlamch('B')**int(log(safmin/eps)/log(la_xlamch('B'))/ &
                     two,KIND=ilp)
           safmx2 = one/safmn2
           if (c == zero) then
              cs = one
              sn = zero
           else if (b == zero) then
              ! swap rows and columns
              cs = zero
              sn = one
              temp = d
              d = a
              a = temp
              b = -c
              c = zero
           else if ((a - d) == zero .and. sign(one,b) /= sign(one,c)) then
              cs = one
              sn = zero
           else
              temp = a - d
              p = half*temp
              bcmax = max(abs(b),abs(c))
              bcmis = min(abs(b),abs(c))*sign(one,b)*sign(one,c)
              scale = max(abs(p),bcmax)
              z = (p/scale)*p + (bcmax/scale)*bcmis
              ! if z is of the order of the machine accuracy, postpone the
              ! decision on the nature of eigenvalues
              if (z >= multpl*eps) then
                 ! real eigenvalues. compute a and d.
                 z = p + sign(sqrt(scale)*sqrt(z),p)
                 a = d + z
                 d = d - (bcmax/z)*bcmis
                 ! compute b and the rotation matrix
                 tau = la_xlapy2(c,z)
                 cs = z/tau
                 sn = c/tau
                 b = b - c
                 c = zero
              else
                 ! complex eigenvalues, or real(almost,KIND=xdp) equal eigenvalues.
                 ! make diagonal elements equal.
                 count = 0
                 sigma = b + c
                 10 continue
                 count = count + 1
                 scale = max(abs(temp),abs(sigma))
                 if (scale >= safmx2) then
                    sigma = sigma*safmn2
                    temp = temp*safmn2
                    if (count <= 20) goto 10
                 end if
                 if (scale <= safmn2) then
                    sigma = sigma*safmx2
                    temp = temp*safmx2
                    if (count <= 20) goto 10
                 end if
                 p = half*temp
                 tau = la_xlapy2(sigma,temp)
                 cs = sqrt(half*(one + abs(sigma)/tau))
                 sn = -(p/(tau*cs))*sign(one,sigma)
                 ! compute [ aa  bb ] = [ a  b ] [ cs -sn ]
                         ! [ cc  dd ]   [ c  d ] [ sn  cs ]
                 aa = a*cs + b*sn
                 bb = -a*sn + b*cs
                 cc = c*cs + d*sn
                 dd = -c*sn + d*cs
                 ! compute [ a  b ] = [ cs  sn ] [ aa  bb ]
                         ! [ c  d ]   [-sn  cs ] [ cc  dd ]
                 a = aa*cs + cc*sn
                 b = bb*cs + dd*sn
                 c = -aa*sn + cc*cs
                 d = -bb*sn + dd*cs
                 temp = half*(a + d)
                 a = temp
                 d = temp
                 if (c /= zero) then
                    if (b /= zero) then
                       if (sign(one,b) == sign(one,c)) then
                          ! real eigenvalues: reduce to upper triangular form
                          sab = sqrt(abs(b))
                          sac = sqrt(abs(c))
                          p = sign(sab*sac,c)
                          tau = one/sqrt(abs(b + c))
                          a = temp + p
                          d = temp - p
                          b = b - c
                          c = zero
                          cs1 = sab*tau
                          sn1 = sac*tau
                          temp = cs*cs1 - sn*sn1
                          sn = cs*sn1 + sn*cs1
                          cs = temp
                       end if
                    else
                       b = -c
                       c = zero
                       temp = cs
                       cs = -sn
                       sn = temp
                    end if
                 end if
              end if
           end if
           ! store eigenvalues in (rt1r,rt1i) and (rt2r,rt2i).
           rt1r = a
           rt2r = d
           if (c == zero) then
              rt1i = zero
              rt2i = zero
           else
              rt1i = sqrt(abs(b))*sqrt(abs(c))
              rt2i = -rt1i
           end if
           return
     end subroutine la_xlanv2
#endif
#ifdef LA_WITH_QP
     !> QLANV2: computes the Schur factorization of a real 2-by-2 nonsymmetric
     !> matrix in standard form:
     !> [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
     !> [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
     !> where either
     !> 1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
     !> 2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
     !> conjugate eigenvalues.

     pure subroutine la_qlanv2(a,b,c,d,rt1r,rt1i,rt2r,rt2i,cs,sn)
        use la_constants_qp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(inout) :: a,b,c,d
           real(qp),intent(out) :: cs,rt1i,rt1r,rt2i,rt2r,sn
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: multpl = 4.0e+0_qp

           ! Local Scalars
           real(qp) :: aa,bb,bcmax,bcmis,cc,cs1,dd,eps,p,sab,sac,scale,sigma,sn1, &
                     tau,temp,z,safmin,safmn2,safmx2
           integer(ilp) :: count
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sign,sqrt
           ! Executable Statements
           safmin = la_qlamch('S')
           eps = la_qlamch('P')
           safmn2 = la_qlamch('B')**int(log(safmin/eps)/log(la_qlamch('B'))/ &
                     two,KIND=ilp)
           safmx2 = one/safmn2
           if (c == zero) then
              cs = one
              sn = zero
           else if (b == zero) then
              ! swap rows and columns
              cs = zero
              sn = one
              temp = d
              d = a
              a = temp
              b = -c
              c = zero
           else if ((a - d) == zero .and. sign(one,b) /= sign(one,c)) then
              cs = one
              sn = zero
           else
              temp = a - d
              p = half*temp
              bcmax = max(abs(b),abs(c))
              bcmis = min(abs(b),abs(c))*sign(one,b)*sign(one,c)
              scale = max(abs(p),bcmax)
              z = (p/scale)*p + (bcmax/scale)*bcmis
              ! if z is of the order of the machine accuracy, postpone the
              ! decision on the nature of eigenvalues
              if (z >= multpl*eps) then
                 ! real eigenvalues. compute a and d.
                 z = p + sign(sqrt(scale)*sqrt(z),p)
                 a = d + z
                 d = d - (bcmax/z)*bcmis
                 ! compute b and the rotation matrix
                 tau = la_qlapy2(c,z)
                 cs = z/tau
                 sn = c/tau
                 b = b - c
                 c = zero
              else
                 ! complex eigenvalues, or real(almost,KIND=qp) equal eigenvalues.
                 ! make diagonal elements equal.
                 count = 0
                 sigma = b + c
                 10 continue
                 count = count + 1
                 scale = max(abs(temp),abs(sigma))
                 if (scale >= safmx2) then
                    sigma = sigma*safmn2
                    temp = temp*safmn2
                    if (count <= 20) goto 10
                 end if
                 if (scale <= safmn2) then
                    sigma = sigma*safmx2
                    temp = temp*safmx2
                    if (count <= 20) goto 10
                 end if
                 p = half*temp
                 tau = la_qlapy2(sigma,temp)
                 cs = sqrt(half*(one + abs(sigma)/tau))
                 sn = -(p/(tau*cs))*sign(one,sigma)
                 ! compute [ aa  bb ] = [ a  b ] [ cs -sn ]
                         ! [ cc  dd ]   [ c  d ] [ sn  cs ]
                 aa = a*cs + b*sn
                 bb = -a*sn + b*cs
                 cc = c*cs + d*sn
                 dd = -c*sn + d*cs
                 ! compute [ a  b ] = [ cs  sn ] [ aa  bb ]
                         ! [ c  d ]   [-sn  cs ] [ cc  dd ]
                 a = aa*cs + cc*sn
                 b = bb*cs + dd*sn
                 c = -aa*sn + cc*cs
                 d = -bb*sn + dd*cs
                 temp = half*(a + d)
                 a = temp
                 d = temp
                 if (c /= zero) then
                    if (b /= zero) then
                       if (sign(one,b) == sign(one,c)) then
                          ! real eigenvalues: reduce to upper triangular form
                          sab = sqrt(abs(b))
                          sac = sqrt(abs(c))
                          p = sign(sab*sac,c)
                          tau = one/sqrt(abs(b + c))
                          a = temp + p
                          d = temp - p
                          b = b - c
                          c = zero
                          cs1 = sab*tau
                          sn1 = sac*tau
                          temp = cs*cs1 - sn*sn1
                          sn = cs*sn1 + sn*cs1
                          cs = temp
                       end if
                    else
                       b = -c
                       c = zero
                       temp = cs
                       cs = -sn
                       sn = temp
                    end if
                 end if
              end if
           end if
           ! store eigenvalues in (rt1r,rt1i) and (rt2r,rt2i).
           rt1r = a
           rt2r = d
           if (c == zero) then
              rt1i = zero
              rt2i = zero
           else
              rt1i = sqrt(abs(b))*sqrt(abs(c))
              rt2i = -rt1i
           end if
           return
     end subroutine la_qlanv2
#endif

     !> SLAEXC: swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
     !> an upper quasi-triangular matrix T by an orthogonal similarity
     !> transformation.
     !> T must be in Schur canonical form, that is, block upper triangular
     !> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
     !> has its diagonal elements equal and its off-diagonal elements of
     !> opposite sign.

     subroutine la_slaexc(wantq,n,t,ldt,q,ldq,j1,n1,n2,work,info)
        use la_constants_sp,only:zero,one,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,ldq,ldt,n,n1,n2
           ! Array Arguments
           real(sp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldd = 4
           integer(ilp),parameter :: ldx = 2

           ! Local Scalars
           integer(ilp) :: ierr,j2,j3,j4,k,nd
           real(sp) :: cs,dnorm,eps,scale,smlnum,sn,t11,t22,t33,tau,tau1,tau2,temp, &
                     thresh,wi1,wi2,wr1,wr2,xnorm
           ! Local Arrays
           real(sp) :: d(ldd,4),u(3),u1(3),u2(3),x(ldx,2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n == 0 .or. n1 == 0 .or. n2 == 0) return
           if (j1 + n1 > n) return
           j2 = j1 + 1
           j3 = j1 + 2
           j4 = j1 + 3
           if (n1 == 1 .and. n2 == 1) then
              ! swap two 1-by-1 blocks.
              t11 = t(j1,j1)
              t22 = t(j2,j2)
              ! determine the transformation to perform the interchange.
              call la_slartg(t(j1,j2),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (j3 <= n) call la_srot(n - j1 - 1,t(j1,j3),ldt,t(j2,j3),ldt,cs,sn)

              call la_srot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
              t(j1,j1) = t22
              t(j2,j2) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_srot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
           else
              ! swapping involves at least one 2-by-2 block.
              ! copy the diagonal block of order n1+n2 to the local array d
              ! and compute its norm.
              nd = n1 + n2
              call la_slacpy('FULL',nd,nd,t(j1,j1),ldt,d,ldd)
              dnorm = la_slange('MAX',nd,nd,d,ldd,work)
              ! compute machine-dependent threshold for test for accepting
              ! swap.
              eps = la_slamch('P')
              smlnum = la_slamch('S')/eps
              thresh = max(ten*eps*dnorm,smlnum)
              ! solve t11*x - x*t22 = scale*t12 for x.
              call la_slasy2(.false.,.false.,-1,n1,n2,d,ldd,d(n1 + 1,n1 + 1),ldd,d(1, &
                         n1 + 1),ldd,scale,x,ldx,xnorm,ierr)
              ! swap the adjacent diagonal blocks.
              k = n1 + n1 + n2 - 3
              go to(10,20,30) k
              10 continue
              ! n1 = 1, n2 = 2: generate elementary reflector h so that:
              ! ( scale, x11, x12 ) h = ( 0, 0, * )
              u(1) = scale
              u(2) = x(1,1)
              u(3) = x(1,2)
              call la_slarfg(3,u(3),u,1,tau)
              u(3) = one
              t11 = t(j1,j1)
              ! perform swap provisionally on diagonal block in d.
              call la_slarfx('L',3,3,u,tau,d,ldd,work)
              call la_slarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(3,3) - t11)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_slarfx('L',3,n - j1 + 1,u,tau,t(j1,j1),ldt,work)
              call la_slarfx('R',j2,3,u,tau,t(1,j1),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j3,j3) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_slarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              20 continue
              ! n1 = 2, n2 = 1: generate elementary reflector h so that:
              ! h (  -x11 ) = ( * )
                ! (  -x21 ) = ( 0 )
                ! ( scale ) = ( 0 )
              u(1) = -x(1,1)
              u(2) = -x(2,1)
              u(3) = scale
              call la_slarfg(3,u(1),u(2),1,tau)
              u(1) = one
              t33 = t(j3,j3)
              ! perform swap provisionally on diagonal block in d.
              call la_slarfx('L',3,3,u,tau,d,ldd,work)
              call la_slarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(2,1)),abs(d(3,1)),abs(d(1,1) - t33)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_slarfx('R',j3,3,u,tau,t(1,j1),ldt,work)
              call la_slarfx('L',3,n - j1,u,tau,t(j1,j2),ldt,work)
              t(j1,j1) = t33
              t(j2,j1) = zero
              t(j3,j1) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_slarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              30 continue
              ! n1 = 2, n2 = 2: generate elementary reflectors h(1) and h(2) so
              ! that:
              ! h(2) h(1) (  -x11  -x12 ) = (  *  * )
                        ! (  -x21  -x22 )   (  0  * )
                        ! ( scale    0  )   (  0  0 )
                        ! (    0  scale )   (  0  0 )
              u1(1) = -x(1,1)
              u1(2) = -x(2,1)
              u1(3) = scale
              call la_slarfg(3,u1(1),u1(2),1,tau1)
              u1(1) = one
              temp = -tau1*(x(1,2) + u1(2)*x(2,2))
              u2(1) = -temp*u1(2) - x(2,2)
              u2(2) = -temp*u1(3)
              u2(3) = scale
              call la_slarfg(3,u2(1),u2(2),1,tau2)
              u2(1) = one
              ! perform swap provisionally on diagonal block in d.
              call la_slarfx('L',3,4,u1,tau1,d,ldd,work)
              call la_slarfx('R',4,3,u1,tau1,d,ldd,work)
              call la_slarfx('L',3,4,u2,tau2,d(2,1),ldd,work)
              call la_slarfx('R',4,3,u2,tau2,d(1,2),ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(4,1)),abs(d(4,2))) &
                        > thresh) go to 50
              ! accept swap: apply transformation to the entire matrix t.
              call la_slarfx('L',3,n - j1 + 1,u1,tau1,t(j1,j1),ldt,work)
              call la_slarfx('R',j4,3,u1,tau1,t(1,j1),ldt,work)
              call la_slarfx('L',3,n - j1 + 1,u2,tau2,t(j2,j1),ldt,work)
              call la_slarfx('R',j4,3,u2,tau2,t(1,j2),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j4,j1) = zero
              t(j4,j2) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_slarfx('R',n,3,u1,tau1,q(1,j1),ldq,work)
                 call la_slarfx('R',n,3,u2,tau2,q(1,j2),ldq,work)
              end if
              40 continue
              if (n2 == 2) then
                 ! standardize new 2-by-2 block t11
                 call la_slanv2(t(j1,j1),t(j1,j2),t(j2,j1),t(j2,j2),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 call la_srot(n - j1 - 1,t(j1,j1 + 2),ldt,t(j2,j1 + 2),ldt,cs,sn)
                 call la_srot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
                 if (wantq) call la_srot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
              if (n1 == 2) then
                 ! standardize new 2-by-2 block t22
                 j3 = j1 + n2
                 j4 = j3 + 1
                 call la_slanv2(t(j3,j3),t(j3,j4),t(j4,j3),t(j4,j4),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 if (j3 + 2 <= n) call la_srot(n - j3 - 1,t(j3,j3 + 2),ldt,t(j4,j3 + 2),ldt,cs, &
                            sn)
                 call la_srot(j3 - 1,t(1,j3),1,t(1,j4),1,cs,sn)
                 if (wantq) call la_srot(n,q(1,j3),1,q(1,j4),1,cs,sn)
              end if
           end if
           return
           ! exit with info = 1 if swap was rejected.
50 info = 1
           return
     end subroutine la_slaexc
     !> DLAEXC: swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
     !> an upper quasi-triangular matrix T by an orthogonal similarity
     !> transformation.
     !> T must be in Schur canonical form, that is, block upper triangular
     !> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
     !> has its diagonal elements equal and its off-diagonal elements of
     !> opposite sign.

     subroutine la_dlaexc(wantq,n,t,ldt,q,ldq,j1,n1,n2,work,info)
        use la_constants_dp,only:zero,one,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,ldq,ldt,n,n1,n2
           ! Array Arguments
           real(dp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldd = 4
           integer(ilp),parameter :: ldx = 2

           ! Local Scalars
           integer(ilp) :: ierr,j2,j3,j4,k,nd
           real(dp) :: cs,dnorm,eps,scale,smlnum,sn,t11,t22,t33,tau,tau1,tau2,temp, &
                     thresh,wi1,wi2,wr1,wr2,xnorm
           ! Local Arrays
           real(dp) :: d(ldd,4),u(3),u1(3),u2(3),x(ldx,2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n == 0 .or. n1 == 0 .or. n2 == 0) return
           if (j1 + n1 > n) return
           j2 = j1 + 1
           j3 = j1 + 2
           j4 = j1 + 3
           if (n1 == 1 .and. n2 == 1) then
              ! swap two 1-by-1 blocks.
              t11 = t(j1,j1)
              t22 = t(j2,j2)
              ! determine the transformation to perform the interchange.
              call la_dlartg(t(j1,j2),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (j3 <= n) call la_drot(n - j1 - 1,t(j1,j3),ldt,t(j2,j3),ldt,cs,sn)

              call la_drot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
              t(j1,j1) = t22
              t(j2,j2) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_drot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
           else
              ! swapping involves at least one 2-by-2 block.
              ! copy the diagonal block of order n1+n2 to the local array d
              ! and compute its norm.
              nd = n1 + n2
              call la_dlacpy('FULL',nd,nd,t(j1,j1),ldt,d,ldd)
              dnorm = la_dlange('MAX',nd,nd,d,ldd,work)
              ! compute machine-dependent threshold for test for accepting
              ! swap.
              eps = la_dlamch('P')
              smlnum = la_dlamch('S')/eps
              thresh = max(ten*eps*dnorm,smlnum)
              ! solve t11*x - x*t22 = scale*t12 for x.
              call la_dlasy2(.false.,.false.,-1,n1,n2,d,ldd,d(n1 + 1,n1 + 1),ldd,d(1, &
                         n1 + 1),ldd,scale,x,ldx,xnorm,ierr)
              ! swap the adjacent diagonal blocks.
              k = n1 + n1 + n2 - 3
              go to(10,20,30) k
              10 continue
              ! n1 = 1, n2 = 2: generate elementary reflector h so that:
              ! ( scale, x11, x12 ) h = ( 0, 0, * )
              u(1) = scale
              u(2) = x(1,1)
              u(3) = x(1,2)
              call la_dlarfg(3,u(3),u,1,tau)
              u(3) = one
              t11 = t(j1,j1)
              ! perform swap provisionally on diagonal block in d.
              call la_dlarfx('L',3,3,u,tau,d,ldd,work)
              call la_dlarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(3,3) - t11)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_dlarfx('L',3,n - j1 + 1,u,tau,t(j1,j1),ldt,work)
              call la_dlarfx('R',j2,3,u,tau,t(1,j1),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j3,j3) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_dlarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              20 continue
              ! n1 = 2, n2 = 1: generate elementary reflector h so that:
              ! h (  -x11 ) = ( * )
                ! (  -x21 ) = ( 0 )
                ! ( scale ) = ( 0 )
              u(1) = -x(1,1)
              u(2) = -x(2,1)
              u(3) = scale
              call la_dlarfg(3,u(1),u(2),1,tau)
              u(1) = one
              t33 = t(j3,j3)
              ! perform swap provisionally on diagonal block in d.
              call la_dlarfx('L',3,3,u,tau,d,ldd,work)
              call la_dlarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(2,1)),abs(d(3,1)),abs(d(1,1) - t33)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_dlarfx('R',j3,3,u,tau,t(1,j1),ldt,work)
              call la_dlarfx('L',3,n - j1,u,tau,t(j1,j2),ldt,work)
              t(j1,j1) = t33
              t(j2,j1) = zero
              t(j3,j1) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_dlarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              30 continue
              ! n1 = 2, n2 = 2: generate elementary reflectors h(1) and h(2) so
              ! that:
              ! h(2) h(1) (  -x11  -x12 ) = (  *  * )
                        ! (  -x21  -x22 )   (  0  * )
                        ! ( scale    0  )   (  0  0 )
                        ! (    0  scale )   (  0  0 )
              u1(1) = -x(1,1)
              u1(2) = -x(2,1)
              u1(3) = scale
              call la_dlarfg(3,u1(1),u1(2),1,tau1)
              u1(1) = one
              temp = -tau1*(x(1,2) + u1(2)*x(2,2))
              u2(1) = -temp*u1(2) - x(2,2)
              u2(2) = -temp*u1(3)
              u2(3) = scale
              call la_dlarfg(3,u2(1),u2(2),1,tau2)
              u2(1) = one
              ! perform swap provisionally on diagonal block in d.
              call la_dlarfx('L',3,4,u1,tau1,d,ldd,work)
              call la_dlarfx('R',4,3,u1,tau1,d,ldd,work)
              call la_dlarfx('L',3,4,u2,tau2,d(2,1),ldd,work)
              call la_dlarfx('R',4,3,u2,tau2,d(1,2),ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(4,1)),abs(d(4,2))) &
                        > thresh) go to 50
              ! accept swap: apply transformation to the entire matrix t.
              call la_dlarfx('L',3,n - j1 + 1,u1,tau1,t(j1,j1),ldt,work)
              call la_dlarfx('R',j4,3,u1,tau1,t(1,j1),ldt,work)
              call la_dlarfx('L',3,n - j1 + 1,u2,tau2,t(j2,j1),ldt,work)
              call la_dlarfx('R',j4,3,u2,tau2,t(1,j2),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j4,j1) = zero
              t(j4,j2) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_dlarfx('R',n,3,u1,tau1,q(1,j1),ldq,work)
                 call la_dlarfx('R',n,3,u2,tau2,q(1,j2),ldq,work)
              end if
              40 continue
              if (n2 == 2) then
                 ! standardize new 2-by-2 block t11
                 call la_dlanv2(t(j1,j1),t(j1,j2),t(j2,j1),t(j2,j2),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 call la_drot(n - j1 - 1,t(j1,j1 + 2),ldt,t(j2,j1 + 2),ldt,cs,sn)
                 call la_drot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
                 if (wantq) call la_drot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
              if (n1 == 2) then
                 ! standardize new 2-by-2 block t22
                 j3 = j1 + n2
                 j4 = j3 + 1
                 call la_dlanv2(t(j3,j3),t(j3,j4),t(j4,j3),t(j4,j4),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 if (j3 + 2 <= n) call la_drot(n - j3 - 1,t(j3,j3 + 2),ldt,t(j4,j3 + 2),ldt,cs, &
                            sn)
                 call la_drot(j3 - 1,t(1,j3),1,t(1,j4),1,cs,sn)
                 if (wantq) call la_drot(n,q(1,j3),1,q(1,j4),1,cs,sn)
              end if
           end if
           return
           ! exit with info = 1 if swap was rejected.
           50 continue
           info = 1
           return
     end subroutine la_dlaexc
#ifdef LA_WITH_XDP
     !> XLAEXC: swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
     !> an upper quasi-triangular matrix T by an orthogonal similarity
     !> transformation.
     !> T must be in Schur canonical form, that is, block upper triangular
     !> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
     !> has its diagonal elements equal and its off-diagonal elements of
     !> opposite sign.

     subroutine la_xlaexc(wantq,n,t,ldt,q,ldq,j1,n1,n2,work,info)
        use la_constants_xdp,only:zero,one,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,ldq,ldt,n,n1,n2
           ! Array Arguments
           real(xdp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldd = 4
           integer(ilp),parameter :: ldx = 2

           ! Local Scalars
           integer(ilp) :: ierr,j2,j3,j4,k,nd
           real(xdp) :: cs,dnorm,eps,scale,smlnum,sn,t11,t22,t33,tau,tau1,tau2,temp, &
                     thresh,wi1,wi2,wr1,wr2,xnorm
           ! Local Arrays
           real(xdp) :: d(ldd,4),u(3),u1(3),u2(3),x(ldx,2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n == 0 .or. n1 == 0 .or. n2 == 0) return
           if (j1 + n1 > n) return
           j2 = j1 + 1
           j3 = j1 + 2
           j4 = j1 + 3
           if (n1 == 1 .and. n2 == 1) then
              ! swap two 1-by-1 blocks.
              t11 = t(j1,j1)
              t22 = t(j2,j2)
              ! determine the transformation to perform the interchange.
              call la_xlartg(t(j1,j2),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (j3 <= n) call la_xrot(n - j1 - 1,t(j1,j3),ldt,t(j2,j3),ldt,cs,sn)

              call la_xrot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
              t(j1,j1) = t22
              t(j2,j2) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_xrot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
           else
              ! swapping involves at least one 2-by-2 block.
              ! copy the diagonal block of order n1+n2 to the local array d
              ! and compute its norm.
              nd = n1 + n2
              call la_xlacpy('FULL',nd,nd,t(j1,j1),ldt,d,ldd)
              dnorm = la_xlange('MAX',nd,nd,d,ldd,work)
              ! compute machine-dependent threshold for test for accepting
              ! swap.
              eps = la_xlamch('P')
              smlnum = la_xlamch('S')/eps
              thresh = max(ten*eps*dnorm,smlnum)
              ! solve t11*x - x*t22 = scale*t12 for x.
              call la_xlasy2(.false.,.false.,-1,n1,n2,d,ldd,d(n1 + 1,n1 + 1),ldd,d(1, &
                         n1 + 1),ldd,scale,x,ldx,xnorm,ierr)
              ! swap the adjacent diagonal blocks.
              k = n1 + n1 + n2 - 3
              go to(10,20,30) k
              10 continue
              ! n1 = 1, n2 = 2: generate elementary reflector h so that:
              ! ( scale, x11, x12 ) h = ( 0, 0, * )
              u(1) = scale
              u(2) = x(1,1)
              u(3) = x(1,2)
              call la_xlarfg(3,u(3),u,1,tau)
              u(3) = one
              t11 = t(j1,j1)
              ! perform swap provisionally on diagonal block in d.
              call la_xlarfx('L',3,3,u,tau,d,ldd,work)
              call la_xlarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(3,3) - t11)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_xlarfx('L',3,n - j1 + 1,u,tau,t(j1,j1),ldt,work)
              call la_xlarfx('R',j2,3,u,tau,t(1,j1),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j3,j3) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_xlarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              20 continue
              ! n1 = 2, n2 = 1: generate elementary reflector h so that:
              ! h (  -x11 ) = ( * )
                ! (  -x21 ) = ( 0 )
                ! ( scale ) = ( 0 )
              u(1) = -x(1,1)
              u(2) = -x(2,1)
              u(3) = scale
              call la_xlarfg(3,u(1),u(2),1,tau)
              u(1) = one
              t33 = t(j3,j3)
              ! perform swap provisionally on diagonal block in d.
              call la_xlarfx('L',3,3,u,tau,d,ldd,work)
              call la_xlarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(2,1)),abs(d(3,1)),abs(d(1,1) - t33)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_xlarfx('R',j3,3,u,tau,t(1,j1),ldt,work)
              call la_xlarfx('L',3,n - j1,u,tau,t(j1,j2),ldt,work)
              t(j1,j1) = t33
              t(j2,j1) = zero
              t(j3,j1) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_xlarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              30 continue
              ! n1 = 2, n2 = 2: generate elementary reflectors h(1) and h(2) so
              ! that:
              ! h(2) h(1) (  -x11  -x12 ) = (  *  * )
                        ! (  -x21  -x22 )   (  0  * )
                        ! ( scale    0  )   (  0  0 )
                        ! (    0  scale )   (  0  0 )
              u1(1) = -x(1,1)
              u1(2) = -x(2,1)
              u1(3) = scale
              call la_xlarfg(3,u1(1),u1(2),1,tau1)
              u1(1) = one
              temp = -tau1*(x(1,2) + u1(2)*x(2,2))
              u2(1) = -temp*u1(2) - x(2,2)
              u2(2) = -temp*u1(3)
              u2(3) = scale
              call la_xlarfg(3,u2(1),u2(2),1,tau2)
              u2(1) = one
              ! perform swap provisionally on diagonal block in d.
              call la_xlarfx('L',3,4,u1,tau1,d,ldd,work)
              call la_xlarfx('R',4,3,u1,tau1,d,ldd,work)
              call la_xlarfx('L',3,4,u2,tau2,d(2,1),ldd,work)
              call la_xlarfx('R',4,3,u2,tau2,d(1,2),ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(4,1)),abs(d(4,2))) &
                        > thresh) go to 50
              ! accept swap: apply transformation to the entire matrix t.
              call la_xlarfx('L',3,n - j1 + 1,u1,tau1,t(j1,j1),ldt,work)
              call la_xlarfx('R',j4,3,u1,tau1,t(1,j1),ldt,work)
              call la_xlarfx('L',3,n - j1 + 1,u2,tau2,t(j2,j1),ldt,work)
              call la_xlarfx('R',j4,3,u2,tau2,t(1,j2),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j4,j1) = zero
              t(j4,j2) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_xlarfx('R',n,3,u1,tau1,q(1,j1),ldq,work)
                 call la_xlarfx('R',n,3,u2,tau2,q(1,j2),ldq,work)
              end if
              40 continue
              if (n2 == 2) then
                 ! standardize new 2-by-2 block t11
                 call la_xlanv2(t(j1,j1),t(j1,j2),t(j2,j1),t(j2,j2),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 call la_xrot(n - j1 - 1,t(j1,j1 + 2),ldt,t(j2,j1 + 2),ldt,cs,sn)
                 call la_xrot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
                 if (wantq) call la_xrot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
              if (n1 == 2) then
                 ! standardize new 2-by-2 block t22
                 j3 = j1 + n2
                 j4 = j3 + 1
                 call la_xlanv2(t(j3,j3),t(j3,j4),t(j4,j3),t(j4,j4),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 if (j3 + 2 <= n) call la_xrot(n - j3 - 1,t(j3,j3 + 2),ldt,t(j4,j3 + 2),ldt,cs, &
                            sn)
                 call la_xrot(j3 - 1,t(1,j3),1,t(1,j4),1,cs,sn)
                 if (wantq) call la_xrot(n,q(1,j3),1,q(1,j4),1,cs,sn)
              end if
           end if
           return
           ! exit with info = 1 if swap was rejected.
           50 continue
           info = 1
           return
     end subroutine la_xlaexc
#endif
#ifdef LA_WITH_QP
     !> QLAEXC: swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
     !> an upper quasi-triangular matrix T by an orthogonal similarity
     !> transformation.
     !> T must be in Schur canonical form, that is, block upper triangular
     !> with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
     !> has its diagonal elements equal and its off-diagonal elements of
     !> opposite sign.

     subroutine la_qlaexc(wantq,n,t,ldt,q,ldq,j1,n1,n2,work,info)
        use la_constants_qp,only:zero,one,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantq
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: j1,ldq,ldt,n,n1,n2
           ! Array Arguments
           real(qp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: ldd = 4
           integer(ilp),parameter :: ldx = 2

           ! Local Scalars
           integer(ilp) :: ierr,j2,j3,j4,k,nd
           real(qp) :: cs,dnorm,eps,scale,smlnum,sn,t11,t22,t33,tau,tau1,tau2,temp, &
                     thresh,wi1,wi2,wr1,wr2,xnorm
           ! Local Arrays
           real(qp) :: d(ldd,4),u(3),u1(3),u2(3),x(ldx,2)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n == 0 .or. n1 == 0 .or. n2 == 0) return
           if (j1 + n1 > n) return
           j2 = j1 + 1
           j3 = j1 + 2
           j4 = j1 + 3
           if (n1 == 1 .and. n2 == 1) then
              ! swap two 1-by-1 blocks.
              t11 = t(j1,j1)
              t22 = t(j2,j2)
              ! determine the transformation to perform the interchange.
              call la_qlartg(t(j1,j2),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (j3 <= n) call la_qrot(n - j1 - 1,t(j1,j3),ldt,t(j2,j3),ldt,cs,sn)

              call la_qrot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
              t(j1,j1) = t22
              t(j2,j2) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_qrot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
           else
              ! swapping involves at least one 2-by-2 block.
              ! copy the diagonal block of order n1+n2 to the local array d
              ! and compute its norm.
              nd = n1 + n2
              call la_qlacpy('FULL',nd,nd,t(j1,j1),ldt,d,ldd)
              dnorm = la_qlange('MAX',nd,nd,d,ldd,work)
              ! compute machine-dependent threshold for test for accepting
              ! swap.
              eps = la_qlamch('P')
              smlnum = la_qlamch('S')/eps
              thresh = max(ten*eps*dnorm,smlnum)
              ! solve t11*x - x*t22 = scale*t12 for x.
              call la_qlasy2(.false.,.false.,-1,n1,n2,d,ldd,d(n1 + 1,n1 + 1),ldd,d(1, &
                         n1 + 1),ldd,scale,x,ldx,xnorm,ierr)
              ! swap the adjacent diagonal blocks.
              k = n1 + n1 + n2 - 3
              go to(10,20,30) k
              10 continue
              ! n1 = 1, n2 = 2: generate elementary reflector h so that:
              ! ( scale, x11, x12 ) h = ( 0, 0, * )
              u(1) = scale
              u(2) = x(1,1)
              u(3) = x(1,2)
              call la_qlarfg(3,u(3),u,1,tau)
              u(3) = one
              t11 = t(j1,j1)
              ! perform swap provisionally on diagonal block in d.
              call la_qlarfx('L',3,3,u,tau,d,ldd,work)
              call la_qlarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(3,3) - t11)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_qlarfx('L',3,n - j1 + 1,u,tau,t(j1,j1),ldt,work)
              call la_qlarfx('R',j2,3,u,tau,t(1,j1),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j3,j3) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_qlarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              20 continue
              ! n1 = 2, n2 = 1: generate elementary reflector h so that:
              ! h (  -x11 ) = ( * )
                ! (  -x21 ) = ( 0 )
                ! ( scale ) = ( 0 )
              u(1) = -x(1,1)
              u(2) = -x(2,1)
              u(3) = scale
              call la_qlarfg(3,u(1),u(2),1,tau)
              u(1) = one
              t33 = t(j3,j3)
              ! perform swap provisionally on diagonal block in d.
              call la_qlarfx('L',3,3,u,tau,d,ldd,work)
              call la_qlarfx('R',3,3,u,tau,d,ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(2,1)),abs(d(3,1)),abs(d(1,1) - t33)) > thresh) go to &
                        50
              ! accept swap: apply transformation to the entire matrix t.
              call la_qlarfx('R',j3,3,u,tau,t(1,j1),ldt,work)
              call la_qlarfx('L',3,n - j1,u,tau,t(j1,j2),ldt,work)
              t(j1,j1) = t33
              t(j2,j1) = zero
              t(j3,j1) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_qlarfx('R',n,3,u,tau,q(1,j1),ldq,work)
              end if
              go to 40
              30 continue
              ! n1 = 2, n2 = 2: generate elementary reflectors h(1) and h(2) so
              ! that:
              ! h(2) h(1) (  -x11  -x12 ) = (  *  * )
                        ! (  -x21  -x22 )   (  0  * )
                        ! ( scale    0  )   (  0  0 )
                        ! (    0  scale )   (  0  0 )
              u1(1) = -x(1,1)
              u1(2) = -x(2,1)
              u1(3) = scale
              call la_qlarfg(3,u1(1),u1(2),1,tau1)
              u1(1) = one
              temp = -tau1*(x(1,2) + u1(2)*x(2,2))
              u2(1) = -temp*u1(2) - x(2,2)
              u2(2) = -temp*u1(3)
              u2(3) = scale
              call la_qlarfg(3,u2(1),u2(2),1,tau2)
              u2(1) = one
              ! perform swap provisionally on diagonal block in d.
              call la_qlarfx('L',3,4,u1,tau1,d,ldd,work)
              call la_qlarfx('R',4,3,u1,tau1,d,ldd,work)
              call la_qlarfx('L',3,4,u2,tau2,d(2,1),ldd,work)
              call la_qlarfx('R',4,3,u2,tau2,d(1,2),ldd,work)
              ! test whether to reject swap.
              if (max(abs(d(3,1)),abs(d(3,2)),abs(d(4,1)),abs(d(4,2))) &
                        > thresh) go to 50
              ! accept swap: apply transformation to the entire matrix t.
              call la_qlarfx('L',3,n - j1 + 1,u1,tau1,t(j1,j1),ldt,work)
              call la_qlarfx('R',j4,3,u1,tau1,t(1,j1),ldt,work)
              call la_qlarfx('L',3,n - j1 + 1,u2,tau2,t(j2,j1),ldt,work)
              call la_qlarfx('R',j4,3,u2,tau2,t(1,j2),ldt,work)
              t(j3,j1) = zero
              t(j3,j2) = zero
              t(j4,j1) = zero
              t(j4,j2) = zero
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_qlarfx('R',n,3,u1,tau1,q(1,j1),ldq,work)
                 call la_qlarfx('R',n,3,u2,tau2,q(1,j2),ldq,work)
              end if
              40 continue
              if (n2 == 2) then
                 ! standardize new 2-by-2 block t11
                 call la_qlanv2(t(j1,j1),t(j1,j2),t(j2,j1),t(j2,j2),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 call la_qrot(n - j1 - 1,t(j1,j1 + 2),ldt,t(j2,j1 + 2),ldt,cs,sn)
                 call la_qrot(j1 - 1,t(1,j1),1,t(1,j2),1,cs,sn)
                 if (wantq) call la_qrot(n,q(1,j1),1,q(1,j2),1,cs,sn)
              end if
              if (n1 == 2) then
                 ! standardize new 2-by-2 block t22
                 j3 = j1 + n2
                 j4 = j3 + 1
                 call la_qlanv2(t(j3,j3),t(j3,j4),t(j4,j3),t(j4,j4),wr1,wi1, &
                           wr2,wi2,cs,sn)
                 if (j3 + 2 <= n) call la_qrot(n - j3 - 1,t(j3,j3 + 2),ldt,t(j4,j3 + 2),ldt,cs, &
                            sn)
                 call la_qrot(j3 - 1,t(1,j3),1,t(1,j4),1,cs,sn)
                 if (wantq) call la_qrot(n,q(1,j3),1,q(1,j4),1,cs,sn)
              end if
           end if
           return
           ! exit with info = 1 if swap was rejected.
           50 continue
           info = 1
           return
     end subroutine la_qlaexc
#endif

     !> STREXC: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
     !> moved to row ILST.
     !> The real Schur form T is reordered by an orthogonal similarity
     !> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
     !> is updated by postmultiplying it with Z.
     !> T must be in Schur canonical form (as returned by SHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_strexc(compq,n,t,ldt,q,ldq,ifst,ilst,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldt,n
           ! Array Arguments
           real(sp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: here,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test the input arguments.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. wantq .and. .not. la_lsame(compq,'N')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('STREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of specified block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ifst > 1) then
              if (t(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (t(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ilst > 1) then
              if (t(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (t(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap block with next one below
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (t(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_slaexc(wantq,n,t,ldt,q,ldq,here,nbf,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (t(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_slaexc(wantq,n,t,ldt,q,ldq,here + 1,1,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_slaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                    here = here + 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_slaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2 by 2 block did split
                       call la_slaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_slaexc(wantq,n,t,ldt,q,ldq,here + 1,1,1,work,info)

                       here = here + 2
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap block with next one above
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_slaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,nbf,work, &
                           info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_slaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,1,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_slaexc(wantq,n,t,ldt,q,ldq,here,nbnext,1,work,info)

                    here = here - 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_slaexc(wantq,n,t,ldt,q,ldq,here - 1,2,1,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2 by 2 block did split
                       call la_slaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_slaexc(wantq,n,t,ldt,q,ldq,here - 1,1,1,work,info)

                       here = here - 2
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           return
     end subroutine la_strexc
     !> DTREXC: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
     !> moved to row ILST.
     !> The real Schur form T is reordered by an orthogonal similarity
     !> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
     !> is updated by postmultiplying it with Z.
     !> T must be in Schur canonical form (as returned by DHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_dtrexc(compq,n,t,ldt,q,ldq,ifst,ilst,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldt,n
           ! Array Arguments
           real(dp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: here,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test the input arguments.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. wantq .and. .not. la_lsame(compq,'N')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of specified block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ifst > 1) then
              if (t(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (t(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ilst > 1) then
              if (t(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (t(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap block with next one below
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (t(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_dlaexc(wantq,n,t,ldt,q,ldq,here,nbf,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (t(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_dlaexc(wantq,n,t,ldt,q,ldq,here + 1,1,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_dlaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                    here = here + 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_dlaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2 by 2 block did split
                       call la_dlaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_dlaexc(wantq,n,t,ldt,q,ldq,here + 1,1,1,work,info)

                       here = here + 2
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap block with next one above
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_dlaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,nbf,work, &
                           info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_dlaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,1,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_dlaexc(wantq,n,t,ldt,q,ldq,here,nbnext,1,work,info)

                    here = here - 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_dlaexc(wantq,n,t,ldt,q,ldq,here - 1,2,1,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2 by 2 block did split
                       call la_dlaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_dlaexc(wantq,n,t,ldt,q,ldq,here - 1,1,1,work,info)

                       here = here - 2
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           return
     end subroutine la_dtrexc
#ifdef LA_WITH_XDP
     !> XTREXC: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
     !> moved to row ILST.
     !> The real Schur form T is reordered by an orthogonal similarity
     !> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
     !> is updated by postmultiplying it with Z.
     !> T must be in Schur canonical form (as returned by XHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_xtrexc(compq,n,t,ldt,q,ldq,ifst,ilst,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldt,n
           ! Array Arguments
           real(xdp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: here,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test the input arguments.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. wantq .and. .not. la_lsame(compq,'N')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('XTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of specified block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ifst > 1) then
              if (t(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (t(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ilst > 1) then
              if (t(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (t(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap block with next one below
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (t(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_xlaexc(wantq,n,t,ldt,q,ldq,here,nbf,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (t(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_xlaexc(wantq,n,t,ldt,q,ldq,here + 1,1,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_xlaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                    here = here + 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_xlaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2 by 2 block did split
                       call la_xlaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_xlaexc(wantq,n,t,ldt,q,ldq,here + 1,1,1,work,info)

                       here = here + 2
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap block with next one above
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_xlaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,nbf,work, &
                           info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_xlaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,1,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_xlaexc(wantq,n,t,ldt,q,ldq,here,nbnext,1,work,info)

                    here = here - 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_xlaexc(wantq,n,t,ldt,q,ldq,here - 1,2,1,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2 by 2 block did split
                       call la_xlaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_xlaexc(wantq,n,t,ldt,q,ldq,here - 1,1,1,work,info)

                       here = here - 2
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           return
     end subroutine la_xtrexc
#endif
#ifdef LA_WITH_QP
     !> QTREXC: reorders the real Schur factorization of a real matrix
     !> A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
     !> moved to row ILST.
     !> The real Schur form T is reordered by an orthogonal similarity
     !> transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
     !> is updated by postmultiplying it with Z.
     !> T must be in Schur canonical form (as returned by QHSEQR), that is,
     !> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
     !> 2-by-2 diagonal block has its diagonal elements equal and its
     !> off-diagonal elements of opposite sign.

     subroutine la_qtrexc(compq,n,t,ldt,q,ldq,ifst,ilst,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(inout) :: ifst,ilst
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldt,n
           ! Array Arguments
           real(qp),intent(inout) :: q(ldq,*),t(ldt,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: here,nbf,nbl,nbnext
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! decode and test the input arguments.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. wantq .and. .not. la_lsame(compq,'N')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           ! determine the first row of specified block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ifst > 1) then
              if (t(ifst,ifst - 1) /= zero) ifst = ifst - 1
           end if
           nbf = 1
           if (ifst < n) then
              if (t(ifst + 1,ifst) /= zero) nbf = 2
           end if
           ! determine the first row of the final block
           ! and find out it is 1 by 1 or 2 by 2.
           if (ilst > 1) then
              if (t(ilst,ilst - 1) /= zero) ilst = ilst - 1
           end if
           nbl = 1
           if (ilst < n) then
              if (t(ilst + 1,ilst) /= zero) nbl = 2
           end if
           if (ifst == ilst) return
           if (ifst < ilst) then
              ! update ilst
              if (nbf == 2 .and. nbl == 1) ilst = ilst - 1
              if (nbf == 1 .and. nbl == 2) ilst = ilst + 1
              here = ifst
              10 continue
              ! swap block with next one below
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here + nbf + 1 <= n) then
                    if (t(here + nbf + 1,here + nbf) /= zero) nbnext = 2
                 end if
                 call la_qlaexc(wantq,n,t,ldt,q,ldq,here,nbf,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here + nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here + 3 <= n) then
                    if (t(here + 3,here + 2) /= zero) nbnext = 2
                 end if
                 call la_qlaexc(wantq,n,t,ldt,q,ldq,here + 1,1,nbnext,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_qlaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                    here = here + 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here + 2,here + 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_qlaexc(wantq,n,t,ldt,q,ldq,here,1,nbnext,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here + 2
                    else
                       ! 2 by 2 block did split
                       call la_qlaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_qlaexc(wantq,n,t,ldt,q,ldq,here + 1,1,1,work,info)

                       here = here + 2
                    end if
                 end if
              end if
              if (here < ilst) go to 10
           else
              here = ifst
              20 continue
              ! swap block with next one above
              if (nbf == 1 .or. nbf == 2) then
                 ! current block either 1 by 1 or 2 by 2
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_qlaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,nbf,work, &
                           info)
                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 here = here - nbnext
                 ! test if 2 by 2 block breaks into two 1 by 1 blocks
                 if (nbf == 2) then
                    if (t(here + 1,here) == zero) nbf = 3
                 end if
              else
                 ! current block consists of two 1 by 1 blocks each of which
                 ! must be swapped individually
                 nbnext = 1
                 if (here >= 3) then
                    if (t(here - 1,here - 2) /= zero) nbnext = 2
                 end if
                 call la_qlaexc(wantq,n,t,ldt,q,ldq,here - nbnext,nbnext,1,work,info)

                 if (info /= 0) then
                    ilst = here
                    return
                 end if
                 if (nbnext == 1) then
                    ! swap two 1 by 1 blocks, no problems possible
                    call la_qlaexc(wantq,n,t,ldt,q,ldq,here,nbnext,1,work,info)

                    here = here - 1
                 else
                    ! recompute nbnext in case 2 by 2 split
                    if (t(here,here - 1) == zero) nbnext = 1
                    if (nbnext == 2) then
                       ! 2 by 2 block did not split
                       call la_qlaexc(wantq,n,t,ldt,q,ldq,here - 1,2,1,work,info)

                       if (info /= 0) then
                          ilst = here
                          return
                       end if
                       here = here - 2
                    else
                       ! 2 by 2 block did split
                       call la_qlaexc(wantq,n,t,ldt,q,ldq,here,1,1,work,info)

                       call la_qlaexc(wantq,n,t,ldt,q,ldq,here - 1,1,1,work,info)

                       here = here - 2
                    end if
                 end if
              end if
              if (here > ilst) go to 20
           end if
           ilst = here
           return
     end subroutine la_qtrexc
#endif

     !> CTREXC: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
     !> is moved to row ILST.
     !> The Schur form T is reordered by a unitary similarity transformation
     !> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
     !> postmultplying it with Z.

     pure subroutine la_ctrexc(compq,n,t,ldt,q,ldq,ifst,ilst,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(in) :: ifst,ilst,ldq,ldt,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(sp),intent(inout) :: q(ldq,*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: k,m1,m2,m3
           real(sp) :: cs
           complex(sp) :: sn,t11,t22,temp
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! decode and test the input parameters.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1 .or. ifst == ilst) return
           if (ifst < ilst) then
              ! move the ifst-th diagonal element forward down the diagonal.
              m1 = 0
              m2 = -1
              m3 = 1
           else
              ! move the ifst-th diagonal element backward up the diagonal.
              m1 = -1
              m2 = 0
              m3 = -1
           end if
           do k = ifst + m1,ilst + m2,m3
              ! interchange the k-th and (k+1)-th diagonal elements.
              t11 = t(k,k)
              t22 = t(k + 1,k + 1)
              ! determine the transformation to perform the interchange.
              call la_clartg(t(k,k + 1),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (k + 2 <= n) call la_crot(n - k - 1,t(k,k + 2),ldt,t(k + 1,k + 2),ldt,cs,sn)

              call la_crot(k - 1,t(1,k),1,t(1,k + 1),1,cs,conjg(sn))
              t(k,k) = t22
              t(k + 1,k + 1) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_crot(n,q(1,k),1,q(1,k + 1),1,cs,conjg(sn))
              end if
           end do
           return
     end subroutine la_ctrexc
     !> ZTREXC: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
     !> is moved to row ILST.
     !> The Schur form T is reordered by a unitary similarity transformation
     !> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
     !> postmultplying it with Z.

     pure subroutine la_ztrexc(compq,n,t,ldt,q,ldq,ifst,ilst,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(in) :: ifst,ilst,ldq,ldt,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(dp),intent(inout) :: q(ldq,*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: k,m1,m2,m3
           real(dp) :: cs
           complex(dp) :: sn,t11,t22,temp
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! decode and test the input parameters.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1 .or. ifst == ilst) return
           if (ifst < ilst) then
              ! move the ifst-th diagonal element forward down the diagonal.
              m1 = 0
              m2 = -1
              m3 = 1
           else
              ! move the ifst-th diagonal element backward up the diagonal.
              m1 = -1
              m2 = 0
              m3 = -1
           end if
           do k = ifst + m1,ilst + m2,m3
              ! interchange the k-th and (k+1)-th diagonal elements.
              t11 = t(k,k)
              t22 = t(k + 1,k + 1)
              ! determine the transformation to perform the interchange.
              call la_zlartg(t(k,k + 1),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (k + 2 <= n) call la_zrot(n - k - 1,t(k,k + 2),ldt,t(k + 1,k + 2),ldt,cs,sn)

              call la_zrot(k - 1,t(1,k),1,t(1,k + 1),1,cs,conjg(sn))
              t(k,k) = t22
              t(k + 1,k + 1) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_zrot(n,q(1,k),1,q(1,k + 1),1,cs,conjg(sn))
              end if
           end do
           return
     end subroutine la_ztrexc
#ifdef LA_WITH_XDP
     !> YTREXC: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
     !> is moved to row ILST.
     !> The Schur form T is reordered by a unitary similarity transformation
     !> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
     !> postmultplying it with Z.

     pure subroutine la_ytrexc(compq,n,t,ldt,q,ldq,ifst,ilst,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(in) :: ifst,ilst,ldq,ldt,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(xdp),intent(inout) :: q(ldq,*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: k,m1,m2,m3
           real(xdp) :: cs
           complex(xdp) :: sn,t11,t22,temp
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! decode and test the input parameters.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('YTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1 .or. ifst == ilst) return
           if (ifst < ilst) then
              ! move the ifst-th diagonal element forward down the diagonal.
              m1 = 0
              m2 = -1
              m3 = 1
           else
              ! move the ifst-th diagonal element backward up the diagonal.
              m1 = -1
              m2 = 0
              m3 = -1
           end if
           do k = ifst + m1,ilst + m2,m3
              ! interchange the k-th and (k+1)-th diagonal elements.
              t11 = t(k,k)
              t22 = t(k + 1,k + 1)
              ! determine the transformation to perform the interchange.
              call la_ylartg(t(k,k + 1),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (k + 2 <= n) call la_yrot(n - k - 1,t(k,k + 2),ldt,t(k + 1,k + 2),ldt,cs,sn)

              call la_yrot(k - 1,t(1,k),1,t(1,k + 1),1,cs,conjg(sn))
              t(k,k) = t22
              t(k + 1,k + 1) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_yrot(n,q(1,k),1,q(1,k + 1),1,cs,conjg(sn))
              end if
           end do
           return
     end subroutine la_ytrexc
#endif
#ifdef LA_WITH_QP
     !> WTREXC: reorders the Schur factorization of a complex matrix
     !> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
     !> is moved to row ILST.
     !> The Schur form T is reordered by a unitary similarity transformation
     !> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
     !> postmultplying it with Z.

     pure subroutine la_wtrexc(compq,n,t,ldt,q,ldq,ifst,ilst,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq
           integer(ilp),intent(in) :: ifst,ilst,ldq,ldt,n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           complex(qp),intent(inout) :: q(ldq,*),t(ldt,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: wantq
           integer(ilp) :: k,m1,m2,m3
           real(qp) :: cs
           complex(qp) :: sn,t11,t22,temp
           ! Intrinsic Functions
           intrinsic :: conjg,max
           ! Executable Statements
           ! decode and test the input parameters.
           info = 0
           wantq = la_lsame(compq,'V')
           if (.not. la_lsame(compq,'N') .and. .not. wantq) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldt < max(1,n)) then
              info = -4
           else if (ldq < 1 .or. (wantq .and. ldq < max(1,n))) then
              info = -6
           else if ((ifst < 1 .or. ifst > n) .and. (n > 0)) then
              info = -7
           else if ((ilst < 1 .or. ilst > n) .and. (n > 0)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WTREXC',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1 .or. ifst == ilst) return
           if (ifst < ilst) then
              ! move the ifst-th diagonal element forward down the diagonal.
              m1 = 0
              m2 = -1
              m3 = 1
           else
              ! move the ifst-th diagonal element backward up the diagonal.
              m1 = -1
              m2 = 0
              m3 = -1
           end if
           do k = ifst + m1,ilst + m2,m3
              ! interchange the k-th and (k+1)-th diagonal elements.
              t11 = t(k,k)
              t22 = t(k + 1,k + 1)
              ! determine the transformation to perform the interchange.
              call la_wlartg(t(k,k + 1),t22 - t11,cs,sn,temp)
              ! apply transformation to the matrix t.
              if (k + 2 <= n) call la_wrot(n - k - 1,t(k,k + 2),ldt,t(k + 1,k + 2),ldt,cs,sn)

              call la_wrot(k - 1,t(1,k),1,t(1,k + 1),1,cs,conjg(sn))
              t(k,k) = t22
              t(k + 1,k + 1) = t11
              if (wantq) then
                 ! accumulate transformation in the matrix q.
                 call la_wrot(n,q(1,k),1,q(1,k + 1),1,cs,conjg(sn))
              end if
           end do
           return
     end subroutine la_wtrexc
#endif

end module la_lapack_eigv_gen_aux
