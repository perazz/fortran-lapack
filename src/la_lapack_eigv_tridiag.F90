!> Symmetric tridiagonal eigenvalues: divide and conquer, rank-one updates, implicit QL and QR
module la_lapack_eigv_tridiag
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level3_gen
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_l2
     use la_lapack_blas_like_l3
     use la_lapack_blas_like_mnorm
     use la_lapack_blas_like_scalar
     use la_lapack_eigv_sym_comp
     use la_lapack_eigv_tridiag2
     use la_lapack_givens_jacobi_rot
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slaed5
     public :: la_slaed6
     public :: la_ssteqr
     public :: la_slaed4
     public :: la_slaed8
     public :: la_slaed9
     public :: la_slaed3
     public :: la_slaed7
     public :: la_slaed2
     public :: la_slaed1
     public :: la_slaed0
     public :: la_dlaed5
     public :: la_dlaed6
     public :: la_dsteqr
     public :: la_dlaed4
     public :: la_dlaed8
     public :: la_dlaed9
     public :: la_dlaed3
     public :: la_dlaed7
     public :: la_dlaed2
     public :: la_dlaed1
     public :: la_dlaed0
#ifdef LA_WITH_XDP
     public :: la_xlaed5
     public :: la_xlaed6
     public :: la_xsteqr
     public :: la_xlaed4
     public :: la_xlaed8
     public :: la_xlaed9
     public :: la_xlaed3
     public :: la_xlaed7
     public :: la_xlaed2
     public :: la_xlaed1
     public :: la_xlaed0
#endif
#ifdef LA_WITH_QP
     public :: la_qlaed5
     public :: la_qlaed6
     public :: la_qsteqr
     public :: la_qlaed4
     public :: la_qlaed8
     public :: la_qlaed9
     public :: la_qlaed3
     public :: la_qlaed7
     public :: la_qlaed2
     public :: la_qlaed1
     public :: la_qlaed0
#endif
     public :: la_claed8
     public :: la_csteqr
     public :: la_claed7
     public :: la_claed0
     public :: la_zlaed8
     public :: la_zsteqr
     public :: la_zlaed7
     public :: la_zlaed0
#ifdef LA_WITH_XDP
     public :: la_ylaed8
     public :: la_ysteqr
     public :: la_ylaed7
     public :: la_ylaed0
#endif
#ifdef LA_WITH_QP
     public :: la_wlaed8
     public :: la_wsteqr
     public :: la_wlaed7
     public :: la_wlaed0
#endif

     contains

     !> This subroutine computes the I-th eigenvalue of a symmetric rank-one
     !> modification of a 2-by-2 diagonal matrix
     !> diag( D )  +  RHO * Z * transpose(Z) .
     !> The diagonal elements in the array D are assumed to satisfy
     !> D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_slaed5(i,d,z,delta,rho,dlam)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(sp),intent(out) :: dlam
           real(sp),intent(in) :: rho
           ! Array Arguments
           real(sp),intent(in) :: d(2),z(2)
           real(sp),intent(out) :: delta(2)
        ! =====================================================================

           ! Local Scalars
           real(sp) :: b,c,del,tau,temp,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           if (i == 1) then
              w = one + two*rho*(z(2)*z(2) - z(1)*z(1))/del
              if (w > zero) then
                 b = del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*del
                 ! b > zero, always
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 dlam = d(1) + tau
                 delta(1) = -z(1)/tau
                 delta(2) = z(2)/(del - tau)
              else
                 b = -del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*del
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 dlam = d(2) + tau
                 delta(1) = -z(1)/(del + tau)
                 delta(2) = -z(2)/tau
              end if
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           else
           ! now i=2
              b = -del + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*del
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              dlam = d(2) + tau
              delta(1) = -z(1)/(del + tau)
              delta(2) = -z(2)/tau
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           end if
           return
     end subroutine la_slaed5
     !> This subroutine computes the I-th eigenvalue of a symmetric rank-one
     !> modification of a 2-by-2 diagonal matrix
     !> diag( D )  +  RHO * Z * transpose(Z) .
     !> The diagonal elements in the array D are assumed to satisfy
     !> D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_dlaed5(i,d,z,delta,rho,dlam)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(dp),intent(out) :: dlam
           real(dp),intent(in) :: rho
           ! Array Arguments
           real(dp),intent(in) :: d(2),z(2)
           real(dp),intent(out) :: delta(2)
        ! =====================================================================

           ! Local Scalars
           real(dp) :: b,c,del,tau,temp,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           if (i == 1) then
              w = one + two*rho*(z(2)*z(2) - z(1)*z(1))/del
              if (w > zero) then
                 b = del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*del
                 ! b > zero, always
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 dlam = d(1) + tau
                 delta(1) = -z(1)/tau
                 delta(2) = z(2)/(del - tau)
              else
                 b = -del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*del
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 dlam = d(2) + tau
                 delta(1) = -z(1)/(del + tau)
                 delta(2) = -z(2)/tau
              end if
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           else
           ! now i=2
              b = -del + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*del
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              dlam = d(2) + tau
              delta(1) = -z(1)/(del + tau)
              delta(2) = -z(2)/tau
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           end if
           return
     end subroutine la_dlaed5
#ifdef LA_WITH_XDP
     !> This subroutine computes the I-th eigenvalue of a symmetric rank-one
     !> modification of a 2-by-2 diagonal matrix
     !> diag( D )  +  RHO * Z * transpose(Z) .
     !> The diagonal elements in the array D are assumed to satisfy
     !> D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_xlaed5(i,d,z,delta,rho,dlam)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(xdp),intent(out) :: dlam
           real(xdp),intent(in) :: rho
           ! Array Arguments
           real(xdp),intent(in) :: d(2),z(2)
           real(xdp),intent(out) :: delta(2)
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: b,c,del,tau,temp,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           if (i == 1) then
              w = one + two*rho*(z(2)*z(2) - z(1)*z(1))/del
              if (w > zero) then
                 b = del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*del
                 ! b > zero, always
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 dlam = d(1) + tau
                 delta(1) = -z(1)/tau
                 delta(2) = z(2)/(del - tau)
              else
                 b = -del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*del
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 dlam = d(2) + tau
                 delta(1) = -z(1)/(del + tau)
                 delta(2) = -z(2)/tau
              end if
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           else
           ! now i=2
              b = -del + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*del
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              dlam = d(2) + tau
              delta(1) = -z(1)/(del + tau)
              delta(2) = -z(2)/tau
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           end if
           return
     end subroutine la_xlaed5
#endif
#ifdef LA_WITH_QP
     !> This subroutine computes the I-th eigenvalue of a symmetric rank-one
     !> modification of a 2-by-2 diagonal matrix
     !> diag( D )  +  RHO * Z * transpose(Z) .
     !> The diagonal elements in the array D are assumed to satisfy
     !> D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_qlaed5(i,d,z,delta,rho,dlam)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(qp),intent(out) :: dlam
           real(qp),intent(in) :: rho
           ! Array Arguments
           real(qp),intent(in) :: d(2),z(2)
           real(qp),intent(out) :: delta(2)
        ! =====================================================================

           ! Local Scalars
           real(qp) :: b,c,del,tau,temp,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           if (i == 1) then
              w = one + two*rho*(z(2)*z(2) - z(1)*z(1))/del
              if (w > zero) then
                 b = del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*del
                 ! b > zero, always
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 dlam = d(1) + tau
                 delta(1) = -z(1)/tau
                 delta(2) = z(2)/(del - tau)
              else
                 b = -del + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*del
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 dlam = d(2) + tau
                 delta(1) = -z(1)/(del + tau)
                 delta(2) = -z(2)/tau
              end if
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           else
           ! now i=2
              b = -del + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*del
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              dlam = d(2) + tau
              delta(1) = -z(1)/(del + tau)
              delta(2) = -z(2)/tau
              temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
              delta(1) = delta(1)/temp
              delta(2) = delta(2)/temp
           end if
           return
     end subroutine la_qlaed5
#endif

     !> SLAED6: computes the positive or negative root (closest to the origin)
     !> of
     !> z(1)        z(2)        z(3)
     !> f(x) =   rho + --------- + ---------- + ---------
     !> d(1)-x      d(2)-x      d(3)-x
     !> It is assumed that
     !> if ORGATI = .true. the root is between d(2) and d(3);
     !> otherwise it is between d(1) and d(2)
     !> This routine will be called by SLAED4 when necessary. In most cases,
     !> the root sought is the smallest in magnitude, though it might not be
     !> in some extremely rare situations.

     pure subroutine la_slaed6(kniter,orgati,rho,d,z,finit,tau,info)
        use la_constants_sp,only:zero,one,two,three,four,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: orgati
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kniter
           real(sp),intent(in) :: finit,rho
           real(sp),intent(out) :: tau
           ! Array Arguments
           real(sp),intent(in) :: d(3),z(3)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40

           ! Local Arrays
           real(sp) :: dscale(3),zscale(3)
           ! Local Scalars
           logical(lk) :: scale
           integer(ilp) :: i,iter,niter
           real(sp) :: a,b,base,c,ddf,df,eps,erretm,eta,f,fc,sclfac,sclinv,small1, &
                     small2,sminv1,sminv2,temp,temp1,temp2,temp3,temp4,lbd,ubd
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           if (orgati) then
              lbd = d(2)
              ubd = d(3)
           else
              lbd = d(1)
              ubd = d(2)
           end if
           if (finit < zero) then
              lbd = zero
           else
              ubd = zero
           end if
           niter = 1
           tau = zero
           if (kniter == 2) then
              if (orgati) then
                 temp = (d(3) - d(2))/two
                 c = rho + z(1)/((d(1) - d(2)) - temp)
                 a = c*(d(2) + d(3)) + z(2) + z(3)
                 b = c*d(2)*d(3) + z(2)*d(3) + z(3)*d(2)
              else
                 temp = (d(1) - d(2))/two
                 c = rho + z(3)/((d(3) - d(2)) - temp)
                 a = c*(d(1) + d(2)) + z(1) + z(2)
                 b = c*d(1)*d(2) + z(1)*d(2) + z(2)*d(1)
              end if
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 tau = b/a
              else if (a <= zero) then
                 tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              if (d(1) == tau .or. d(2) == tau .or. d(3) == tau) then
                 tau = zero
              else
                 temp = finit + tau*z(1)/(d(1)*(d(1) - tau)) + tau*z(2)/(d(2)*(d(2) - tau)) &
                            + tau*z(3)/(d(3)*(d(3) - tau))
                 if (temp <= zero) then
                    lbd = tau
                 else
                    ubd = tau
                 end if
                 if (abs(finit) <= abs(temp)) tau = zero
              end if
           end if
           ! get machine parameters for possible scaling to avoid overflow
           ! modified by sven: parameters small1, sminv1, small2,
           ! sminv2, eps are not saved anymore between one call to the
           ! others but recomputed at each call
           eps = la_slamch('EPSILON')
           base = la_slamch('BASE')
           small1 = base**(int(log(la_slamch('SAFMIN'))/log(base)/three,KIND=ilp))

           sminv1 = one/small1
           small2 = small1*small1
           sminv2 = sminv1*sminv1
           ! determine if scaling of inputs necessary to avoid overflow
           ! when computing 1/temp**3
           if (orgati) then
              temp = min(abs(d(2) - tau),abs(d(3) - tau))
           else
              temp = min(abs(d(1) - tau),abs(d(2) - tau))
           end if
           scale = .false.
           if (temp <= small1) then
              scale = .true.
              if (temp <= small2) then
              ! scale up by power of radix nearest 1/safmin**(2/3)
                 sclfac = sminv2
                 sclinv = small2
              else
              ! scale up by power of radix nearest 1/safmin**(1/3)
                 sclfac = sminv1
                 sclinv = small1
              end if
              ! scaling up safe because d, z, tau scaled elsewhere to be o(1)
              do i = 1,3
                 dscale(i) = d(i)*sclfac
                 zscale(i) = z(i)*sclfac
              end do
              tau = tau*sclfac
              lbd = lbd*sclfac
              ubd = ubd*sclfac
           else
              ! copy d and z to dscale and zscale
              do i = 1,3
                 dscale(i) = d(i)
                 zscale(i) = z(i)
              end do
           end if
           fc = zero
           df = zero
           ddf = zero
           do i = 1,3
              temp = one/(dscale(i) - tau)
              temp1 = zscale(i)*temp
              temp2 = temp1*temp
              temp3 = temp2*temp
              fc = fc + temp1/dscale(i)
              df = df + temp2
              ddf = ddf + temp3
           end do
           f = finit + tau*fc
           if (abs(f) <= zero) go to 60
           if (f <= zero) then
              lbd = tau
           else
              ubd = tau
           end if
              ! iteration begins -- use gragg-thornton-warner cubic convergent
                                  ! scheme
           ! it is not hard to see that
                 ! 1) iterations will go up monotonically
                    ! if finit < 0;
                 ! 2) iterations will go down monotonically
                    ! if finit > 0.
           iter = niter + 1
           loop_50: do niter = iter,maxit
              if (orgati) then
                 temp1 = dscale(2) - tau
                 temp2 = dscale(3) - tau
              else
                 temp1 = dscale(1) - tau
                 temp2 = dscale(2) - tau
              end if
              a = (temp1 + temp2)*f - temp1*temp2*df
              b = temp1*temp2*f
              c = f - (temp1 + temp2)*df + temp1*temp2*ddf
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 eta = b/a
              else if (a <= zero) then
                 eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (f*eta >= zero) then
                 eta = -f/df
              end if
              tau = tau + eta
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              fc = zero
              erretm = zero
              df = zero
              ddf = zero
              do i = 1,3
                 if ((dscale(i) - tau) /= zero) then
                    temp = one/(dscale(i) - tau)
                    temp1 = zscale(i)*temp
                    temp2 = temp1*temp
                    temp3 = temp2*temp
                    temp4 = temp1/dscale(i)
                    fc = fc + temp4
                    erretm = erretm + abs(temp4)
                    df = df + temp2
                    ddf = ddf + temp3
                 else
                    go to 60
                 end if
              end do
              f = finit + tau*fc
              erretm = eight*(abs(finit) + abs(tau)*erretm) + abs(tau)*df
              if ((abs(f) <= four*eps*erretm) .or. ((ubd - lbd) <= four*eps*abs(tau))) go to &
                        60
              if (f <= zero) then
                 lbd = tau
              else
                 ubd = tau
              end if
           end do loop_50
           info = 1
           60 continue
           ! undo scaling
           if (scale) tau = tau*sclinv
           return
     end subroutine la_slaed6
     !> DLAED6: computes the positive or negative root (closest to the origin)
     !> of
     !> z(1)        z(2)        z(3)
     !> f(x) =   rho + --------- + ---------- + ---------
     !> d(1)-x      d(2)-x      d(3)-x
     !> It is assumed that
     !> if ORGATI = .true. the root is between d(2) and d(3);
     !> otherwise it is between d(1) and d(2)
     !> This routine will be called by DLAED4 when necessary. In most cases,
     !> the root sought is the smallest in magnitude, though it might not be
     !> in some extremely rare situations.

     pure subroutine la_dlaed6(kniter,orgati,rho,d,z,finit,tau,info)
        use la_constants_dp,only:zero,one,two,three,four,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: orgati
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kniter
           real(dp),intent(in) :: finit,rho
           real(dp),intent(out) :: tau
           ! Array Arguments
           real(dp),intent(in) :: d(3),z(3)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40

           ! Local Arrays
           real(dp) :: dscale(3),zscale(3)
           ! Local Scalars
           logical(lk) :: scale
           integer(ilp) :: i,iter,niter
           real(dp) :: a,b,base,c,ddf,df,eps,erretm,eta,f,fc,sclfac,sclinv,small1, &
                     small2,sminv1,sminv2,temp,temp1,temp2,temp3,temp4,lbd,ubd
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           if (orgati) then
              lbd = d(2)
              ubd = d(3)
           else
              lbd = d(1)
              ubd = d(2)
           end if
           if (finit < zero) then
              lbd = zero
           else
              ubd = zero
           end if
           niter = 1
           tau = zero
           if (kniter == 2) then
              if (orgati) then
                 temp = (d(3) - d(2))/two
                 c = rho + z(1)/((d(1) - d(2)) - temp)
                 a = c*(d(2) + d(3)) + z(2) + z(3)
                 b = c*d(2)*d(3) + z(2)*d(3) + z(3)*d(2)
              else
                 temp = (d(1) - d(2))/two
                 c = rho + z(3)/((d(3) - d(2)) - temp)
                 a = c*(d(1) + d(2)) + z(1) + z(2)
                 b = c*d(1)*d(2) + z(1)*d(2) + z(2)*d(1)
              end if
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 tau = b/a
              else if (a <= zero) then
                 tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              if (d(1) == tau .or. d(2) == tau .or. d(3) == tau) then
                 tau = zero
              else
                 temp = finit + tau*z(1)/(d(1)*(d(1) - tau)) + tau*z(2)/(d(2)*(d(2) - tau)) &
                            + tau*z(3)/(d(3)*(d(3) - tau))
                 if (temp <= zero) then
                    lbd = tau
                 else
                    ubd = tau
                 end if
                 if (abs(finit) <= abs(temp)) tau = zero
              end if
           end if
           ! get machine parameters for possible scaling to avoid overflow
           ! modified by sven: parameters small1, sminv1, small2,
           ! sminv2, eps are not saved anymore between one call to the
           ! others but recomputed at each call
           eps = la_dlamch('EPSILON')
           base = la_dlamch('BASE')
           small1 = base**(int(log(la_dlamch('SAFMIN'))/log(base)/three,KIND=ilp))

           sminv1 = one/small1
           small2 = small1*small1
           sminv2 = sminv1*sminv1
           ! determine if scaling of inputs necessary to avoid overflow
           ! when computing 1/temp**3
           if (orgati) then
              temp = min(abs(d(2) - tau),abs(d(3) - tau))
           else
              temp = min(abs(d(1) - tau),abs(d(2) - tau))
           end if
           scale = .false.
           if (temp <= small1) then
              scale = .true.
              if (temp <= small2) then
              ! scale up by power of radix nearest 1/safmin**(2/3)
                 sclfac = sminv2
                 sclinv = small2
              else
              ! scale up by power of radix nearest 1/safmin**(1/3)
                 sclfac = sminv1
                 sclinv = small1
              end if
              ! scaling up safe because d, z, tau scaled elsewhere to be o(1)
              do i = 1,3
                 dscale(i) = d(i)*sclfac
                 zscale(i) = z(i)*sclfac
              end do
              tau = tau*sclfac
              lbd = lbd*sclfac
              ubd = ubd*sclfac
           else
              ! copy d and z to dscale and zscale
              do i = 1,3
                 dscale(i) = d(i)
                 zscale(i) = z(i)
              end do
           end if
           fc = zero
           df = zero
           ddf = zero
           do i = 1,3
              temp = one/(dscale(i) - tau)
              temp1 = zscale(i)*temp
              temp2 = temp1*temp
              temp3 = temp2*temp
              fc = fc + temp1/dscale(i)
              df = df + temp2
              ddf = ddf + temp3
           end do
           f = finit + tau*fc
           if (abs(f) <= zero) go to 60
           if (f <= zero) then
              lbd = tau
           else
              ubd = tau
           end if
              ! iteration begins -- use gragg-thornton-warner cubic convergent
                                  ! scheme
           ! it is not hard to see that
                 ! 1) iterations will go up monotonically
                    ! if finit < 0;
                 ! 2) iterations will go down monotonically
                    ! if finit > 0.
           iter = niter + 1
           loop_50: do niter = iter,maxit
              if (orgati) then
                 temp1 = dscale(2) - tau
                 temp2 = dscale(3) - tau
              else
                 temp1 = dscale(1) - tau
                 temp2 = dscale(2) - tau
              end if
              a = (temp1 + temp2)*f - temp1*temp2*df
              b = temp1*temp2*f
              c = f - (temp1 + temp2)*df + temp1*temp2*ddf
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 eta = b/a
              else if (a <= zero) then
                 eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (f*eta >= zero) then
                 eta = -f/df
              end if
              tau = tau + eta
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              fc = zero
              erretm = zero
              df = zero
              ddf = zero
              do i = 1,3
                 if ((dscale(i) - tau) /= zero) then
                    temp = one/(dscale(i) - tau)
                    temp1 = zscale(i)*temp
                    temp2 = temp1*temp
                    temp3 = temp2*temp
                    temp4 = temp1/dscale(i)
                    fc = fc + temp4
                    erretm = erretm + abs(temp4)
                    df = df + temp2
                    ddf = ddf + temp3
                 else
                    go to 60
                 end if
              end do
              f = finit + tau*fc
              erretm = eight*(abs(finit) + abs(tau)*erretm) + abs(tau)*df
              if ((abs(f) <= four*eps*erretm) .or. ((ubd - lbd) <= four*eps*abs(tau))) go to &
                        60
              if (f <= zero) then
                 lbd = tau
              else
                 ubd = tau
              end if
           end do loop_50
           info = 1
           60 continue
           ! undo scaling
           if (scale) tau = tau*sclinv
           return
     end subroutine la_dlaed6
#ifdef LA_WITH_XDP
     !> XLAED6: computes the positive or negative root (closest to the origin)
     !> of
     !> z(1)        z(2)        z(3)
     !> f(x) =   rho + --------- + ---------- + ---------
     !> d(1)-x      d(2)-x      d(3)-x
     !> It is assumed that
     !> if ORGATI = .true. the root is between d(2) and d(3);
     !> otherwise it is between d(1) and d(2)
     !> This routine will be called by XLAED4 when necessary. In most cases,
     !> the root sought is the smallest in magnitude, though it might not be
     !> in some extremely rare situations.

     pure subroutine la_xlaed6(kniter,orgati,rho,d,z,finit,tau,info)
        use la_constants_xdp,only:zero,one,two,three,four,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: orgati
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kniter
           real(xdp),intent(in) :: finit,rho
           real(xdp),intent(out) :: tau
           ! Array Arguments
           real(xdp),intent(in) :: d(3),z(3)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40

           ! Local Arrays
           real(xdp) :: dscale(3),zscale(3)
           ! Local Scalars
           logical(lk) :: scale
           integer(ilp) :: i,iter,niter
           real(xdp) :: a,b,base,c,ddf,df,eps,erretm,eta,f,fc,sclfac,sclinv,small1, &
                     small2,sminv1,sminv2,temp,temp1,temp2,temp3,temp4,lbd,ubd
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           if (orgati) then
              lbd = d(2)
              ubd = d(3)
           else
              lbd = d(1)
              ubd = d(2)
           end if
           if (finit < zero) then
              lbd = zero
           else
              ubd = zero
           end if
           niter = 1
           tau = zero
           if (kniter == 2) then
              if (orgati) then
                 temp = (d(3) - d(2))/two
                 c = rho + z(1)/((d(1) - d(2)) - temp)
                 a = c*(d(2) + d(3)) + z(2) + z(3)
                 b = c*d(2)*d(3) + z(2)*d(3) + z(3)*d(2)
              else
                 temp = (d(1) - d(2))/two
                 c = rho + z(3)/((d(3) - d(2)) - temp)
                 a = c*(d(1) + d(2)) + z(1) + z(2)
                 b = c*d(1)*d(2) + z(1)*d(2) + z(2)*d(1)
              end if
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 tau = b/a
              else if (a <= zero) then
                 tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              if (d(1) == tau .or. d(2) == tau .or. d(3) == tau) then
                 tau = zero
              else
                 temp = finit + tau*z(1)/(d(1)*(d(1) - tau)) + tau*z(2)/(d(2)*(d(2) - tau)) &
                            + tau*z(3)/(d(3)*(d(3) - tau))
                 if (temp <= zero) then
                    lbd = tau
                 else
                    ubd = tau
                 end if
                 if (abs(finit) <= abs(temp)) tau = zero
              end if
           end if
           ! get machine parameters for possible scaling to avoid overflow
           ! modified by sven: parameters small1, sminv1, small2,
           ! sminv2, eps are not saved anymore between one call to the
           ! others but recomputed at each call
           eps = la_xlamch('EPSILON')
           base = la_xlamch('BASE')
           small1 = base**(int(log(la_xlamch('SAFMIN'))/log(base)/three,KIND=ilp))

           sminv1 = one/small1
           small2 = small1*small1
           sminv2 = sminv1*sminv1
           ! determine if scaling of inputs necessary to avoid overflow
           ! when computing 1/temp**3
           if (orgati) then
              temp = min(abs(d(2) - tau),abs(d(3) - tau))
           else
              temp = min(abs(d(1) - tau),abs(d(2) - tau))
           end if
           scale = .false.
           if (temp <= small1) then
              scale = .true.
              if (temp <= small2) then
              ! scale up by power of radix nearest 1/safmin**(2/3)
                 sclfac = sminv2
                 sclinv = small2
              else
              ! scale up by power of radix nearest 1/safmin**(1/3)
                 sclfac = sminv1
                 sclinv = small1
              end if
              ! scaling up safe because d, z, tau scaled elsewhere to be o(1)
              do i = 1,3
                 dscale(i) = d(i)*sclfac
                 zscale(i) = z(i)*sclfac
              end do
              tau = tau*sclfac
              lbd = lbd*sclfac
              ubd = ubd*sclfac
           else
              ! copy d and z to dscale and zscale
              do i = 1,3
                 dscale(i) = d(i)
                 zscale(i) = z(i)
              end do
           end if
           fc = zero
           df = zero
           ddf = zero
           do i = 1,3
              temp = one/(dscale(i) - tau)
              temp1 = zscale(i)*temp
              temp2 = temp1*temp
              temp3 = temp2*temp
              fc = fc + temp1/dscale(i)
              df = df + temp2
              ddf = ddf + temp3
           end do
           f = finit + tau*fc
           if (abs(f) <= zero) go to 60
           if (f <= zero) then
              lbd = tau
           else
              ubd = tau
           end if
              ! iteration begins -- use gragg-thornton-warner cubic convergent
                                  ! scheme
           ! it is not hard to see that
                 ! 1) iterations will go up monotonically
                    ! if finit < 0;
                 ! 2) iterations will go down monotonically
                    ! if finit > 0.
           iter = niter + 1
           loop_50: do niter = iter,maxit
              if (orgati) then
                 temp1 = dscale(2) - tau
                 temp2 = dscale(3) - tau
              else
                 temp1 = dscale(1) - tau
                 temp2 = dscale(2) - tau
              end if
              a = (temp1 + temp2)*f - temp1*temp2*df
              b = temp1*temp2*f
              c = f - (temp1 + temp2)*df + temp1*temp2*ddf
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 eta = b/a
              else if (a <= zero) then
                 eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (f*eta >= zero) then
                 eta = -f/df
              end if
              tau = tau + eta
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              fc = zero
              erretm = zero
              df = zero
              ddf = zero
              do i = 1,3
                 if ((dscale(i) - tau) /= zero) then
                    temp = one/(dscale(i) - tau)
                    temp1 = zscale(i)*temp
                    temp2 = temp1*temp
                    temp3 = temp2*temp
                    temp4 = temp1/dscale(i)
                    fc = fc + temp4
                    erretm = erretm + abs(temp4)
                    df = df + temp2
                    ddf = ddf + temp3
                 else
                    go to 60
                 end if
              end do
              f = finit + tau*fc
              erretm = eight*(abs(finit) + abs(tau)*erretm) + abs(tau)*df
              if ((abs(f) <= four*eps*erretm) .or. ((ubd - lbd) <= four*eps*abs(tau))) go to &
                        60
              if (f <= zero) then
                 lbd = tau
              else
                 ubd = tau
              end if
           end do loop_50
           info = 1
           60 continue
           ! undo scaling
           if (scale) tau = tau*sclinv
           return
     end subroutine la_xlaed6
#endif
#ifdef LA_WITH_QP
     !> QLAED6: computes the positive or negative root (closest to the origin)
     !> of
     !> z(1)        z(2)        z(3)
     !> f(x) =   rho + --------- + ---------- + ---------
     !> d(1)-x      d(2)-x      d(3)-x
     !> It is assumed that
     !> if ORGATI = .true. the root is between d(2) and d(3);
     !> otherwise it is between d(1) and d(2)
     !> This routine will be called by QLAED4 when necessary. In most cases,
     !> the root sought is the smallest in magnitude, though it might not be
     !> in some extremely rare situations.

     pure subroutine la_qlaed6(kniter,orgati,rho,d,z,finit,tau,info)
        use la_constants_qp,only:zero,one,two,three,four,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: orgati
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kniter
           real(qp),intent(in) :: finit,rho
           real(qp),intent(out) :: tau
           ! Array Arguments
           real(qp),intent(in) :: d(3),z(3)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40

           ! Local Arrays
           real(qp) :: dscale(3),zscale(3)
           ! Local Scalars
           logical(lk) :: scale
           integer(ilp) :: i,iter,niter
           real(qp) :: a,b,base,c,ddf,df,eps,erretm,eta,f,fc,sclfac,sclinv,small1, &
                     small2,sminv1,sminv2,temp,temp1,temp2,temp3,temp4,lbd,ubd
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           if (orgati) then
              lbd = d(2)
              ubd = d(3)
           else
              lbd = d(1)
              ubd = d(2)
           end if
           if (finit < zero) then
              lbd = zero
           else
              ubd = zero
           end if
           niter = 1
           tau = zero
           if (kniter == 2) then
              if (orgati) then
                 temp = (d(3) - d(2))/two
                 c = rho + z(1)/((d(1) - d(2)) - temp)
                 a = c*(d(2) + d(3)) + z(2) + z(3)
                 b = c*d(2)*d(3) + z(2)*d(3) + z(3)*d(2)
              else
                 temp = (d(1) - d(2))/two
                 c = rho + z(3)/((d(3) - d(2)) - temp)
                 a = c*(d(1) + d(2)) + z(1) + z(2)
                 b = c*d(1)*d(2) + z(1)*d(2) + z(2)*d(1)
              end if
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 tau = b/a
              else if (a <= zero) then
                 tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              if (d(1) == tau .or. d(2) == tau .or. d(3) == tau) then
                 tau = zero
              else
                 temp = finit + tau*z(1)/(d(1)*(d(1) - tau)) + tau*z(2)/(d(2)*(d(2) - tau)) &
                            + tau*z(3)/(d(3)*(d(3) - tau))
                 if (temp <= zero) then
                    lbd = tau
                 else
                    ubd = tau
                 end if
                 if (abs(finit) <= abs(temp)) tau = zero
              end if
           end if
           ! get machine parameters for possible scaling to avoid overflow
           ! modified by sven: parameters small1, sminv1, small2,
           ! sminv2, eps are not saved anymore between one call to the
           ! others but recomputed at each call
           eps = la_qlamch('EPSILON')
           base = la_qlamch('BASE')
           small1 = base**(int(log(la_qlamch('SAFMIN'))/log(base)/three,KIND=ilp))

           sminv1 = one/small1
           small2 = small1*small1
           sminv2 = sminv1*sminv1
           ! determine if scaling of inputs necessary to avoid overflow
           ! when computing 1/temp**3
           if (orgati) then
              temp = min(abs(d(2) - tau),abs(d(3) - tau))
           else
              temp = min(abs(d(1) - tau),abs(d(2) - tau))
           end if
           scale = .false.
           if (temp <= small1) then
              scale = .true.
              if (temp <= small2) then
              ! scale up by power of radix nearest 1/safmin**(2/3)
                 sclfac = sminv2
                 sclinv = small2
              else
              ! scale up by power of radix nearest 1/safmin**(1/3)
                 sclfac = sminv1
                 sclinv = small1
              end if
              ! scaling up safe because d, z, tau scaled elsewhere to be o(1)
              do i = 1,3
                 dscale(i) = d(i)*sclfac
                 zscale(i) = z(i)*sclfac
              end do
              tau = tau*sclfac
              lbd = lbd*sclfac
              ubd = ubd*sclfac
           else
              ! copy d and z to dscale and zscale
              do i = 1,3
                 dscale(i) = d(i)
                 zscale(i) = z(i)
              end do
           end if
           fc = zero
           df = zero
           ddf = zero
           do i = 1,3
              temp = one/(dscale(i) - tau)
              temp1 = zscale(i)*temp
              temp2 = temp1*temp
              temp3 = temp2*temp
              fc = fc + temp1/dscale(i)
              df = df + temp2
              ddf = ddf + temp3
           end do
           f = finit + tau*fc
           if (abs(f) <= zero) go to 60
           if (f <= zero) then
              lbd = tau
           else
              ubd = tau
           end if
              ! iteration begins -- use gragg-thornton-warner cubic convergent
                                  ! scheme
           ! it is not hard to see that
                 ! 1) iterations will go up monotonically
                    ! if finit < 0;
                 ! 2) iterations will go down monotonically
                    ! if finit > 0.
           iter = niter + 1
           loop_50: do niter = iter,maxit
              if (orgati) then
                 temp1 = dscale(2) - tau
                 temp2 = dscale(3) - tau
              else
                 temp1 = dscale(1) - tau
                 temp2 = dscale(2) - tau
              end if
              a = (temp1 + temp2)*f - temp1*temp2*df
              b = temp1*temp2*f
              c = f - (temp1 + temp2)*df + temp1*temp2*ddf
              temp = max(abs(a),abs(b),abs(c))
              a = a/temp
              b = b/temp
              c = c/temp
              if (c == zero) then
                 eta = b/a
              else if (a <= zero) then
                 eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
              end if
              if (f*eta >= zero) then
                 eta = -f/df
              end if
              tau = tau + eta
              if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/two
              fc = zero
              erretm = zero
              df = zero
              ddf = zero
              do i = 1,3
                 if ((dscale(i) - tau) /= zero) then
                    temp = one/(dscale(i) - tau)
                    temp1 = zscale(i)*temp
                    temp2 = temp1*temp
                    temp3 = temp2*temp
                    temp4 = temp1/dscale(i)
                    fc = fc + temp4
                    erretm = erretm + abs(temp4)
                    df = df + temp2
                    ddf = ddf + temp3
                 else
                    go to 60
                 end if
              end do
              f = finit + tau*fc
              erretm = eight*(abs(finit) + abs(tau)*erretm) + abs(tau)*df
              if ((abs(f) <= four*eps*erretm) .or. ((ubd - lbd) <= four*eps*abs(tau))) go to &
                        60
              if (f <= zero) then
                 lbd = tau
              else
                 ubd = tau
              end if
           end do loop_50
           info = 1
           60 continue
           ! undo scaling
           if (scale) tau = tau*sclinv
           return
     end subroutine la_qlaed6
#endif

     !> SSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band symmetric matrix can also be found
     !> if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to
     !> tridiagonal form.

     pure subroutine la_ssteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_sp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(sp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(sp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = one
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_slamch('E')
           eps2 = eps**2
           safmin = la_slamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_slaset('FULL',n,n,zero,one,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_slanst('M',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_slascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_slascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_slascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_slascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_slae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_slaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_slasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_slae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_slapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_slartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_slasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_slae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_slaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_slasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_slae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_slapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_slartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_slasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_slascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_slascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_slascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_slascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 190
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_slasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_sswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           190 continue
           return
     end subroutine la_ssteqr
     !> DSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band symmetric matrix can also be found
     !> if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
     !> tridiagonal form.

     pure subroutine la_dsteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_dp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(dp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(dp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = one
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_dlamch('E')
           eps2 = eps**2
           safmin = la_dlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_dlaset('FULL',n,n,zero,one,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_dlanst('M',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_dlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_dlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_dlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_dlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_dlae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_dlaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_dlasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_dlae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_dlapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_dlartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_dlasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_dlae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_dlaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_dlasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_dlae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_dlapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_dlartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_dlasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_dlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_dlascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_dlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_dlascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 190
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_dlasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_dswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           190 continue
           return
     end subroutine la_dsteqr
#ifdef LA_WITH_XDP
     !> XSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band symmetric matrix can also be found
     !> if XSYTRD or XSPTRD or XSBTRD has been used to reduce this matrix to
     !> tridiagonal form.

     pure subroutine la_xsteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_xdp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(xdp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(xdp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('XSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = one
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_xlamch('E')
           eps2 = eps**2
           safmin = la_xlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_xlaset('FULL',n,n,zero,one,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_xlanst('M',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_xlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_xlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_xlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_xlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_xlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_xlaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_xlasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_xlae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_xlapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_xlartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_xlasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_xlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_xlaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_xlasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_xlae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_xlapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_xlartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_xlasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_xlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_xlascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_xlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_xlascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 190
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_xlasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_xswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           190 continue
           return
     end subroutine la_xsteqr
#endif
#ifdef LA_WITH_QP
     !> QSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band symmetric matrix can also be found
     !> if QSYTRD or QSPTRD or QSBTRD has been used to reduce this matrix to
     !> tridiagonal form.

     pure subroutine la_qsteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_qp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(qp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(qp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = one
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_qlamch('E')
           eps2 = eps**2
           safmin = la_qlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_qlaset('FULL',n,n,zero,one,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_qlanst('M',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_qlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_qlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_qlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_qlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_qlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_qlaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_qlasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_qlae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_qlapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_qlartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_qlasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_qlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_qlaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_qlasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_qlae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_qlapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_qlartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_qlasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_qlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_qlascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_qlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_qlascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 190
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_qlasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_qswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           190 continue
           return
     end subroutine la_qsteqr
#endif

     !> This subroutine computes the I-th updated eigenvalue of a symmetric
     !> rank-one modification to a diagonal matrix whose elements are
     !> given in the array d, and that
     !> D(i) < D(j)  for  i < j
     !> and that RHO > 0.  This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D )  +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_slaed4(n,i,d,z,delta,rho,dlam,info)
        use la_constants_sp,only:zero,one,two,three,four,eight,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(sp),intent(out) :: dlam
           real(sp),intent(in) :: rho
           ! Array Arguments
           real(sp),intent(in) :: d(*),z(*)
           real(sp),intent(out) :: delta(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(sp) :: a,b,c,del,dltlb,dltub,dphi,dpsi,dw,eps,erretm,eta,midpt,phi, &
                     prew,psi,rhoinv,tau,temp,temp1,w
           ! Local Arrays
           real(sp) :: zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
               ! presumably, i=1 upon entry
              dlam = d(1) + rho*z(1)*z(1)
              delta(1) = one
              return
           end if
           if (n == 2) then
              call la_slaed5(i,d,z,delta,rho,dlam)
              return
           end if
           ! compute machine epsilon
           eps = la_slamch('EPSILON')
           rhoinv = one/rho
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              midpt = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/delta(ii) + z(n)*z(n)/delta(n)
              if (w <= zero) then
                 temp = z(n - 1)*z(n - 1)/(d(n) - d(n - 1) + rho) + z(n)*z(n)/rho
                 if (c <= temp) then
                    tau = rho
                 else
                    del = d(n) - d(n - 1)
                    a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*del
                    if (a < zero) then
                       tau = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                 end if
                 ! it can be proved that
                     ! d(n)+rho/2 <= lambda(n) < d(n)+tau <= d(n)+rho
                 dltlb = midpt
                 dltub = rho
              else
                 del = d(n) - d(n - 1)
                 a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*del
                 if (a < zero) then
                    tau = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 ! it can be proved that
                     ! d(n) < d(n)+tau < lambda(n) < d(n)+rho/2
                 dltlb = zero
                 dltub = midpt
              end if
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 dlam = d(i) + tau
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              c = w - delta(n - 1)*dpsi - delta(n)*dphi
              a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
              b = delta(n - 1)*delta(n)*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                ! eta = b/a
                 ! eta = rho - tau
                 eta = dltub - tau
              else if (a >= zero) then
                 eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta > zero) eta = -w/(dpsi + dphi)
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              tau = tau + eta
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    dlam = d(i) + tau
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 c = w - delta(n - 1)*dpsi - delta(n)*dphi
                 a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
                 b = delta(n - 1)*delta(n)*w
                 if (a >= zero) then
                    eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta > zero) eta = -w/(dpsi + dphi)
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 temp = z(n)/delta(n)
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              dlam = d(i) + tau
              go to 250
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              del = d(ip1) - d(i)
              midpt = del/two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/delta(i) + z(ip1)*z(ip1)/delta(ip1)
              if (w > zero) then
                 ! d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 a = c*del + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*del
                 if (a > zero) then
                    tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 dltlb = zero
                 dltub = midpt
              else
                 ! (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 a = c*del - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*del
                 if (a < zero) then
                    tau = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 dltlb = -midpt
                 dltub = zero
              end if
              if (orgati) then
                 do j = 1,n
                    delta(j) = (d(j) - d(i)) - tau
                 end do
              else
                 do j = 1,n
                    delta(j) = (d(j) - d(ip1)) - tau
                 end do
              end if
              if (orgati) then
                 ii = i
              else
                 ii = i + 1
              end if
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              w = rhoinv + phi + psi
              ! w is the value of the secular function with
              ! its ii-th element removed.
              swtch3 = .false.
              if (orgati) then
                 if (w < zero) swtch3 = .true.
              else
                 if (w > zero) swtch3 = .true.
              end if
              if (ii == 1 .or. ii == n) swtch3 = .false.
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                        *dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 if (orgati) then
                    dlam = d(i) + tau
                 else
                    dlam = d(ip1) + tau
                 end if
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 if (orgati) then
                    c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                              **2
                 else
                    c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                              **2
                 end if
                 a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                 b = delta(i)*delta(ip1)*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                       end if
                    end if
                    eta = b/a
                 else if (a <= zero) then
                    eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 end if
              else
                 ! interpolation using three most relevant poles
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/delta(iim1)
                    temp1 = temp1*temp1
                    c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                 else
                    temp1 = z(iip1)/delta(iip1)
                    temp1 = temp1*temp1
                    c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1))*temp1
                    zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 call la_slaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                 if (info /= 0) go to 250
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              prew = w
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau + eta) &
                        *dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              tau = tau + eta
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_240: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    if (orgati) then
                       dlam = d(i) + tau
                    else
                       dlam = d(ip1) + tau
                    end if
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                                    **2
                       else
                          c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                                    **2
                       end if
                    else
                       temp = z(ii)/delta(ii)
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - delta(i)*dpsi - delta(ip1)*dphi
                    end if
                    a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                    b = delta(i)*delta(ip1)*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)

                             else
                                a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                             end if
                          else
                             a = delta(i)*delta(i)*dpsi + delta(ip1)*delta(ip1) &
                                       *dphi
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 else
                    ! interpolation using three most relevant poles
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - delta(iim1)*dpsi - delta(iip1)*dphi
                       zz(1) = delta(iim1)*delta(iim1)*dpsi
                       zz(3) = delta(iip1)*delta(iip1)*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/delta(iim1)
                          temp1 = temp1*temp1
                          c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1)) &
                                    *temp1
                          zz(1) = z(iim1)*z(iim1)
                          zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                       else
                          temp1 = z(iip1)/delta(iip1)
                          temp1 = temp1*temp1
                          c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1)) &
                                    *temp1
                          zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    call la_slaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                    if (info /= 0) go to 250
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 prew = w
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/delta(j)
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 temp = z(ii)/delta(ii)
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                           *dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_240
              ! return with info = 1, niter = maxit and not converged
              info = 1
              if (orgati) then
                 dlam = d(i) + tau
              else
                 dlam = d(ip1) + tau
              end if
           end if
           250 continue
           return
     end subroutine la_slaed4
     !> This subroutine computes the I-th updated eigenvalue of a symmetric
     !> rank-one modification to a diagonal matrix whose elements are
     !> given in the array d, and that
     !> D(i) < D(j)  for  i < j
     !> and that RHO > 0.  This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D )  +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_dlaed4(n,i,d,z,delta,rho,dlam,info)
        use la_constants_dp,only:zero,one,two,three,four,eight,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(dp),intent(out) :: dlam
           real(dp),intent(in) :: rho
           ! Array Arguments
           real(dp),intent(in) :: d(*),z(*)
           real(dp),intent(out) :: delta(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(dp) :: a,b,c,del,dltlb,dltub,dphi,dpsi,dw,eps,erretm,eta,midpt,phi, &
                     prew,psi,rhoinv,tau,temp,temp1,w
           ! Local Arrays
           real(dp) :: zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
               ! presumably, i=1 upon entry
              dlam = d(1) + rho*z(1)*z(1)
              delta(1) = one
              return
           end if
           if (n == 2) then
              call la_dlaed5(i,d,z,delta,rho,dlam)
              return
           end if
           ! compute machine epsilon
           eps = la_dlamch('EPSILON')
           rhoinv = one/rho
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              midpt = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/delta(ii) + z(n)*z(n)/delta(n)
              if (w <= zero) then
                 temp = z(n - 1)*z(n - 1)/(d(n) - d(n - 1) + rho) + z(n)*z(n)/rho
                 if (c <= temp) then
                    tau = rho
                 else
                    del = d(n) - d(n - 1)
                    a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*del
                    if (a < zero) then
                       tau = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                 end if
                 ! it can be proved that
                     ! d(n)+rho/2 <= lambda(n) < d(n)+tau <= d(n)+rho
                 dltlb = midpt
                 dltub = rho
              else
                 del = d(n) - d(n - 1)
                 a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*del
                 if (a < zero) then
                    tau = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 ! it can be proved that
                     ! d(n) < d(n)+tau < lambda(n) < d(n)+rho/2
                 dltlb = zero
                 dltub = midpt
              end if
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 dlam = d(i) + tau
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              c = w - delta(n - 1)*dpsi - delta(n)*dphi
              a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
              b = delta(n - 1)*delta(n)*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                ! eta = b/a
                 ! eta = rho - tau
                 eta = dltub - tau
              else if (a >= zero) then
                 eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta > zero) eta = -w/(dpsi + dphi)
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              tau = tau + eta
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    dlam = d(i) + tau
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 c = w - delta(n - 1)*dpsi - delta(n)*dphi
                 a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
                 b = delta(n - 1)*delta(n)*w
                 if (a >= zero) then
                    eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta > zero) eta = -w/(dpsi + dphi)
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 temp = z(n)/delta(n)
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              dlam = d(i) + tau
              go to 250
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              del = d(ip1) - d(i)
              midpt = del/two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/delta(i) + z(ip1)*z(ip1)/delta(ip1)
              if (w > zero) then
                 ! d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 a = c*del + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*del
                 if (a > zero) then
                    tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 dltlb = zero
                 dltub = midpt
              else
                 ! (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 a = c*del - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*del
                 if (a < zero) then
                    tau = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 dltlb = -midpt
                 dltub = zero
              end if
              if (orgati) then
                 do j = 1,n
                    delta(j) = (d(j) - d(i)) - tau
                 end do
              else
                 do j = 1,n
                    delta(j) = (d(j) - d(ip1)) - tau
                 end do
              end if
              if (orgati) then
                 ii = i
              else
                 ii = i + 1
              end if
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              w = rhoinv + phi + psi
              ! w is the value of the secular function with
              ! its ii-th element removed.
              swtch3 = .false.
              if (orgati) then
                 if (w < zero) swtch3 = .true.
              else
                 if (w > zero) swtch3 = .true.
              end if
              if (ii == 1 .or. ii == n) swtch3 = .false.
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                        *dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 if (orgati) then
                    dlam = d(i) + tau
                 else
                    dlam = d(ip1) + tau
                 end if
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 if (orgati) then
                    c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                              **2
                 else
                    c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                              **2
                 end if
                 a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                 b = delta(i)*delta(ip1)*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                       end if
                    end if
                    eta = b/a
                 else if (a <= zero) then
                    eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 end if
              else
                 ! interpolation using three most relevant poles
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/delta(iim1)
                    temp1 = temp1*temp1
                    c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                 else
                    temp1 = z(iip1)/delta(iip1)
                    temp1 = temp1*temp1
                    c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1))*temp1
                    zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 call la_dlaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                 if (info /= 0) go to 250
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              prew = w
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau + eta) &
                        *dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              tau = tau + eta
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_240: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    if (orgati) then
                       dlam = d(i) + tau
                    else
                       dlam = d(ip1) + tau
                    end if
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                                    **2
                       else
                          c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                                    **2
                       end if
                    else
                       temp = z(ii)/delta(ii)
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - delta(i)*dpsi - delta(ip1)*dphi
                    end if
                    a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                    b = delta(i)*delta(ip1)*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)

                             else
                                a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                             end if
                          else
                             a = delta(i)*delta(i)*dpsi + delta(ip1)*delta(ip1) &
                                       *dphi
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 else
                    ! interpolation using three most relevant poles
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - delta(iim1)*dpsi - delta(iip1)*dphi
                       zz(1) = delta(iim1)*delta(iim1)*dpsi
                       zz(3) = delta(iip1)*delta(iip1)*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/delta(iim1)
                          temp1 = temp1*temp1
                          c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1)) &
                                    *temp1
                          zz(1) = z(iim1)*z(iim1)
                          zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                       else
                          temp1 = z(iip1)/delta(iip1)
                          temp1 = temp1*temp1
                          c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1)) &
                                    *temp1
                          zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    call la_dlaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                    if (info /= 0) go to 250
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 prew = w
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/delta(j)
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 temp = z(ii)/delta(ii)
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                           *dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_240
              ! return with info = 1, niter = maxit and not converged
              info = 1
              if (orgati) then
                 dlam = d(i) + tau
              else
                 dlam = d(ip1) + tau
              end if
           end if
           250 continue
           return
     end subroutine la_dlaed4
#ifdef LA_WITH_XDP
     !> This subroutine computes the I-th updated eigenvalue of a symmetric
     !> rank-one modification to a diagonal matrix whose elements are
     !> given in the array d, and that
     !> D(i) < D(j)  for  i < j
     !> and that RHO > 0.  This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D )  +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_xlaed4(n,i,d,z,delta,rho,dlam,info)
        use la_constants_xdp,only:zero,one,two,three,four,eight,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(xdp),intent(out) :: dlam
           real(xdp),intent(in) :: rho
           ! Array Arguments
           real(xdp),intent(in) :: d(*),z(*)
           real(xdp),intent(out) :: delta(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(xdp) :: a,b,c,del,dltlb,dltub,dphi,dpsi,dw,eps,erretm,eta,midpt,phi, &
                     prew,psi,rhoinv,tau,temp,temp1,w
           ! Local Arrays
           real(xdp) :: zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
               ! presumably, i=1 upon entry
              dlam = d(1) + rho*z(1)*z(1)
              delta(1) = one
              return
           end if
           if (n == 2) then
              call la_xlaed5(i,d,z,delta,rho,dlam)
              return
           end if
           ! compute machine epsilon
           eps = la_xlamch('EPSILON')
           rhoinv = one/rho
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              midpt = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/delta(ii) + z(n)*z(n)/delta(n)
              if (w <= zero) then
                 temp = z(n - 1)*z(n - 1)/(d(n) - d(n - 1) + rho) + z(n)*z(n)/rho
                 if (c <= temp) then
                    tau = rho
                 else
                    del = d(n) - d(n - 1)
                    a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*del
                    if (a < zero) then
                       tau = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                 end if
                 ! it can be proved that
                     ! d(n)+rho/2 <= lambda(n) < d(n)+tau <= d(n)+rho
                 dltlb = midpt
                 dltub = rho
              else
                 del = d(n) - d(n - 1)
                 a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*del
                 if (a < zero) then
                    tau = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 ! it can be proved that
                     ! d(n) < d(n)+tau < lambda(n) < d(n)+rho/2
                 dltlb = zero
                 dltub = midpt
              end if
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 dlam = d(i) + tau
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              c = w - delta(n - 1)*dpsi - delta(n)*dphi
              a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
              b = delta(n - 1)*delta(n)*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                ! eta = b/a
                 ! eta = rho - tau
                 eta = dltub - tau
              else if (a >= zero) then
                 eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta > zero) eta = -w/(dpsi + dphi)
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              tau = tau + eta
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    dlam = d(i) + tau
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 c = w - delta(n - 1)*dpsi - delta(n)*dphi
                 a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
                 b = delta(n - 1)*delta(n)*w
                 if (a >= zero) then
                    eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta > zero) eta = -w/(dpsi + dphi)
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 temp = z(n)/delta(n)
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              dlam = d(i) + tau
              go to 250
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              del = d(ip1) - d(i)
              midpt = del/two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/delta(i) + z(ip1)*z(ip1)/delta(ip1)
              if (w > zero) then
                 ! d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 a = c*del + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*del
                 if (a > zero) then
                    tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 dltlb = zero
                 dltub = midpt
              else
                 ! (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 a = c*del - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*del
                 if (a < zero) then
                    tau = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 dltlb = -midpt
                 dltub = zero
              end if
              if (orgati) then
                 do j = 1,n
                    delta(j) = (d(j) - d(i)) - tau
                 end do
              else
                 do j = 1,n
                    delta(j) = (d(j) - d(ip1)) - tau
                 end do
              end if
              if (orgati) then
                 ii = i
              else
                 ii = i + 1
              end if
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              w = rhoinv + phi + psi
              ! w is the value of the secular function with
              ! its ii-th element removed.
              swtch3 = .false.
              if (orgati) then
                 if (w < zero) swtch3 = .true.
              else
                 if (w > zero) swtch3 = .true.
              end if
              if (ii == 1 .or. ii == n) swtch3 = .false.
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                        *dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 if (orgati) then
                    dlam = d(i) + tau
                 else
                    dlam = d(ip1) + tau
                 end if
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 if (orgati) then
                    c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                              **2
                 else
                    c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                              **2
                 end if
                 a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                 b = delta(i)*delta(ip1)*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                       end if
                    end if
                    eta = b/a
                 else if (a <= zero) then
                    eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 end if
              else
                 ! interpolation using three most relevant poles
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/delta(iim1)
                    temp1 = temp1*temp1
                    c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                 else
                    temp1 = z(iip1)/delta(iip1)
                    temp1 = temp1*temp1
                    c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1))*temp1
                    zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 call la_xlaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                 if (info /= 0) go to 250
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              prew = w
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau + eta) &
                        *dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              tau = tau + eta
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_240: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    if (orgati) then
                       dlam = d(i) + tau
                    else
                       dlam = d(ip1) + tau
                    end if
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                                    **2
                       else
                          c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                                    **2
                       end if
                    else
                       temp = z(ii)/delta(ii)
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - delta(i)*dpsi - delta(ip1)*dphi
                    end if
                    a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                    b = delta(i)*delta(ip1)*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)

                             else
                                a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                             end if
                          else
                             a = delta(i)*delta(i)*dpsi + delta(ip1)*delta(ip1) &
                                       *dphi
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 else
                    ! interpolation using three most relevant poles
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - delta(iim1)*dpsi - delta(iip1)*dphi
                       zz(1) = delta(iim1)*delta(iim1)*dpsi
                       zz(3) = delta(iip1)*delta(iip1)*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/delta(iim1)
                          temp1 = temp1*temp1
                          c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1)) &
                                    *temp1
                          zz(1) = z(iim1)*z(iim1)
                          zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                       else
                          temp1 = z(iip1)/delta(iip1)
                          temp1 = temp1*temp1
                          c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1)) &
                                    *temp1
                          zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    call la_xlaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                    if (info /= 0) go to 250
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 prew = w
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/delta(j)
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 temp = z(ii)/delta(ii)
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                           *dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_240
              ! return with info = 1, niter = maxit and not converged
              info = 1
              if (orgati) then
                 dlam = d(i) + tau
              else
                 dlam = d(ip1) + tau
              end if
           end if
           250 continue
           return
     end subroutine la_xlaed4
#endif
#ifdef LA_WITH_QP
     !> This subroutine computes the I-th updated eigenvalue of a symmetric
     !> rank-one modification to a diagonal matrix whose elements are
     !> given in the array d, and that
     !> D(i) < D(j)  for  i < j
     !> and that RHO > 0.  This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D )  +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_qlaed4(n,i,d,z,delta,rho,dlam,info)
        use la_constants_qp,only:zero,one,two,three,four,eight,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(qp),intent(out) :: dlam
           real(qp),intent(in) :: rho
           ! Array Arguments
           real(qp),intent(in) :: d(*),z(*)
           real(qp),intent(out) :: delta(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(qp) :: a,b,c,del,dltlb,dltub,dphi,dpsi,dw,eps,erretm,eta,midpt,phi, &
                     prew,psi,rhoinv,tau,temp,temp1,w
           ! Local Arrays
           real(qp) :: zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
               ! presumably, i=1 upon entry
              dlam = d(1) + rho*z(1)*z(1)
              delta(1) = one
              return
           end if
           if (n == 2) then
              call la_qlaed5(i,d,z,delta,rho,dlam)
              return
           end if
           ! compute machine epsilon
           eps = la_qlamch('EPSILON')
           rhoinv = one/rho
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              midpt = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/delta(ii) + z(n)*z(n)/delta(n)
              if (w <= zero) then
                 temp = z(n - 1)*z(n - 1)/(d(n) - d(n - 1) + rho) + z(n)*z(n)/rho
                 if (c <= temp) then
                    tau = rho
                 else
                    del = d(n) - d(n - 1)
                    a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*del
                    if (a < zero) then
                       tau = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                 end if
                 ! it can be proved that
                     ! d(n)+rho/2 <= lambda(n) < d(n)+tau <= d(n)+rho
                 dltlb = midpt
                 dltub = rho
              else
                 del = d(n) - d(n - 1)
                 a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*del
                 if (a < zero) then
                    tau = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 ! it can be proved that
                     ! d(n) < d(n)+tau < lambda(n) < d(n)+rho/2
                 dltlb = zero
                 dltub = midpt
              end if
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 dlam = d(i) + tau
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              c = w - delta(n - 1)*dpsi - delta(n)*dphi
              a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
              b = delta(n - 1)*delta(n)*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                ! eta = b/a
                 ! eta = rho - tau
                 eta = dltub - tau
              else if (a >= zero) then
                 eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
              else
                 eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta > zero) eta = -w/(dpsi + dphi)
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              tau = tau + eta
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/delta(n)
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    dlam = d(i) + tau
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 c = w - delta(n - 1)*dpsi - delta(n)*dphi
                 a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
                 b = delta(n - 1)*delta(n)*w
                 if (a >= zero) then
                    eta = (a + sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a - sqrt(abs(a*a - four*b*c)))
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta > zero) eta = -w/(dpsi + dphi)
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 temp = z(n)/delta(n)
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)

                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              dlam = d(i) + tau
              go to 250
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              del = d(ip1) - d(i)
              midpt = del/two
              do j = 1,n
                 delta(j) = (d(j) - d(i)) - midpt
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/delta(j)
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/delta(j)
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/delta(i) + z(ip1)*z(ip1)/delta(ip1)
              if (w > zero) then
                 ! d(i)< the ith eigenvalue < (d(i)+d(i+1))/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 a = c*del + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*del
                 if (a > zero) then
                    tau = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 dltlb = zero
                 dltub = midpt
              else
                 ! (d(i)+d(i+1))/2 <= the ith eigenvalue < d(i+1)
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 a = c*del - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*del
                 if (a < zero) then
                    tau = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 dltlb = -midpt
                 dltub = zero
              end if
              if (orgati) then
                 do j = 1,n
                    delta(j) = (d(j) - d(i)) - tau
                 end do
              else
                 do j = 1,n
                    delta(j) = (d(j) - d(ip1)) - tau
                 end do
              end if
              if (orgati) then
                 ii = i
              else
                 ii = i + 1
              end if
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              w = rhoinv + phi + psi
              ! w is the value of the secular function with
              ! its ii-th element removed.
              swtch3 = .false.
              if (orgati) then
                 if (w < zero) swtch3 = .true.
              else
                 if (w > zero) swtch3 = .true.
              end if
              if (ii == 1 .or. ii == n) swtch3 = .false.
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                        *dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 if (orgati) then
                    dlam = d(i) + tau
                 else
                    dlam = d(ip1) + tau
                 end if
                 go to 250
              end if
              if (w <= zero) then
                 dltlb = max(dltlb,tau)
              else
                 dltub = min(dltub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 if (orgati) then
                    c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                              **2
                 else
                    c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                              **2
                 end if
                 a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                 b = delta(i)*delta(ip1)*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                       end if
                    end if
                    eta = b/a
                 else if (a <= zero) then
                    eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 else
                    eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 end if
              else
                 ! interpolation using three most relevant poles
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/delta(iim1)
                    temp1 = temp1*temp1
                    c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                 else
                    temp1 = z(iip1)/delta(iip1)
                    temp1 = temp1*temp1
                    c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1))*temp1
                    zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 call la_qlaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                 if (info /= 0) go to 250
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              temp = tau + eta
              if (temp > dltub .or. temp < dltlb) then
                 if (w < zero) then
                    eta = (dltub - tau)/two
                 else
                    eta = (dltlb - tau)/two
                 end if
              end if
              prew = w
              do j = 1,n
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/delta(j)
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/delta(j)
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              temp = z(ii)/delta(ii)
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau + eta) &
                        *dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              tau = tau + eta
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_240: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    if (orgati) then
                       dlam = d(i) + tau
                    else
                       dlam = d(ip1) + tau
                    end if
                    go to 250
                 end if
                 if (w <= zero) then
                    dltlb = max(dltlb,tau)
                 else
                    dltub = min(dltub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i)) &
                                    **2
                       else
                          c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1)) &
                                    **2
                       end if
                    else
                       temp = z(ii)/delta(ii)
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - delta(i)*dpsi - delta(ip1)*dphi
                    end if
                    a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                    b = delta(i)*delta(ip1)*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)

                             else
                                a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                             end if
                          else
                             a = delta(i)*delta(i)*dpsi + delta(ip1)*delta(ip1) &
                                       *dphi
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 else
                    ! interpolation using three most relevant poles
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - delta(iim1)*dpsi - delta(iip1)*dphi
                       zz(1) = delta(iim1)*delta(iim1)*dpsi
                       zz(3) = delta(iip1)*delta(iip1)*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/delta(iim1)
                          temp1 = temp1*temp1
                          c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1)) &
                                    *temp1
                          zz(1) = z(iim1)*z(iim1)
                          zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                       else
                          temp1 = z(iip1)/delta(iip1)
                          temp1 = temp1*temp1
                          c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1)) &
                                    *temp1
                          zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    call la_qlaed6(niter,orgati,c,delta(iim1),zz,w,eta,info)
                    if (info /= 0) go to 250
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 temp = tau + eta
                 if (temp > dltub .or. temp < dltlb) then
                    if (w < zero) then
                       eta = (dltub - tau)/two
                    else
                       eta = (dltlb - tau)/two
                    end if
                 end if
                 do j = 1,n
                    delta(j) = delta(j) - eta
                 end do
                 tau = tau + eta
                 prew = w
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/delta(j)
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/delta(j)
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 temp = z(ii)/delta(ii)
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp) + abs(tau) &
                           *dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_240
              ! return with info = 1, niter = maxit and not converged
              info = 1
              if (orgati) then
                 dlam = d(i) + tau
              else
                 dlam = d(ip1) + tau
              end if
           end if
           250 continue
           return
     end subroutine la_qlaed4
#endif

     !> SLAED8: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_slaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,z,dlamda, &
               q2,ldq2,w,perm,givptr,givcol,givnum,indxp,indx,info)
        use la_constants_sp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,icompq,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(sp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(sp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(sp),intent(out) :: dlamda(*),givnum(2,*),q2(ldq2,*),w(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: mone = -1.0_sp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(sp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -3
           else if (icompq == 1 .and. qsiz < n) then
              info = -4
           else if (ldq < max(1,n)) then
              info = -7
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -10
           else if (ldq2 < max(1,n)) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('SLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_sscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_sscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_slamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_isamax(n,z,1)
           jmax = la_isamax(n,d,1)
           eps = la_slamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              if (icompq == 0) then
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                 end do
              else
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                    call la_scopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
                 end do
                 call la_slacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              end if
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 110
              else
                 jlam = j
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           if (j > n) go to 100
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_slapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 if (icompq == 1) then
                    call la_srot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j &
                              ))),1,c,s)
                 end if
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           110 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           if (icompq == 0) then
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
              end do
           else
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
                 call la_scopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
           end if
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              if (icompq == 0) then
                 call la_scopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              else
                 call la_scopy(n - k,dlamda(k + 1),1,d(k + 1),1)
                 call la_slacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
              end if
           end if
           return
     end subroutine la_slaed8
     !> DLAED8: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_dlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,z,dlamda, &
               q2,ldq2,w,perm,givptr,givcol,givnum,indxp,indx,info)
        use la_constants_dp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,icompq,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(dp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(dp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(dp),intent(out) :: dlamda(*),givnum(2,*),q2(ldq2,*),w(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: mone = -1.0_dp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(dp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -3
           else if (icompq == 1 .and. qsiz < n) then
              info = -4
           else if (ldq < max(1,n)) then
              info = -7
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -10
           else if (ldq2 < max(1,n)) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('DLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_dscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_dscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_dlamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_idamax(n,z,1)
           jmax = la_idamax(n,d,1)
           eps = la_dlamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              if (icompq == 0) then
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                 end do
              else
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                    call la_dcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
                 end do
                 call la_dlacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              end if
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 110
              else
                 jlam = j
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           if (j > n) go to 100
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_dlapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 if (icompq == 1) then
                    call la_drot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j &
                              ))),1,c,s)
                 end if
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           110 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           if (icompq == 0) then
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
              end do
           else
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
                 call la_dcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
           end if
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              if (icompq == 0) then
                 call la_dcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              else
                 call la_dcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
                 call la_dlacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
              end if
           end if
           return
     end subroutine la_dlaed8
#ifdef LA_WITH_XDP
     !> XLAED8: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_xlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,z,dlamda, &
               q2,ldq2,w,perm,givptr,givcol,givnum,indxp,indx,info)
        use la_constants_xdp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,icompq,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(xdp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(xdp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(xdp),intent(out) :: dlamda(*),givnum(2,*),q2(ldq2,*),w(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: mone = -1.0_xdp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(xdp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -3
           else if (icompq == 1 .and. qsiz < n) then
              info = -4
           else if (ldq < max(1,n)) then
              info = -7
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -10
           else if (ldq2 < max(1,n)) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('XLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_xscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_xscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_xlamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_ixamax(n,z,1)
           jmax = la_ixamax(n,d,1)
           eps = la_xlamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              if (icompq == 0) then
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                 end do
              else
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                    call la_xcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
                 end do
                 call la_xlacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              end if
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 110
              else
                 jlam = j
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           if (j > n) go to 100
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_xlapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 if (icompq == 1) then
                    call la_xrot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j &
                              ))),1,c,s)
                 end if
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           110 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           if (icompq == 0) then
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
              end do
           else
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
                 call la_xcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
           end if
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              if (icompq == 0) then
                 call la_xcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              else
                 call la_xcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
                 call la_xlacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
              end if
           end if
           return
     end subroutine la_xlaed8
#endif
#ifdef LA_WITH_QP
     !> QLAED8: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_qlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,z,dlamda, &
               q2,ldq2,w,perm,givptr,givcol,givnum,indxp,indx,info)
        use la_constants_qp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,icompq,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(qp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(qp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(qp),intent(out) :: dlamda(*),givnum(2,*),q2(ldq2,*),w(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: mone = -1.0_qp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(qp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -3
           else if (icompq == 1 .and. qsiz < n) then
              info = -4
           else if (ldq < max(1,n)) then
              info = -7
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -10
           else if (ldq2 < max(1,n)) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('QLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_qscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_qscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_qlamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_iqamax(n,z,1)
           jmax = la_iqamax(n,d,1)
           eps = la_qlamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              if (icompq == 0) then
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                 end do
              else
                 do j = 1,n
                    perm(j) = indxq(indx(j))
                    call la_qcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
                 end do
                 call la_qlacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              end if
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 110
              else
                 jlam = j
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           if (j > n) go to 100
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_qlapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 if (icompq == 1) then
                    call la_qrot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j &
                              ))),1,c,s)
                 end if
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           110 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           if (icompq == 0) then
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
              end do
           else
              do j = 1,n
                 jp = indxp(j)
                 dlamda(j) = d(jp)
                 perm(j) = indxq(indx(jp))
                 call la_qcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
           end if
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              if (icompq == 0) then
                 call la_qcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              else
                 call la_qcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
                 call la_qlacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
              end if
           end if
           return
     end subroutine la_qlaed8
#endif

     !> SLAED9: finds the roots of the secular equation, as defined by the
     !> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
     !> appropriate calls to SLAED4 and then stores the new matrix of
     !> eigenvectors for use in calculating the next level of Z vectors.

     pure subroutine la_slaed9(k,kstart,kstop,n,d,q,ldq,rho,dlamda,w,s,lds,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,kstart,kstop,ldq,lds,n
           real(sp),intent(in) :: rho
           ! Array Arguments
           real(sp),intent(out) :: d(*),q(ldq,*),s(lds,*)
           real(sp),intent(inout) :: dlamda(*),w(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (kstart < 1 .or. kstart > max(1,k)) then
              info = -2
           else if (max(1,kstop) < kstart .or. kstop > max(1,k)) then
              info = -3
           else if (n < k) then
              info = -4
           else if (ldq < max(1,k)) then
              info = -7
           else if (lds < max(1,k)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('SLAED9',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,n
              dlamda(i) = la_slamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = kstart,kstop
              call la_slaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1 .or. k == 2) then
              do i = 1,k
                 do j = 1,k
                    s(j,i) = q(j,i)
                 end do
              end do
              go to 120
           end if
           ! compute updated w.
           call la_scopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_scopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i,1))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 q(i,j) = w(i)/q(i,j)
              end do
              temp = la_snrm2(k,q(1,j),1)
              do i = 1,k
                 s(i,j) = q(i,j)/temp
              end do
           end do
           120 continue
           return
     end subroutine la_slaed9
     !> DLAED9: finds the roots of the secular equation, as defined by the
     !> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
     !> appropriate calls to DLAED4 and then stores the new matrix of
     !> eigenvectors for use in calculating the next level of Z vectors.

     pure subroutine la_dlaed9(k,kstart,kstop,n,d,q,ldq,rho,dlamda,w,s,lds,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,kstart,kstop,ldq,lds,n
           real(dp),intent(in) :: rho
           ! Array Arguments
           real(dp),intent(out) :: d(*),q(ldq,*),s(lds,*)
           real(dp),intent(inout) :: dlamda(*),w(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (kstart < 1 .or. kstart > max(1,k)) then
              info = -2
           else if (max(1,kstop) < kstart .or. kstop > max(1,k)) then
              info = -3
           else if (n < k) then
              info = -4
           else if (ldq < max(1,k)) then
              info = -7
           else if (lds < max(1,k)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('DLAED9',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,n
              dlamda(i) = la_dlamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = kstart,kstop
              call la_dlaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1 .or. k == 2) then
              do i = 1,k
                 do j = 1,k
                    s(j,i) = q(j,i)
                 end do
              end do
              go to 120
           end if
           ! compute updated w.
           call la_dcopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_dcopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i,1))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 q(i,j) = w(i)/q(i,j)
              end do
              temp = la_dnrm2(k,q(1,j),1)
              do i = 1,k
                 s(i,j) = q(i,j)/temp
              end do
           end do
           120 continue
           return
     end subroutine la_dlaed9
#ifdef LA_WITH_XDP
     !> XLAED9: finds the roots of the secular equation, as defined by the
     !> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
     !> appropriate calls to XLAED4 and then stores the new matrix of
     !> eigenvectors for use in calculating the next level of Z vectors.

     pure subroutine la_xlaed9(k,kstart,kstop,n,d,q,ldq,rho,dlamda,w,s,lds,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,kstart,kstop,ldq,lds,n
           real(xdp),intent(in) :: rho
           ! Array Arguments
           real(xdp),intent(out) :: d(*),q(ldq,*),s(lds,*)
           real(xdp),intent(inout) :: dlamda(*),w(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (kstart < 1 .or. kstart > max(1,k)) then
              info = -2
           else if (max(1,kstop) < kstart .or. kstop > max(1,k)) then
              info = -3
           else if (n < k) then
              info = -4
           else if (ldq < max(1,k)) then
              info = -7
           else if (lds < max(1,k)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('XLAED9',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,n
              dlamda(i) = la_xlamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = kstart,kstop
              call la_xlaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1 .or. k == 2) then
              do i = 1,k
                 do j = 1,k
                    s(j,i) = q(j,i)
                 end do
              end do
              go to 120
           end if
           ! compute updated w.
           call la_xcopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_xcopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i,1))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 q(i,j) = w(i)/q(i,j)
              end do
              temp = la_xnrm2(k,q(1,j),1)
              do i = 1,k
                 s(i,j) = q(i,j)/temp
              end do
           end do
           120 continue
           return
     end subroutine la_xlaed9
#endif
#ifdef LA_WITH_QP
     !> QLAED9: finds the roots of the secular equation, as defined by the
     !> values in D, Z, and RHO, between KSTART and KSTOP.  It makes the
     !> appropriate calls to QLAED4 and then stores the new matrix of
     !> eigenvectors for use in calculating the next level of Z vectors.

     pure subroutine la_qlaed9(k,kstart,kstop,n,d,q,ldq,rho,dlamda,w,s,lds,info)

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,kstart,kstop,ldq,lds,n
           real(qp),intent(in) :: rho
           ! Array Arguments
           real(qp),intent(out) :: d(*),q(ldq,*),s(lds,*)
           real(qp),intent(inout) :: dlamda(*),w(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (kstart < 1 .or. kstart > max(1,k)) then
              info = -2
           else if (max(1,kstop) < kstart .or. kstop > max(1,k)) then
              info = -3
           else if (n < k) then
              info = -4
           else if (ldq < max(1,k)) then
              info = -7
           else if (lds < max(1,k)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('QLAED9',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,n
              dlamda(i) = la_qlamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = kstart,kstop
              call la_qlaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1 .or. k == 2) then
              do i = 1,k
                 do j = 1,k
                    s(j,i) = q(j,i)
                 end do
              end do
              go to 120
           end if
           ! compute updated w.
           call la_qcopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_qcopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i,1))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 q(i,j) = w(i)/q(i,j)
              end do
              temp = la_qnrm2(k,q(1,j),1)
              do i = 1,k
                 s(i,j) = q(i,j)/temp
              end do
           end do
           120 continue
           return
     end subroutine la_qlaed9
#endif

     !> SLAED3: finds the roots of the secular equation, as defined by the
     !> values in D, W, and RHO, between 1 and K.  It makes the
     !> appropriate calls to SLAED4 and then updates the eigenvectors by
     !> multiplying the matrix of eigenvectors of the pair of eigensystems
     !> being combined by the matrix of eigenvectors of the K-by-K system
     !> which is solved here.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_slaed3(k,n,n1,d,q,ldq,rho,dlamda,q2,indx,ctot,w,s,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,n,n1
           real(sp),intent(in) :: rho
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),indx(*)
           real(sp),intent(out) :: d(*),q(ldq,*),s(*)
           real(sp),intent(inout) :: dlamda(*),w(*)
           real(sp),intent(in) :: q2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,iq2,j,n12,n2,n23
           real(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (n < k) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SLAED3',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dlamda(i) = la_slamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = 1,k
              call la_slaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1) go to 110
           if (k == 2) then
              do j = 1,k
                 w(1) = q(1,j)
                 w(2) = q(2,j)
                 ii = indx(1)
                 q(1,j) = w(ii)
                 ii = indx(2)
                 q(2,j) = w(ii)
              end do
              go to 110
           end if
           ! compute updated w.
           call la_scopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_scopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 s(i) = w(i)/q(i,j)
              end do
              temp = la_snrm2(k,s,1)
              do i = 1,k
                 ii = indx(i)
                 q(i,j) = s(ii)/temp
              end do
           end do
           ! compute the updated eigenvectors.
           110 continue
           n2 = n - n1
           n12 = ctot(1) + ctot(2)
           n23 = ctot(2) + ctot(3)
           call la_slacpy('A',n23,k,q(ctot(1) + 1,1),ldq,s,n23)
           iq2 = n1*n12 + 1
           if (n23 /= 0) then
              call la_sgemm('N','N',n2,k,n23,one,q2(iq2),n2,s,n23,zero,q(n1 + 1, &
                        1),ldq)
           else
              call la_slaset('A',n2,k,zero,zero,q(n1 + 1,1),ldq)
           end if
           call la_slacpy('A',n12,k,q,ldq,s,n12)
           if (n12 /= 0) then
              call la_sgemm('N','N',n1,k,n12,one,q2,n1,s,n12,zero,q,ldq)
           else
              call la_slaset('A',n1,k,zero,zero,q(1,1),ldq)
           end if
           120 continue
           return
     end subroutine la_slaed3
     !> DLAED3: finds the roots of the secular equation, as defined by the
     !> values in D, W, and RHO, between 1 and K.  It makes the
     !> appropriate calls to DLAED4 and then updates the eigenvectors by
     !> multiplying the matrix of eigenvectors of the pair of eigensystems
     !> being combined by the matrix of eigenvectors of the K-by-K system
     !> which is solved here.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_dlaed3(k,n,n1,d,q,ldq,rho,dlamda,q2,indx,ctot,w,s,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,n,n1
           real(dp),intent(in) :: rho
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),indx(*)
           real(dp),intent(out) :: d(*),q(ldq,*),s(*)
           real(dp),intent(inout) :: dlamda(*),w(*)
           real(dp),intent(in) :: q2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,iq2,j,n12,n2,n23
           real(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (n < k) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DLAED3',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dlamda(i) = la_dlamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = 1,k
              call la_dlaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1) go to 110
           if (k == 2) then
              do j = 1,k
                 w(1) = q(1,j)
                 w(2) = q(2,j)
                 ii = indx(1)
                 q(1,j) = w(ii)
                 ii = indx(2)
                 q(2,j) = w(ii)
              end do
              go to 110
           end if
           ! compute updated w.
           call la_dcopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_dcopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 s(i) = w(i)/q(i,j)
              end do
              temp = la_dnrm2(k,s,1)
              do i = 1,k
                 ii = indx(i)
                 q(i,j) = s(ii)/temp
              end do
           end do
           ! compute the updated eigenvectors.
           110 continue
           n2 = n - n1
           n12 = ctot(1) + ctot(2)
           n23 = ctot(2) + ctot(3)
           call la_dlacpy('A',n23,k,q(ctot(1) + 1,1),ldq,s,n23)
           iq2 = n1*n12 + 1
           if (n23 /= 0) then
              call la_dgemm('N','N',n2,k,n23,one,q2(iq2),n2,s,n23,zero,q(n1 + 1, &
                        1),ldq)
           else
              call la_dlaset('A',n2,k,zero,zero,q(n1 + 1,1),ldq)
           end if
           call la_dlacpy('A',n12,k,q,ldq,s,n12)
           if (n12 /= 0) then
              call la_dgemm('N','N',n1,k,n12,one,q2,n1,s,n12,zero,q,ldq)
           else
              call la_dlaset('A',n1,k,zero,zero,q(1,1),ldq)
           end if
           120 continue
           return
     end subroutine la_dlaed3
#ifdef LA_WITH_XDP
     !> XLAED3: finds the roots of the secular equation, as defined by the
     !> values in D, W, and RHO, between 1 and K.  It makes the
     !> appropriate calls to XLAED4 and then updates the eigenvectors by
     !> multiplying the matrix of eigenvectors of the pair of eigensystems
     !> being combined by the matrix of eigenvectors of the K-by-K system
     !> which is solved here.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_xlaed3(k,n,n1,d,q,ldq,rho,dlamda,q2,indx,ctot,w,s,info)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,n,n1
           real(xdp),intent(in) :: rho
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),indx(*)
           real(xdp),intent(out) :: d(*),q(ldq,*),s(*)
           real(xdp),intent(inout) :: dlamda(*),w(*)
           real(xdp),intent(in) :: q2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,iq2,j,n12,n2,n23
           real(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (n < k) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('XLAED3',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dlamda(i) = la_xlamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = 1,k
              call la_xlaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1) go to 110
           if (k == 2) then
              do j = 1,k
                 w(1) = q(1,j)
                 w(2) = q(2,j)
                 ii = indx(1)
                 q(1,j) = w(ii)
                 ii = indx(2)
                 q(2,j) = w(ii)
              end do
              go to 110
           end if
           ! compute updated w.
           call la_xcopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_xcopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 s(i) = w(i)/q(i,j)
              end do
              temp = la_xnrm2(k,s,1)
              do i = 1,k
                 ii = indx(i)
                 q(i,j) = s(ii)/temp
              end do
           end do
           ! compute the updated eigenvectors.
           110 continue
           n2 = n - n1
           n12 = ctot(1) + ctot(2)
           n23 = ctot(2) + ctot(3)
           call la_xlacpy('A',n23,k,q(ctot(1) + 1,1),ldq,s,n23)
           iq2 = n1*n12 + 1
           if (n23 /= 0) then
              call la_xgemm('N','N',n2,k,n23,one,q2(iq2),n2,s,n23,zero,q(n1 + 1, &
                        1),ldq)
           else
              call la_xlaset('A',n2,k,zero,zero,q(n1 + 1,1),ldq)
           end if
           call la_xlacpy('A',n12,k,q,ldq,s,n12)
           if (n12 /= 0) then
              call la_xgemm('N','N',n1,k,n12,one,q2,n1,s,n12,zero,q,ldq)
           else
              call la_xlaset('A',n1,k,zero,zero,q(1,1),ldq)
           end if
           120 continue
           return
     end subroutine la_xlaed3
#endif
#ifdef LA_WITH_QP
     !> QLAED3: finds the roots of the secular equation, as defined by the
     !> values in D, W, and RHO, between 1 and K.  It makes the
     !> appropriate calls to QLAED4 and then updates the eigenvectors by
     !> multiplying the matrix of eigenvectors of the pair of eigensystems
     !> being combined by the matrix of eigenvectors of the K-by-K system
     !> which is solved here.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_qlaed3(k,n,n1,d,q,ldq,rho,dlamda,q2,indx,ctot,w,s,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,n,n1
           real(qp),intent(in) :: rho
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),indx(*)
           real(qp),intent(out) :: d(*),q(ldq,*),s(*)
           real(qp),intent(inout) :: dlamda(*),w(*)
           real(qp),intent(in) :: q2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,ii,iq2,j,n12,n2,n23
           real(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (k < 0) then
              info = -1
           else if (n < k) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QLAED3',-info)
              return
           end if
           ! quick return if possible
           if (k == 0) return
           ! modify values dlamda(i) to make sure all dlamda(i)-dlamda(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dlamda(i) by 2*dlamda(i)-dlamda(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dlamda(i) if it is 1; this makes the subsequent
           ! subtractions dlamda(i)-dlamda(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dlamda(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dlamda(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dlamda(i) = la_qlamc3(dlamda(i),dlamda(i)) - dlamda(i)
           end do
           do j = 1,k
              call la_qlaed4(k,j,dlamda,w,q(1,j),rho,d(j),info)
              ! if the zero finder fails, the computation is terminated.
              if (info /= 0) go to 120
           end do
           if (k == 1) go to 110
           if (k == 2) then
              do j = 1,k
                 w(1) = q(1,j)
                 w(2) = q(2,j)
                 ii = indx(1)
                 q(1,j) = w(ii)
                 ii = indx(2)
                 q(2,j) = w(ii)
              end do
              go to 110
           end if
           ! compute updated w.
           call la_qcopy(k,w,1,s,1)
           ! initialize w(i) = q(i,i)
           call la_qcopy(k,q,ldq + 1,w,1)
           do j = 1,k
              do i = 1,j - 1
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
              do i = j + 1,k
                 w(i) = w(i)*(q(i,j)/(dlamda(i) - dlamda(j)))
              end do
           end do
           do i = 1,k
              w(i) = sign(sqrt(-w(i)),s(i))
           end do
           ! compute eigenvectors of the modified rank-1 modification.
           do j = 1,k
              do i = 1,k
                 s(i) = w(i)/q(i,j)
              end do
              temp = la_qnrm2(k,s,1)
              do i = 1,k
                 ii = indx(i)
                 q(i,j) = s(ii)/temp
              end do
           end do
           ! compute the updated eigenvectors.
           110 continue
           n2 = n - n1
           n12 = ctot(1) + ctot(2)
           n23 = ctot(2) + ctot(3)
           call la_qlacpy('A',n23,k,q(ctot(1) + 1,1),ldq,s,n23)
           iq2 = n1*n12 + 1
           if (n23 /= 0) then
              call la_qgemm('N','N',n2,k,n23,one,q2(iq2),n2,s,n23,zero,q(n1 + 1, &
                        1),ldq)
           else
              call la_qlaset('A',n2,k,zero,zero,q(n1 + 1,1),ldq)
           end if
           call la_qlacpy('A',n12,k,q,ldq,s,n12)
           if (n12 /= 0) then
              call la_qgemm('N','N',n1,k,n12,one,q2,n1,s,n12,zero,q,ldq)
           else
              call la_qlaset('A',n1,k,zero,zero,q(1,1),ldq)
           end if
           120 continue
           return
     end subroutine la_qlaed3
#endif

     !> SLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense symmetric matrix
     !> that has been reduced to tridiagonal form.  SLAED1 handles
     !> the case in which all eigenvalues and eigenvectors of a symmetric
     !> tridiagonal matrix are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**Tu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine SLAED8.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine SLAED4 (as called by SLAED9).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_slaed7(icompq,n,qsiz,tlvls,curlvl,curpbm,d,q,ldq,indxq,rho, &
               cutpnt,qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,icompq,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(sp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(sp),intent(inout) :: d(*),givnum(2,*),q(ldq,*),qstore(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,ldq2, &
                     n1,n2,ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (icompq == 1 .and. qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('SLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_slaed8 and la_slaed9.
           if (icompq == 1) then
              ldq2 = qsiz
           else
              ldq2 = n
           end if
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           is = iq2 + n*ldq2
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_slaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,work(iz),work(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_slaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,work(iz),work( &
            idlmda),work(iq2),ldq2,work(iw),perm(prmptr(curr)),givptr(curr + 1), &
            givcol(1,givptr(curr)),givnum(1,givptr(curr)),iwork(indxp),iwork(indx), &
                       info)
           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_slaed9(k,1,k,n,d,work(is),k,rho,work(idlmda),work(iw), &
                        qstore(qptr(curr)),k,info)
              if (info /= 0) go to 30
              if (icompq == 1) then
                 call la_sgemm('N','N',qsiz,k,k,one,work(iq2),ldq2,qstore(qptr( &
                           curr)),k,zero,q,ldq)
              end if
              qptr(curr + 1) = qptr(curr) + k**2
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_slamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           30 continue
           return
     end subroutine la_slaed7
     !> DLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense symmetric matrix
     !> that has been reduced to tridiagonal form.  DLAED1 handles
     !> the case in which all eigenvalues and eigenvectors of a symmetric
     !> tridiagonal matrix are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**Tu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine DLAED8.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine DLAED4 (as called by DLAED9).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_dlaed7(icompq,n,qsiz,tlvls,curlvl,curpbm,d,q,ldq,indxq,rho, &
               cutpnt,qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,icompq,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(dp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(dp),intent(inout) :: d(*),givnum(2,*),q(ldq,*),qstore(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,ldq2, &
                     n1,n2,ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (icompq == 1 .and. qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('DLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_dlaed8 and la_dlaed9.
           if (icompq == 1) then
              ldq2 = qsiz
           else
              ldq2 = n
           end if
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           is = iq2 + n*ldq2
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_dlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,work(iz),work(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_dlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,work(iz),work( &
            idlmda),work(iq2),ldq2,work(iw),perm(prmptr(curr)),givptr(curr + 1), &
            givcol(1,givptr(curr)),givnum(1,givptr(curr)),iwork(indxp),iwork(indx), &
                       info)
           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_dlaed9(k,1,k,n,d,work(is),k,rho,work(idlmda),work(iw), &
                        qstore(qptr(curr)),k,info)
              if (info /= 0) go to 30
              if (icompq == 1) then
                 call la_dgemm('N','N',qsiz,k,k,one,work(iq2),ldq2,qstore(qptr( &
                           curr)),k,zero,q,ldq)
              end if
              qptr(curr + 1) = qptr(curr) + k**2
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_dlamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           30 continue
           return
     end subroutine la_dlaed7
#ifdef LA_WITH_XDP
     !> XLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense symmetric matrix
     !> that has been reduced to tridiagonal form.  XLAED1 handles
     !> the case in which all eigenvalues and eigenvectors of a symmetric
     !> tridiagonal matrix are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**Tu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine XLAED8.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine XLAED4 (as called by XLAED9).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_xlaed7(icompq,n,qsiz,tlvls,curlvl,curpbm,d,q,ldq,indxq,rho, &
               cutpnt,qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,iwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,icompq,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(xdp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(xdp),intent(inout) :: d(*),givnum(2,*),q(ldq,*),qstore(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,ldq2, &
                     n1,n2,ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (icompq == 1 .and. qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('XLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_xlaed8 and la_xlaed9.
           if (icompq == 1) then
              ldq2 = qsiz
           else
              ldq2 = n
           end if
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           is = iq2 + n*ldq2
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_xlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,work(iz),work(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_xlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,work(iz),work( &
            idlmda),work(iq2),ldq2,work(iw),perm(prmptr(curr)),givptr(curr + 1), &
            givcol(1,givptr(curr)),givnum(1,givptr(curr)),iwork(indxp),iwork(indx), &
                       info)
           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_xlaed9(k,1,k,n,d,work(is),k,rho,work(idlmda),work(iw), &
                        qstore(qptr(curr)),k,info)
              if (info /= 0) go to 30
              if (icompq == 1) then
                 call la_xgemm('N','N',qsiz,k,k,one,work(iq2),ldq2,qstore(qptr( &
                           curr)),k,zero,q,ldq)
              end if
              qptr(curr + 1) = qptr(curr) + k**2
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_xlamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           30 continue
           return
     end subroutine la_xlaed7
#endif
#ifdef LA_WITH_QP
     !> QLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense symmetric matrix
     !> that has been reduced to tridiagonal form.  QLAED1 handles
     !> the case in which all eigenvalues and eigenvectors of a symmetric
     !> tridiagonal matrix are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**Tu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine QLAED8.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine QLAED4 (as called by QLAED9).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_qlaed7(icompq,n,qsiz,tlvls,curlvl,curpbm,d,q,ldq,indxq,rho, &
               cutpnt,qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,icompq,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(qp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(qp),intent(inout) :: d(*),givnum(2,*),q(ldq,*),qstore(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,ldq2, &
                     n1,n2,ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 1) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (icompq == 1 .and. qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('QLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_qlaed8 and la_qlaed9.
           if (icompq == 1) then
              ldq2 = qsiz
           else
              ldq2 = n
           end if
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           is = iq2 + n*ldq2
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_qlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,work(iz),work(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_qlaed8(icompq,k,n,qsiz,d,q,ldq,indxq,rho,cutpnt,work(iz),work( &
            idlmda),work(iq2),ldq2,work(iw),perm(prmptr(curr)),givptr(curr + 1), &
            givcol(1,givptr(curr)),givnum(1,givptr(curr)),iwork(indxp),iwork(indx), &
                       info)
           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_qlaed9(k,1,k,n,d,work(is),k,rho,work(idlmda),work(iw), &
                        qstore(qptr(curr)),k,info)
              if (info /= 0) go to 30
              if (icompq == 1) then
                 call la_qgemm('N','N',qsiz,k,k,one,work(iq2),ldq2,qstore(qptr( &
                           curr)),k,zero,q,ldq)
              end if
              qptr(curr + 1) = qptr(curr) + k**2
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_qlamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           30 continue
           return
     end subroutine la_qlaed7
#endif

     !> SLAED2: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_slaed2(k,n,n1,d,q,ldq,indxq,rho,z,dlamda,w,q2,indx,indxc, &
                indxp,coltyp,info)
        use la_constants_sp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldq,n,n1
           real(sp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),indx(*),indxc(*),indxp(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(sp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(sp),intent(out) :: dlamda(*),q2(*),w(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: mone = -1.0_sp

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,imax,iq1,iq2,j,jmax,js,k2,n1p1,n2,nj,pj
           real(sp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (min(1, (n/2)) > n1 .or. (n/2) < n1) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('SLAED2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_sscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1.  since z is the concatenation of
           ! two normalized vectors, norm2(z) = sqrt(2).
           t = one/sqrt(two)
           call la_sscal(n,t,z,1)
           ! rho = abs( norm(z)**2 * rho )
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = n1p1,n
              indxq(i) = indxq(i) + n1
           end do
           ! re-integrate the deflated parts from the last pass
           do i = 1,n
              dlamda(i) = d(indxq(i))
           end do
           call la_slamrg(n1,n2,dlamda,1,1,indxc)
           do i = 1,n
              indx(i) = indxq(indxc(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_isamax(n,z,1)
           jmax = la_isamax(n,d,1)
           eps = la_slamch('EPSILON')
           tol = eight*eps*max(abs(d(jmax)),abs(z(imax)))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              iq2 = 1
              do j = 1,n
                 i = indx(j)
                 call la_scopy(n,q(1,i),1,q2(iq2),1)
                 dlamda(j) = d(i)
                 iq2 = iq2 + n
              end do
              call la_slacpy('A',n,n,q2,n,q,ldq)
              call la_scopy(n,dlamda,1,d,1)
              go to 190
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           do i = 1,n1
              coltyp(i) = 1
           end do
           do i = n1p1,n
              coltyp(i) = 3
           end do
           k = 0
           k2 = n + 1
           do j = 1,n
              nj = indx(j)
              if (rho*abs(z(nj)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 coltyp(nj) = 4
                 indxp(k2) = nj
                 if (j == n) go to 100
              else
                 pj = nj
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           nj = indx(j)
           if (j > n) go to 100
           if (rho*abs(z(nj)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              coltyp(nj) = 4
              indxp(k2) = nj
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(pj)
              c = z(nj)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_slapy2(c,s)
              t = d(nj) - d(pj)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(nj) = tau
                 z(pj) = zero
                 if (coltyp(nj) /= coltyp(pj)) coltyp(nj) = 2
                 coltyp(pj) = 4
                 call la_srot(n,q(1,pj),1,q(1,nj),1,c,s)
                 t = d(pj)*c**2 + d(nj)*s**2
                 d(nj) = d(pj)*s**2 + d(nj)*c**2
                 d(pj) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(pj) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = pj
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = pj
                    end if
                 else
                    indxp(k2 + i - 1) = pj
                 end if
                 pj = nj
              else
                 k = k + 1
                 dlamda(k) = d(pj)
                 w(k) = z(pj)
                 indxp(k) = pj
                 pj = nj
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           dlamda(k) = d(pj)
           w(k) = z(pj)
           indxp(k) = pj
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four uniform groups (although one or more of these groups may be
           ! empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 1,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 1
           psm(2) = 1 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           k = n - ctot(4)
           ! fill out the indxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's.
           do j = 1,n
              js = indxp(j)
              ct = coltyp(js)
              indx(psm(ct)) = js
              indxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           i = 1
           iq1 = 1
           iq2 = 1 + (ctot(1) + ctot(2))*n1
           do j = 1,ctot(1)
              js = indx(i)
              call la_scopy(n1,q(1,js),1,q2(iq1),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
           end do
           do j = 1,ctot(2)
              js = indx(i)
              call la_scopy(n1,q(1,js),1,q2(iq1),1)
              call la_scopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
              iq2 = iq2 + n2
           end do
           do j = 1,ctot(3)
              js = indx(i)
              call la_scopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq2 = iq2 + n2
           end do
           iq1 = iq2
           do j = 1,ctot(4)
              js = indx(i)
              call la_scopy(n,q(1,js),1,q2(iq2),1)
              iq2 = iq2 + n
              z(i) = d(js)
              i = i + 1
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_slacpy('A',n,ctot(4),q2(iq1),n,q(1,k + 1),ldq)
              call la_scopy(n - k,z(k + 1),1,d(k + 1),1)
           end if
           ! copy ctot into coltyp for referencing in la_slaed3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           190 continue
           return
     end subroutine la_slaed2
     !> DLAED2: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_dlaed2(k,n,n1,d,q,ldq,indxq,rho,z,dlamda,w,q2,indx,indxc, &
                indxp,coltyp,info)
        use la_constants_dp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldq,n,n1
           real(dp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),indx(*),indxc(*),indxp(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(dp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(dp),intent(out) :: dlamda(*),q2(*),w(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: mone = -1.0_dp

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,imax,iq1,iq2,j,jmax,js,k2,n1p1,n2,nj,pj
           real(dp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (min(1, (n/2)) > n1 .or. (n/2) < n1) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('DLAED2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_dscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1.  since z is the concatenation of
           ! two normalized vectors, norm2(z) = sqrt(2).
           t = one/sqrt(two)
           call la_dscal(n,t,z,1)
           ! rho = abs( norm(z)**2 * rho )
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = n1p1,n
              indxq(i) = indxq(i) + n1
           end do
           ! re-integrate the deflated parts from the last pass
           do i = 1,n
              dlamda(i) = d(indxq(i))
           end do
           call la_dlamrg(n1,n2,dlamda,1,1,indxc)
           do i = 1,n
              indx(i) = indxq(indxc(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_idamax(n,z,1)
           jmax = la_idamax(n,d,1)
           eps = la_dlamch('EPSILON')
           tol = eight*eps*max(abs(d(jmax)),abs(z(imax)))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              iq2 = 1
              do j = 1,n
                 i = indx(j)
                 call la_dcopy(n,q(1,i),1,q2(iq2),1)
                 dlamda(j) = d(i)
                 iq2 = iq2 + n
              end do
              call la_dlacpy('A',n,n,q2,n,q,ldq)
              call la_dcopy(n,dlamda,1,d,1)
              go to 190
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           do i = 1,n1
              coltyp(i) = 1
           end do
           do i = n1p1,n
              coltyp(i) = 3
           end do
           k = 0
           k2 = n + 1
           do j = 1,n
              nj = indx(j)
              if (rho*abs(z(nj)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 coltyp(nj) = 4
                 indxp(k2) = nj
                 if (j == n) go to 100
              else
                 pj = nj
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           nj = indx(j)
           if (j > n) go to 100
           if (rho*abs(z(nj)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              coltyp(nj) = 4
              indxp(k2) = nj
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(pj)
              c = z(nj)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_dlapy2(c,s)
              t = d(nj) - d(pj)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(nj) = tau
                 z(pj) = zero
                 if (coltyp(nj) /= coltyp(pj)) coltyp(nj) = 2
                 coltyp(pj) = 4
                 call la_drot(n,q(1,pj),1,q(1,nj),1,c,s)
                 t = d(pj)*c**2 + d(nj)*s**2
                 d(nj) = d(pj)*s**2 + d(nj)*c**2
                 d(pj) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(pj) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = pj
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = pj
                    end if
                 else
                    indxp(k2 + i - 1) = pj
                 end if
                 pj = nj
              else
                 k = k + 1
                 dlamda(k) = d(pj)
                 w(k) = z(pj)
                 indxp(k) = pj
                 pj = nj
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           dlamda(k) = d(pj)
           w(k) = z(pj)
           indxp(k) = pj
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four uniform groups (although one or more of these groups may be
           ! empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 1,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 1
           psm(2) = 1 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           k = n - ctot(4)
           ! fill out the indxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's.
           do j = 1,n
              js = indxp(j)
              ct = coltyp(js)
              indx(psm(ct)) = js
              indxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           i = 1
           iq1 = 1
           iq2 = 1 + (ctot(1) + ctot(2))*n1
           do j = 1,ctot(1)
              js = indx(i)
              call la_dcopy(n1,q(1,js),1,q2(iq1),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
           end do
           do j = 1,ctot(2)
              js = indx(i)
              call la_dcopy(n1,q(1,js),1,q2(iq1),1)
              call la_dcopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
              iq2 = iq2 + n2
           end do
           do j = 1,ctot(3)
              js = indx(i)
              call la_dcopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq2 = iq2 + n2
           end do
           iq1 = iq2
           do j = 1,ctot(4)
              js = indx(i)
              call la_dcopy(n,q(1,js),1,q2(iq2),1)
              iq2 = iq2 + n
              z(i) = d(js)
              i = i + 1
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_dlacpy('A',n,ctot(4),q2(iq1),n,q(1,k + 1),ldq)
              call la_dcopy(n - k,z(k + 1),1,d(k + 1),1)
           end if
           ! copy ctot into coltyp for referencing in la_dlaed3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           190 continue
           return
     end subroutine la_dlaed2
#ifdef LA_WITH_XDP
     !> XLAED2: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_xlaed2(k,n,n1,d,q,ldq,indxq,rho,z,dlamda,w,q2,indx,indxc, &
                indxp,coltyp,info)
        use la_constants_xdp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldq,n,n1
           real(xdp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),indx(*),indxc(*),indxp(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(xdp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(xdp),intent(out) :: dlamda(*),q2(*),w(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: mone = -1.0_xdp

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,imax,iq1,iq2,j,jmax,js,k2,n1p1,n2,nj,pj
           real(xdp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (min(1, (n/2)) > n1 .or. (n/2) < n1) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('XLAED2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_xscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1.  since z is the concatenation of
           ! two normalized vectors, norm2(z) = sqrt(2).
           t = one/sqrt(two)
           call la_xscal(n,t,z,1)
           ! rho = abs( norm(z)**2 * rho )
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = n1p1,n
              indxq(i) = indxq(i) + n1
           end do
           ! re-integrate the deflated parts from the last pass
           do i = 1,n
              dlamda(i) = d(indxq(i))
           end do
           call la_xlamrg(n1,n2,dlamda,1,1,indxc)
           do i = 1,n
              indx(i) = indxq(indxc(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_ixamax(n,z,1)
           jmax = la_ixamax(n,d,1)
           eps = la_xlamch('EPSILON')
           tol = eight*eps*max(abs(d(jmax)),abs(z(imax)))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              iq2 = 1
              do j = 1,n
                 i = indx(j)
                 call la_xcopy(n,q(1,i),1,q2(iq2),1)
                 dlamda(j) = d(i)
                 iq2 = iq2 + n
              end do
              call la_xlacpy('A',n,n,q2,n,q,ldq)
              call la_xcopy(n,dlamda,1,d,1)
              go to 190
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           do i = 1,n1
              coltyp(i) = 1
           end do
           do i = n1p1,n
              coltyp(i) = 3
           end do
           k = 0
           k2 = n + 1
           do j = 1,n
              nj = indx(j)
              if (rho*abs(z(nj)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 coltyp(nj) = 4
                 indxp(k2) = nj
                 if (j == n) go to 100
              else
                 pj = nj
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           nj = indx(j)
           if (j > n) go to 100
           if (rho*abs(z(nj)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              coltyp(nj) = 4
              indxp(k2) = nj
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(pj)
              c = z(nj)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_xlapy2(c,s)
              t = d(nj) - d(pj)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(nj) = tau
                 z(pj) = zero
                 if (coltyp(nj) /= coltyp(pj)) coltyp(nj) = 2
                 coltyp(pj) = 4
                 call la_xrot(n,q(1,pj),1,q(1,nj),1,c,s)
                 t = d(pj)*c**2 + d(nj)*s**2
                 d(nj) = d(pj)*s**2 + d(nj)*c**2
                 d(pj) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(pj) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = pj
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = pj
                    end if
                 else
                    indxp(k2 + i - 1) = pj
                 end if
                 pj = nj
              else
                 k = k + 1
                 dlamda(k) = d(pj)
                 w(k) = z(pj)
                 indxp(k) = pj
                 pj = nj
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           dlamda(k) = d(pj)
           w(k) = z(pj)
           indxp(k) = pj
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four uniform groups (although one or more of these groups may be
           ! empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 1,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 1
           psm(2) = 1 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           k = n - ctot(4)
           ! fill out the indxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's.
           do j = 1,n
              js = indxp(j)
              ct = coltyp(js)
              indx(psm(ct)) = js
              indxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           i = 1
           iq1 = 1
           iq2 = 1 + (ctot(1) + ctot(2))*n1
           do j = 1,ctot(1)
              js = indx(i)
              call la_xcopy(n1,q(1,js),1,q2(iq1),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
           end do
           do j = 1,ctot(2)
              js = indx(i)
              call la_xcopy(n1,q(1,js),1,q2(iq1),1)
              call la_xcopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
              iq2 = iq2 + n2
           end do
           do j = 1,ctot(3)
              js = indx(i)
              call la_xcopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq2 = iq2 + n2
           end do
           iq1 = iq2
           do j = 1,ctot(4)
              js = indx(i)
              call la_xcopy(n,q(1,js),1,q2(iq2),1)
              iq2 = iq2 + n
              z(i) = d(js)
              i = i + 1
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_xlacpy('A',n,ctot(4),q2(iq1),n,q(1,k + 1),ldq)
              call la_xcopy(n - k,z(k + 1),1,d(k + 1),1)
           end if
           ! copy ctot into coltyp for referencing in la_xlaed3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           190 continue
           return
     end subroutine la_xlaed2
#endif
#ifdef LA_WITH_QP
     !> QLAED2: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_qlaed2(k,n,n1,d,q,ldq,indxq,rho,z,dlamda,w,q2,indx,indxc, &
                indxp,coltyp,info)
        use la_constants_qp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldq,n,n1
           real(qp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),indx(*),indxc(*),indxp(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(qp),intent(inout) :: d(*),q(ldq,*),z(*)
           real(qp),intent(out) :: dlamda(*),q2(*),w(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: mone = -1.0_qp

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,imax,iq1,iq2,j,jmax,js,k2,n1p1,n2,nj,pj
           real(qp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (min(1, (n/2)) > n1 .or. (n/2) < n1) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('QLAED2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_qscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1.  since z is the concatenation of
           ! two normalized vectors, norm2(z) = sqrt(2).
           t = one/sqrt(two)
           call la_qscal(n,t,z,1)
           ! rho = abs( norm(z)**2 * rho )
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = n1p1,n
              indxq(i) = indxq(i) + n1
           end do
           ! re-integrate the deflated parts from the last pass
           do i = 1,n
              dlamda(i) = d(indxq(i))
           end do
           call la_qlamrg(n1,n2,dlamda,1,1,indxc)
           do i = 1,n
              indx(i) = indxq(indxc(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_iqamax(n,z,1)
           jmax = la_iqamax(n,d,1)
           eps = la_qlamch('EPSILON')
           tol = eight*eps*max(abs(d(jmax)),abs(z(imax)))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              iq2 = 1
              do j = 1,n
                 i = indx(j)
                 call la_qcopy(n,q(1,i),1,q2(iq2),1)
                 dlamda(j) = d(i)
                 iq2 = iq2 + n
              end do
              call la_qlacpy('A',n,n,q2,n,q,ldq)
              call la_qcopy(n,dlamda,1,d,1)
              go to 190
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           do i = 1,n1
              coltyp(i) = 1
           end do
           do i = n1p1,n
              coltyp(i) = 3
           end do
           k = 0
           k2 = n + 1
           do j = 1,n
              nj = indx(j)
              if (rho*abs(z(nj)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 coltyp(nj) = 4
                 indxp(k2) = nj
                 if (j == n) go to 100
              else
                 pj = nj
                 go to 80
              end if
           end do
           80 continue
           j = j + 1
           nj = indx(j)
           if (j > n) go to 100
           if (rho*abs(z(nj)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              coltyp(nj) = 4
              indxp(k2) = nj
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(pj)
              c = z(nj)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_qlapy2(c,s)
              t = d(nj) - d(pj)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(nj) = tau
                 z(pj) = zero
                 if (coltyp(nj) /= coltyp(pj)) coltyp(nj) = 2
                 coltyp(pj) = 4
                 call la_qrot(n,q(1,pj),1,q(1,nj),1,c,s)
                 t = d(pj)*c**2 + d(nj)*s**2
                 d(nj) = d(pj)*s**2 + d(nj)*c**2
                 d(pj) = t
                 k2 = k2 - 1
                 i = 1
                 90 continue
                 if (k2 + i <= n) then
                    if (d(pj) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = pj
                       i = i + 1
                       go to 90
                    else
                       indxp(k2 + i - 1) = pj
                    end if
                 else
                    indxp(k2 + i - 1) = pj
                 end if
                 pj = nj
              else
                 k = k + 1
                 dlamda(k) = d(pj)
                 w(k) = z(pj)
                 indxp(k) = pj
                 pj = nj
              end if
           end if
           go to 80
           100 continue
           ! record the last eigenvalue.
           k = k + 1
           dlamda(k) = d(pj)
           w(k) = z(pj)
           indxp(k) = pj
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four uniform groups (although one or more of these groups may be
           ! empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 1,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 1
           psm(2) = 1 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           k = n - ctot(4)
           ! fill out the indxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's.
           do j = 1,n
              js = indxp(j)
              ct = coltyp(js)
              indx(psm(ct)) = js
              indxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           i = 1
           iq1 = 1
           iq2 = 1 + (ctot(1) + ctot(2))*n1
           do j = 1,ctot(1)
              js = indx(i)
              call la_qcopy(n1,q(1,js),1,q2(iq1),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
           end do
           do j = 1,ctot(2)
              js = indx(i)
              call la_qcopy(n1,q(1,js),1,q2(iq1),1)
              call la_qcopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq1 = iq1 + n1
              iq2 = iq2 + n2
           end do
           do j = 1,ctot(3)
              js = indx(i)
              call la_qcopy(n2,q(n1 + 1,js),1,q2(iq2),1)
              z(i) = d(js)
              i = i + 1
              iq2 = iq2 + n2
           end do
           iq1 = iq2
           do j = 1,ctot(4)
              js = indx(i)
              call la_qcopy(n,q(1,js),1,q2(iq2),1)
              iq2 = iq2 + n
              z(i) = d(js)
              i = i + 1
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_qlacpy('A',n,ctot(4),q2(iq1),n,q(1,k + 1),ldq)
              call la_qcopy(n - k,z(k + 1),1,d(k + 1),1)
           end if
           ! copy ctot into coltyp for referencing in la_qlaed3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           190 continue
           return
     end subroutine la_qlaed2
#endif

     !> SLAED1: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix.  This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and eigenvectors of a tridiagonal matrix.  SLAED7 handles
     !> the case in which eigenvalues only or eigenvalues and eigenvectors
     !> of a full symmetric matrix (which was reduced to tridiagonal form)
     !> are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**T*u, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine SLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine SLAED4 (as called by SLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_slaed1(n,d,q,ldq,indxq,rho,cutpnt,work,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,n
           integer(ilp),intent(out) :: info
           real(sp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: indxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),q(ldq,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,cpp1,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,n1, &
                     n2
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (ldq < max(1,n)) then
              info = -4
           else if (min(1,n/2) > cutpnt .or. (n/2) < cutpnt) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SLAED1',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are integer pointers which indicate
           ! the portion of the workspace
           ! used by a particular array in la_slaed2 and la_slaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           call la_scopy(cutpnt,q(cutpnt,1),ldq,work(iz),1)
           cpp1 = cutpnt + 1
           call la_scopy(n - cutpnt,q(cpp1,cpp1),ldq,work(iz + cutpnt),1)
           ! deflate eigenvalues.
           call la_slaed2(k,n,cutpnt,d,q,ldq,indxq,rho,work(iz),work(idlmda), &
           work(iw),work(iq2),iwork(indx),iwork(indxc),iwork(indxp),iwork(coltyp), &
                     info)
           if (info /= 0) go to 20
           ! solve secular equation.
           if (k /= 0) then
              is = (iwork(coltyp) + iwork(coltyp + 1))*cutpnt + (iwork(coltyp + 1) + iwork( &
                        coltyp + 2))*(n - cutpnt) + iq2
              call la_slaed3(k,n,cutpnt,d,q,ldq,rho,work(idlmda),work(iq2),iwork( &
                         indxc),iwork(coltyp),work(iw),work(is),info)
              if (info /= 0) go to 20
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_slamrg(n1,n2,d,1,-1,indxq)
           else
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           20 continue
           return
     end subroutine la_slaed1
     !> DLAED1: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix.  This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and eigenvectors of a tridiagonal matrix.  DLAED7 handles
     !> the case in which eigenvalues only or eigenvalues and eigenvectors
     !> of a full symmetric matrix (which was reduced to tridiagonal form)
     !> are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**T*u, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine DLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine DLAED4 (as called by DLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_dlaed1(n,d,q,ldq,indxq,rho,cutpnt,work,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,n
           integer(ilp),intent(out) :: info
           real(dp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: indxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),q(ldq,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,n1,n2, &
                     zpp1
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (ldq < max(1,n)) then
              info = -4
           else if (min(1,n/2) > cutpnt .or. (n/2) < cutpnt) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DLAED1',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are integer pointers which indicate
           ! the portion of the workspace
           ! used by a particular array in la_dlaed2 and la_dlaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           call la_dcopy(cutpnt,q(cutpnt,1),ldq,work(iz),1)
           zpp1 = cutpnt + 1
           call la_dcopy(n - cutpnt,q(zpp1,zpp1),ldq,work(iz + cutpnt),1)
           ! deflate eigenvalues.
           call la_dlaed2(k,n,cutpnt,d,q,ldq,indxq,rho,work(iz),work(idlmda), &
           work(iw),work(iq2),iwork(indx),iwork(indxc),iwork(indxp),iwork(coltyp), &
                     info)
           if (info /= 0) go to 20
           ! solve secular equation.
           if (k /= 0) then
              is = (iwork(coltyp) + iwork(coltyp + 1))*cutpnt + (iwork(coltyp + 1) + iwork( &
                        coltyp + 2))*(n - cutpnt) + iq2
              call la_dlaed3(k,n,cutpnt,d,q,ldq,rho,work(idlmda),work(iq2),iwork( &
                         indxc),iwork(coltyp),work(iw),work(is),info)
              if (info /= 0) go to 20
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_dlamrg(n1,n2,d,1,-1,indxq)
           else
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           20 continue
           return
     end subroutine la_dlaed1
#ifdef LA_WITH_XDP
     !> XLAED1: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix.  This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and eigenvectors of a tridiagonal matrix.  XLAED7 handles
     !> the case in which eigenvalues only or eigenvalues and eigenvectors
     !> of a full symmetric matrix (which was reduced to tridiagonal form)
     !> are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**T*u, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine XLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine XLAED4 (as called by XLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_xlaed1(n,d,q,ldq,indxq,rho,cutpnt,work,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,n
           integer(ilp),intent(out) :: info
           real(xdp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: indxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: d(*),q(ldq,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,n1,n2, &
                     ypp1
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (ldq < max(1,n)) then
              info = -4
           else if (min(1,n/2) > cutpnt .or. (n/2) < cutpnt) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('XLAED1',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are integer pointers which indicate
           ! the portion of the workspace
           ! used by a particular array in la_xlaed2 and la_xlaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           call la_xcopy(cutpnt,q(cutpnt,1),ldq,work(iz),1)
           ypp1 = cutpnt + 1
           call la_xcopy(n - cutpnt,q(ypp1,ypp1),ldq,work(iz + cutpnt),1)
           ! deflate eigenvalues.
           call la_xlaed2(k,n,cutpnt,d,q,ldq,indxq,rho,work(iz),work(idlmda), &
           work(iw),work(iq2),iwork(indx),iwork(indxc),iwork(indxp),iwork(coltyp), &
                     info)
           if (info /= 0) go to 20
           ! solve secular equation.
           if (k /= 0) then
              is = (iwork(coltyp) + iwork(coltyp + 1))*cutpnt + (iwork(coltyp + 1) + iwork( &
                        coltyp + 2))*(n - cutpnt) + iq2
              call la_xlaed3(k,n,cutpnt,d,q,ldq,rho,work(idlmda),work(iq2),iwork( &
                         indxc),iwork(coltyp),work(iw),work(is),info)
              if (info /= 0) go to 20
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_xlamrg(n1,n2,d,1,-1,indxq)
           else
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           20 continue
           return
     end subroutine la_xlaed1
#endif
#ifdef LA_WITH_QP
     !> QLAED1: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix.  This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and eigenvectors of a tridiagonal matrix.  QLAED7 handles
     !> the case in which eigenvalues only or eigenvalues and eigenvectors
     !> of a full symmetric matrix (which was reduced to tridiagonal form)
     !> are desired.
     !> T = Q(in) ( D(in) + RHO * Z*Z**T ) Q**T(in) = Q(out) * D(out) * Q**T(out)
     !> where Z = Q**T*u, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine QLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine QLAED4 (as called by QLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_qlaed1(n,d,q,ldq,indxq,rho,cutpnt,work,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,n
           integer(ilp),intent(out) :: info
           real(qp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: indxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),q(ldq,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,i,idlmda,indx,indxc,indxp,iq2,is,iw,iz,k,n1,n2, &
                     wpp1
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if (ldq < max(1,n)) then
              info = -4
           else if (min(1,n/2) > cutpnt .or. (n/2) < cutpnt) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QLAED1',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are integer pointers which indicate
           ! the portion of the workspace
           ! used by a particular array in la_qlaed2 and la_qlaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq2 = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           call la_qcopy(cutpnt,q(cutpnt,1),ldq,work(iz),1)
           wpp1 = cutpnt + 1
           call la_qcopy(n - cutpnt,q(wpp1,wpp1),ldq,work(iz + cutpnt),1)
           ! deflate eigenvalues.
           call la_qlaed2(k,n,cutpnt,d,q,ldq,indxq,rho,work(iz),work(idlmda), &
           work(iw),work(iq2),iwork(indx),iwork(indxc),iwork(indxp),iwork(coltyp), &
                     info)
           if (info /= 0) go to 20
           ! solve secular equation.
           if (k /= 0) then
              is = (iwork(coltyp) + iwork(coltyp + 1))*cutpnt + (iwork(coltyp + 1) + iwork( &
                        coltyp + 2))*(n - cutpnt) + iq2
              call la_qlaed3(k,n,cutpnt,d,q,ldq,rho,work(idlmda),work(iq2),iwork( &
                         indxc),iwork(coltyp),work(iw),work(is),info)
              if (info /= 0) go to 20
           ! prepare the indxq sorting permutation.
              n1 = k
              n2 = n - k
              call la_qlamrg(n1,n2,d,1,-1,indxq)
           else
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           20 continue
           return
     end subroutine la_qlaed1
#endif

     !> SLAED0: computes all eigenvalues and corresponding eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.

     pure subroutine la_slaed0(icompq,qsiz,n,d,e,q,ldq,qstore,ldqs,work,iwork,info &
               )
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldq,ldqs,n,qsiz
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),e(*),q(ldq,*)
           real(sp),intent(out) :: qstore(ldqs,*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat,subpbs, &
                     tlvls
           real(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,real
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 2) then
              info = -1
           else if ((icompq == 1) .and. (qsiz < max(0,n))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -7
           else if (ldqs < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'SLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           if (icompq /= 2) then
              ! set up workspaces for eigenvalues only/accumulate new vectors
              ! routine
              temp = log(real(n,KIND=sp))/log(two)
              lgn = int(temp,KIND=ilp)
              if (2**lgn < n) lgn = lgn + 1
              if (2**lgn < n) lgn = lgn + 1
              iprmpt = indxq + n + 1
              iperm = iprmpt + n*lgn
              iqptr = iperm + n*lgn
              igivpt = iqptr + n + 2
              igivcl = igivpt + n*lgn
              igivnm = 1
              iq = igivnm + 2*n*lgn
              iwrem = iq + n**2 + 1
              ! initialize pointers
              do i = 0,subpbs
                 iwork(iprmpt + i) = 1
                 iwork(igivpt + i) = 1
              end do
              iwork(iqptr) = 1
           end if
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           loop_70: do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              if (icompq == 2) then
                 call la_ssteqr('I',matsiz,d(submat),e(submat),q(submat,submat), &
                           ldq,work,info)
                 if (info /= 0) go to 130
              else
                 call la_ssteqr('I',matsiz,d(submat),e(submat),work(iq - 1 + iwork( &
                           iqptr + curr)),matsiz,work,info)
                 if (info /= 0) go to 130
                 if (icompq == 1) then
                    call la_sgemm('N','N',qsiz,matsiz,matsiz,one,q(1,submat),ldq, &
                    work(iq - 1 + iwork(iqptr + curr)),matsiz,zero,qstore(1,submat),ldqs)

                 end if
                 iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
                 curr = curr + 1
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do loop_70
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              loop_90: do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.
           ! la_slaed1 is used only for the full eigensystem of a tridiagonal
           ! matrix.
           ! la_slaed7 handles the cases in which eigenvalues only or eigenvalues
           ! and eigenvectors of a full symmetric matrix (which was reduced to
           ! tridiagonal form) are desired.
                 if (icompq == 2) then
                    call la_slaed1(matsiz,d(submat),q(submat,submat),ldq,iwork( &
                    indxq + submat),e(submat + msd2 - 1),msd2,work,iwork(subpbs + 1),info)

                 else
                    call la_slaed7(icompq,matsiz,qsiz,tlvls,curlvl,curprb,d(submat), &
                    qstore(1,submat),ldqs,iwork(indxq + submat),e(submat + msd2 - 1),msd2, &
                    work(iq),iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt), &
                    iwork(igivcl),work(igivnm),work(iwrem),iwork(subpbs + 1),info)

                 end if
                 if (info /= 0) go to 130
                 iwork(i/2 + 1) = iwork(i + 2)
              end do loop_90
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           if (icompq == 1) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_scopy(qsiz,qstore(1,j),1,q(1,i),1)
              end do
              call la_scopy(n,work,1,d,1)
           else if (icompq == 2) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_scopy(n,q(1,j),1,work(n*i + 1),1)
              end do
              call la_scopy(n,work,1,d,1)
              call la_slacpy('A',n,n,work(n + 1),n,q,ldq)
           else
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
              end do
              call la_scopy(n,work,1,d,1)
           end if
           go to 140
           130 continue
           info = submat*(n + 1) + submat + matsiz - 1
           140 continue
           return
     end subroutine la_slaed0
     !> DLAED0: computes all eigenvalues and corresponding eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.

     pure subroutine la_dlaed0(icompq,qsiz,n,d,e,q,ldq,qstore,ldqs,work,iwork,info &
               )
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldq,ldqs,n,qsiz
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),e(*),q(ldq,*)
           real(dp),intent(out) :: qstore(ldqs,*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat,subpbs, &
                     tlvls
           real(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 2) then
              info = -1
           else if ((icompq == 1) .and. (qsiz < max(0,n))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -7
           else if (ldqs < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'DLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           if (icompq /= 2) then
              ! set up workspaces for eigenvalues only/accumulate new vectors
              ! routine
              temp = log(real(n,KIND=dp))/log(two)
              lgn = int(temp,KIND=ilp)
              if (2**lgn < n) lgn = lgn + 1
              if (2**lgn < n) lgn = lgn + 1
              iprmpt = indxq + n + 1
              iperm = iprmpt + n*lgn
              iqptr = iperm + n*lgn
              igivpt = iqptr + n + 2
              igivcl = igivpt + n*lgn
              igivnm = 1
              iq = igivnm + 2*n*lgn
              iwrem = iq + n**2 + 1
              ! initialize pointers
              do i = 0,subpbs
                 iwork(iprmpt + i) = 1
                 iwork(igivpt + i) = 1
              end do
              iwork(iqptr) = 1
           end if
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           loop_70: do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              if (icompq == 2) then
                 call la_dsteqr('I',matsiz,d(submat),e(submat),q(submat,submat), &
                           ldq,work,info)
                 if (info /= 0) go to 130
              else
                 call la_dsteqr('I',matsiz,d(submat),e(submat),work(iq - 1 + iwork( &
                           iqptr + curr)),matsiz,work,info)
                 if (info /= 0) go to 130
                 if (icompq == 1) then
                    call la_dgemm('N','N',qsiz,matsiz,matsiz,one,q(1,submat),ldq, &
                    work(iq - 1 + iwork(iqptr + curr)),matsiz,zero,qstore(1,submat),ldqs)

                 end if
                 iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
                 curr = curr + 1
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do loop_70
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              loop_90: do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.
           ! la_dlaed1 is used only for the full eigensystem of a tridiagonal
           ! matrix.
           ! la_dlaed7 handles the cases in which eigenvalues only or eigenvalues
           ! and eigenvectors of a full symmetric matrix (which was reduced to
           ! tridiagonal form) are desired.
                 if (icompq == 2) then
                    call la_dlaed1(matsiz,d(submat),q(submat,submat),ldq,iwork( &
                    indxq + submat),e(submat + msd2 - 1),msd2,work,iwork(subpbs + 1),info)

                 else
                    call la_dlaed7(icompq,matsiz,qsiz,tlvls,curlvl,curprb,d(submat), &
                    qstore(1,submat),ldqs,iwork(indxq + submat),e(submat + msd2 - 1),msd2, &
                    work(iq),iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt), &
                    iwork(igivcl),work(igivnm),work(iwrem),iwork(subpbs + 1),info)

                 end if
                 if (info /= 0) go to 130
                 iwork(i/2 + 1) = iwork(i + 2)
              end do loop_90
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           if (icompq == 1) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_dcopy(qsiz,qstore(1,j),1,q(1,i),1)
              end do
              call la_dcopy(n,work,1,d,1)
           else if (icompq == 2) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_dcopy(n,q(1,j),1,work(n*i + 1),1)
              end do
              call la_dcopy(n,work,1,d,1)
              call la_dlacpy('A',n,n,work(n + 1),n,q,ldq)
           else
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
              end do
              call la_dcopy(n,work,1,d,1)
           end if
           go to 140
           130 continue
           info = submat*(n + 1) + submat + matsiz - 1
           140 continue
           return
     end subroutine la_dlaed0
#ifdef LA_WITH_XDP
     !> XLAED0: computes all eigenvalues and corresponding eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.

     pure subroutine la_xlaed0(icompq,qsiz,n,d,e,q,ldq,qstore,ldqs,work,iwork,info &
               )
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldq,ldqs,n,qsiz
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: d(*),e(*),q(ldq,*)
           real(xdp),intent(out) :: qstore(ldqs,*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat,subpbs, &
                     tlvls
           real(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 2) then
              info = -1
           else if ((icompq == 1) .and. (qsiz < max(0,n))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -7
           else if (ldqs < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'XLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           if (icompq /= 2) then
              ! set up workspaces for eigenvalues only/accumulate new vectors
              ! routine
              temp = log(real(n,KIND=xdp))/log(two)
              lgn = int(temp,KIND=ilp)
              if (2**lgn < n) lgn = lgn + 1
              if (2**lgn < n) lgn = lgn + 1
              iprmpt = indxq + n + 1
              iperm = iprmpt + n*lgn
              iqptr = iperm + n*lgn
              igivpt = iqptr + n + 2
              igivcl = igivpt + n*lgn
              igivnm = 1
              iq = igivnm + 2*n*lgn
              iwrem = iq + n**2 + 1
              ! initialize pointers
              do i = 0,subpbs
                 iwork(iprmpt + i) = 1
                 iwork(igivpt + i) = 1
              end do
              iwork(iqptr) = 1
           end if
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           loop_70: do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              if (icompq == 2) then
                 call la_xsteqr('I',matsiz,d(submat),e(submat),q(submat,submat), &
                           ldq,work,info)
                 if (info /= 0) go to 130
              else
                 call la_xsteqr('I',matsiz,d(submat),e(submat),work(iq - 1 + iwork( &
                           iqptr + curr)),matsiz,work,info)
                 if (info /= 0) go to 130
                 if (icompq == 1) then
                    call la_xgemm('N','N',qsiz,matsiz,matsiz,one,q(1,submat),ldq, &
                    work(iq - 1 + iwork(iqptr + curr)),matsiz,zero,qstore(1,submat),ldqs)

                 end if
                 iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
                 curr = curr + 1
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do loop_70
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              loop_90: do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.
           ! la_xlaed1 is used only for the full eigensystem of a tridiagonal
           ! matrix.
           ! la_xlaed7 handles the cases in which eigenvalues only or eigenvalues
           ! and eigenvectors of a full symmetric matrix (which was reduced to
           ! tridiagonal form) are desired.
                 if (icompq == 2) then
                    call la_xlaed1(matsiz,d(submat),q(submat,submat),ldq,iwork( &
                    indxq + submat),e(submat + msd2 - 1),msd2,work,iwork(subpbs + 1),info)

                 else
                    call la_xlaed7(icompq,matsiz,qsiz,tlvls,curlvl,curprb,d(submat), &
                    qstore(1,submat),ldqs,iwork(indxq + submat),e(submat + msd2 - 1),msd2, &
                    work(iq),iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt), &
                    iwork(igivcl),work(igivnm),work(iwrem),iwork(subpbs + 1),info)

                 end if
                 if (info /= 0) go to 130
                 iwork(i/2 + 1) = iwork(i + 2)
              end do loop_90
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           if (icompq == 1) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_xcopy(qsiz,qstore(1,j),1,q(1,i),1)
              end do
              call la_xcopy(n,work,1,d,1)
           else if (icompq == 2) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_xcopy(n,q(1,j),1,work(n*i + 1),1)
              end do
              call la_xcopy(n,work,1,d,1)
              call la_xlacpy('A',n,n,work(n + 1),n,q,ldq)
           else
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
              end do
              call la_xcopy(n,work,1,d,1)
           end if
           go to 140
           130 continue
           info = submat*(n + 1) + submat + matsiz - 1
           140 continue
           return
     end subroutine la_xlaed0
#endif
#ifdef LA_WITH_QP
     !> QLAED0: computes all eigenvalues and corresponding eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.

     pure subroutine la_qlaed0(icompq,qsiz,n,d,e,q,ldq,qstore,ldqs,work,iwork,info &
               )
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldq,ldqs,n,qsiz
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),e(*),q(ldq,*)
           real(qp),intent(out) :: qstore(ldqs,*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat,subpbs, &
                     tlvls
           real(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (icompq < 0 .or. icompq > 2) then
              info = -1
           else if ((icompq == 1) .and. (qsiz < max(0,n))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -7
           else if (ldqs < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'QLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           if (icompq /= 2) then
              ! set up workspaces for eigenvalues only/accumulate new vectors
              ! routine
              temp = log(real(n,KIND=qp))/log(two)
              lgn = int(temp,KIND=ilp)
              if (2**lgn < n) lgn = lgn + 1
              if (2**lgn < n) lgn = lgn + 1
              iprmpt = indxq + n + 1
              iperm = iprmpt + n*lgn
              iqptr = iperm + n*lgn
              igivpt = iqptr + n + 2
              igivcl = igivpt + n*lgn
              igivnm = 1
              iq = igivnm + 2*n*lgn
              iwrem = iq + n**2 + 1
              ! initialize pointers
              do i = 0,subpbs
                 iwork(iprmpt + i) = 1
                 iwork(igivpt + i) = 1
              end do
              iwork(iqptr) = 1
           end if
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           loop_70: do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              if (icompq == 2) then
                 call la_qsteqr('I',matsiz,d(submat),e(submat),q(submat,submat), &
                           ldq,work,info)
                 if (info /= 0) go to 130
              else
                 call la_qsteqr('I',matsiz,d(submat),e(submat),work(iq - 1 + iwork( &
                           iqptr + curr)),matsiz,work,info)
                 if (info /= 0) go to 130
                 if (icompq == 1) then
                    call la_qgemm('N','N',qsiz,matsiz,matsiz,one,q(1,submat),ldq, &
                    work(iq - 1 + iwork(iqptr + curr)),matsiz,zero,qstore(1,submat),ldqs)

                 end if
                 iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
                 curr = curr + 1
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do loop_70
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              loop_90: do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.
           ! la_qlaed1 is used only for the full eigensystem of a tridiagonal
           ! matrix.
           ! la_qlaed7 handles the cases in which eigenvalues only or eigenvalues
           ! and eigenvectors of a full symmetric matrix (which was reduced to
           ! tridiagonal form) are desired.
                 if (icompq == 2) then
                    call la_qlaed1(matsiz,d(submat),q(submat,submat),ldq,iwork( &
                    indxq + submat),e(submat + msd2 - 1),msd2,work,iwork(subpbs + 1),info)

                 else
                    call la_qlaed7(icompq,matsiz,qsiz,tlvls,curlvl,curprb,d(submat), &
                    qstore(1,submat),ldqs,iwork(indxq + submat),e(submat + msd2 - 1),msd2, &
                    work(iq),iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt), &
                    iwork(igivcl),work(igivnm),work(iwrem),iwork(subpbs + 1),info)

                 end if
                 if (info /= 0) go to 130
                 iwork(i/2 + 1) = iwork(i + 2)
              end do loop_90
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           if (icompq == 1) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_qcopy(qsiz,qstore(1,j),1,q(1,i),1)
              end do
              call la_qcopy(n,work,1,d,1)
           else if (icompq == 2) then
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
                 call la_qcopy(n,q(1,j),1,work(n*i + 1),1)
              end do
              call la_qcopy(n,work,1,d,1)
              call la_qlacpy('A',n,n,work(n + 1),n,q,ldq)
           else
              do i = 1,n
                 j = iwork(indxq + i)
                 work(i) = d(j)
              end do
              call la_qcopy(n,work,1,d,1)
           end if
           go to 140
           130 continue
           info = submat*(n + 1) + submat + matsiz - 1
           140 continue
           return
     end subroutine la_qlaed0
#endif

     !> CLAED8: merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_claed8(k,n,qsiz,q,ldq,d,rho,cutpnt,z,dlamda,q2,ldq2,w, &
               indxp,indx,indxq,perm,givptr,givcol,givnum,info)
        use la_constants_sp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(sp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(sp),intent(inout) :: d(*),z(*)
           real(sp),intent(out) :: dlamda(*),givnum(2,*),w(*)
           complex(sp),intent(inout) :: q(ldq,*)
           complex(sp),intent(out) :: q2(ldq2,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: mone = -1.0_sp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(sp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -5
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -8
           else if (ldq2 < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('CLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_sscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_sscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_slamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_isamax(n,z,1)
           jmax = la_isamax(n,d,1)
           eps = la_slamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! -- except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              do j = 1,n
                 perm(j) = indxq(indx(j))
                 call la_ccopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
              call la_clacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 100
              else
                 jlam = j
                 go to 70
              end if
           end do
           70 continue
           j = j + 1
           if (j > n) go to 90
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_slapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 call la_csrot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j) &
                           )),1,c,s)
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 80 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 80
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 70
           90 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           100 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           do j = 1,n
              jp = indxp(j)
              dlamda(j) = d(jp)
              perm(j) = indxq(indx(jp))
              call la_ccopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_scopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              call la_clacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
           end if
           return
     end subroutine la_claed8
     !> ZLAED8 merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_zlaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,z,dlamda,q2,ldq2,w, &
               indxp,indx,indxq,perm,givptr,givcol,givnum,info)
        use la_constants_dp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(dp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(dp),intent(inout) :: d(*),z(*)
           real(dp),intent(out) :: dlamda(*),givnum(2,*),w(*)
           complex(dp),intent(inout) :: q(ldq,*)
           complex(dp),intent(out) :: q2(ldq2,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: mone = -1.0_dp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(dp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -5
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -8
           else if (ldq2 < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('ZLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_dscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_dscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_dlamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_idamax(n,z,1)
           jmax = la_idamax(n,d,1)
           eps = la_dlamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! -- except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              do j = 1,n
                 perm(j) = indxq(indx(j))
                 call la_zcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
              call la_zlacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 100
              else
                 jlam = j
                 go to 70
              end if
           end do
           70 continue
           j = j + 1
           if (j > n) go to 90
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_dlapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 call la_zdrot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j) &
                           )),1,c,s)
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 80 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 80
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 70
           90 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           100 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           do j = 1,n
              jp = indxp(j)
              dlamda(j) = d(jp)
              perm(j) = indxq(indx(jp))
              call la_zcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_dcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              call la_zlacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
           end if
           return
     end subroutine la_zlaed8
#ifdef LA_WITH_XDP
     !> YLAED8 merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_ylaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,z,dlamda,q2,ldq2,w, &
               indxp,indx,indxq,perm,givptr,givcol,givnum,info)
        use la_constants_xdp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(xdp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(xdp),intent(inout) :: d(*),z(*)
           real(xdp),intent(out) :: dlamda(*),givnum(2,*),w(*)
           complex(xdp),intent(inout) :: q(ldq,*)
           complex(xdp),intent(out) :: q2(ldq2,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: mone = -1.0_xdp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(xdp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -5
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -8
           else if (ldq2 < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('YLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_xscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_xscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_xlamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_ixamax(n,z,1)
           jmax = la_ixamax(n,d,1)
           eps = la_xlamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! -- except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              do j = 1,n
                 perm(j) = indxq(indx(j))
                 call la_ycopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
              call la_ylacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 100
              else
                 jlam = j
                 go to 70
              end if
           end do
           70 continue
           j = j + 1
           if (j > n) go to 90
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_xlapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 call la_yxrot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j) &
                           )),1,c,s)
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 80 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 80
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 70
           90 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           100 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           do j = 1,n
              jp = indxp(j)
              dlamda(j) = d(jp)
              perm(j) = indxq(indx(jp))
              call la_ycopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_xcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              call la_ylacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
           end if
           return
     end subroutine la_ylaed8
#endif
#ifdef LA_WITH_QP
     !> WLAED8 merges the two sets of eigenvalues together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> eigenvalues are close together or if there is a tiny element in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.

     pure subroutine la_wlaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,z,dlamda,q2,ldq2,w, &
               indxp,indx,indxq,perm,givptr,givcol,givnum,info)
        use la_constants_qp,only:zero,one,two,eight
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: cutpnt,ldq,ldq2,n,qsiz
           integer(ilp),intent(out) :: givptr,info,k
           real(qp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(2,*),indx(*),indxp(*),perm(*)
           integer(ilp),intent(inout) :: indxq(*)
           real(qp),intent(inout) :: d(*),z(*)
           real(qp),intent(out) :: dlamda(*),givnum(2,*),w(*)
           complex(qp),intent(inout) :: q(ldq,*)
           complex(qp),intent(out) :: q2(ldq2,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: mone = -1.0_qp

           ! Local Scalars
           integer(ilp) :: i,imax,j,jlam,jmax,jp,k2,n1,n1p1,n2
           real(qp) :: c,eps,s,t,tau,tol
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -5
           else if (cutpnt < min(1,n) .or. cutpnt > n) then
              info = -8
           else if (ldq2 < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('WLAED8',-info)
              return
           end if
           ! need to initialize givptr to o here in case of quick exit
           ! to prevent an unspecified code behavior (usually sigfault)
           ! when iwork array on entry to *stedc is not zeroed
           ! (or at least some iwork entries which used in *laed7 for givptr).
           givptr = 0
           ! quick return if possible
           if (n == 0) return
           n1 = cutpnt
           n2 = n - n1
           n1p1 = n1 + 1
           if (rho < zero) then
              call la_qscal(n2,mone,z(n1p1),1)
           end if
           ! normalize z so that norm(z) = 1
           t = one/sqrt(two)
           do j = 1,n
              indx(j) = j
           end do
           call la_qscal(n,t,z,1)
           rho = abs(two*rho)
           ! sort the eigenvalues into increasing order
           do i = cutpnt + 1,n
              indxq(i) = indxq(i) + cutpnt
           end do
           do i = 1,n
              dlamda(i) = d(indxq(i))
              w(i) = z(indxq(i))
           end do
           i = 1
           j = cutpnt + 1
           call la_qlamrg(n1,n2,dlamda,1,1,indx)
           do i = 1,n
              d(i) = dlamda(indx(i))
              z(i) = w(indx(i))
           end do
           ! calculate the allowable deflation tolerance
           imax = la_iqamax(n,z,1)
           jmax = la_iqamax(n,d,1)
           eps = la_qlamch('EPSILON')
           tol = eight*eps*abs(d(jmax))
           ! if the rank-1 modifier is small enough, no more needs to be done
           ! -- except to reorganize q so that its columns correspond with the
           ! elements in d.
           if (rho*abs(z(imax)) <= tol) then
              k = 0
              do j = 1,n
                 perm(j) = indxq(indx(j))
                 call la_wcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
              end do
              call la_wlacpy('A',qsiz,n,q2(1,1),ldq2,q(1,1),ldq)
              return
           end if
           ! if there are multiple eigenvalues then the problem deflates.  here
           ! the number of equal eigenvalues are found.  as each equal
           ! eigenvalue is found, an elementary reflector is computed to rotate
           ! the corresponding eigensubspace so that the corresponding
           ! components of z are zero in this new basis.
           k = 0
           k2 = n + 1
           do j = 1,n
              if (rho*abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 indxp(k2) = j
                 if (j == n) go to 100
              else
                 jlam = j
                 go to 70
              end if
           end do
           70 continue
           j = j + 1
           if (j > n) go to 90
           if (rho*abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              indxp(k2) = j
           else
              ! check if eigenvalues are close enough to allow deflation.
              s = z(jlam)
              c = z(j)
              ! find sqrt(a**2+b**2) without overflow or
              ! destructive underflow.
              tau = la_qlapy2(c,s)
              t = d(j) - d(jlam)
              c = c/tau
              s = -s/tau
              if (abs(t*c*s) <= tol) then
                 ! deflation is possible.
                 z(j) = tau
                 z(jlam) = zero
                 ! record the appropriate givens rotation
                 givptr = givptr + 1
                 givcol(1,givptr) = indxq(indx(jlam))
                 givcol(2,givptr) = indxq(indx(j))
                 givnum(1,givptr) = c
                 givnum(2,givptr) = s
                 call la_wqrot(qsiz,q(1,indxq(indx(jlam))),1,q(1,indxq(indx(j) &
                           )),1,c,s)
                 t = d(jlam)*c*c + d(j)*s*s
                 d(j) = d(jlam)*s*s + d(j)*c*c
                 d(jlam) = t
                 k2 = k2 - 1
                 i = 1
                 80 continue
                 if (k2 + i <= n) then
                    if (d(jlam) < d(indxp(k2 + i))) then
                       indxp(k2 + i - 1) = indxp(k2 + i)
                       indxp(k2 + i) = jlam
                       i = i + 1
                       go to 80
                    else
                       indxp(k2 + i - 1) = jlam
                    end if
                 else
                    indxp(k2 + i - 1) = jlam
                 end if
                 jlam = j
              else
                 k = k + 1
                 w(k) = z(jlam)
                 dlamda(k) = d(jlam)
                 indxp(k) = jlam
                 jlam = j
              end if
           end if
           go to 70
           90 continue
           ! record the last eigenvalue.
           k = k + 1
           w(k) = z(jlam)
           dlamda(k) = d(jlam)
           indxp(k) = jlam
           100 continue
           ! sort the eigenvalues and corresponding eigenvectors into dlamda
           ! and q2 respectively.  the eigenvalues/vectors which were not
           ! deflated go into the first k slots of dlamda and q2 respectively,
           ! while those which were deflated go into the last n - k slots.
           do j = 1,n
              jp = indxp(j)
              dlamda(j) = d(jp)
              perm(j) = indxq(indx(jp))
              call la_wcopy(qsiz,q(1,perm(j)),1,q2(1,j),1)
           end do
           ! the deflated eigenvalues and their corresponding vectors go back
           ! into the last n - k slots of d and q respectively.
           if (k < n) then
              call la_qcopy(n - k,dlamda(k + 1),1,d(k + 1),1)
              call la_wlacpy('A',qsiz,n - k,q2(1,k + 1),ldq2,q(1,k + 1),ldq)
           end if
           return
     end subroutine la_wlaed8
#endif

     !> CSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
     !> matrix to tridiagonal form.

     pure subroutine la_csteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_sp,only:zero,one,two,three,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(sp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = cone
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_slamch('E')
           eps2 = eps**2
           safmin = la_slamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_claset('FULL',n,n,czero,cone,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_slanst('I',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_slascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_slascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_slascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_slascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_slae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_slaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_clasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_slae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_slapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_slartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_clasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_slae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_slaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_clasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_slae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_slapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_slartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_clasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_slascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_slascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_slascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_slascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot == nmaxit) then
              do i = 1,n - 1
                 if (e(i) /= zero) info = info + 1
              end do
              return
           end if
           go to 10
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_slasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_cswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           return
     end subroutine la_csteqr
     !> ZSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
     !> matrix to tridiagonal form.

     pure subroutine la_zsteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_dp,only:zero,one,two,three,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(dp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = cone
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_dlamch('E')
           eps2 = eps**2
           safmin = la_dlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_zlaset('FULL',n,n,czero,cone,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_dlanst('I',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_dlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_dlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_dlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_dlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_dlae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_dlaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_zlasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_dlae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_dlapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_dlartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_zlasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_dlae2 or la_slaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_dlaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_zlasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_dlae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_dlapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_dlartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_zlasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_dlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_dlascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_dlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_dlascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot == nmaxit) then
              do i = 1,n - 1
                 if (e(i) /= zero) info = info + 1
              end do
              return
           end if
           go to 10
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_dlasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_zswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           return
     end subroutine la_zsteqr
#ifdef LA_WITH_XDP
     !> YSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if YHETRD or YHPTRD or YHBTRD has been used to reduce this
     !> matrix to tridiagonal form.

     pure subroutine la_ysteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_xdp,only:zero,one,two,three,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(xdp),intent(inout) :: d(*),e(*)
           real(xdp),intent(out) :: work(*)
           complex(xdp),intent(inout) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(xdp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('YSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = cone
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_xlamch('E')
           eps2 = eps**2
           safmin = la_xlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_ylaset('FULL',n,n,czero,cone,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_xlanst('I',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_xlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_xlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_xlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_xlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_xlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_xlaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_ylasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_xlae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_xlapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_xlartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_ylasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_xlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_xlaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_ylasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_xlae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_xlapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_xlartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_ylasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_xlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_xlascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_xlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_xlascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot == nmaxit) then
              do i = 1,n - 1
                 if (e(i) /= zero) info = info + 1
              end do
              return
           end if
           go to 10
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_xlasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_yswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           return
     end subroutine la_ysteqr
#endif
#ifdef LA_WITH_QP
     !> WSTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the implicit QL or QR method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if WHETRD or WHPTRD or WHBTRD has been used to reduce this
     !> matrix to tridiagonal form.

     pure subroutine la_wsteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_qp,only:zero,one,two,three,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,icompz,ii,iscale,j,jtot,k,l,l1,lend,lendm1,lendp1,lendsv, &
                      lm1,lsv,m,mm,mm1,nm1,nmaxit
           real(qp) :: anorm,b,c,eps,eps2,f,g,p,r,rt1,rt2,s,safmax,safmin,ssfmax, &
                     ssfmin,tst
           ! Intrinsic Functions
           intrinsic :: abs,max,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (la_lsame(compz,'N')) then
              icompz = 0
           else if (la_lsame(compz,'V')) then
              icompz = 1
           else if (la_lsame(compz,'I')) then
              icompz = 2
           else
              icompz = -1
           end if
           if (icompz < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if ((ldz < 1) .or. (icompz > 0 .and. ldz < max(1,n))) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WSTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz == 2) z(1,1) = cone
              return
           end if
           ! determine the unit roundoff and over/underflow thresholds.
           eps = la_qlamch('E')
           eps2 = eps**2
           safmin = la_qlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues and eigenvectors of the tridiagonal
           ! matrix.
           if (icompz == 2) call la_wlaset('FULL',n,n,czero,cone,z,ldz)
           nmaxit = n*maxit
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           nm1 = n - 1
           10 continue
           if (l1 > n) go to 160
           if (l1 > 1) e(l1 - 1) = zero
           if (l1 <= nm1) then
              do m = l1,nm1
                 tst = abs(e(m))
                 if (tst == zero) go to 30
                 if (tst <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) then
                    e(m) = zero
                    go to 30
                 end if
              end do
           end if
           m = n
           30 continue
           l = l1
           lsv = l
           lend = m
           lendsv = lend
           l1 = m + 1
           if (lend == l) go to 10
           ! scale submatrix in rows and columns l to lend
           anorm = la_qlanst('I',lend - l + 1,d(l),e(l))
           iscale = 0
           if (anorm == zero) go to 10
           if (anorm > ssfmax) then
              iscale = 1
              call la_qlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_qlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_qlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_qlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend > l) then
              ! ql iteration
              ! look for small subdiagonal element.
              40 continue
              if (l /= lend) then
                 lendm1 = lend - 1
                 do m = l,lendm1
                    tst = abs(e(m))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m + 1)) + safmin) go to 60
                 end do
              end if
              m = lend
              60 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 80
              ! if remaining matrix is 2-by-2, use la_qlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l + 1) then
                 if (icompz > 0) then
                    call la_qlaev2(d(l),e(l),d(l + 1),rt1,rt2,c,s)
                    work(l) = c
                    work(n - 1 + l) = s
                    call la_wlasr('R','V','B',n,2,work(l),work(n - 1 + l),z(1,l), &
                              ldz)
                 else
                    call la_qlae2(d(l),e(l),d(l + 1),rt1,rt2)
                 end if
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 40
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l + 1) - p)/(two*e(l))
              r = la_qlapy2(g,one)
              g = d(m) - p + (e(l)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              mm1 = m - 1
              do i = mm1,l,-1
                 f = s*e(i)
                 b = c*e(i)
                 call la_qlartg(g,f,c,s,r)
                 if (i /= m - 1) e(i + 1) = r
                 g = d(i + 1) - p
                 r = (d(i) - g)*s + two*c*b
                 p = s*r
                 d(i + 1) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = -s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = m - l + 1
                 call la_wlasr('R','V','B',n,mm,work(l),work(n - 1 + l),z(1,l),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(l) = g
              go to 40
              ! eigenvalue found.
              80 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 40
              go to 140
           else
              ! qr iteration
              ! look for small superdiagonal element.
              90 continue
              if (l /= lend) then
                 lendp1 = lend + 1
                 do m = l,lendp1,-1
                    tst = abs(e(m - 1))**2
                    if (tst <= (eps2*abs(d(m)))*abs(d(m - 1)) + safmin) go to 110
                 end do
              end if
              m = lend
              110 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 130
              ! if remaining matrix is 2-by-2, use la_qlae2 or la_dlaev2
              ! to compute its eigensystem.
              if (m == l - 1) then
                 if (icompz > 0) then
                    call la_qlaev2(d(l - 1),e(l - 1),d(l),rt1,rt2,c,s)
                    work(m) = c
                    work(n - 1 + m) = s
                    call la_wlasr('R','V','F',n,2,work(m),work(n - 1 + m),z(1,l - 1), &
                              ldz)
                 else
                    call la_qlae2(d(l - 1),e(l - 1),d(l),rt1,rt2)
                 end if
                 d(l - 1) = rt1
                 d(l) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 90
                 go to 140
              end if
              if (jtot == nmaxit) go to 140
              jtot = jtot + 1
              ! form shift.
              g = (d(l - 1) - p)/(two*e(l - 1))
              r = la_qlapy2(g,one)
              g = d(m) - p + (e(l - 1)/(g + sign(r,g)))
              s = one
              c = one
              p = zero
              ! inner loop
              lm1 = l - 1
              do i = m,lm1
                 f = s*e(i)
                 b = c*e(i)
                 call la_qlartg(g,f,c,s,r)
                 if (i /= m) e(i - 1) = r
                 g = d(i) - p
                 r = (d(i + 1) - g)*s + two*c*b
                 p = s*r
                 d(i) = g + p
                 g = c*r - b
                 ! if eigenvectors are desired, then save rotations.
                 if (icompz > 0) then
                    work(i) = c
                    work(n - 1 + i) = s
                 end if
              end do
              ! if eigenvectors are desired, then apply saved rotations.
              if (icompz > 0) then
                 mm = l - m + 1
                 call la_wlasr('R','V','F',n,mm,work(m),work(n - 1 + m),z(1,m),ldz &
                           )
              end if
              d(l) = d(l) - p
              e(lm1) = g
              go to 90
              ! eigenvalue found.
              130 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 90
              go to 140
           end if
           ! undo scaling if necessary
           140 continue
           if (iscale == 1) then
              call la_qlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_qlascl('G',0,0,ssfmax,anorm,lendsv - lsv,1,e(lsv),n,info)

           else if (iscale == 2) then
              call la_qlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv),n,info)

              call la_qlascl('G',0,0,ssfmin,anorm,lendsv - lsv,1,e(lsv),n,info)

           end if
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot == nmaxit) then
              do i = 1,n - 1
                 if (e(i) /= zero) info = info + 1
              end do
              return
           end if
           go to 10
           ! order eigenvalues and eigenvectors.
           160 continue
           if (icompz == 0) then
              ! use quick sort
              call la_qlasrt('I',n,d,info)
           else
              ! use selection sort to minimize swaps of eigenvectors
              do ii = 2,n
                 i = ii - 1
                 k = i
                 p = d(i)
                 do j = ii,n
                    if (d(j) < p) then
                       k = j
                       p = d(j)
                    end if
                 end do
                 if (k /= i) then
                    d(k) = d(i)
                    d(i) = p
                    call la_wswap(n,z(1,i),1,z(1,k),1)
                 end if
              end do
           end if
           return
     end subroutine la_wsteqr
#endif

     !> CLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense or banded
     !> Hermitian matrix that has been reduced to tridiagonal form.
     !> T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)
     !> where Z = Q**Hu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine SLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine SLAED4 (as called by SLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_claed7(n,cutpnt,qsiz,tlvls,curlvl,curpbm,d,q,ldq,rho,indxq, &
               qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,rwork,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(sp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(sp),intent(inout) :: d(*),givnum(2,*),qstore(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: q(ldq,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq,iw,iz,k,n1,n2, &
                     ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq<0 .or. icompq>1 ) then
              ! info = -1
           ! else if( n<0 ) then
           if (n < 0) then
              info = -1
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('CLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_slaed2 and la_slaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_slaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,rwork(iz),rwork(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_claed8(k,n,qsiz,q,ldq,d,rho,cutpnt,rwork(iz),rwork(idlmda), &
           work,qsiz,rwork(iw),iwork(indxp),iwork(indx),indxq,perm(prmptr(curr)), &
           givptr(curr + 1),givcol(1,givptr(curr)),givnum(1,givptr(curr)),info)

           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_slaed9(k,1,k,n,d,rwork(iq),k,rho,rwork(idlmda),rwork(iw), &
                        qstore(qptr(curr)),k,info)
              call la_clacrm(qsiz,k,work,qsiz,qstore(qptr(curr)),k,q,ldq,rwork( &
                        iq))
              qptr(curr + 1) = qptr(curr) + k**2
              if (info /= 0) then
                 return
              end if
           ! prepare the indxq sorting premutation.
              n1 = k
              n2 = n - k
              call la_slamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           return
     end subroutine la_claed7
     !> ZLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense or banded
     !> Hermitian matrix that has been reduced to tridiagonal form.
     !> T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)
     !> where Z = Q**Hu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine DLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine DLAED4 (as called by SLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_zlaed7(n,cutpnt,qsiz,tlvls,curlvl,curpbm,d,q,ldq,rho,indxq, &
               qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,rwork,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(dp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(dp),intent(inout) :: d(*),givnum(2,*),qstore(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: q(ldq,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq,iw,iz,k,n1,n2, &
                     ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq<0 .or. icompq>1 ) then
              ! info = -1
           ! else if( n<0 ) then
           if (n < 0) then
              info = -1
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_dlaed2 and la_slaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_dlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,rwork(iz),rwork(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_zlaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,rwork(iz),rwork(idlmda), &
           work,qsiz,rwork(iw),iwork(indxp),iwork(indx),indxq,perm(prmptr(curr)), &
           givptr(curr + 1),givcol(1,givptr(curr)),givnum(1,givptr(curr)),info)

           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_dlaed9(k,1,k,n,d,rwork(iq),k,rho,rwork(idlmda),rwork(iw), &
                        qstore(qptr(curr)),k,info)
              call la_zlacrm(qsiz,k,work,qsiz,qstore(qptr(curr)),k,q,ldq,rwork( &
                        iq))
              qptr(curr + 1) = qptr(curr) + k**2
              if (info /= 0) then
                 return
              end if
           ! prepare the indxq sorting premutation.
              n1 = k
              n2 = n - k
              call la_dlamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           return
     end subroutine la_zlaed7
#ifdef LA_WITH_XDP
     !> YLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense or banded
     !> Hermitian matrix that has been reduced to tridiagonal form.
     !> T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)
     !> where Z = Q**Hu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine XLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine XLAED4 (as called by DLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_ylaed7(n,cutpnt,qsiz,tlvls,curlvl,curpbm,d,q,ldq,rho,indxq, &
               qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,rwork,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(xdp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(xdp),intent(inout) :: d(*),givnum(2,*),qstore(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(inout) :: q(ldq,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq,iw,iz,k,n1,n2, &
                     ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq<0 .or. icompq>1 ) then
              ! info = -1
           ! else if( n<0 ) then
           if (n < 0) then
              info = -1
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('YLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_xlaed2 and la_dlaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_xlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,rwork(iz),rwork(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_ylaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,rwork(iz),rwork(idlmda), &
           work,qsiz,rwork(iw),iwork(indxp),iwork(indx),indxq,perm(prmptr(curr)), &
           givptr(curr + 1),givcol(1,givptr(curr)),givnum(1,givptr(curr)),info)

           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_xlaed9(k,1,k,n,d,rwork(iq),k,rho,rwork(idlmda),rwork(iw), &
                        qstore(qptr(curr)),k,info)
              call la_ylacrm(qsiz,k,work,qsiz,qstore(qptr(curr)),k,q,ldq,rwork( &
                        iq))
              qptr(curr + 1) = qptr(curr) + k**2
              if (info /= 0) then
                 return
              end if
           ! prepare the indxq sorting premutation.
              n1 = k
              n2 = n - k
              call la_xlamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           return
     end subroutine la_ylaed7
#endif
#ifdef LA_WITH_QP
     !> WLAED7: computes the updated eigensystem of a diagonal
     !> matrix after modification by a rank-one symmetric matrix. This
     !> routine is used only for the eigenproblem which requires all
     !> eigenvalues and optionally eigenvectors of a dense or banded
     !> Hermitian matrix that has been reduced to tridiagonal form.
     !> T = Q(in) ( D(in) + RHO * Z*Z**H ) Q**H(in) = Q(out) * D(out) * Q**H(out)
     !> where Z = Q**Hu, u is a vector of length N with ones in the
     !> CUTPNT and CUTPNT + 1 th elements and zeros elsewhere.
     !> The eigenvectors of the original matrix are stored in Q, and the
     !> eigenvalues are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple eigenvalues or if there is a zero in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine QLAED2.
     !> The second stage consists of calculating the updated
     !> eigenvalues. This is done by finding the roots of the secular
     !> equation via the routine QLAED4 (as called by DLAED3).
     !> This routine also calculates the eigenvectors of the current
     !> problem.
     !> The final stage consists of computing the updated eigenvectors
     !> directly using the updated eigenvalues.  The eigenvectors for
     !> the current problem are multiplied with the eigenvectors from
     !> the overall problem.

     pure subroutine la_wlaed7(n,cutpnt,qsiz,tlvls,curlvl,curpbm,d,q,ldq,rho,indxq, &
               qstore,qptr,prmptr,perm,givptr,givcol,givnum,work,rwork,iwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,cutpnt,ldq,n,qsiz,tlvls
           integer(ilp),intent(out) :: info
           real(qp),intent(inout) :: rho
           ! Array Arguments
           integer(ilp),intent(inout) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)

           integer(ilp),intent(out) :: indxq(*),iwork(*)
           real(qp),intent(inout) :: d(*),givnum(2,*),qstore(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: q(ldq,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: coltyp,curr,i,idlmda,indx,indxc,indxp,iq,iw,iz,k,n1,n2, &
                     ptr
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq<0 .or. icompq>1 ) then
              ! info = -1
           ! else if( n<0 ) then
           if (n < 0) then
              info = -1
           else if (min(1,n) > cutpnt .or. n < cutpnt) then
              info = -2
           else if (qsiz < n) then
              info = -3
           else if (ldq < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WLAED7',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_qlaed2 and la_dlaed3.
           iz = 1
           idlmda = iz + n
           iw = idlmda + n
           iq = iw + n
           indx = 1
           indxc = indx + n
           coltyp = indxc + n
           indxp = coltyp + n
           ! form the z-vector which consists of the last row of q_1 and the
           ! first row of q_2.
           ptr = 1 + 2**tlvls
           do i = 1,curlvl - 1
              ptr = ptr + 2**(tlvls - i)
           end do
           curr = ptr + curpbm
           call la_qlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                     qstore,qptr,rwork(iz),rwork(iz + n),info)
           ! when solving the final problem, we no longer need the stored data,
           ! so we will overwrite the data from this level onto the previously
           ! used storage space.
           if (curlvl == tlvls) then
              qptr(curr) = 1
              prmptr(curr) = 1
              givptr(curr) = 1
           end if
           ! sort and deflate eigenvalues.
           call la_wlaed8(k,n,qsiz,q,ldq,d,rho,cutpnt,rwork(iz),rwork(idlmda), &
           work,qsiz,rwork(iw),iwork(indxp),iwork(indx),indxq,perm(prmptr(curr)), &
           givptr(curr + 1),givcol(1,givptr(curr)),givnum(1,givptr(curr)),info)

           prmptr(curr + 1) = prmptr(curr) + n
           givptr(curr + 1) = givptr(curr + 1) + givptr(curr)
           ! solve secular equation.
           if (k /= 0) then
              call la_qlaed9(k,1,k,n,d,rwork(iq),k,rho,rwork(idlmda),rwork(iw), &
                        qstore(qptr(curr)),k,info)
              call la_wlacrm(qsiz,k,work,qsiz,qstore(qptr(curr)),k,q,ldq,rwork( &
                        iq))
              qptr(curr + 1) = qptr(curr) + k**2
              if (info /= 0) then
                 return
              end if
           ! prepare the indxq sorting premutation.
              n1 = k
              n2 = n - k
              call la_qlamrg(n1,n2,d,1,-1,indxq)
           else
              qptr(curr + 1) = qptr(curr)
              do i = 1,n
                 indxq(i) = i
              end do
           end if
           return
     end subroutine la_wlaed7
#endif

     !> Using the divide and conquer method, CLAED0: computes all eigenvalues
     !> of a symmetric tridiagonal matrix which is one diagonal block of
     !> those from reducing a dense or band Hermitian matrix and
     !> corresponding eigenvectors of the dense or band matrix.

     pure subroutine la_claed0(qsiz,n,d,e,q,ldq,qstore,ldqs,rwork,iwork,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldqs,n,qsiz
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: q(ldq,*)
           complex(sp),intent(out) :: qstore(ldqs,*)
        ! =====================================================================
        ! warning:      n could be as big as qsiz!

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,ll,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat, &
                     subpbs,tlvls
           real(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,real
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq < 0 .or. icompq > 2 ) then
              ! info = -1
           ! else if( ( icompq == 1 ) .and. ( qsiz < max( 0, n ) ) )
          ! $        then
           if (qsiz < max(0,n)) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (ldqs < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'CLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           ! set up workspaces for eigenvalues only/accumulate new vectors
           ! routine
           temp = log(real(n,KIND=sp))/log(two)
           lgn = int(temp,KIND=ilp)
           if (2**lgn < n) lgn = lgn + 1
           if (2**lgn < n) lgn = lgn + 1
           iprmpt = indxq + n + 1
           iperm = iprmpt + n*lgn
           iqptr = iperm + n*lgn
           igivpt = iqptr + n + 2
           igivcl = igivpt + n*lgn
           igivnm = 1
           iq = igivnm + 2*n*lgn
           iwrem = iq + n**2 + 1
           ! initialize pointers
           do i = 0,subpbs
              iwork(iprmpt + i) = 1
              iwork(igivpt + i) = 1
           end do
           iwork(iqptr) = 1
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              ll = iq - 1 + iwork(iqptr + curr)
              call la_ssteqr('I',matsiz,d(submat),e(submat),rwork(ll),matsiz, &
                        rwork,info)
              call la_clacrm(qsiz,matsiz,q(1,submat),ldq,rwork(ll),matsiz,qstore( &
                        1,submat),ldqs,rwork(iwrem))
              iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
              curr = curr + 1
              if (info > 0) then
                 info = submat*(n + 1) + submat + matsiz - 1
                 return
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.  la_claed7 handles the case
           ! when the eigenvectors of a full or band hermitian matrix (which
           ! was reduced to tridiagonal form) are desired.
           ! i am free to use q as a valuable working space until loop 150.
                 call la_claed7(matsiz,msd2,qsiz,tlvls,curlvl,curprb,d(submat), &
                 qstore(1,submat),ldqs,e(submat + msd2 - 1),iwork(indxq + submat),rwork(iq), &
                 iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt),iwork(igivcl), &
                           rwork(igivnm),q(1,submat),rwork(iwrem),iwork(subpbs + 1),info)
                 if (info > 0) then
                    info = submat*(n + 1) + submat + matsiz - 1
                    return
                 end if
                 iwork(i/2 + 1) = iwork(i + 2)
              end do
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           do i = 1,n
              j = iwork(indxq + i)
              rwork(i) = d(j)
              call la_ccopy(qsiz,qstore(1,j),1,q(1,i),1)
           end do
           call la_scopy(n,rwork,1,d,1)
           return
     end subroutine la_claed0
     !> Using the divide and conquer method, ZLAED0: computes all eigenvalues
     !> of a symmetric tridiagonal matrix which is one diagonal block of
     !> those from reducing a dense or band Hermitian matrix and
     !> corresponding eigenvectors of the dense or band matrix.

     pure subroutine la_zlaed0(qsiz,n,d,e,q,ldq,qstore,ldqs,rwork,iwork,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldqs,n,qsiz
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: q(ldq,*)
           complex(dp),intent(out) :: qstore(ldqs,*)
        ! =====================================================================
        ! warning:      n could be as big as qsiz!

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,ll,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat, &
                     subpbs,tlvls
           real(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq < 0 .or. icompq > 2 ) then
              ! info = -1
           ! else if( ( icompq == 1 ) .and. ( qsiz < max( 0, n ) ) )
          ! $        then
           if (qsiz < max(0,n)) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (ldqs < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'ZLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           ! set up workspaces for eigenvalues only/accumulate new vectors
           ! routine
           temp = log(real(n,KIND=dp))/log(two)
           lgn = int(temp,KIND=ilp)
           if (2**lgn < n) lgn = lgn + 1
           if (2**lgn < n) lgn = lgn + 1
           iprmpt = indxq + n + 1
           iperm = iprmpt + n*lgn
           iqptr = iperm + n*lgn
           igivpt = iqptr + n + 2
           igivcl = igivpt + n*lgn
           igivnm = 1
           iq = igivnm + 2*n*lgn
           iwrem = iq + n**2 + 1
           ! initialize pointers
           do i = 0,subpbs
              iwork(iprmpt + i) = 1
              iwork(igivpt + i) = 1
           end do
           iwork(iqptr) = 1
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              ll = iq - 1 + iwork(iqptr + curr)
              call la_dsteqr('I',matsiz,d(submat),e(submat),rwork(ll),matsiz, &
                        rwork,info)
              call la_zlacrm(qsiz,matsiz,q(1,submat),ldq,rwork(ll),matsiz,qstore( &
                        1,submat),ldqs,rwork(iwrem))
              iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
              curr = curr + 1
              if (info > 0) then
                 info = submat*(n + 1) + submat + matsiz - 1
                 return
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.  la_zlaed7 handles the case
           ! when the eigenvectors of a full or band hermitian matrix (which
           ! was reduced to tridiagonal form) are desired.
           ! i am free to use q as a valuable working space until loop 150.
                 call la_zlaed7(matsiz,msd2,qsiz,tlvls,curlvl,curprb,d(submat), &
                 qstore(1,submat),ldqs,e(submat + msd2 - 1),iwork(indxq + submat),rwork(iq), &
                 iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt),iwork(igivcl), &
                           rwork(igivnm),q(1,submat),rwork(iwrem),iwork(subpbs + 1),info)
                 if (info > 0) then
                    info = submat*(n + 1) + submat + matsiz - 1
                    return
                 end if
                 iwork(i/2 + 1) = iwork(i + 2)
              end do
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           do i = 1,n
              j = iwork(indxq + i)
              rwork(i) = d(j)
              call la_zcopy(qsiz,qstore(1,j),1,q(1,i),1)
           end do
           call la_dcopy(n,rwork,1,d,1)
           return
     end subroutine la_zlaed0
#ifdef LA_WITH_XDP
     !> Using the divide and conquer method, YLAED0: computes all eigenvalues
     !> of a symmetric tridiagonal matrix which is one diagonal block of
     !> those from reducing a dense or band Hermitian matrix and
     !> corresponding eigenvectors of the dense or band matrix.

     pure subroutine la_ylaed0(qsiz,n,d,e,q,ldq,qstore,ldqs,rwork,iwork,info)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldqs,n,qsiz
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: d(*),e(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(inout) :: q(ldq,*)
           complex(xdp),intent(out) :: qstore(ldqs,*)
        ! =====================================================================
        ! warning:      n could be as big as qsiz!

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,ll,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat, &
                     subpbs,tlvls
           real(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq < 0 .or. icompq > 2 ) then
              ! info = -1
           ! else if( ( icompq == 1 ) .and. ( qsiz < max( 0, n ) ) )
          ! $        then
           if (qsiz < max(0,n)) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (ldqs < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('YLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'YLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           ! set up workspaces for eigenvalues only/accumulate new vectors
           ! routine
           temp = log(real(n,KIND=xdp))/log(two)
           lgn = int(temp,KIND=ilp)
           if (2**lgn < n) lgn = lgn + 1
           if (2**lgn < n) lgn = lgn + 1
           iprmpt = indxq + n + 1
           iperm = iprmpt + n*lgn
           iqptr = iperm + n*lgn
           igivpt = iqptr + n + 2
           igivcl = igivpt + n*lgn
           igivnm = 1
           iq = igivnm + 2*n*lgn
           iwrem = iq + n**2 + 1
           ! initialize pointers
           do i = 0,subpbs
              iwork(iprmpt + i) = 1
              iwork(igivpt + i) = 1
           end do
           iwork(iqptr) = 1
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              ll = iq - 1 + iwork(iqptr + curr)
              call la_xsteqr('I',matsiz,d(submat),e(submat),rwork(ll),matsiz, &
                        rwork,info)
              call la_ylacrm(qsiz,matsiz,q(1,submat),ldq,rwork(ll),matsiz,qstore( &
                        1,submat),ldqs,rwork(iwrem))
              iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
              curr = curr + 1
              if (info > 0) then
                 info = submat*(n + 1) + submat + matsiz - 1
                 return
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.  la_ylaed7 handles the case
           ! when the eigenvectors of a full or band hermitian matrix (which
           ! was reduced to tridiagonal form) are desired.
           ! i am free to use q as a valuable working space until loop 150.
                 call la_ylaed7(matsiz,msd2,qsiz,tlvls,curlvl,curprb,d(submat), &
                 qstore(1,submat),ldqs,e(submat + msd2 - 1),iwork(indxq + submat),rwork(iq), &
                 iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt),iwork(igivcl), &
                           rwork(igivnm),q(1,submat),rwork(iwrem),iwork(subpbs + 1),info)
                 if (info > 0) then
                    info = submat*(n + 1) + submat + matsiz - 1
                    return
                 end if
                 iwork(i/2 + 1) = iwork(i + 2)
              end do
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           do i = 1,n
              j = iwork(indxq + i)
              rwork(i) = d(j)
              call la_ycopy(qsiz,qstore(1,j),1,q(1,i),1)
           end do
           call la_xcopy(n,rwork,1,d,1)
           return
     end subroutine la_ylaed0
#endif
#ifdef LA_WITH_QP
     !> Using the divide and conquer method, WLAED0: computes all eigenvalues
     !> of a symmetric tridiagonal matrix which is one diagonal block of
     !> those from reducing a dense or band Hermitian matrix and
     !> corresponding eigenvectors of the dense or band matrix.

     pure subroutine la_wlaed0(qsiz,n,d,e,q,ldq,qstore,ldqs,rwork,iwork,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldq,ldqs,n,qsiz
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: q(ldq,*)
           complex(qp),intent(out) :: qstore(ldqs,*)
        ! =====================================================================
        ! warning:      n could be as big as qsiz!

           ! Local Scalars
           integer(ilp) :: curlvl,curprb,curr,i,igivcl,igivnm,igivpt,indxq,iperm,iprmpt, &
           iq,iqptr,iwrem,j,k,lgn,ll,matsiz,msd2,smlsiz,smm1,spm1,spm2,submat, &
                     subpbs,tlvls
           real(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! if( icompq < 0 .or. icompq > 2 ) then
              ! info = -1
           ! else if( ( icompq == 1 ) .and. ( qsiz < max( 0, n ) ) )
          ! $        then
           if (qsiz < max(0,n)) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldq < max(1,n)) then
              info = -6
           else if (ldqs < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WLAED0',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'WLAED0',' ',0,0,0,0)
           ! determine the size and placement of the submatrices, and save in
           ! the leading elements of iwork.
           iwork(1) = n
           subpbs = 1
           tlvls = 0
           10 continue
           if (iwork(subpbs) > smlsiz) then
              do j = subpbs,1,-1
                 iwork(2*j) = (iwork(j) + 1)/2
                 iwork(2*j - 1) = iwork(j)/2
              end do
              tlvls = tlvls + 1
              subpbs = 2*subpbs
              go to 10
           end if
           do j = 2,subpbs
              iwork(j) = iwork(j) + iwork(j - 1)
           end do
           ! divide the matrix into subpbs submatrices of size at most smlsiz+1
           ! using rank-1 modifications (cuts).
           spm1 = subpbs - 1
           do i = 1,spm1
              submat = iwork(i) + 1
              smm1 = submat - 1
              d(smm1) = d(smm1) - abs(e(smm1))
              d(submat) = d(submat) - abs(e(smm1))
           end do
           indxq = 4*n + 3
           ! set up workspaces for eigenvalues only/accumulate new vectors
           ! routine
           temp = log(real(n,KIND=qp))/log(two)
           lgn = int(temp,KIND=ilp)
           if (2**lgn < n) lgn = lgn + 1
           if (2**lgn < n) lgn = lgn + 1
           iprmpt = indxq + n + 1
           iperm = iprmpt + n*lgn
           iqptr = iperm + n*lgn
           igivpt = iqptr + n + 2
           igivcl = igivpt + n*lgn
           igivnm = 1
           iq = igivnm + 2*n*lgn
           iwrem = iq + n**2 + 1
           ! initialize pointers
           do i = 0,subpbs
              iwork(iprmpt + i) = 1
              iwork(igivpt + i) = 1
           end do
           iwork(iqptr) = 1
           ! solve each submatrix eigenproblem at the bottom of the divide and
           ! conquer tree.
           curr = 0
           do i = 0,spm1
              if (i == 0) then
                 submat = 1
                 matsiz = iwork(1)
              else
                 submat = iwork(i) + 1
                 matsiz = iwork(i + 1) - iwork(i)
              end if
              ll = iq - 1 + iwork(iqptr + curr)
              call la_qsteqr('I',matsiz,d(submat),e(submat),rwork(ll),matsiz, &
                        rwork,info)
              call la_wlacrm(qsiz,matsiz,q(1,submat),ldq,rwork(ll),matsiz,qstore( &
                        1,submat),ldqs,rwork(iwrem))
              iwork(iqptr + curr + 1) = iwork(iqptr + curr) + matsiz**2
              curr = curr + 1
              if (info > 0) then
                 info = submat*(n + 1) + submat + matsiz - 1
                 return
              end if
              k = 1
              do j = submat,iwork(i + 1)
                 iwork(indxq + j) = k
                 k = k + 1
              end do
           end do
           ! successively merge eigensystems of adjacent submatrices
           ! into eigensystem for the corresponding larger matrix.
           ! while ( subpbs > 1 )
           curlvl = 1
           80 continue
           if (subpbs > 1) then
              spm2 = subpbs - 2
              do i = 0,spm2,2
                 if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                    curprb = 0
                 else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                    curprb = curprb + 1
                 end if
           ! merge lower order eigensystems (of size msd2 and matsiz - msd2)
           ! into an eigensystem of size matsiz.  la_wlaed7 handles the case
           ! when the eigenvectors of a full or band hermitian matrix (which
           ! was reduced to tridiagonal form) are desired.
           ! i am free to use q as a valuable working space until loop 150.
                 call la_wlaed7(matsiz,msd2,qsiz,tlvls,curlvl,curprb,d(submat), &
                 qstore(1,submat),ldqs,e(submat + msd2 - 1),iwork(indxq + submat),rwork(iq), &
                 iwork(iqptr),iwork(iprmpt),iwork(iperm),iwork(igivpt),iwork(igivcl), &
                           rwork(igivnm),q(1,submat),rwork(iwrem),iwork(subpbs + 1),info)
                 if (info > 0) then
                    info = submat*(n + 1) + submat + matsiz - 1
                    return
                 end if
                 iwork(i/2 + 1) = iwork(i + 2)
              end do
              subpbs = subpbs/2
              curlvl = curlvl + 1
              go to 80
           end if
           ! end while
           ! re-merge the eigenvalues/vectors which were deflated at the final
           ! merge step.
           do i = 1,n
              j = iwork(indxq + i)
              rwork(i) = d(j)
              call la_wcopy(qsiz,qstore(1,j),1,q(1,i),1)
           end do
           call la_qcopy(n,rwork,1,d,1)
           return
     end subroutine la_wlaed0
#endif

end module la_lapack_eigv_tridiag
