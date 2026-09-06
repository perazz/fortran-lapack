!> Bidiagonal singular values by divide and conquer, with its secular-equation and merge kernels
module la_lapack_eigv_svd_bidiag_dc
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level3_gen
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l2
     use la_lapack_blas_like_mnorm
     use la_lapack_blas_like_scalar
     use la_lapack_eigv_tridiag
     use la_lapack_eigv_tridiag2
     use la_lapack_givens_jacobi_rot
     use la_lapack_svd_bidiag_qr
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slasd5
     public :: la_slasdt
     public :: la_slasd4
     public :: la_slasd7
     public :: la_slasd8
     public :: la_slasd3
     public :: la_slasd6
     public :: la_slasd2
     public :: la_slasd1
     public :: la_sbdsdc
     public :: la_slasd0
     public :: la_slasda
     public :: la_slasdq
     public :: la_dlasd5
     public :: la_dlasdt
     public :: la_dlasd4
     public :: la_dlasd7
     public :: la_dlasd8
     public :: la_dlasd3
     public :: la_dlasd6
     public :: la_dlasd2
     public :: la_dlasd1
     public :: la_dbdsdc
     public :: la_dlasd0
     public :: la_dlasda
     public :: la_dlasdq
#ifdef LA_WITH_XDP
     public :: la_xlasd5
     public :: la_xlasdt
     public :: la_xlasd4
     public :: la_xlasd7
     public :: la_xlasd8
     public :: la_xlasd3
     public :: la_xlasd6
     public :: la_xlasd2
     public :: la_xlasd1
     public :: la_xbdsdc
     public :: la_xlasd0
     public :: la_xlasda
     public :: la_xlasdq
#endif
#ifdef LA_WITH_QP
     public :: la_qlasd5
     public :: la_qlasdt
     public :: la_qlasd4
     public :: la_qlasd7
     public :: la_qlasd8
     public :: la_qlasd3
     public :: la_qlasd6
     public :: la_qlasd2
     public :: la_qlasd1
     public :: la_qbdsdc
     public :: la_qlasd0
     public :: la_qlasda
     public :: la_qlasdq
#endif

     contains

     !> This subroutine computes the square root of the I-th eigenvalue
     !> of a positive symmetric rank-one modification of a 2-by-2 diagonal
     !> matrix
     !> diag( D ) * diag( D ) +  RHO * Z * transpose(Z) .
     !> The diagonal entries in the array D are assumed to satisfy
     !> 0 <= D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_slasd5(i,d,z,delta,rho,dsigma,work)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(sp),intent(out) :: dsigma
           real(sp),intent(in) :: rho
           ! Array Arguments
           real(sp),intent(in) :: d(2),z(2)
           real(sp),intent(out) :: delta(2),work(2)
        ! =====================================================================

           ! Local Scalars
           real(sp) :: b,c,del,delsq,tau,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           delsq = del*(d(2) + d(1))
           if (i == 1) then
              w = one + four*rho*(z(2)*z(2)/(d(1) + three*d(2)) - z(1)*z(1)/( &
                        three*d(1) + d(2)))/del
              if (w > zero) then
                 b = delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*delsq
                 ! b > zero, always
                 ! the following tau is dsigma * dsigma - d( 1 ) * d( 1 )
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 ! the following tau is dsigma - d( 1 )
                 tau = tau/(d(1) + sqrt(d(1)*d(1) + tau))
                 dsigma = d(1) + tau
                 delta(1) = -tau
                 delta(2) = del - tau
                 work(1) = two*d(1) + tau
                 work(2) = (d(1) + tau) + d(2)
                 ! delta( 1 ) = -z( 1 ) / tau
                 ! delta( 2 ) = z( 2 ) / ( del-tau )
              else
                 b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*delsq
                 ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 ! the following tau is dsigma - d( 2 )
                 tau = tau/(d(2) + sqrt(abs(d(2)*d(2) + tau)))
                 dsigma = d(2) + tau
                 delta(1) = -(del + tau)
                 delta(2) = -tau
                 work(1) = d(1) + tau + d(2)
                 work(2) = two*d(2) + tau
                 ! delta( 1 ) = -z( 1 ) / ( del+tau )
                 ! delta( 2 ) = -z( 2 ) / tau
              end if
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           else
              ! now i=2
              b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*delsq
              ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              ! the following tau is dsigma - d( 2 )
              tau = tau/(d(2) + sqrt(d(2)*d(2) + tau))
              dsigma = d(2) + tau
              delta(1) = -(del + tau)
              delta(2) = -tau
              work(1) = d(1) + tau + d(2)
              work(2) = two*d(2) + tau
              ! delta( 1 ) = -z( 1 ) / ( del+tau )
              ! delta( 2 ) = -z( 2 ) / tau
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           end if
           return
     end subroutine la_slasd5
     !> This subroutine computes the square root of the I-th eigenvalue
     !> of a positive symmetric rank-one modification of a 2-by-2 diagonal
     !> matrix
     !> diag( D ) * diag( D ) +  RHO * Z * transpose(Z) .
     !> The diagonal entries in the array D are assumed to satisfy
     !> 0 <= D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_dlasd5(i,d,z,delta,rho,dsigma,work)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(dp),intent(out) :: dsigma
           real(dp),intent(in) :: rho
           ! Array Arguments
           real(dp),intent(in) :: d(2),z(2)
           real(dp),intent(out) :: delta(2),work(2)
        ! =====================================================================

           ! Local Scalars
           real(dp) :: b,c,del,delsq,tau,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           delsq = del*(d(2) + d(1))
           if (i == 1) then
              w = one + four*rho*(z(2)*z(2)/(d(1) + three*d(2)) - z(1)*z(1)/( &
                        three*d(1) + d(2)))/del
              if (w > zero) then
                 b = delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*delsq
                 ! b > zero, always
                 ! the following tau is dsigma * dsigma - d( 1 ) * d( 1 )
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 ! the following tau is dsigma - d( 1 )
                 tau = tau/(d(1) + sqrt(d(1)*d(1) + tau))
                 dsigma = d(1) + tau
                 delta(1) = -tau
                 delta(2) = del - tau
                 work(1) = two*d(1) + tau
                 work(2) = (d(1) + tau) + d(2)
                 ! delta( 1 ) = -z( 1 ) / tau
                 ! delta( 2 ) = z( 2 ) / ( del-tau )
              else
                 b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*delsq
                 ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 ! the following tau is dsigma - d( 2 )
                 tau = tau/(d(2) + sqrt(abs(d(2)*d(2) + tau)))
                 dsigma = d(2) + tau
                 delta(1) = -(del + tau)
                 delta(2) = -tau
                 work(1) = d(1) + tau + d(2)
                 work(2) = two*d(2) + tau
                 ! delta( 1 ) = -z( 1 ) / ( del+tau )
                 ! delta( 2 ) = -z( 2 ) / tau
              end if
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           else
              ! now i=2
              b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*delsq
              ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              ! the following tau is dsigma - d( 2 )
              tau = tau/(d(2) + sqrt(d(2)*d(2) + tau))
              dsigma = d(2) + tau
              delta(1) = -(del + tau)
              delta(2) = -tau
              work(1) = d(1) + tau + d(2)
              work(2) = two*d(2) + tau
              ! delta( 1 ) = -z( 1 ) / ( del+tau )
              ! delta( 2 ) = -z( 2 ) / tau
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           end if
           return
     end subroutine la_dlasd5
#ifdef LA_WITH_XDP
     !> This subroutine computes the square root of the I-th eigenvalue
     !> of a positive symmetric rank-one modification of a 2-by-2 diagonal
     !> matrix
     !> diag( D ) * diag( D ) +  RHO * Z * transpose(Z) .
     !> The diagonal entries in the array D are assumed to satisfy
     !> 0 <= D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_xlasd5(i,d,z,delta,rho,dsigma,work)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(xdp),intent(out) :: dsigma
           real(xdp),intent(in) :: rho
           ! Array Arguments
           real(xdp),intent(in) :: d(2),z(2)
           real(xdp),intent(out) :: delta(2),work(2)
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: b,c,del,delsq,tau,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           delsq = del*(d(2) + d(1))
           if (i == 1) then
              w = one + four*rho*(z(2)*z(2)/(d(1) + three*d(2)) - z(1)*z(1)/( &
                        three*d(1) + d(2)))/del
              if (w > zero) then
                 b = delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*delsq
                 ! b > zero, always
                 ! the following tau is dsigma * dsigma - d( 1 ) * d( 1 )
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 ! the following tau is dsigma - d( 1 )
                 tau = tau/(d(1) + sqrt(d(1)*d(1) + tau))
                 dsigma = d(1) + tau
                 delta(1) = -tau
                 delta(2) = del - tau
                 work(1) = two*d(1) + tau
                 work(2) = (d(1) + tau) + d(2)
                 ! delta( 1 ) = -z( 1 ) / tau
                 ! delta( 2 ) = z( 2 ) / ( del-tau )
              else
                 b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*delsq
                 ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 ! the following tau is dsigma - d( 2 )
                 tau = tau/(d(2) + sqrt(abs(d(2)*d(2) + tau)))
                 dsigma = d(2) + tau
                 delta(1) = -(del + tau)
                 delta(2) = -tau
                 work(1) = d(1) + tau + d(2)
                 work(2) = two*d(2) + tau
                 ! delta( 1 ) = -z( 1 ) / ( del+tau )
                 ! delta( 2 ) = -z( 2 ) / tau
              end if
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           else
              ! now i=2
              b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*delsq
              ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              ! the following tau is dsigma - d( 2 )
              tau = tau/(d(2) + sqrt(d(2)*d(2) + tau))
              dsigma = d(2) + tau
              delta(1) = -(del + tau)
              delta(2) = -tau
              work(1) = d(1) + tau + d(2)
              work(2) = two*d(2) + tau
              ! delta( 1 ) = -z( 1 ) / ( del+tau )
              ! delta( 2 ) = -z( 2 ) / tau
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           end if
           return
     end subroutine la_xlasd5
#endif
#ifdef LA_WITH_QP
     !> This subroutine computes the square root of the I-th eigenvalue
     !> of a positive symmetric rank-one modification of a 2-by-2 diagonal
     !> matrix
     !> diag( D ) * diag( D ) +  RHO * Z * transpose(Z) .
     !> The diagonal entries in the array D are assumed to satisfy
     !> 0 <= D(i) < D(j)  for  i < j .
     !> We also assume RHO > 0 and that the Euclidean norm of the vector
     !> Z is one.

     pure subroutine la_qlasd5(i,d,z,delta,rho,dsigma,work)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i
           real(qp),intent(out) :: dsigma
           real(qp),intent(in) :: rho
           ! Array Arguments
           real(qp),intent(in) :: d(2),z(2)
           real(qp),intent(out) :: delta(2),work(2)
        ! =====================================================================

           ! Local Scalars
           real(qp) :: b,c,del,delsq,tau,w
           ! Intrinsic Functions
           intrinsic :: abs,sqrt
           ! Executable Statements
           del = d(2) - d(1)
           delsq = del*(d(2) + d(1))
           if (i == 1) then
              w = one + four*rho*(z(2)*z(2)/(d(1) + three*d(2)) - z(1)*z(1)/( &
                        three*d(1) + d(2)))/del
              if (w > zero) then
                 b = delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(1)*z(1)*delsq
                 ! b > zero, always
                 ! the following tau is dsigma * dsigma - d( 1 ) * d( 1 )
                 tau = two*c/(b + sqrt(abs(b*b - four*c)))
                 ! the following tau is dsigma - d( 1 )
                 tau = tau/(d(1) + sqrt(d(1)*d(1) + tau))
                 dsigma = d(1) + tau
                 delta(1) = -tau
                 delta(2) = del - tau
                 work(1) = two*d(1) + tau
                 work(2) = (d(1) + tau) + d(2)
                 ! delta( 1 ) = -z( 1 ) / tau
                 ! delta( 2 ) = z( 2 ) / ( del-tau )
              else
                 b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
                 c = rho*z(2)*z(2)*delsq
                 ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
                 if (b > zero) then
                    tau = -two*c/(b + sqrt(b*b + four*c))
                 else
                    tau = (b - sqrt(b*b + four*c))/two
                 end if
                 ! the following tau is dsigma - d( 2 )
                 tau = tau/(d(2) + sqrt(abs(d(2)*d(2) + tau)))
                 dsigma = d(2) + tau
                 delta(1) = -(del + tau)
                 delta(2) = -tau
                 work(1) = d(1) + tau + d(2)
                 work(2) = two*d(2) + tau
                 ! delta( 1 ) = -z( 1 ) / ( del+tau )
                 ! delta( 2 ) = -z( 2 ) / tau
              end if
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           else
              ! now i=2
              b = -delsq + rho*(z(1)*z(1) + z(2)*z(2))
              c = rho*z(2)*z(2)*delsq
              ! the following tau is dsigma * dsigma - d( 2 ) * d( 2 )
              if (b > zero) then
                 tau = (b + sqrt(b*b + four*c))/two
              else
                 tau = two*c/(-b + sqrt(b*b + four*c))
              end if
              ! the following tau is dsigma - d( 2 )
              tau = tau/(d(2) + sqrt(d(2)*d(2) + tau))
              dsigma = d(2) + tau
              delta(1) = -(del + tau)
              delta(2) = -tau
              work(1) = d(1) + tau + d(2)
              work(2) = two*d(2) + tau
              ! delta( 1 ) = -z( 1 ) / ( del+tau )
              ! delta( 2 ) = -z( 2 ) / tau
              ! temp = sqrt( delta( 1 )*delta( 1 )+delta( 2 )*delta( 2 ) )
              ! delta( 1 ) = delta( 1 ) / temp
              ! delta( 2 ) = delta( 2 ) / temp
           end if
           return
     end subroutine la_qlasd5
#endif

     !> SLASDT: creates a tree of subproblems for bidiagonal divide and
     !> conquer.

     pure subroutine la_slasdt(n,lvl,nd,inode,ndiml,ndimr,msub)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: lvl,nd
           integer(ilp),intent(in) :: msub,n
           ! Array Arguments
           integer(ilp),intent(out) :: inode(*),ndiml(*),ndimr(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,il,ir,llst,maxn,ncrnt,nlvl
           real(sp) :: temp
           ! Intrinsic Functions
           intrinsic :: int,log,max,real
           ! Executable Statements
           ! find the number of levels on the tree.
           maxn = max(1,n)
           temp = log(real(maxn,KIND=sp)/real(msub + 1,KIND=sp))/log(two)
           lvl = int(temp,KIND=ilp) + 1
           i = n/2
           inode(1) = i + 1
           ndiml(1) = i
           ndimr(1) = n - i - 1
           il = 0
           ir = 1
           llst = 1
           do nlvl = 1,lvl - 1
              ! constructing the tree at (nlvl+1)-st level. the number of
              ! nodes created on this level is llst * 2.
              do i = 0,llst - 1
                 il = il + 2
                 ir = ir + 2
                 ncrnt = llst + i
                 ndiml(il) = ndiml(ncrnt)/2
                 ndimr(il) = ndiml(ncrnt) - ndiml(il) - 1
                 inode(il) = inode(ncrnt) - ndimr(il) - 1
                 ndiml(ir) = ndimr(ncrnt)/2
                 ndimr(ir) = ndimr(ncrnt) - ndiml(ir) - 1
                 inode(ir) = inode(ncrnt) + ndiml(ir) + 1
              end do
              llst = llst*2
           end do
           nd = llst*2 - 1
           return
     end subroutine la_slasdt
     !> DLASDT: creates a tree of subproblems for bidiagonal divide and
     !> conquer.

     pure subroutine la_dlasdt(n,lvl,nd,inode,ndiml,ndimr,msub)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: lvl,nd
           integer(ilp),intent(in) :: msub,n
           ! Array Arguments
           integer(ilp),intent(out) :: inode(*),ndiml(*),ndimr(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,il,ir,llst,maxn,ncrnt,nlvl
           real(dp) :: temp
           ! Intrinsic Functions
           intrinsic :: real,int,log,max
           ! Executable Statements
           ! find the number of levels on the tree.
           maxn = max(1,n)
           temp = log(real(maxn,KIND=dp)/real(msub + 1,KIND=dp))/log(two)
           lvl = int(temp,KIND=ilp) + 1
           i = n/2
           inode(1) = i + 1
           ndiml(1) = i
           ndimr(1) = n - i - 1
           il = 0
           ir = 1
           llst = 1
           do nlvl = 1,lvl - 1
              ! constructing the tree at (nlvl+1)-st level. the number of
              ! nodes created on this level is llst * 2.
              do i = 0,llst - 1
                 il = il + 2
                 ir = ir + 2
                 ncrnt = llst + i
                 ndiml(il) = ndiml(ncrnt)/2
                 ndimr(il) = ndiml(ncrnt) - ndiml(il) - 1
                 inode(il) = inode(ncrnt) - ndimr(il) - 1
                 ndiml(ir) = ndimr(ncrnt)/2
                 ndimr(ir) = ndimr(ncrnt) - ndiml(ir) - 1
                 inode(ir) = inode(ncrnt) + ndiml(ir) + 1
              end do
              llst = llst*2
           end do
           nd = llst*2 - 1
           return
     end subroutine la_dlasdt
#ifdef LA_WITH_XDP
     !> XLASDT: creates a tree of subproblems for bidiagonal divide and
     !> conquer.

     pure subroutine la_xlasdt(n,lvl,nd,inode,ndiml,ndimr,msub)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: lvl,nd
           integer(ilp),intent(in) :: msub,n
           ! Array Arguments
           integer(ilp),intent(out) :: inode(*),ndiml(*),ndimr(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,il,ir,llst,maxn,ncrnt,nlvl
           real(xdp) :: temp
           ! Intrinsic Functions
           intrinsic :: real,int,log,max
           ! Executable Statements
           ! find the number of levels on the tree.
           maxn = max(1,n)
           temp = log(real(maxn,KIND=xdp)/real(msub + 1,KIND=xdp))/log(two)
           lvl = int(temp,KIND=ilp) + 1
           i = n/2
           inode(1) = i + 1
           ndiml(1) = i
           ndimr(1) = n - i - 1
           il = 0
           ir = 1
           llst = 1
           do nlvl = 1,lvl - 1
              ! constructing the tree at (nlvl+1)-st level. the number of
              ! nodes created on this level is llst * 2.
              do i = 0,llst - 1
                 il = il + 2
                 ir = ir + 2
                 ncrnt = llst + i
                 ndiml(il) = ndiml(ncrnt)/2
                 ndimr(il) = ndiml(ncrnt) - ndiml(il) - 1
                 inode(il) = inode(ncrnt) - ndimr(il) - 1
                 ndiml(ir) = ndimr(ncrnt)/2
                 ndimr(ir) = ndimr(ncrnt) - ndiml(ir) - 1
                 inode(ir) = inode(ncrnt) + ndiml(ir) + 1
              end do
              llst = llst*2
           end do
           nd = llst*2 - 1
           return
     end subroutine la_xlasdt
#endif
#ifdef LA_WITH_QP
     !> QLASDT: creates a tree of subproblems for bidiagonal divide and
     !> conquer.

     pure subroutine la_qlasdt(n,lvl,nd,inode,ndiml,ndimr,msub)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: lvl,nd
           integer(ilp),intent(in) :: msub,n
           ! Array Arguments
           integer(ilp),intent(out) :: inode(*),ndiml(*),ndimr(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,il,ir,llst,maxn,ncrnt,nlvl
           real(qp) :: temp
           ! Intrinsic Functions
           intrinsic :: real,int,log,max
           ! Executable Statements
           ! find the number of levels on the tree.
           maxn = max(1,n)
           temp = log(real(maxn,KIND=qp)/real(msub + 1,KIND=qp))/log(two)
           lvl = int(temp,KIND=ilp) + 1
           i = n/2
           inode(1) = i + 1
           ndiml(1) = i
           ndimr(1) = n - i - 1
           il = 0
           ir = 1
           llst = 1
           do nlvl = 1,lvl - 1
              ! constructing the tree at (nlvl+1)-st level. the number of
              ! nodes created on this level is llst * 2.
              do i = 0,llst - 1
                 il = il + 2
                 ir = ir + 2
                 ncrnt = llst + i
                 ndiml(il) = ndiml(ncrnt)/2
                 ndimr(il) = ndiml(ncrnt) - ndiml(il) - 1
                 inode(il) = inode(ncrnt) - ndimr(il) - 1
                 ndiml(ir) = ndimr(ncrnt)/2
                 ndimr(ir) = ndimr(ncrnt) - ndiml(ir) - 1
                 inode(ir) = inode(ncrnt) + ndiml(ir) + 1
              end do
              llst = llst*2
           end do
           nd = llst*2 - 1
           return
     end subroutine la_qlasdt
#endif

     !> This subroutine computes the square root of the I-th updated
     !> eigenvalue of a positive symmetric rank-one modification to
     !> a positive diagonal matrix whose entries are given as the squares
     !> of the corresponding entries in the array d, and that
     !> 0 <= D(i) < D(j)  for  i < j
     !> and that RHO > 0. This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D ) * diag( D ) +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_slasd4(n,i,d,z,delta,rho,sigma,work,info)
        use la_constants_sp,only:zero,one,two,three,four,eight,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(sp),intent(in) :: rho
           real(sp),intent(out) :: sigma
           ! Array Arguments
           real(sp),intent(in) :: d(*),z(*)
           real(sp),intent(out) :: delta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 400

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3,geomavg
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(sp) :: a,b,c,delsq,delsq2,sq2,dphi,dpsi,dtiim,dtiip,dtipsq,dtisq, &
           dtnsq,dtnsq1,dw,eps,erretm,eta,phi,prew,psi,rhoinv,sglb,sgub,tau,tau2, &
                     temp,temp1,temp2,w
           ! Local Arrays
           real(sp) :: dd(3),zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
              ! presumably, i=1 upon entry
              sigma = sqrt(d(1)*d(1) + rho*z(1)*z(1))
              delta(1) = one
              work(1) = one
              return
           end if
           if (n == 2) then
              call la_slasd5(i,d,z,delta,rho,sigma,work)
              return
           end if
           ! compute machine epsilon
           eps = la_slamch('EPSILON')
           rhoinv = one/rho
           tau2 = zero
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              temp = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              temp1 = temp/(d(n) + sqrt(d(n)*d(n) + temp))
              do j = 1,n
                 work(j) = d(j) + d(n) + temp1
                 delta(j) = (d(j) - d(n)) - temp1
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/(delta(j)*work(j))
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/(delta(ii)*work(ii)) + z(n)*z(n)/(delta(n) &
                        *work(n))
              if (w <= zero) then
                 temp1 = sqrt(d(n)*d(n) + rho)
                 temp = z(n - 1)*z(n - 1)/((d(n - 1) + temp1)*(d(n) - d(n - 1) + rho/(d(n) + &
                           temp1))) + z(n)*z(n)/rho
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (c <= temp) then
                    tau = rho
                 else
                    delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                    a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*delsq
                    if (a < zero) then
                       tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                    tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 end if
                 ! it can be proved that
                     ! d(n)^2+rho/2 <= sigma_n^2 < d(n)^2+tau2 <= d(n)^2+rho
              else
                 delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                 a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*delsq
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (a < zero) then
                    tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 ! it can be proved that
                 ! d(n)^2 < d(n)^2+tau2 < sigma(n)^2 < d(n)^2+rho/2
              end if
              ! the following tau is to approximate sigma_n - d( n )
               ! tau = tau2 / ( d( n )+sqrt( d( n )*d( n )+tau2 ) )
              sigma = d(n) + tau
              do j = 1,n
                 delta(j) = (d(j) - d(n)) - tau
                 work(j) = d(j) + d(n) + tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(delta(j)*work(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/(delta(n)*work(n))
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              ! calculate the new step
              niter = niter + 1
              dtnsq1 = work(n - 1)*delta(n - 1)
              dtnsq = work(n)*delta(n)
              c = w - dtnsq1*dpsi - dtnsq*dphi
              a = (dtnsq + dtnsq1)*w - dtnsq*dtnsq1*(dpsi + dphi)
              b = dtnsq*dtnsq1*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                 eta = rho - sigma*sigma
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
              temp = eta - dtnsq
              if (temp > rho) eta = rho + dtnsq
              eta = eta/(sigma + sqrt(eta + sigma*sigma))
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 delta(j) = delta(j) - eta
                 work(j) = work(j) + eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              tau2 = work(n)*delta(n)
              temp = z(n)/tau2
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    go to 240
                 end if
                 ! calculate the new step
                 dtnsq1 = work(n - 1)*delta(n - 1)
                 dtnsq = work(n)*delta(n)
                 c = w - dtnsq1*dpsi - dtnsq*dphi
                 a = (dtnsq + dtnsq1)*w - dtnsq1*dtnsq*(dpsi + dphi)
                 b = dtnsq1*dtnsq*w
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
                 temp = eta - dtnsq
                 if (temp <= zero) eta = eta/two
                 eta = eta/(sigma + sqrt(eta + sigma*sigma))
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    delta(j) = delta(j) - eta
                    work(j) = work(j) + eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 tau2 = work(n)*delta(n)
                 temp = z(n)/tau2
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $             + abs( tau2 )*( dpsi+dphi )
                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              go to 240
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              delsq = (d(ip1) - d(i))*(d(ip1) + d(i))
              delsq2 = delsq/two
              sq2 = sqrt((d(i)*d(i) + d(ip1)*d(ip1))/two)
              temp = delsq2/(d(i) + sq2)
              do j = 1,n
                 work(j) = d(j) + d(i) + temp
                 delta(j) = (d(j) - d(i)) - temp
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/(work(j)*delta(j))
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/(work(j)*delta(j))
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/(work(i)*delta(i)) + z(ip1)*z(ip1)/(work(ip1) &
                        *delta(ip1))
              geomavg = .false.
              if (w > zero) then
                 ! d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 ii = i
                 sglb = zero
                 sgub = delsq2/(d(i) + sq2)
                 a = c*delsq + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*delsq
                 if (a > zero) then
                    tau2 = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau2 = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( i )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( i ).
                 tau = tau2/(d(i) + sqrt(d(i)*d(i) + tau2))
                 temp = sqrt(eps)
                 if ((d(i) <= temp*d(ip1)) .and. (abs(z(i)) <= temp) .and. (d(i) > zero)) then
                    tau = min(ten*d(i),sgub)
                    geomavg = .true.
                 end if
              else
                 ! (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 ii = ip1
                 sglb = -delsq2/(d(ii) + sq2)
                 sgub = zero
                 a = c*delsq - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*delsq
                 if (a < zero) then
                    tau2 = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau2 = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( ip1 )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( ip1 ).
                 tau = tau2/(d(ip1) + sqrt(abs(d(ip1)*d(ip1) + tau2)))
              end if
              sigma = d(ii) + tau
              do j = 1,n
                 work(j) = d(j) + d(ii) + tau
                 delta(j) = (d(j) - d(ii)) - tau
              end do
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
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
              temp = z(ii)/(work(ii)*delta(ii))
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              if (w <= zero) then
                 sglb = max(sglb,tau)
              else
                 sgub = min(sgub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 dtipsq = work(ip1)*delta(ip1)
                 dtisq = work(i)*delta(i)
                 if (orgati) then
                    c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                 else
                    c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                 end if
                 a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                 b = dtipsq*dtisq*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
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
                 dtiim = work(iim1)*delta(iim1)
                 dtiip = work(iip1)*delta(iip1)
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/dtiim
                    temp1 = temp1*temp1
                    c = (temp - dtiip*(dpsi + dphi)) - (d(iim1) - d(iip1))*(d(iim1) + d( &
                              iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    if (dpsi < temp1) then
                       zz(3) = dtiip*dtiip*dphi
                    else
                       zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                    end if
                 else
                    temp1 = z(iip1)/dtiip
                    temp1 = temp1*temp1
                    c = (temp - dtiim*(dpsi + dphi)) - (d(iip1) - d(iim1))*(d(iim1) + d( &
                              iip1))*temp1
                    if (dphi < temp1) then
                       zz(1) = dtiim*dtiim*dpsi
                    else
                       zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                    end if
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 dd(1) = dtiim
                 dd(2) = delta(ii)*work(ii)
                 dd(3) = dtiip
                 call la_slaed6(niter,orgati,c,dd,zz,w,eta,info)
                 if (info /= 0) then
                    ! if info is not 0, i.e., la_slaed6 failed, switch back
                    ! to 2 pole interpolation.
                    swtch3 = .false.
                    info = 0
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (orgati) then
                       c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                    else
                       c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (orgati) then
                             a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                          else
                             a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 end if
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              eta = eta/(sigma + sqrt(sigma*sigma + eta))
              temp = tau + eta
              if (temp > sgub .or. temp < sglb) then
                 if (w < zero) then
                    eta = (sgub - tau)/two
                 else
                    eta = (sglb - tau)/two
                 end if
                 if (geomavg) then
                    if (w < zero) then
                       if (tau > zero) then
                          eta = sqrt(sgub*tau) - tau
                       end if
                    else
                       if (sglb > zero) then
                          eta = sqrt(sglb*tau) - tau
                       end if
                    end if
                 end if
              end if
              prew = w
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 work(j) = work(j) + eta
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              tau2 = work(ii)*delta(ii)
              temp = z(ii)/tau2
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              ! main loop to update the values of the array   delta and work
              iter = niter + 1
              loop_230: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
           ! $          .or. (sgub-sglb)<=eight*abs(sgub+sglb) ) then
                    go to 240
                 end if
                 if (w <= zero) then
                    sglb = max(sglb,tau)
                 else
                    sgub = min(sgub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                       else
                          c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                       end if
                    else
                       temp = z(ii)/(work(ii)*delta(ii))
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - dtisq*dpsi - dtipsq*dphi
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                             else
                                a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                             end if
                          else
                             a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
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
                    dtiim = work(iim1)*delta(iim1)
                    dtiip = work(iip1)*delta(iip1)
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - dtiim*dpsi - dtiip*dphi
                       zz(1) = dtiim*dtiim*dpsi
                       zz(3) = dtiip*dtiip*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/dtiim
                          temp1 = temp1*temp1
                          temp2 = (d(iim1) - d(iip1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiip*(dpsi + dphi) - temp2
                          zz(1) = z(iim1)*z(iim1)
                          if (dpsi < temp1) then
                             zz(3) = dtiip*dtiip*dphi
                          else
                             zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                          end if
                       else
                          temp1 = z(iip1)/dtiip
                          temp1 = temp1*temp1
                          temp2 = (d(iip1) - d(iim1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiim*(dpsi + dphi) - temp2
                          if (dphi < temp1) then
                             zz(1) = dtiim*dtiim*dpsi
                          else
                             zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                          end if
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    dd(1) = dtiim
                    dd(2) = delta(ii)*work(ii)
                    dd(3) = dtiip
                    call la_slaed6(niter,orgati,c,dd,zz,w,eta,info)
                    if (info /= 0) then
                       ! if info is not 0, i.e., la_slaed6 failed, switch
                       ! back to two pole interpolation
                       swtch3 = .false.
                       info = 0
                       dtipsq = work(ip1)*delta(ip1)
                       dtisq = work(i)*delta(i)
                       if (.not. swtch) then
                          if (orgati) then
                             c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                          else
                             c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                          end if
                       else
                          temp = z(ii)/(work(ii)*delta(ii))
                          if (orgati) then
                             dpsi = dpsi + temp*temp
                          else
                             dphi = dphi + temp*temp
                          end if
                          c = w - dtisq*dpsi - dtipsq*dphi
                       end if
                       a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                       b = dtipsq*dtisq*w
                       if (c == zero) then
                          if (a == zero) then
                             if (.not. swtch) then
                                if (orgati) then
                                   a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                                else
                                   a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                                end if
                             else
                                a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
                             end if
                          end if
                          eta = b/a
                       else if (a <= zero) then
                          eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                       else
                          eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                       end if
                    end if
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 eta = eta/(sigma + sqrt(sigma*sigma + eta))
                 temp = tau + eta
                 if (temp > sgub .or. temp < sglb) then
                    if (w < zero) then
                       eta = (sgub - tau)/two
                    else
                       eta = (sglb - tau)/two
                    end if
                    if (geomavg) then
                       if (w < zero) then
                          if (tau > zero) then
                             eta = sqrt(sgub*tau) - tau
                          end if
                       else
                          if (sglb > zero) then
                             eta = sqrt(sglb*tau) - tau
                          end if
                       end if
                    end if
                 end if
                 prew = w
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    work(j) = work(j) + eta
                    delta(j) = delta(j) - eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/(work(j)*delta(j))
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 tau2 = work(ii)*delta(ii)
                 temp = z(ii)/tau2
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $             + abs( tau2 )*dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_230
              ! return with info = 1, niter = maxit and not converged
              info = 1
           end if
           240 continue
           return
     end subroutine la_slasd4
     !> This subroutine computes the square root of the I-th updated
     !> eigenvalue of a positive symmetric rank-one modification to
     !> a positive diagonal matrix whose entries are given as the squares
     !> of the corresponding entries in the array d, and that
     !> 0 <= D(i) < D(j)  for  i < j
     !> and that RHO > 0. This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D ) * diag( D ) +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_dlasd4(n,i,d,z,delta,rho,sigma,work,info)
        use la_constants_dp,only:zero,one,two,three,four,eight,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(dp),intent(in) :: rho
           real(dp),intent(out) :: sigma
           ! Array Arguments
           real(dp),intent(in) :: d(*),z(*)
           real(dp),intent(out) :: delta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 400

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3,geomavg
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(dp) :: a,b,c,delsq,delsq2,sq2,dphi,dpsi,dtiim,dtiip,dtipsq,dtisq, &
           dtnsq,dtnsq1,dw,eps,erretm,eta,phi,prew,psi,rhoinv,sglb,sgub,tau,tau2, &
                     temp,temp1,temp2,w
           ! Local Arrays
           real(dp) :: dd(3),zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
              ! presumably, i=1 upon entry
              sigma = sqrt(d(1)*d(1) + rho*z(1)*z(1))
              delta(1) = one
              work(1) = one
              return
           end if
           if (n == 2) then
              call la_dlasd5(i,d,z,delta,rho,sigma,work)
              return
           end if
           ! compute machine epsilon
           eps = la_dlamch('EPSILON')
           rhoinv = one/rho
           tau2 = zero
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              temp = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              temp1 = temp/(d(n) + sqrt(d(n)*d(n) + temp))
              do j = 1,n
                 work(j) = d(j) + d(n) + temp1
                 delta(j) = (d(j) - d(n)) - temp1
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/(delta(j)*work(j))
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/(delta(ii)*work(ii)) + z(n)*z(n)/(delta(n) &
                        *work(n))
              if (w <= zero) then
                 temp1 = sqrt(d(n)*d(n) + rho)
                 temp = z(n - 1)*z(n - 1)/((d(n - 1) + temp1)*(d(n) - d(n - 1) + rho/(d(n) + &
                           temp1))) + z(n)*z(n)/rho
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (c <= temp) then
                    tau = rho
                 else
                    delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                    a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*delsq
                    if (a < zero) then
                       tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                    tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 end if
                 ! it can be proved that
                     ! d(n)^2+rho/2 <= sigma_n^2 < d(n)^2+tau2 <= d(n)^2+rho
              else
                 delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                 a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*delsq
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (a < zero) then
                    tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 ! it can be proved that
                 ! d(n)^2 < d(n)^2+tau2 < sigma(n)^2 < d(n)^2+rho/2
              end if
              ! the following tau is to approximate sigma_n - d( n )
               ! tau = tau2 / ( d( n )+sqrt( d( n )*d( n )+tau2 ) )
              sigma = d(n) + tau
              do j = 1,n
                 delta(j) = (d(j) - d(n)) - tau
                 work(j) = d(j) + d(n) + tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(delta(j)*work(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/(delta(n)*work(n))
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              ! calculate the new step
              niter = niter + 1
              dtnsq1 = work(n - 1)*delta(n - 1)
              dtnsq = work(n)*delta(n)
              c = w - dtnsq1*dpsi - dtnsq*dphi
              a = (dtnsq + dtnsq1)*w - dtnsq*dtnsq1*(dpsi + dphi)
              b = dtnsq*dtnsq1*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                 eta = rho - sigma*sigma
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
              temp = eta - dtnsq
              if (temp > rho) eta = rho + dtnsq
              eta = eta/(sigma + sqrt(eta + sigma*sigma))
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 delta(j) = delta(j) - eta
                 work(j) = work(j) + eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              tau2 = work(n)*delta(n)
              temp = z(n)/tau2
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    go to 240
                 end if
                 ! calculate the new step
                 dtnsq1 = work(n - 1)*delta(n - 1)
                 dtnsq = work(n)*delta(n)
                 c = w - dtnsq1*dpsi - dtnsq*dphi
                 a = (dtnsq + dtnsq1)*w - dtnsq1*dtnsq*(dpsi + dphi)
                 b = dtnsq1*dtnsq*w
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
                 temp = eta - dtnsq
                 if (temp <= zero) eta = eta/two
                 eta = eta/(sigma + sqrt(eta + sigma*sigma))
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    delta(j) = delta(j) - eta
                    work(j) = work(j) + eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 tau2 = work(n)*delta(n)
                 temp = z(n)/tau2
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $             + abs( tau2 )*( dpsi+dphi )
                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              go to 240
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              delsq = (d(ip1) - d(i))*(d(ip1) + d(i))
              delsq2 = delsq/two
              sq2 = sqrt((d(i)*d(i) + d(ip1)*d(ip1))/two)
              temp = delsq2/(d(i) + sq2)
              do j = 1,n
                 work(j) = d(j) + d(i) + temp
                 delta(j) = (d(j) - d(i)) - temp
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/(work(j)*delta(j))
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/(work(j)*delta(j))
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/(work(i)*delta(i)) + z(ip1)*z(ip1)/(work(ip1) &
                        *delta(ip1))
              geomavg = .false.
              if (w > zero) then
                 ! d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 ii = i
                 sglb = zero
                 sgub = delsq2/(d(i) + sq2)
                 a = c*delsq + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*delsq
                 if (a > zero) then
                    tau2 = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau2 = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( i )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( i ).
                 tau = tau2/(d(i) + sqrt(d(i)*d(i) + tau2))
                 temp = sqrt(eps)
                 if ((d(i) <= temp*d(ip1)) .and. (abs(z(i)) <= temp) .and. (d(i) > zero)) then
                    tau = min(ten*d(i),sgub)
                    geomavg = .true.
                 end if
              else
                 ! (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 ii = ip1
                 sglb = -delsq2/(d(ii) + sq2)
                 sgub = zero
                 a = c*delsq - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*delsq
                 if (a < zero) then
                    tau2 = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau2 = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( ip1 )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( ip1 ).
                 tau = tau2/(d(ip1) + sqrt(abs(d(ip1)*d(ip1) + tau2)))
              end if
              sigma = d(ii) + tau
              do j = 1,n
                 work(j) = d(j) + d(ii) + tau
                 delta(j) = (d(j) - d(ii)) - tau
              end do
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
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
              temp = z(ii)/(work(ii)*delta(ii))
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              if (w <= zero) then
                 sglb = max(sglb,tau)
              else
                 sgub = min(sgub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 dtipsq = work(ip1)*delta(ip1)
                 dtisq = work(i)*delta(i)
                 if (orgati) then
                    c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                 else
                    c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                 end if
                 a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                 b = dtipsq*dtisq*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
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
                 dtiim = work(iim1)*delta(iim1)
                 dtiip = work(iip1)*delta(iip1)
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/dtiim
                    temp1 = temp1*temp1
                    c = (temp - dtiip*(dpsi + dphi)) - (d(iim1) - d(iip1))*(d(iim1) + d( &
                              iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    if (dpsi < temp1) then
                       zz(3) = dtiip*dtiip*dphi
                    else
                       zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                    end if
                 else
                    temp1 = z(iip1)/dtiip
                    temp1 = temp1*temp1
                    c = (temp - dtiim*(dpsi + dphi)) - (d(iip1) - d(iim1))*(d(iim1) + d( &
                              iip1))*temp1
                    if (dphi < temp1) then
                       zz(1) = dtiim*dtiim*dpsi
                    else
                       zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                    end if
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 dd(1) = dtiim
                 dd(2) = delta(ii)*work(ii)
                 dd(3) = dtiip
                 call la_dlaed6(niter,orgati,c,dd,zz,w,eta,info)
                 if (info /= 0) then
                    ! if info is not 0, i.e., la_dlaed6 failed, switch back
                    ! to 2 pole interpolation.
                    swtch3 = .false.
                    info = 0
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (orgati) then
                       c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                    else
                       c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (orgati) then
                             a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                          else
                             a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 end if
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              eta = eta/(sigma + sqrt(sigma*sigma + eta))
              temp = tau + eta
              if (temp > sgub .or. temp < sglb) then
                 if (w < zero) then
                    eta = (sgub - tau)/two
                 else
                    eta = (sglb - tau)/two
                 end if
                 if (geomavg) then
                    if (w < zero) then
                       if (tau > zero) then
                          eta = sqrt(sgub*tau) - tau
                       end if
                    else
                       if (sglb > zero) then
                          eta = sqrt(sglb*tau) - tau
                       end if
                    end if
                 end if
              end if
              prew = w
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 work(j) = work(j) + eta
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              tau2 = work(ii)*delta(ii)
              temp = z(ii)/tau2
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              ! main loop to update the values of the array   delta and work
              iter = niter + 1
              loop_230: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
           ! $          .or. (sgub-sglb)<=eight*abs(sgub+sglb) ) then
                    go to 240
                 end if
                 if (w <= zero) then
                    sglb = max(sglb,tau)
                 else
                    sgub = min(sgub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                       else
                          c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                       end if
                    else
                       temp = z(ii)/(work(ii)*delta(ii))
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - dtisq*dpsi - dtipsq*dphi
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                             else
                                a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                             end if
                          else
                             a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
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
                    dtiim = work(iim1)*delta(iim1)
                    dtiip = work(iip1)*delta(iip1)
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - dtiim*dpsi - dtiip*dphi
                       zz(1) = dtiim*dtiim*dpsi
                       zz(3) = dtiip*dtiip*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/dtiim
                          temp1 = temp1*temp1
                          temp2 = (d(iim1) - d(iip1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiip*(dpsi + dphi) - temp2
                          zz(1) = z(iim1)*z(iim1)
                          if (dpsi < temp1) then
                             zz(3) = dtiip*dtiip*dphi
                          else
                             zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                          end if
                       else
                          temp1 = z(iip1)/dtiip
                          temp1 = temp1*temp1
                          temp2 = (d(iip1) - d(iim1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiim*(dpsi + dphi) - temp2
                          if (dphi < temp1) then
                             zz(1) = dtiim*dtiim*dpsi
                          else
                             zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                          end if
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    dd(1) = dtiim
                    dd(2) = delta(ii)*work(ii)
                    dd(3) = dtiip
                    call la_dlaed6(niter,orgati,c,dd,zz,w,eta,info)
                    if (info /= 0) then
                       ! if info is not 0, i.e., la_dlaed6 failed, switch
                       ! back to two pole interpolation
                       swtch3 = .false.
                       info = 0
                       dtipsq = work(ip1)*delta(ip1)
                       dtisq = work(i)*delta(i)
                       if (.not. swtch) then
                          if (orgati) then
                             c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                          else
                             c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                          end if
                       else
                          temp = z(ii)/(work(ii)*delta(ii))
                          if (orgati) then
                             dpsi = dpsi + temp*temp
                          else
                             dphi = dphi + temp*temp
                          end if
                          c = w - dtisq*dpsi - dtipsq*dphi
                       end if
                       a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                       b = dtipsq*dtisq*w
                       if (c == zero) then
                          if (a == zero) then
                             if (.not. swtch) then
                                if (orgati) then
                                   a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                                else
                                   a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                                end if
                             else
                                a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
                             end if
                          end if
                          eta = b/a
                       else if (a <= zero) then
                          eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                       else
                          eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                       end if
                    end if
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 eta = eta/(sigma + sqrt(sigma*sigma + eta))
                 temp = tau + eta
                 if (temp > sgub .or. temp < sglb) then
                    if (w < zero) then
                       eta = (sgub - tau)/two
                    else
                       eta = (sglb - tau)/two
                    end if
                    if (geomavg) then
                       if (w < zero) then
                          if (tau > zero) then
                             eta = sqrt(sgub*tau) - tau
                          end if
                       else
                          if (sglb > zero) then
                             eta = sqrt(sglb*tau) - tau
                          end if
                       end if
                    end if
                 end if
                 prew = w
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    work(j) = work(j) + eta
                    delta(j) = delta(j) - eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/(work(j)*delta(j))
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 tau2 = work(ii)*delta(ii)
                 temp = z(ii)/tau2
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $             + abs( tau2 )*dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_230
              ! return with info = 1, niter = maxit and not converged
              info = 1
           end if
           240 continue
           return
     end subroutine la_dlasd4
#ifdef LA_WITH_XDP
     !> This subroutine computes the square root of the I-th updated
     !> eigenvalue of a positive symmetric rank-one modification to
     !> a positive diagonal matrix whose entries are given as the squares
     !> of the corresponding entries in the array d, and that
     !> 0 <= D(i) < D(j)  for  i < j
     !> and that RHO > 0. This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D ) * diag( D ) +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_xlasd4(n,i,d,z,delta,rho,sigma,work,info)
        use la_constants_xdp,only:zero,one,two,three,four,eight,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(xdp),intent(in) :: rho
           real(xdp),intent(out) :: sigma
           ! Array Arguments
           real(xdp),intent(in) :: d(*),z(*)
           real(xdp),intent(out) :: delta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 400

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3,geomavg
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(xdp) :: a,b,c,delsq,delsq2,sq2,dphi,dpsi,dtiim,dtiip,dtipsq,dtisq, &
           dtnsq,dtnsq1,dw,eps,erretm,eta,phi,prew,psi,rhoinv,sglb,sgub,tau,tau2, &
                     temp,temp1,temp2,w
           ! Local Arrays
           real(xdp) :: dd(3),zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
              ! presumably, i=1 upon entry
              sigma = sqrt(d(1)*d(1) + rho*z(1)*z(1))
              delta(1) = one
              work(1) = one
              return
           end if
           if (n == 2) then
              call la_xlasd5(i,d,z,delta,rho,sigma,work)
              return
           end if
           ! compute machine epsilon
           eps = la_xlamch('EPSILON')
           rhoinv = one/rho
           tau2 = zero
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              temp = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              temp1 = temp/(d(n) + sqrt(d(n)*d(n) + temp))
              do j = 1,n
                 work(j) = d(j) + d(n) + temp1
                 delta(j) = (d(j) - d(n)) - temp1
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/(delta(j)*work(j))
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/(delta(ii)*work(ii)) + z(n)*z(n)/(delta(n) &
                        *work(n))
              if (w <= zero) then
                 temp1 = sqrt(d(n)*d(n) + rho)
                 temp = z(n - 1)*z(n - 1)/((d(n - 1) + temp1)*(d(n) - d(n - 1) + rho/(d(n) + &
                           temp1))) + z(n)*z(n)/rho
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (c <= temp) then
                    tau = rho
                 else
                    delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                    a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*delsq
                    if (a < zero) then
                       tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                    tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 end if
                 ! it can be proved that
                     ! d(n)^2+rho/2 <= sigma_n^2 < d(n)^2+tau2 <= d(n)^2+rho
              else
                 delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                 a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*delsq
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (a < zero) then
                    tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 ! it can be proved that
                 ! d(n)^2 < d(n)^2+tau2 < sigma(n)^2 < d(n)^2+rho/2
              end if
              ! the following tau is to approximate sigma_n - d( n )
               ! tau = tau2 / ( d( n )+sqrt( d( n )*d( n )+tau2 ) )
              sigma = d(n) + tau
              do j = 1,n
                 delta(j) = (d(j) - d(n)) - tau
                 work(j) = d(j) + d(n) + tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(delta(j)*work(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/(delta(n)*work(n))
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              ! calculate the new step
              niter = niter + 1
              dtnsq1 = work(n - 1)*delta(n - 1)
              dtnsq = work(n)*delta(n)
              c = w - dtnsq1*dpsi - dtnsq*dphi
              a = (dtnsq + dtnsq1)*w - dtnsq*dtnsq1*(dpsi + dphi)
              b = dtnsq*dtnsq1*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                 eta = rho - sigma*sigma
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
              temp = eta - dtnsq
              if (temp > rho) eta = rho + dtnsq
              eta = eta/(sigma + sqrt(eta + sigma*sigma))
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 delta(j) = delta(j) - eta
                 work(j) = work(j) + eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              tau2 = work(n)*delta(n)
              temp = z(n)/tau2
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    go to 240
                 end if
                 ! calculate the new step
                 dtnsq1 = work(n - 1)*delta(n - 1)
                 dtnsq = work(n)*delta(n)
                 c = w - dtnsq1*dpsi - dtnsq*dphi
                 a = (dtnsq + dtnsq1)*w - dtnsq1*dtnsq*(dpsi + dphi)
                 b = dtnsq1*dtnsq*w
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
                 temp = eta - dtnsq
                 if (temp <= zero) eta = eta/two
                 eta = eta/(sigma + sqrt(eta + sigma*sigma))
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    delta(j) = delta(j) - eta
                    work(j) = work(j) + eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 tau2 = work(n)*delta(n)
                 temp = z(n)/tau2
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $             + abs( tau2 )*( dpsi+dphi )
                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              go to 240
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              delsq = (d(ip1) - d(i))*(d(ip1) + d(i))
              delsq2 = delsq/two
              sq2 = sqrt((d(i)*d(i) + d(ip1)*d(ip1))/two)
              temp = delsq2/(d(i) + sq2)
              do j = 1,n
                 work(j) = d(j) + d(i) + temp
                 delta(j) = (d(j) - d(i)) - temp
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/(work(j)*delta(j))
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/(work(j)*delta(j))
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/(work(i)*delta(i)) + z(ip1)*z(ip1)/(work(ip1) &
                        *delta(ip1))
              geomavg = .false.
              if (w > zero) then
                 ! d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 ii = i
                 sglb = zero
                 sgub = delsq2/(d(i) + sq2)
                 a = c*delsq + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*delsq
                 if (a > zero) then
                    tau2 = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau2 = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( i )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( i ).
                 tau = tau2/(d(i) + sqrt(d(i)*d(i) + tau2))
                 temp = sqrt(eps)
                 if ((d(i) <= temp*d(ip1)) .and. (abs(z(i)) <= temp) .and. (d(i) > zero)) then
                    tau = min(ten*d(i),sgub)
                    geomavg = .true.
                 end if
              else
                 ! (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 ii = ip1
                 sglb = -delsq2/(d(ii) + sq2)
                 sgub = zero
                 a = c*delsq - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*delsq
                 if (a < zero) then
                    tau2 = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau2 = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( ip1 )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( ip1 ).
                 tau = tau2/(d(ip1) + sqrt(abs(d(ip1)*d(ip1) + tau2)))
              end if
              sigma = d(ii) + tau
              do j = 1,n
                 work(j) = d(j) + d(ii) + tau
                 delta(j) = (d(j) - d(ii)) - tau
              end do
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
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
              temp = z(ii)/(work(ii)*delta(ii))
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              if (w <= zero) then
                 sglb = max(sglb,tau)
              else
                 sgub = min(sgub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 dtipsq = work(ip1)*delta(ip1)
                 dtisq = work(i)*delta(i)
                 if (orgati) then
                    c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                 else
                    c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                 end if
                 a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                 b = dtipsq*dtisq*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
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
                 dtiim = work(iim1)*delta(iim1)
                 dtiip = work(iip1)*delta(iip1)
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/dtiim
                    temp1 = temp1*temp1
                    c = (temp - dtiip*(dpsi + dphi)) - (d(iim1) - d(iip1))*(d(iim1) + d( &
                              iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    if (dpsi < temp1) then
                       zz(3) = dtiip*dtiip*dphi
                    else
                       zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                    end if
                 else
                    temp1 = z(iip1)/dtiip
                    temp1 = temp1*temp1
                    c = (temp - dtiim*(dpsi + dphi)) - (d(iip1) - d(iim1))*(d(iim1) + d( &
                              iip1))*temp1
                    if (dphi < temp1) then
                       zz(1) = dtiim*dtiim*dpsi
                    else
                       zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                    end if
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 dd(1) = dtiim
                 dd(2) = delta(ii)*work(ii)
                 dd(3) = dtiip
                 call la_xlaed6(niter,orgati,c,dd,zz,w,eta,info)
                 if (info /= 0) then
                    ! if info is not 0, i.e., la_xlaed6 failed, switch back
                    ! to 2 pole interpolation.
                    swtch3 = .false.
                    info = 0
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (orgati) then
                       c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                    else
                       c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (orgati) then
                             a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                          else
                             a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 end if
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              eta = eta/(sigma + sqrt(sigma*sigma + eta))
              temp = tau + eta
              if (temp > sgub .or. temp < sglb) then
                 if (w < zero) then
                    eta = (sgub - tau)/two
                 else
                    eta = (sglb - tau)/two
                 end if
                 if (geomavg) then
                    if (w < zero) then
                       if (tau > zero) then
                          eta = sqrt(sgub*tau) - tau
                       end if
                    else
                       if (sglb > zero) then
                          eta = sqrt(sglb*tau) - tau
                       end if
                    end if
                 end if
              end if
              prew = w
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 work(j) = work(j) + eta
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              tau2 = work(ii)*delta(ii)
              temp = z(ii)/tau2
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              ! main loop to update the values of the array   delta and work
              iter = niter + 1
              loop_230: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
           ! $          .or. (sgub-sglb)<=eight*abs(sgub+sglb) ) then
                    go to 240
                 end if
                 if (w <= zero) then
                    sglb = max(sglb,tau)
                 else
                    sgub = min(sgub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                       else
                          c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                       end if
                    else
                       temp = z(ii)/(work(ii)*delta(ii))
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - dtisq*dpsi - dtipsq*dphi
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                             else
                                a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                             end if
                          else
                             a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
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
                    dtiim = work(iim1)*delta(iim1)
                    dtiip = work(iip1)*delta(iip1)
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - dtiim*dpsi - dtiip*dphi
                       zz(1) = dtiim*dtiim*dpsi
                       zz(3) = dtiip*dtiip*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/dtiim
                          temp1 = temp1*temp1
                          temp2 = (d(iim1) - d(iip1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiip*(dpsi + dphi) - temp2
                          zz(1) = z(iim1)*z(iim1)
                          if (dpsi < temp1) then
                             zz(3) = dtiip*dtiip*dphi
                          else
                             zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                          end if
                       else
                          temp1 = z(iip1)/dtiip
                          temp1 = temp1*temp1
                          temp2 = (d(iip1) - d(iim1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiim*(dpsi + dphi) - temp2
                          if (dphi < temp1) then
                             zz(1) = dtiim*dtiim*dpsi
                          else
                             zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                          end if
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    dd(1) = dtiim
                    dd(2) = delta(ii)*work(ii)
                    dd(3) = dtiip
                    call la_xlaed6(niter,orgati,c,dd,zz,w,eta,info)
                    if (info /= 0) then
                       ! if info is not 0, i.e., la_xlaed6 failed, switch
                       ! back to two pole interpolation
                       swtch3 = .false.
                       info = 0
                       dtipsq = work(ip1)*delta(ip1)
                       dtisq = work(i)*delta(i)
                       if (.not. swtch) then
                          if (orgati) then
                             c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                          else
                             c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                          end if
                       else
                          temp = z(ii)/(work(ii)*delta(ii))
                          if (orgati) then
                             dpsi = dpsi + temp*temp
                          else
                             dphi = dphi + temp*temp
                          end if
                          c = w - dtisq*dpsi - dtipsq*dphi
                       end if
                       a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                       b = dtipsq*dtisq*w
                       if (c == zero) then
                          if (a == zero) then
                             if (.not. swtch) then
                                if (orgati) then
                                   a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                                else
                                   a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                                end if
                             else
                                a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
                             end if
                          end if
                          eta = b/a
                       else if (a <= zero) then
                          eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                       else
                          eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                       end if
                    end if
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 eta = eta/(sigma + sqrt(sigma*sigma + eta))
                 temp = tau + eta
                 if (temp > sgub .or. temp < sglb) then
                    if (w < zero) then
                       eta = (sgub - tau)/two
                    else
                       eta = (sglb - tau)/two
                    end if
                    if (geomavg) then
                       if (w < zero) then
                          if (tau > zero) then
                             eta = sqrt(sgub*tau) - tau
                          end if
                       else
                          if (sglb > zero) then
                             eta = sqrt(sglb*tau) - tau
                          end if
                       end if
                    end if
                 end if
                 prew = w
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    work(j) = work(j) + eta
                    delta(j) = delta(j) - eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/(work(j)*delta(j))
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 tau2 = work(ii)*delta(ii)
                 temp = z(ii)/tau2
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $             + abs( tau2 )*dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_230
              ! return with info = 1, niter = maxit and not converged
              info = 1
           end if
           240 continue
           return
     end subroutine la_xlasd4
#endif
#ifdef LA_WITH_QP
     !> This subroutine computes the square root of the I-th updated
     !> eigenvalue of a positive symmetric rank-one modification to
     !> a positive diagonal matrix whose entries are given as the squares
     !> of the corresponding entries in the array d, and that
     !> 0 <= D(i) < D(j)  for  i < j
     !> and that RHO > 0. This is arranged by the calling routine, and is
     !> no loss in generality.  The rank-one modified system is thus
     !> diag( D ) * diag( D ) +  RHO * Z * Z_transpose.
     !> where we assume the Euclidean norm of Z is 1.
     !> The method consists of approximating the rational functions in the
     !> secular equation by simpler interpolating rational functions.

     pure subroutine la_qlasd4(n,i,d,z,delta,rho,sigma,work,info)
        use la_constants_qp,only:zero,one,two,three,four,eight,ten
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i,n
           integer(ilp),intent(out) :: info
           real(qp),intent(in) :: rho
           real(qp),intent(out) :: sigma
           ! Array Arguments
           real(qp),intent(in) :: d(*),z(*)
           real(qp),intent(out) :: delta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 400

           ! Local Scalars
           logical(lk) :: orgati,swtch,swtch3,geomavg
           integer(ilp) :: ii,iim1,iip1,ip1,iter,j,niter
           real(qp) :: a,b,c,delsq,delsq2,sq2,dphi,dpsi,dtiim,dtiip,dtipsq,dtisq, &
           dtnsq,dtnsq1,dw,eps,erretm,eta,phi,prew,psi,rhoinv,sglb,sgub,tau,tau2, &
                     temp,temp1,temp2,w
           ! Local Arrays
           real(qp) :: dd(3),zz(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           ! since this routine is called in an inner loop, we do no argument
           ! checking.
           ! quick return for n=1 and 2.
           info = 0
           if (n == 1) then
              ! presumably, i=1 upon entry
              sigma = sqrt(d(1)*d(1) + rho*z(1)*z(1))
              delta(1) = one
              work(1) = one
              return
           end if
           if (n == 2) then
              call la_qlasd5(i,d,z,delta,rho,sigma,work)
              return
           end if
           ! compute machine epsilon
           eps = la_qlamch('EPSILON')
           rhoinv = one/rho
           tau2 = zero
           ! the case i = n
           if (i == n) then
              ! initialize some basic variables
              ii = n - 1
              niter = 1
              ! calculate initial guess
              temp = rho/two
              ! if ||z||_2 is not one, then temp should be set to
              ! rho * ||z||_2^2 / two
              temp1 = temp/(d(n) + sqrt(d(n)*d(n) + temp))
              do j = 1,n
                 work(j) = d(j) + d(n) + temp1
                 delta(j) = (d(j) - d(n)) - temp1
              end do
              psi = zero
              do j = 1,n - 2
                 psi = psi + z(j)*z(j)/(delta(j)*work(j))
              end do
              c = rhoinv + psi
              w = c + z(ii)*z(ii)/(delta(ii)*work(ii)) + z(n)*z(n)/(delta(n) &
                        *work(n))
              if (w <= zero) then
                 temp1 = sqrt(d(n)*d(n) + rho)
                 temp = z(n - 1)*z(n - 1)/((d(n - 1) + temp1)*(d(n) - d(n - 1) + rho/(d(n) + &
                           temp1))) + z(n)*z(n)/rho
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (c <= temp) then
                    tau = rho
                 else
                    delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                    a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*delsq
                    if (a < zero) then
                       tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                    else
                       tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                    end if
                    tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 end if
                 ! it can be proved that
                     ! d(n)^2+rho/2 <= sigma_n^2 < d(n)^2+tau2 <= d(n)^2+rho
              else
                 delsq = (d(n) - d(n - 1))*(d(n) + d(n - 1))
                 a = -c*delsq + z(n - 1)*z(n - 1) + z(n)*z(n)
                 b = z(n)*z(n)*delsq
                 ! the following tau2 is to approximate
                 ! sigma_n^2 - d( n )*d( n )
                 if (a < zero) then
                    tau2 = two*b/(sqrt(a*a + four*b*c) - a)
                 else
                    tau2 = (a + sqrt(a*a + four*b*c))/(two*c)
                 end if
                 tau = tau2/(d(n) + sqrt(d(n)*d(n) + tau2))
                 ! it can be proved that
                 ! d(n)^2 < d(n)^2+tau2 < sigma(n)^2 < d(n)^2+rho/2
              end if
              ! the following tau is to approximate sigma_n - d( n )
               ! tau = tau2 / ( d( n )+sqrt( d( n )*d( n )+tau2 ) )
              sigma = d(n) + tau
              do j = 1,n
                 delta(j) = (d(j) - d(n)) - tau
                 work(j) = d(j) + d(n) + tau
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(delta(j)*work(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              temp = z(n)/(delta(n)*work(n))
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              ! calculate the new step
              niter = niter + 1
              dtnsq1 = work(n - 1)*delta(n - 1)
              dtnsq = work(n)*delta(n)
              c = w - dtnsq1*dpsi - dtnsq*dphi
              a = (dtnsq + dtnsq1)*w - dtnsq*dtnsq1*(dpsi + dphi)
              b = dtnsq*dtnsq1*w
              if (c < zero) c = abs(c)
              if (c == zero) then
                 eta = rho - sigma*sigma
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
              temp = eta - dtnsq
              if (temp > rho) eta = rho + dtnsq
              eta = eta/(sigma + sqrt(eta + sigma*sigma))
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 delta(j) = delta(j) - eta
                 work(j) = work(j) + eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,ii
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              tau2 = work(n)*delta(n)
              temp = z(n)/tau2
              phi = z(n)*temp
              dphi = temp*temp
              erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $          + abs( tau2 )*( dpsi+dphi )
              w = rhoinv + phi + psi
              ! main loop to update the values of the array   delta
              iter = niter + 1
              loop_90: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
                    go to 240
                 end if
                 ! calculate the new step
                 dtnsq1 = work(n - 1)*delta(n - 1)
                 dtnsq = work(n)*delta(n)
                 c = w - dtnsq1*dpsi - dtnsq*dphi
                 a = (dtnsq + dtnsq1)*w - dtnsq1*dtnsq*(dpsi + dphi)
                 b = dtnsq1*dtnsq*w
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
                 temp = eta - dtnsq
                 if (temp <= zero) eta = eta/two
                 eta = eta/(sigma + sqrt(eta + sigma*sigma))
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    delta(j) = delta(j) - eta
                    work(j) = work(j) + eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,ii
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 tau2 = work(n)*delta(n)
                 temp = z(n)/tau2
                 phi = z(n)*temp
                 dphi = temp*temp
                 erretm = eight*(-phi - psi) + erretm - phi + rhoinv
          ! $             + abs( tau2 )*( dpsi+dphi )
                 w = rhoinv + phi + psi
              end do loop_90
              ! return with info = 1, niter = maxit and not converged
              info = 1
              go to 240
              ! end for the case i = n
           else
              ! the case for i < n
              niter = 1
              ip1 = i + 1
              ! calculate initial guess
              delsq = (d(ip1) - d(i))*(d(ip1) + d(i))
              delsq2 = delsq/two
              sq2 = sqrt((d(i)*d(i) + d(ip1)*d(ip1))/two)
              temp = delsq2/(d(i) + sq2)
              do j = 1,n
                 work(j) = d(j) + d(i) + temp
                 delta(j) = (d(j) - d(i)) - temp
              end do
              psi = zero
              do j = 1,i - 1
                 psi = psi + z(j)*z(j)/(work(j)*delta(j))
              end do
              phi = zero
              do j = n,i + 2,-1
                 phi = phi + z(j)*z(j)/(work(j)*delta(j))
              end do
              c = rhoinv + psi + phi
              w = c + z(i)*z(i)/(work(i)*delta(i)) + z(ip1)*z(ip1)/(work(ip1) &
                        *delta(ip1))
              geomavg = .false.
              if (w > zero) then
                 ! d(i)^2 < the ith sigma^2 < (d(i)^2+d(i+1)^2)/2
                 ! we choose d(i) as origin.
                 orgati = .true.
                 ii = i
                 sglb = zero
                 sgub = delsq2/(d(i) + sq2)
                 a = c*delsq + z(i)*z(i) + z(ip1)*z(ip1)
                 b = z(i)*z(i)*delsq
                 if (a > zero) then
                    tau2 = two*b/(a + sqrt(abs(a*a - four*b*c)))
                 else
                    tau2 = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( i )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( i ).
                 tau = tau2/(d(i) + sqrt(d(i)*d(i) + tau2))
                 temp = sqrt(eps)
                 if ((d(i) <= temp*d(ip1)) .and. (abs(z(i)) <= temp) .and. (d(i) > zero)) then
                    tau = min(ten*d(i),sgub)
                    geomavg = .true.
                 end if
              else
                 ! (d(i)^2+d(i+1)^2)/2 <= the ith sigma^2 < d(i+1)^2/2
                 ! we choose d(i+1) as origin.
                 orgati = .false.
                 ii = ip1
                 sglb = -delsq2/(d(ii) + sq2)
                 sgub = zero
                 a = c*delsq - z(i)*z(i) - z(ip1)*z(ip1)
                 b = z(ip1)*z(ip1)*delsq
                 if (a < zero) then
                    tau2 = two*b/(a - sqrt(abs(a*a + four*b*c)))
                 else
                    tau2 = -(a + sqrt(abs(a*a + four*b*c)))/(two*c)
                 end if
                 ! tau2 now is an estimation of sigma^2 - d( ip1 )^2. the
                 ! following, however, is the corresponding estimation of
                 ! sigma - d( ip1 ).
                 tau = tau2/(d(ip1) + sqrt(abs(d(ip1)*d(ip1) + tau2)))
              end if
              sigma = d(ii) + tau
              do j = 1,n
                 work(j) = d(j) + d(ii) + tau
                 delta(j) = (d(j) - d(ii)) - tau
              end do
              iim1 = ii - 1
              iip1 = ii + 1
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
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
              temp = z(ii)/(work(ii)*delta(ii))
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = w + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              ! test for convergence
              if (abs(w) <= eps*erretm) then
                 go to 240
              end if
              if (w <= zero) then
                 sglb = max(sglb,tau)
              else
                 sgub = min(sgub,tau)
              end if
              ! calculate the new step
              niter = niter + 1
              if (.not. swtch3) then
                 dtipsq = work(ip1)*delta(ip1)
                 dtisq = work(i)*delta(i)
                 if (orgati) then
                    c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                 else
                    c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                 end if
                 a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                 b = dtipsq*dtisq*w
                 if (c == zero) then
                    if (a == zero) then
                       if (orgati) then
                          a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                       else
                          a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
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
                 dtiim = work(iim1)*delta(iim1)
                 dtiip = work(iip1)*delta(iip1)
                 temp = rhoinv + psi + phi
                 if (orgati) then
                    temp1 = z(iim1)/dtiim
                    temp1 = temp1*temp1
                    c = (temp - dtiip*(dpsi + dphi)) - (d(iim1) - d(iip1))*(d(iim1) + d( &
                              iip1))*temp1
                    zz(1) = z(iim1)*z(iim1)
                    if (dpsi < temp1) then
                       zz(3) = dtiip*dtiip*dphi
                    else
                       zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                    end if
                 else
                    temp1 = z(iip1)/dtiip
                    temp1 = temp1*temp1
                    c = (temp - dtiim*(dpsi + dphi)) - (d(iip1) - d(iim1))*(d(iim1) + d( &
                              iip1))*temp1
                    if (dphi < temp1) then
                       zz(1) = dtiim*dtiim*dpsi
                    else
                       zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                    end if
                    zz(3) = z(iip1)*z(iip1)
                 end if
                 zz(2) = z(ii)*z(ii)
                 dd(1) = dtiim
                 dd(2) = delta(ii)*work(ii)
                 dd(3) = dtiip
                 call la_qlaed6(niter,orgati,c,dd,zz,w,eta,info)
                 if (info /= 0) then
                    ! if info is not 0, i.e., la_qlaed6 failed, switch back
                    ! to 2 pole interpolation.
                    swtch3 = .false.
                    info = 0
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (orgati) then
                       c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                    else
                       c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (orgati) then
                             a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                          else
                             a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                          end if
                       end if
                       eta = b/a
                    else if (a <= zero) then
                       eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                    else
                       eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                    end if
                 end if
              end if
              ! note, eta should be positive if w is negative, and
              ! eta should be negative otherwise. however,
              ! if for some reason caused by roundoff, eta*w > 0,
              ! we simply use one newton step instead. this way
              ! will guarantee eta*w < 0.
              if (w*eta >= zero) eta = -w/dw
              eta = eta/(sigma + sqrt(sigma*sigma + eta))
              temp = tau + eta
              if (temp > sgub .or. temp < sglb) then
                 if (w < zero) then
                    eta = (sgub - tau)/two
                 else
                    eta = (sglb - tau)/two
                 end if
                 if (geomavg) then
                    if (w < zero) then
                       if (tau > zero) then
                          eta = sqrt(sgub*tau) - tau
                       end if
                    else
                       if (sglb > zero) then
                          eta = sqrt(sglb*tau) - tau
                       end if
                    end if
                 end if
              end if
              prew = w
              tau = tau + eta
              sigma = sigma + eta
              do j = 1,n
                 work(j) = work(j) + eta
                 delta(j) = delta(j) - eta
              end do
              ! evaluate psi and the derivative dpsi
              dpsi = zero
              psi = zero
              erretm = zero
              do j = 1,iim1
                 temp = z(j)/(work(j)*delta(j))
                 psi = psi + z(j)*temp
                 dpsi = dpsi + temp*temp
                 erretm = erretm + psi
              end do
              erretm = abs(erretm)
              ! evaluate phi and the derivative dphi
              dphi = zero
              phi = zero
              do j = n,iip1,-1
                 temp = z(j)/(work(j)*delta(j))
                 phi = phi + z(j)*temp
                 dphi = dphi + temp*temp
                 erretm = erretm + phi
              end do
              tau2 = work(ii)*delta(ii)
              temp = z(ii)/tau2
              dw = dpsi + dphi + temp*temp
              temp = z(ii)*temp
              w = rhoinv + phi + psi + temp
              erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $          + abs( tau2 )*dw
              swtch = .false.
              if (orgati) then
                 if (-w > abs(prew)/ten) swtch = .true.
              else
                 if (w > abs(prew)/ten) swtch = .true.
              end if
              ! main loop to update the values of the array   delta and work
              iter = niter + 1
              loop_230: do niter = iter,maxit
                 ! test for convergence
                 if (abs(w) <= eps*erretm) then
           ! $          .or. (sgub-sglb)<=eight*abs(sgub+sglb) ) then
                    go to 240
                 end if
                 if (w <= zero) then
                    sglb = max(sglb,tau)
                 else
                    sgub = min(sgub,tau)
                 end if
                 ! calculate the new step
                 if (.not. swtch3) then
                    dtipsq = work(ip1)*delta(ip1)
                    dtisq = work(i)*delta(i)
                    if (.not. swtch) then
                       if (orgati) then
                          c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                       else
                          c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                       end if
                    else
                       temp = z(ii)/(work(ii)*delta(ii))
                       if (orgati) then
                          dpsi = dpsi + temp*temp
                       else
                          dphi = dphi + temp*temp
                       end if
                       c = w - dtisq*dpsi - dtipsq*dphi
                    end if
                    a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                    b = dtipsq*dtisq*w
                    if (c == zero) then
                       if (a == zero) then
                          if (.not. swtch) then
                             if (orgati) then
                                a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                             else
                                a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                             end if
                          else
                             a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
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
                    dtiim = work(iim1)*delta(iim1)
                    dtiip = work(iip1)*delta(iip1)
                    temp = rhoinv + psi + phi
                    if (swtch) then
                       c = temp - dtiim*dpsi - dtiip*dphi
                       zz(1) = dtiim*dtiim*dpsi
                       zz(3) = dtiip*dtiip*dphi
                    else
                       if (orgati) then
                          temp1 = z(iim1)/dtiim
                          temp1 = temp1*temp1
                          temp2 = (d(iim1) - d(iip1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiip*(dpsi + dphi) - temp2
                          zz(1) = z(iim1)*z(iim1)
                          if (dpsi < temp1) then
                             zz(3) = dtiip*dtiip*dphi
                          else
                             zz(3) = dtiip*dtiip*((dpsi - temp1) + dphi)
                          end if
                       else
                          temp1 = z(iip1)/dtiip
                          temp1 = temp1*temp1
                          temp2 = (d(iip1) - d(iim1))*(d(iim1) + d(iip1))*temp1
                          c = temp - dtiim*(dpsi + dphi) - temp2
                          if (dphi < temp1) then
                             zz(1) = dtiim*dtiim*dpsi
                          else
                             zz(1) = dtiim*dtiim*(dpsi + (dphi - temp1))
                          end if
                          zz(3) = z(iip1)*z(iip1)
                       end if
                    end if
                    dd(1) = dtiim
                    dd(2) = delta(ii)*work(ii)
                    dd(3) = dtiip
                    call la_qlaed6(niter,orgati,c,dd,zz,w,eta,info)
                    if (info /= 0) then
                       ! if info is not 0, i.e., la_qlaed6 failed, switch
                       ! back to two pole interpolation
                       swtch3 = .false.
                       info = 0
                       dtipsq = work(ip1)*delta(ip1)
                       dtisq = work(i)*delta(i)
                       if (.not. swtch) then
                          if (orgati) then
                             c = w - dtipsq*dw + delsq*(z(i)/dtisq)**2
                          else
                             c = w - dtisq*dw - delsq*(z(ip1)/dtipsq)**2
                          end if
                       else
                          temp = z(ii)/(work(ii)*delta(ii))
                          if (orgati) then
                             dpsi = dpsi + temp*temp
                          else
                             dphi = dphi + temp*temp
                          end if
                          c = w - dtisq*dpsi - dtipsq*dphi
                       end if
                       a = (dtipsq + dtisq)*w - dtipsq*dtisq*dw
                       b = dtipsq*dtisq*w
                       if (c == zero) then
                          if (a == zero) then
                             if (.not. swtch) then
                                if (orgati) then
                                   a = z(i)*z(i) + dtipsq*dtipsq*(dpsi + dphi)
                                else
                                   a = z(ip1)*z(ip1) + dtisq*dtisq*(dpsi + dphi)
                                end if
                             else
                                a = dtisq*dtisq*dpsi + dtipsq*dtipsq*dphi
                             end if
                          end if
                          eta = b/a
                       else if (a <= zero) then
                          eta = (a - sqrt(abs(a*a - four*b*c)))/(two*c)
                       else
                          eta = two*b/(a + sqrt(abs(a*a - four*b*c)))
                       end if
                    end if
                 end if
                 ! note, eta should be positive if w is negative, and
                 ! eta should be negative otherwise. however,
                 ! if for some reason caused by roundoff, eta*w > 0,
                 ! we simply use one newton step instead. this way
                 ! will guarantee eta*w < 0.
                 if (w*eta >= zero) eta = -w/dw
                 eta = eta/(sigma + sqrt(sigma*sigma + eta))
                 temp = tau + eta
                 if (temp > sgub .or. temp < sglb) then
                    if (w < zero) then
                       eta = (sgub - tau)/two
                    else
                       eta = (sglb - tau)/two
                    end if
                    if (geomavg) then
                       if (w < zero) then
                          if (tau > zero) then
                             eta = sqrt(sgub*tau) - tau
                          end if
                       else
                          if (sglb > zero) then
                             eta = sqrt(sglb*tau) - tau
                          end if
                       end if
                    end if
                 end if
                 prew = w
                 tau = tau + eta
                 sigma = sigma + eta
                 do j = 1,n
                    work(j) = work(j) + eta
                    delta(j) = delta(j) - eta
                 end do
                 ! evaluate psi and the derivative dpsi
                 dpsi = zero
                 psi = zero
                 erretm = zero
                 do j = 1,iim1
                    temp = z(j)/(work(j)*delta(j))
                    psi = psi + z(j)*temp
                    dpsi = dpsi + temp*temp
                    erretm = erretm + psi
                 end do
                 erretm = abs(erretm)
                 ! evaluate phi and the derivative dphi
                 dphi = zero
                 phi = zero
                 do j = n,iip1,-1
                    temp = z(j)/(work(j)*delta(j))
                    phi = phi + z(j)*temp
                    dphi = dphi + temp*temp
                    erretm = erretm + phi
                 end do
                 tau2 = work(ii)*delta(ii)
                 temp = z(ii)/tau2
                 dw = dpsi + dphi + temp*temp
                 temp = z(ii)*temp
                 w = rhoinv + phi + psi + temp
                 erretm = eight*(phi - psi) + erretm + two*rhoinv + three*abs(temp)
          ! $             + abs( tau2 )*dw
                 if (w*prew > zero .and. abs(w) > abs(prew)/ten) swtch = .not. swtch
              end do loop_230
              ! return with info = 1, niter = maxit and not converged
              info = 1
           end if
           240 continue
           return
     end subroutine la_qlasd4
#endif

     !> SLASD7: merges the two sets of singular values together into a single
     !> sorted set. Then it tries to deflate the size of the problem. There
     !> are two ways in which deflation can occur:  when two or more singular
     !> values are close together or if there is a tiny entry in the Z
     !> vector. For each such occurrence the order of the related
     !> secular equation problem is reduced by one.
     !> SLASD7 is called from SLASD6.

     pure subroutine la_slasd7(icompq,nl,nr,sqre,k,d,z,zw,vf,vfw,vl,vlw,alpha, &
     beta,dsigma,idx,idxp,idxq,perm,givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
        use la_constants_sp,only:zero,one,two,eight

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(sp),intent(in) :: alpha,beta
           real(sp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),idx(*),idxp(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(sp),intent(inout) :: d(*),vf(*),vl(*)
           real(sp),intent(out) :: dsigma(*),givnum(ldgnum,*),vfw(*),vlw(*),z(*),zw(*)

        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(sp) :: eps,hlftol,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -22
           else if (ldgnum < n) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('SLASD7',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           if (icompq == 1) then
              givptr = 0
           end if
           ! generate the first part of the vector z and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vl(nlp1)
           vl(nlp1) = zero
           tau = vf(nlp1)
           do i = nl,1,-1
              z(i + 1) = alpha*vl(i)
              vl(i) = zero
              vf(i + 1) = vf(i)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           vf(1) = tau
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vf(i)
              vf(i) = zero
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and zw are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              zw(i) = z(idxq(i))
              vfw(i) = vf(idxq(i))
              vlw(i) = vl(idxq(i))
           end do
           call la_slamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = zw(idxi)
              vf(i) = vfw(idxi)
              vl(i) = vlw(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_slamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 if (j == n) go to 100
              else
                 jprev = j
                 go to 70
              end if
           end do
           70 continue
           j = jprev
           80 continue
           j = j + 1
           if (j > n) go to 90
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_slapy2(c,s)
                 z(j) = tau
                 z(jprev) = zero
                 c = c/tau
                 s = -s/tau
                 ! record the appropriate givens rotation
                 if (icompq == 1) then
                    givptr = givptr + 1
                    idxjp = idxq(idx(jprev) + 1)
                    idxj = idxq(idx(j) + 1)
                    if (idxjp <= nlp1) then
                       idxjp = idxjp - 1
                    end if
                    if (idxj <= nlp1) then
                       idxj = idxj - 1
                    end if
                    givcol(givptr,2) = idxjp
                    givcol(givptr,1) = idxj
                    givnum(givptr,2) = c
                    givnum(givptr,1) = s
                 end if
                 call la_srot(1,vf(jprev),1,vf(j),1,c,s)
                 call la_srot(1,vl(jprev),1,vl(j),1,c,s)
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 zw(k) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 80
           90 continue
           ! record the last singular value.
           k = k + 1
           zw(k) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           100 continue
           ! sort the singular values into dsigma. the singular values which
           ! were not deflated go into the first k slots of dsigma, except
           ! that dsigma(1) is treated separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              vfw(j) = vf(jp)
              vlw(j) = vl(jp)
           end do
           if (icompq == 1) then
              do j = 2,n
                 jp = idxp(j)
                 perm(j) = idxq(idx(jp) + 1)
                 if (perm(j) <= nlp1) then
                    perm(j) = perm(j) - 1
                 end if
              end do
           end if
           ! the deflated singular values go back into the last n - k slots of
           ! d.
           call la_scopy(n - k,dsigma(k + 1),1,d(k + 1),1)
           ! determine dsigma(1), dsigma(2), z(1), vf(1), vl(1), vf(m), and
           ! vl(m).
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_slapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = -z(m)/z(1)
              end if
              call la_srot(1,vf(m),1,vf(1),1,c,s)
              call la_srot(1,vl(m),1,vl(1),1,c,s)
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! restore z, vf, and vl.
           call la_scopy(k - 1,zw(2),1,z(2),1)
           call la_scopy(n - 1,vfw(2),1,vf(2),1)
           call la_scopy(n - 1,vlw(2),1,vl(2),1)
           return
     end subroutine la_slasd7
     !> DLASD7: merges the two sets of singular values together into a single
     !> sorted set. Then it tries to deflate the size of the problem. There
     !> are two ways in which deflation can occur:  when two or more singular
     !> values are close together or if there is a tiny entry in the Z
     !> vector. For each such occurrence the order of the related
     !> secular equation problem is reduced by one.
     !> DLASD7 is called from DLASD6.

     pure subroutine la_dlasd7(icompq,nl,nr,sqre,k,d,z,zw,vf,vfw,vl,vlw,alpha, &
     beta,dsigma,idx,idxp,idxq,perm,givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
        use la_constants_dp,only:zero,one,two,eight

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(dp),intent(in) :: alpha,beta
           real(dp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),idx(*),idxp(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(dp),intent(inout) :: d(*),vf(*),vl(*)
           real(dp),intent(out) :: dsigma(*),givnum(ldgnum,*),vfw(*),vlw(*),z(*),zw(*)

        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(dp) :: eps,hlftol,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -22
           else if (ldgnum < n) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('DLASD7',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           if (icompq == 1) then
              givptr = 0
           end if
           ! generate the first part of the vector z and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vl(nlp1)
           vl(nlp1) = zero
           tau = vf(nlp1)
           do i = nl,1,-1
              z(i + 1) = alpha*vl(i)
              vl(i) = zero
              vf(i + 1) = vf(i)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           vf(1) = tau
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vf(i)
              vf(i) = zero
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and zw are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              zw(i) = z(idxq(i))
              vfw(i) = vf(idxq(i))
              vlw(i) = vl(idxq(i))
           end do
           call la_dlamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = zw(idxi)
              vf(i) = vfw(idxi)
              vl(i) = vlw(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_dlamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 if (j == n) go to 100
              else
                 jprev = j
                 go to 70
              end if
           end do
           70 continue
           j = jprev
           80 continue
           j = j + 1
           if (j > n) go to 90
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_dlapy2(c,s)
                 z(j) = tau
                 z(jprev) = zero
                 c = c/tau
                 s = -s/tau
                 ! record the appropriate givens rotation
                 if (icompq == 1) then
                    givptr = givptr + 1
                    idxjp = idxq(idx(jprev) + 1)
                    idxj = idxq(idx(j) + 1)
                    if (idxjp <= nlp1) then
                       idxjp = idxjp - 1
                    end if
                    if (idxj <= nlp1) then
                       idxj = idxj - 1
                    end if
                    givcol(givptr,2) = idxjp
                    givcol(givptr,1) = idxj
                    givnum(givptr,2) = c
                    givnum(givptr,1) = s
                 end if
                 call la_drot(1,vf(jprev),1,vf(j),1,c,s)
                 call la_drot(1,vl(jprev),1,vl(j),1,c,s)
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 zw(k) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 80
           90 continue
           ! record the last singular value.
           k = k + 1
           zw(k) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           100 continue
           ! sort the singular values into dsigma. the singular values which
           ! were not deflated go into the first k slots of dsigma, except
           ! that dsigma(1) is treated separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              vfw(j) = vf(jp)
              vlw(j) = vl(jp)
           end do
           if (icompq == 1) then
              do j = 2,n
                 jp = idxp(j)
                 perm(j) = idxq(idx(jp) + 1)
                 if (perm(j) <= nlp1) then
                    perm(j) = perm(j) - 1
                 end if
              end do
           end if
           ! the deflated singular values go back into the last n - k slots of
           ! d.
           call la_dcopy(n - k,dsigma(k + 1),1,d(k + 1),1)
           ! determine dsigma(1), dsigma(2), z(1), vf(1), vl(1), vf(m), and
           ! vl(m).
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_dlapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = -z(m)/z(1)
              end if
              call la_drot(1,vf(m),1,vf(1),1,c,s)
              call la_drot(1,vl(m),1,vl(1),1,c,s)
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! restore z, vf, and vl.
           call la_dcopy(k - 1,zw(2),1,z(2),1)
           call la_dcopy(n - 1,vfw(2),1,vf(2),1)
           call la_dcopy(n - 1,vlw(2),1,vl(2),1)
           return
     end subroutine la_dlasd7
#ifdef LA_WITH_XDP
     !> XLASD7: merges the two sets of singular values together into a single
     !> sorted set. Then it tries to deflate the size of the problem. There
     !> are two ways in which deflation can occur:  when two or more singular
     !> values are close together or if there is a tiny entry in the Z
     !> vector. For each such occurrence the order of the related
     !> secular equation problem is reduced by one.
     !> XLASD7 is called from XLASD6.

     pure subroutine la_xlasd7(icompq,nl,nr,sqre,k,d,z,zw,vf,vfw,vl,vlw,alpha, &
     beta,dsigma,idx,idxp,idxq,perm,givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
        use la_constants_xdp,only:zero,one,two,eight

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(xdp),intent(in) :: alpha,beta
           real(xdp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),idx(*),idxp(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(xdp),intent(inout) :: d(*),vf(*),vl(*)
           real(xdp),intent(out) :: dsigma(*),givnum(ldgnum,*),vfw(*),vlw(*),z(*),zw(*)

        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(xdp) :: eps,hlftol,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -22
           else if (ldgnum < n) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('XLASD7',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           if (icompq == 1) then
              givptr = 0
           end if
           ! generate the first part of the vector z and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vl(nlp1)
           vl(nlp1) = zero
           tau = vf(nlp1)
           do i = nl,1,-1
              z(i + 1) = alpha*vl(i)
              vl(i) = zero
              vf(i + 1) = vf(i)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           vf(1) = tau
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vf(i)
              vf(i) = zero
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and zw are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              zw(i) = z(idxq(i))
              vfw(i) = vf(idxq(i))
              vlw(i) = vl(idxq(i))
           end do
           call la_xlamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = zw(idxi)
              vf(i) = vfw(idxi)
              vl(i) = vlw(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_xlamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 if (j == n) go to 100
              else
                 jprev = j
                 go to 70
              end if
           end do
           70 continue
           j = jprev
           80 continue
           j = j + 1
           if (j > n) go to 90
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_xlapy2(c,s)
                 z(j) = tau
                 z(jprev) = zero
                 c = c/tau
                 s = -s/tau
                 ! record the appropriate givens rotation
                 if (icompq == 1) then
                    givptr = givptr + 1
                    idxjp = idxq(idx(jprev) + 1)
                    idxj = idxq(idx(j) + 1)
                    if (idxjp <= nlp1) then
                       idxjp = idxjp - 1
                    end if
                    if (idxj <= nlp1) then
                       idxj = idxj - 1
                    end if
                    givcol(givptr,2) = idxjp
                    givcol(givptr,1) = idxj
                    givnum(givptr,2) = c
                    givnum(givptr,1) = s
                 end if
                 call la_xrot(1,vf(jprev),1,vf(j),1,c,s)
                 call la_xrot(1,vl(jprev),1,vl(j),1,c,s)
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 zw(k) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 80
           90 continue
           ! record the last singular value.
           k = k + 1
           zw(k) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           100 continue
           ! sort the singular values into dsigma. the singular values which
           ! were not deflated go into the first k slots of dsigma, except
           ! that dsigma(1) is treated separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              vfw(j) = vf(jp)
              vlw(j) = vl(jp)
           end do
           if (icompq == 1) then
              do j = 2,n
                 jp = idxp(j)
                 perm(j) = idxq(idx(jp) + 1)
                 if (perm(j) <= nlp1) then
                    perm(j) = perm(j) - 1
                 end if
              end do
           end if
           ! the deflated singular values go back into the last n - k slots of
           ! d.
           call la_xcopy(n - k,dsigma(k + 1),1,d(k + 1),1)
           ! determine dsigma(1), dsigma(2), z(1), vf(1), vl(1), vf(m), and
           ! vl(m).
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_xlapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = -z(m)/z(1)
              end if
              call la_xrot(1,vf(m),1,vf(1),1,c,s)
              call la_xrot(1,vl(m),1,vl(1),1,c,s)
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! restore z, vf, and vl.
           call la_xcopy(k - 1,zw(2),1,z(2),1)
           call la_xcopy(n - 1,vfw(2),1,vf(2),1)
           call la_xcopy(n - 1,vlw(2),1,vl(2),1)
           return
     end subroutine la_xlasd7
#endif
#ifdef LA_WITH_QP
     !> QLASD7: merges the two sets of singular values together into a single
     !> sorted set. Then it tries to deflate the size of the problem. There
     !> are two ways in which deflation can occur:  when two or more singular
     !> values are close together or if there is a tiny entry in the Z
     !> vector. For each such occurrence the order of the related
     !> secular equation problem is reduced by one.
     !> QLASD7 is called from QLASD6.

     pure subroutine la_qlasd7(icompq,nl,nr,sqre,k,d,z,zw,vf,vfw,vl,vlw,alpha, &
     beta,dsigma,idx,idxp,idxq,perm,givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
        use la_constants_qp,only:zero,one,two,eight

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(qp),intent(in) :: alpha,beta
           real(qp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),idx(*),idxp(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(qp),intent(inout) :: d(*),vf(*),vl(*)
           real(qp),intent(out) :: dsigma(*),givnum(ldgnum,*),vfw(*),vlw(*),z(*),zw(*)

        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(qp) :: eps,hlftol,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -22
           else if (ldgnum < n) then
              info = -24
           end if
           if (info /= 0) then
              call la_xerbla('QLASD7',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           if (icompq == 1) then
              givptr = 0
           end if
           ! generate the first part of the vector z and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vl(nlp1)
           vl(nlp1) = zero
           tau = vf(nlp1)
           do i = nl,1,-1
              z(i + 1) = alpha*vl(i)
              vl(i) = zero
              vf(i + 1) = vf(i)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           vf(1) = tau
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vf(i)
              vf(i) = zero
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and zw are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              zw(i) = z(idxq(i))
              vfw(i) = vf(idxq(i))
              vlw(i) = vl(idxq(i))
           end do
           call la_qlamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = zw(idxi)
              vf(i) = vfw(idxi)
              vl(i) = vlw(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_qlamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 if (j == n) go to 100
              else
                 jprev = j
                 go to 70
              end if
           end do
           70 continue
           j = jprev
           80 continue
           j = j + 1
           if (j > n) go to 90
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_qlapy2(c,s)
                 z(j) = tau
                 z(jprev) = zero
                 c = c/tau
                 s = -s/tau
                 ! record the appropriate givens rotation
                 if (icompq == 1) then
                    givptr = givptr + 1
                    idxjp = idxq(idx(jprev) + 1)
                    idxj = idxq(idx(j) + 1)
                    if (idxjp <= nlp1) then
                       idxjp = idxjp - 1
                    end if
                    if (idxj <= nlp1) then
                       idxj = idxj - 1
                    end if
                    givcol(givptr,2) = idxjp
                    givcol(givptr,1) = idxj
                    givnum(givptr,2) = c
                    givnum(givptr,1) = s
                 end if
                 call la_qrot(1,vf(jprev),1,vf(j),1,c,s)
                 call la_qrot(1,vl(jprev),1,vl(j),1,c,s)
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 zw(k) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 80
           90 continue
           ! record the last singular value.
           k = k + 1
           zw(k) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           100 continue
           ! sort the singular values into dsigma. the singular values which
           ! were not deflated go into the first k slots of dsigma, except
           ! that dsigma(1) is treated separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              vfw(j) = vf(jp)
              vlw(j) = vl(jp)
           end do
           if (icompq == 1) then
              do j = 2,n
                 jp = idxp(j)
                 perm(j) = idxq(idx(jp) + 1)
                 if (perm(j) <= nlp1) then
                    perm(j) = perm(j) - 1
                 end if
              end do
           end if
           ! the deflated singular values go back into the last n - k slots of
           ! d.
           call la_qcopy(n - k,dsigma(k + 1),1,d(k + 1),1)
           ! determine dsigma(1), dsigma(2), z(1), vf(1), vl(1), vf(m), and
           ! vl(m).
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_qlapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = -z(m)/z(1)
              end if
              call la_qrot(1,vf(m),1,vf(1),1,c,s)
              call la_qrot(1,vl(m),1,vl(1),1,c,s)
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! restore z, vf, and vl.
           call la_qcopy(k - 1,zw(2),1,z(2),1)
           call la_qcopy(n - 1,vfw(2),1,vf(2),1)
           call la_qcopy(n - 1,vlw(2),1,vl(2),1)
           return
     end subroutine la_qlasd7
#endif

     !> SLASD8: finds the square roots of the roots of the secular equation,
     !> as defined by the values in DSIGMA and Z. It makes the appropriate
     !> calls to SLASD4, and stores, for each  element in D, the distance
     !> to its two nearest poles (elements in DSIGMA). It also updates
     !> the arrays VF and VL, the first and last components of all the
     !> right singular vectors of the original bidiagonal matrix.
     !> SLASD8 is called from SLASD6.

     pure subroutine la_slasd8(icompq,k,d,z,vf,vl,difl,difr,lddifr,dsigma,work, &
               info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,k,lddifr
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(out) :: d(*),difl(*),difr(lddifr,*),work(*)
           real(sp),intent(inout) :: dsigma(*),vf(*),vl(*),z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iwk1,iwk2,iwk2i,iwk3,iwk3i,j
           real(sp) :: diflj,difrj,dj,dsigj,dsigjp,rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (k < 1) then
              info = -2
           else if (lddifr < k) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SLASD8',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              difl(1) = d(1)
              if (icompq == 1) then
                 difl(2) = one
                 difr(1,2) = one
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_slamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! book keeping.
           iwk1 = 1
           iwk2 = iwk1 + k
           iwk3 = iwk2 + k
           iwk2i = iwk2 - 1
           iwk3i = iwk3 - 1
           ! normalize z.
           rho = la_snrm2(k,z,1)
           call la_slascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! initialize work(iwk3).
           call la_slaset('A',k,1,one,one,work(iwk3),k)
           ! compute the updated singular values, the arrays difl, difr,
           ! and the updated z.
           do j = 1,k
              call la_slasd4(k,j,dsigma,z,work(iwk1),rho,d(j),work(iwk2),info)

              ! if the root finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
              work(iwk3i + j) = work(iwk3i + j)*work(j)*work(iwk2i + j)
              difl(j) = -work(j)
              difr(j,1) = -work(j + 1)
              do i = 1,j - 1
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
              do i = j + 1,k
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = sign(sqrt(abs(work(iwk3i + i))),z(i))
           end do
           ! update vf and vl.
           do j = 1,k
              diflj = difl(j)
              dj = d(j)
              dsigj = -dsigma(j)
              if (j < k) then
                 difrj = -difr(j,1)
                 dsigjp = -dsigma(j + 1)
              end if
              work(j) = -z(j)/diflj/(dsigma(j) + dj)
              do i = 1,j - 1
                 work(i) = z(i)/(la_slamc3(dsigma(i),dsigj) - diflj)/(dsigma(i) &
                           + dj)
              end do
              do i = j + 1,k
                 work(i) = z(i)/(la_slamc3(dsigma(i),dsigjp) + difrj)/(dsigma(i &
                           ) + dj)
              end do
              temp = la_snrm2(k,work,1)
              work(iwk2i + j) = la_sdot(k,work,1,vf,1)/temp
              work(iwk3i + j) = la_sdot(k,work,1,vl,1)/temp
              if (icompq == 1) then
                 difr(j,2) = temp
              end if
           end do
           call la_scopy(k,work(iwk2),1,vf,1)
           call la_scopy(k,work(iwk3),1,vl,1)
           return
     end subroutine la_slasd8
     !> DLASD8: finds the square roots of the roots of the secular equation,
     !> as defined by the values in DSIGMA and Z. It makes the appropriate
     !> calls to DLASD4, and stores, for each  element in D, the distance
     !> to its two nearest poles (elements in DSIGMA). It also updates
     !> the arrays VF and VL, the first and last components of all the
     !> right singular vectors of the original bidiagonal matrix.
     !> DLASD8 is called from DLASD6.

     pure subroutine la_dlasd8(icompq,k,d,z,vf,vl,difl,difr,lddifr,dsigma,work, &
               info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,k,lddifr
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(out) :: d(*),difl(*),difr(lddifr,*),work(*)
           real(dp),intent(inout) :: dsigma(*),vf(*),vl(*),z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iwk1,iwk2,iwk2i,iwk3,iwk3i,j
           real(dp) :: diflj,difrj,dj,dsigj,dsigjp,rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (k < 1) then
              info = -2
           else if (lddifr < k) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DLASD8',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              difl(1) = d(1)
              if (icompq == 1) then
                 difl(2) = one
                 difr(1,2) = one
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_dlamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! book keeping.
           iwk1 = 1
           iwk2 = iwk1 + k
           iwk3 = iwk2 + k
           iwk2i = iwk2 - 1
           iwk3i = iwk3 - 1
           ! normalize z.
           rho = la_dnrm2(k,z,1)
           call la_dlascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! initialize work(iwk3).
           call la_dlaset('A',k,1,one,one,work(iwk3),k)
           ! compute the updated singular values, the arrays difl, difr,
           ! and the updated z.
           do j = 1,k
              call la_dlasd4(k,j,dsigma,z,work(iwk1),rho,d(j),work(iwk2),info)

              ! if the root finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
              work(iwk3i + j) = work(iwk3i + j)*work(j)*work(iwk2i + j)
              difl(j) = -work(j)
              difr(j,1) = -work(j + 1)
              do i = 1,j - 1
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
              do i = j + 1,k
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = sign(sqrt(abs(work(iwk3i + i))),z(i))
           end do
           ! update vf and vl.
           do j = 1,k
              diflj = difl(j)
              dj = d(j)
              dsigj = -dsigma(j)
              if (j < k) then
                 difrj = -difr(j,1)
                 dsigjp = -dsigma(j + 1)
              end if
              work(j) = -z(j)/diflj/(dsigma(j) + dj)
              do i = 1,j - 1
                 work(i) = z(i)/(la_dlamc3(dsigma(i),dsigj) - diflj)/(dsigma(i) &
                           + dj)
              end do
              do i = j + 1,k
                 work(i) = z(i)/(la_dlamc3(dsigma(i),dsigjp) + difrj)/(dsigma(i &
                           ) + dj)
              end do
              temp = la_dnrm2(k,work,1)
              work(iwk2i + j) = la_ddot(k,work,1,vf,1)/temp
              work(iwk3i + j) = la_ddot(k,work,1,vl,1)/temp
              if (icompq == 1) then
                 difr(j,2) = temp
              end if
           end do
           call la_dcopy(k,work(iwk2),1,vf,1)
           call la_dcopy(k,work(iwk3),1,vl,1)
           return
     end subroutine la_dlasd8
#ifdef LA_WITH_XDP
     !> XLASD8: finds the square roots of the roots of the secular equation,
     !> as defined by the values in DSIGMA and Z. It makes the appropriate
     !> calls to XLASD4, and stores, for each  element in D, the distance
     !> to its two nearest poles (elements in DSIGMA). It also updates
     !> the arrays VF and VL, the first and last components of all the
     !> right singular vectors of the original bidiagonal matrix.
     !> XLASD8 is called from XLASD6.

     pure subroutine la_xlasd8(icompq,k,d,z,vf,vl,difl,difr,lddifr,dsigma,work, &
               info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,k,lddifr
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(xdp),intent(out) :: d(*),difl(*),difr(lddifr,*),work(*)
           real(xdp),intent(inout) :: dsigma(*),vf(*),vl(*),z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iwk1,iwk2,iwk2i,iwk3,iwk3i,j
           real(xdp) :: diflj,difrj,dj,dsigj,dsigjp,rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (k < 1) then
              info = -2
           else if (lddifr < k) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XLASD8',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              difl(1) = d(1)
              if (icompq == 1) then
                 difl(2) = one
                 difr(1,2) = one
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_xlamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! book keeping.
           iwk1 = 1
           iwk2 = iwk1 + k
           iwk3 = iwk2 + k
           iwk2i = iwk2 - 1
           iwk3i = iwk3 - 1
           ! normalize z.
           rho = la_xnrm2(k,z,1)
           call la_xlascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! initialize work(iwk3).
           call la_xlaset('A',k,1,one,one,work(iwk3),k)
           ! compute the updated singular values, the arrays difl, difr,
           ! and the updated z.
           do j = 1,k
              call la_xlasd4(k,j,dsigma,z,work(iwk1),rho,d(j),work(iwk2),info)

              ! if the root finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
              work(iwk3i + j) = work(iwk3i + j)*work(j)*work(iwk2i + j)
              difl(j) = -work(j)
              difr(j,1) = -work(j + 1)
              do i = 1,j - 1
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
              do i = j + 1,k
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = sign(sqrt(abs(work(iwk3i + i))),z(i))
           end do
           ! update vf and vl.
           do j = 1,k
              diflj = difl(j)
              dj = d(j)
              dsigj = -dsigma(j)
              if (j < k) then
                 difrj = -difr(j,1)
                 dsigjp = -dsigma(j + 1)
              end if
              work(j) = -z(j)/diflj/(dsigma(j) + dj)
              do i = 1,j - 1
                 work(i) = z(i)/(la_xlamc3(dsigma(i),dsigj) - diflj)/(dsigma(i) &
                           + dj)
              end do
              do i = j + 1,k
                 work(i) = z(i)/(la_xlamc3(dsigma(i),dsigjp) + difrj)/(dsigma(i &
                           ) + dj)
              end do
              temp = la_xnrm2(k,work,1)
              work(iwk2i + j) = la_xdot(k,work,1,vf,1)/temp
              work(iwk3i + j) = la_xdot(k,work,1,vl,1)/temp
              if (icompq == 1) then
                 difr(j,2) = temp
              end if
           end do
           call la_xcopy(k,work(iwk2),1,vf,1)
           call la_xcopy(k,work(iwk3),1,vl,1)
           return
     end subroutine la_xlasd8
#endif
#ifdef LA_WITH_QP
     !> QLASD8: finds the square roots of the roots of the secular equation,
     !> as defined by the values in DSIGMA and Z. It makes the appropriate
     !> calls to QLASD4, and stores, for each  element in D, the distance
     !> to its two nearest poles (elements in DSIGMA). It also updates
     !> the arrays VF and VL, the first and last components of all the
     !> right singular vectors of the original bidiagonal matrix.
     !> QLASD8 is called from QLASD6.

     pure subroutine la_qlasd8(icompq,k,d,z,vf,vl,difl,difr,lddifr,dsigma,work, &
               info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,k,lddifr
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(out) :: d(*),difl(*),difr(lddifr,*),work(*)
           real(qp),intent(inout) :: dsigma(*),vf(*),vl(*),z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iwk1,iwk2,iwk2i,iwk3,iwk3i,j
           real(qp) :: diflj,difrj,dj,dsigj,dsigjp,rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (k < 1) then
              info = -2
           else if (lddifr < k) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QLASD8',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              difl(1) = d(1)
              if (icompq == 1) then
                 difl(2) = one
                 difr(1,2) = one
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dlambda(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_qlamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! book keeping.
           iwk1 = 1
           iwk2 = iwk1 + k
           iwk3 = iwk2 + k
           iwk2i = iwk2 - 1
           iwk3i = iwk3 - 1
           ! normalize z.
           rho = la_qnrm2(k,z,1)
           call la_qlascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! initialize work(iwk3).
           call la_qlaset('A',k,1,one,one,work(iwk3),k)
           ! compute the updated singular values, the arrays difl, difr,
           ! and the updated z.
           do j = 1,k
              call la_qlasd4(k,j,dsigma,z,work(iwk1),rho,d(j),work(iwk2),info)

              ! if the root finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
              work(iwk3i + j) = work(iwk3i + j)*work(j)*work(iwk2i + j)
              difl(j) = -work(j)
              difr(j,1) = -work(j + 1)
              do i = 1,j - 1
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
              do i = j + 1,k
                 work(iwk3i + i) = work(iwk3i + i)*work(i)*work(iwk2i + i)/(dsigma(i) - &
                           dsigma(j))/(dsigma(i) + dsigma(j))
              end do
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = sign(sqrt(abs(work(iwk3i + i))),z(i))
           end do
           ! update vf and vl.
           do j = 1,k
              diflj = difl(j)
              dj = d(j)
              dsigj = -dsigma(j)
              if (j < k) then
                 difrj = -difr(j,1)
                 dsigjp = -dsigma(j + 1)
              end if
              work(j) = -z(j)/diflj/(dsigma(j) + dj)
              do i = 1,j - 1
                 work(i) = z(i)/(la_qlamc3(dsigma(i),dsigj) - diflj)/(dsigma(i) &
                           + dj)
              end do
              do i = j + 1,k
                 work(i) = z(i)/(la_qlamc3(dsigma(i),dsigjp) + difrj)/(dsigma(i &
                           ) + dj)
              end do
              temp = la_qnrm2(k,work,1)
              work(iwk2i + j) = la_qdot(k,work,1,vf,1)/temp
              work(iwk3i + j) = la_qdot(k,work,1,vl,1)/temp
              if (icompq == 1) then
                 difr(j,2) = temp
              end if
           end do
           call la_qcopy(k,work(iwk2),1,vf,1)
           call la_qcopy(k,work(iwk3),1,vl,1)
           return
     end subroutine la_qlasd8
#endif

     !> SLASD3: finds all the square roots of the roots of the secular
     !> equation, as defined by the values in D and Z.  It makes the
     !> appropriate calls to SLASD4 and then updates the singular
     !> vectors by matrix multiplication.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.
     !> SLASD3 is called from SLASD1.

     pure subroutine la_slasd3(nl,nr,sqre,k,d,q,ldq,dsigma,u,ldu,u2,ldu2,vt,ldvt, &
                vt2,ldvt2,idxc,ctot,z,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),idxc(*)
           real(sp),intent(out) :: d(*),q(ldq,*),u(ldu,*),vt(ldvt,*)
           real(sp),intent(inout) :: dsigma(*),vt2(ldvt2,*),z(*)
           real(sp),intent(in) :: u2(ldu2,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: ctemp,i,j,jc,ktemp,m,n,nlp1,nlp2,nrp1
           real(sp) :: rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           nlp1 = nl + 1
           nlp2 = nl + 2
           if ((k < 1) .or. (k > n)) then
              info = -4
           else if (ldq < k) then
              info = -7
           else if (ldu < n) then
              info = -10
           else if (ldu2 < n) then
              info = -12
           else if (ldvt < m) then
              info = -14
           else if (ldvt2 < m) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('SLASD3',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              call la_scopy(m,vt2(1,1),ldvt2,vt(1,1),ldvt)
              if (z(1) > zero) then
                 call la_scopy(n,u2(1,1),1,u(1,1),1)
              else
                 do i = 1,n
                    u(i,1) = -u2(i,1)
                 end do
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dsigma(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_slamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! keep a copy of z.
           call la_scopy(k,z,1,q,1)
           ! normalize z.
           rho = la_snrm2(k,z,1)
           call la_slascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! find the new singular values.
           do j = 1,k
              call la_slasd4(k,j,dsigma,z,u(1,j),rho,d(j),vt(1,j),info)

              ! if the zero finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = u(i,k)*vt(i,k)
              do j = 1,i - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j))/(dsigma(i &
                           ) + dsigma(j)))
              end do
              do j = i,k - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j + 1))/(dsigma( &
                           i) + dsigma(j + 1)))
              end do
              z(i) = sign(sqrt(abs(z(i))),q(i,1))
           end do
           ! compute left singular vectors of the modified diagonal matrix,
           ! and store related information for the right singular vectors.
           do i = 1,k
              vt(1,i) = z(1)/u(1,i)/vt(1,i)
              u(1,i) = negone
              do j = 2,k
                 vt(j,i) = z(j)/u(j,i)/vt(j,i)
                 u(j,i) = dsigma(j)*vt(j,i)
              end do
              temp = la_snrm2(k,u(1,i),1)
              q(1,i) = u(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(j,i) = u(jc,i)/temp
              end do
           end do
           ! update the left singular vector matrix.
           if (k == 2) then
              call la_sgemm('N','N',n,k,k,one,u2,ldu2,q,ldq,zero,u,ldu)
              go to 100
           end if
           if (ctot(1) > 0) then
              call la_sgemm('N','N',nl,k,ctot(1),one,u2(1,2),ldu2,q(2,1),ldq, &
                         zero,u(1,1),ldu)
              if (ctot(3) > 0) then
                 ktemp = 2 + ctot(1) + ctot(2)
                 call la_sgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q( &
                           ktemp,1),ldq,one,u(1,1),ldu)
              end if
           else if (ctot(3) > 0) then
              ktemp = 2 + ctot(1) + ctot(2)
              call la_sgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q(ktemp, &
                        1),ldq,zero,u(1,1),ldu)
           else
              call la_slacpy('F',nl,k,u2,ldu2,u,ldu)
           end if
           call la_scopy(k,q(1,1),ldq,u(nlp1,1),ldu)
           ktemp = 2 + ctot(1)
           ctemp = ctot(2) + ctot(3)
           call la_sgemm('N','N',nr,k,ctemp,one,u2(nlp2,ktemp),ldu2,q(ktemp,1), &
                     ldq,zero,u(nlp2,1),ldu)
           ! generate the right singular vectors.
           100 continue
           do i = 1,k
              temp = la_snrm2(k,vt(1,i),1)
              q(i,1) = vt(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(i,j) = vt(jc,i)/temp
              end do
           end do
           ! update the right singular vector matrix.
           if (k == 2) then
              call la_sgemm('N','N',k,m,k,one,q,ldq,vt2,ldvt2,zero,vt,ldvt)

              return
           end if
           ktemp = 1 + ctot(1)
           call la_sgemm('N','N',k,nlp1,ktemp,one,q(1,1),ldq,vt2(1,1),ldvt2, &
                     zero,vt(1,1),ldvt)
           ktemp = 2 + ctot(1) + ctot(2)
           if (ktemp <= ldvt2) call la_sgemm('N','N',k,nlp1,ctot(3),one,q(1,ktemp), &
                     ldq,vt2(ktemp,1),ldvt2,one,vt(1,1),ldvt)
           ktemp = ctot(1) + 1
           nrp1 = nr + sqre
           if (ktemp > 1) then
              do i = 1,k
                 q(i,ktemp) = q(i,1)
              end do
              do i = nlp2,m
                 vt2(ktemp,i) = vt2(1,i)
              end do
           end if
           ctemp = 1 + ctot(2) + ctot(3)
           call la_sgemm('N','N',k,nrp1,ctemp,one,q(1,ktemp),ldq,vt2(ktemp,nlp2) &
                     ,ldvt2,zero,vt(1,nlp2),ldvt)
           return
     end subroutine la_slasd3
     !> DLASD3: finds all the square roots of the roots of the secular
     !> equation, as defined by the values in D and Z.  It makes the
     !> appropriate calls to DLASD4 and then updates the singular
     !> vectors by matrix multiplication.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.
     !> DLASD3 is called from DLASD1.

     pure subroutine la_dlasd3(nl,nr,sqre,k,d,q,ldq,dsigma,u,ldu,u2,ldu2,vt,ldvt, &
                vt2,ldvt2,idxc,ctot,z,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),idxc(*)
           real(dp),intent(out) :: d(*),q(ldq,*),u(ldu,*),vt(ldvt,*)
           real(dp),intent(inout) :: dsigma(*),vt2(ldvt2,*),z(*)
           real(dp),intent(in) :: u2(ldu2,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: ctemp,i,j,jc,ktemp,m,n,nlp1,nlp2,nrp1
           real(dp) :: rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           nlp1 = nl + 1
           nlp2 = nl + 2
           if ((k < 1) .or. (k > n)) then
              info = -4
           else if (ldq < k) then
              info = -7
           else if (ldu < n) then
              info = -10
           else if (ldu2 < n) then
              info = -12
           else if (ldvt < m) then
              info = -14
           else if (ldvt2 < m) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('DLASD3',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              call la_dcopy(m,vt2(1,1),ldvt2,vt(1,1),ldvt)
              if (z(1) > zero) then
                 call la_dcopy(n,u2(1,1),1,u(1,1),1)
              else
                 do i = 1,n
                    u(i,1) = -u2(i,1)
                 end do
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dsigma(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_dlamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! keep a copy of z.
           call la_dcopy(k,z,1,q,1)
           ! normalize z.
           rho = la_dnrm2(k,z,1)
           call la_dlascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! find the new singular values.
           do j = 1,k
              call la_dlasd4(k,j,dsigma,z,u(1,j),rho,d(j),vt(1,j),info)

              ! if the zero finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = u(i,k)*vt(i,k)
              do j = 1,i - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j))/(dsigma(i &
                           ) + dsigma(j)))
              end do
              do j = i,k - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j + 1))/(dsigma( &
                           i) + dsigma(j + 1)))
              end do
              z(i) = sign(sqrt(abs(z(i))),q(i,1))
           end do
           ! compute left singular vectors of the modified diagonal matrix,
           ! and store related information for the right singular vectors.
           do i = 1,k
              vt(1,i) = z(1)/u(1,i)/vt(1,i)
              u(1,i) = negone
              do j = 2,k
                 vt(j,i) = z(j)/u(j,i)/vt(j,i)
                 u(j,i) = dsigma(j)*vt(j,i)
              end do
              temp = la_dnrm2(k,u(1,i),1)
              q(1,i) = u(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(j,i) = u(jc,i)/temp
              end do
           end do
           ! update the left singular vector matrix.
           if (k == 2) then
              call la_dgemm('N','N',n,k,k,one,u2,ldu2,q,ldq,zero,u,ldu)
              go to 100
           end if
           if (ctot(1) > 0) then
              call la_dgemm('N','N',nl,k,ctot(1),one,u2(1,2),ldu2,q(2,1),ldq, &
                         zero,u(1,1),ldu)
              if (ctot(3) > 0) then
                 ktemp = 2 + ctot(1) + ctot(2)
                 call la_dgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q( &
                           ktemp,1),ldq,one,u(1,1),ldu)
              end if
           else if (ctot(3) > 0) then
              ktemp = 2 + ctot(1) + ctot(2)
              call la_dgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q(ktemp, &
                        1),ldq,zero,u(1,1),ldu)
           else
              call la_dlacpy('F',nl,k,u2,ldu2,u,ldu)
           end if
           call la_dcopy(k,q(1,1),ldq,u(nlp1,1),ldu)
           ktemp = 2 + ctot(1)
           ctemp = ctot(2) + ctot(3)
           call la_dgemm('N','N',nr,k,ctemp,one,u2(nlp2,ktemp),ldu2,q(ktemp,1), &
                     ldq,zero,u(nlp2,1),ldu)
           ! generate the right singular vectors.
           100 continue
           do i = 1,k
              temp = la_dnrm2(k,vt(1,i),1)
              q(i,1) = vt(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(i,j) = vt(jc,i)/temp
              end do
           end do
           ! update the right singular vector matrix.
           if (k == 2) then
              call la_dgemm('N','N',k,m,k,one,q,ldq,vt2,ldvt2,zero,vt,ldvt)

              return
           end if
           ktemp = 1 + ctot(1)
           call la_dgemm('N','N',k,nlp1,ktemp,one,q(1,1),ldq,vt2(1,1),ldvt2, &
                     zero,vt(1,1),ldvt)
           ktemp = 2 + ctot(1) + ctot(2)
           if (ktemp <= ldvt2) call la_dgemm('N','N',k,nlp1,ctot(3),one,q(1,ktemp), &
                     ldq,vt2(ktemp,1),ldvt2,one,vt(1,1),ldvt)
           ktemp = ctot(1) + 1
           nrp1 = nr + sqre
           if (ktemp > 1) then
              do i = 1,k
                 q(i,ktemp) = q(i,1)
              end do
              do i = nlp2,m
                 vt2(ktemp,i) = vt2(1,i)
              end do
           end if
           ctemp = 1 + ctot(2) + ctot(3)
           call la_dgemm('N','N',k,nrp1,ctemp,one,q(1,ktemp),ldq,vt2(ktemp,nlp2) &
                     ,ldvt2,zero,vt(1,nlp2),ldvt)
           return
     end subroutine la_dlasd3
#ifdef LA_WITH_XDP
     !> XLASD3: finds all the square roots of the roots of the secular
     !> equation, as defined by the values in D and Z.  It makes the
     !> appropriate calls to XLASD4 and then updates the singular
     !> vectors by matrix multiplication.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.
     !> XLASD3 is called from XLASD1.

     pure subroutine la_xlasd3(nl,nr,sqre,k,d,q,ldq,dsigma,u,ldu,u2,ldu2,vt,ldvt, &
                vt2,ldvt2,idxc,ctot,z,info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),idxc(*)
           real(xdp),intent(out) :: d(*),q(ldq,*),u(ldu,*),vt(ldvt,*)
           real(xdp),intent(inout) :: dsigma(*),vt2(ldvt2,*),z(*)
           real(xdp),intent(in) :: u2(ldu2,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: ctemp,i,j,jc,ktemp,m,n,nlp1,nlp2,nrp1
           real(xdp) :: rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           nlp1 = nl + 1
           nlp2 = nl + 2
           if ((k < 1) .or. (k > n)) then
              info = -4
           else if (ldq < k) then
              info = -7
           else if (ldu < n) then
              info = -10
           else if (ldu2 < n) then
              info = -12
           else if (ldvt < m) then
              info = -14
           else if (ldvt2 < m) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('XLASD3',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              call la_xcopy(m,vt2(1,1),ldvt2,vt(1,1),ldvt)
              if (z(1) > zero) then
                 call la_xcopy(n,u2(1,1),1,u(1,1),1)
              else
                 do i = 1,n
                    u(i,1) = -u2(i,1)
                 end do
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dsigma(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_xlamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! keep a copy of z.
           call la_xcopy(k,z,1,q,1)
           ! normalize z.
           rho = la_xnrm2(k,z,1)
           call la_xlascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! find the new singular values.
           do j = 1,k
              call la_xlasd4(k,j,dsigma,z,u(1,j),rho,d(j),vt(1,j),info)

              ! if the zero finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = u(i,k)*vt(i,k)
              do j = 1,i - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j))/(dsigma(i &
                           ) + dsigma(j)))
              end do
              do j = i,k - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j + 1))/(dsigma( &
                           i) + dsigma(j + 1)))
              end do
              z(i) = sign(sqrt(abs(z(i))),q(i,1))
           end do
           ! compute left singular vectors of the modified diagonal matrix,
           ! and store related information for the right singular vectors.
           do i = 1,k
              vt(1,i) = z(1)/u(1,i)/vt(1,i)
              u(1,i) = negone
              do j = 2,k
                 vt(j,i) = z(j)/u(j,i)/vt(j,i)
                 u(j,i) = dsigma(j)*vt(j,i)
              end do
              temp = la_xnrm2(k,u(1,i),1)
              q(1,i) = u(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(j,i) = u(jc,i)/temp
              end do
           end do
           ! update the left singular vector matrix.
           if (k == 2) then
              call la_xgemm('N','N',n,k,k,one,u2,ldu2,q,ldq,zero,u,ldu)
              go to 100
           end if
           if (ctot(1) > 0) then
              call la_xgemm('N','N',nl,k,ctot(1),one,u2(1,2),ldu2,q(2,1),ldq, &
                         zero,u(1,1),ldu)
              if (ctot(3) > 0) then
                 ktemp = 2 + ctot(1) + ctot(2)
                 call la_xgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q( &
                           ktemp,1),ldq,one,u(1,1),ldu)
              end if
           else if (ctot(3) > 0) then
              ktemp = 2 + ctot(1) + ctot(2)
              call la_xgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q(ktemp, &
                        1),ldq,zero,u(1,1),ldu)
           else
              call la_xlacpy('F',nl,k,u2,ldu2,u,ldu)
           end if
           call la_xcopy(k,q(1,1),ldq,u(nlp1,1),ldu)
           ktemp = 2 + ctot(1)
           ctemp = ctot(2) + ctot(3)
           call la_xgemm('N','N',nr,k,ctemp,one,u2(nlp2,ktemp),ldu2,q(ktemp,1), &
                     ldq,zero,u(nlp2,1),ldu)
           ! generate the right singular vectors.
           100 continue
           do i = 1,k
              temp = la_xnrm2(k,vt(1,i),1)
              q(i,1) = vt(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(i,j) = vt(jc,i)/temp
              end do
           end do
           ! update the right singular vector matrix.
           if (k == 2) then
              call la_xgemm('N','N',k,m,k,one,q,ldq,vt2,ldvt2,zero,vt,ldvt)

              return
           end if
           ktemp = 1 + ctot(1)
           call la_xgemm('N','N',k,nlp1,ktemp,one,q(1,1),ldq,vt2(1,1),ldvt2, &
                     zero,vt(1,1),ldvt)
           ktemp = 2 + ctot(1) + ctot(2)
           if (ktemp <= ldvt2) call la_xgemm('N','N',k,nlp1,ctot(3),one,q(1,ktemp), &
                     ldq,vt2(ktemp,1),ldvt2,one,vt(1,1),ldvt)
           ktemp = ctot(1) + 1
           nrp1 = nr + sqre
           if (ktemp > 1) then
              do i = 1,k
                 q(i,ktemp) = q(i,1)
              end do
              do i = nlp2,m
                 vt2(ktemp,i) = vt2(1,i)
              end do
           end if
           ctemp = 1 + ctot(2) + ctot(3)
           call la_xgemm('N','N',k,nrp1,ctemp,one,q(1,ktemp),ldq,vt2(ktemp,nlp2) &
                     ,ldvt2,zero,vt(1,nlp2),ldvt)
           return
     end subroutine la_xlasd3
#endif
#ifdef LA_WITH_QP
     !> QLASD3: finds all the square roots of the roots of the secular
     !> equation, as defined by the values in D and Z.  It makes the
     !> appropriate calls to QLASD4 and then updates the singular
     !> vectors by matrix multiplication.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.
     !> QLASD3 is called from QLASD1.

     pure subroutine la_qlasd3(nl,nr,sqre,k,d,q,ldq,dsigma,u,ldu,u2,ldu2,vt,ldvt, &
                vt2,ldvt2,idxc,ctot,z,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,ldq,ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           ! Array Arguments
           integer(ilp),intent(in) :: ctot(*),idxc(*)
           real(qp),intent(out) :: d(*),q(ldq,*),u(ldu,*),vt(ldvt,*)
           real(qp),intent(inout) :: dsigma(*),vt2(ldvt2,*),z(*)
           real(qp),intent(in) :: u2(ldu2,*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: ctemp,i,j,jc,ktemp,m,n,nlp1,nlp2,nrp1
           real(qp) :: rho,temp
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           nlp1 = nl + 1
           nlp2 = nl + 2
           if ((k < 1) .or. (k > n)) then
              info = -4
           else if (ldq < k) then
              info = -7
           else if (ldu < n) then
              info = -10
           else if (ldu2 < n) then
              info = -12
           else if (ldvt < m) then
              info = -14
           else if (ldvt2 < m) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('QLASD3',-info)
              return
           end if
           ! quick return if possible
           if (k == 1) then
              d(1) = abs(z(1))
              call la_qcopy(m,vt2(1,1),ldvt2,vt(1,1),ldvt)
              if (z(1) > zero) then
                 call la_qcopy(n,u2(1,1),1,u(1,1),1)
              else
                 do i = 1,n
                    u(i,1) = -u2(i,1)
                 end do
              end if
              return
           end if
           ! modify values dsigma(i) to make sure all dsigma(i)-dsigma(j) can
           ! be computed with high relative accuracy (barring over/underflow).
           ! this is a problem on machines without a guard digit in
           ! add/subtract (cray xmp, cray ymp, cray c 90 and cray 2).
           ! the following code replaces dsigma(i) by 2*dsigma(i)-dsigma(i),
           ! which on any of these machines zeros out the bottommost
           ! bit of dsigma(i) if it is 1; this makes the subsequent
           ! subtractions dsigma(i)-dsigma(j) unproblematic when cancellation
           ! occurs. on binary machines with a guard digit (almost all
           ! machines) it does not change dsigma(i) at all. on hexadecimal
           ! and decimal machines with a guard digit, it slightly
           ! changes the bottommost bits of dsigma(i). it does not account
           ! for hexadecimal or decimal machines without guard digits
           ! (we know of none). we use a subroutine call to compute
           ! 2*dsigma(i) to prevent optimizing compilers from eliminating
           ! this code.
           do i = 1,k
              dsigma(i) = la_qlamc3(dsigma(i),dsigma(i)) - dsigma(i)
           end do
           ! keep a copy of z.
           call la_qcopy(k,z,1,q,1)
           ! normalize z.
           rho = la_qnrm2(k,z,1)
           call la_qlascl('G',0,0,rho,one,k,1,z,k,info)
           rho = rho*rho
           ! find the new singular values.
           do j = 1,k
              call la_qlasd4(k,j,dsigma,z,u(1,j),rho,d(j),vt(1,j),info)

              ! if the zero finder fails, report the convergence failure.
              if (info /= 0) then
                 return
              end if
           end do
           ! compute updated z.
           do i = 1,k
              z(i) = u(i,k)*vt(i,k)
              do j = 1,i - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j))/(dsigma(i &
                           ) + dsigma(j)))
              end do
              do j = i,k - 1
                 z(i) = z(i)*(u(i,j)*vt(i,j)/(dsigma(i) - dsigma(j + 1))/(dsigma( &
                           i) + dsigma(j + 1)))
              end do
              z(i) = sign(sqrt(abs(z(i))),q(i,1))
           end do
           ! compute left singular vectors of the modified diagonal matrix,
           ! and store related information for the right singular vectors.
           do i = 1,k
              vt(1,i) = z(1)/u(1,i)/vt(1,i)
              u(1,i) = negone
              do j = 2,k
                 vt(j,i) = z(j)/u(j,i)/vt(j,i)
                 u(j,i) = dsigma(j)*vt(j,i)
              end do
              temp = la_qnrm2(k,u(1,i),1)
              q(1,i) = u(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(j,i) = u(jc,i)/temp
              end do
           end do
           ! update the left singular vector matrix.
           if (k == 2) then
              call la_qgemm('N','N',n,k,k,one,u2,ldu2,q,ldq,zero,u,ldu)
              go to 100
           end if
           if (ctot(1) > 0) then
              call la_qgemm('N','N',nl,k,ctot(1),one,u2(1,2),ldu2,q(2,1),ldq, &
                         zero,u(1,1),ldu)
              if (ctot(3) > 0) then
                 ktemp = 2 + ctot(1) + ctot(2)
                 call la_qgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q( &
                           ktemp,1),ldq,one,u(1,1),ldu)
              end if
           else if (ctot(3) > 0) then
              ktemp = 2 + ctot(1) + ctot(2)
              call la_qgemm('N','N',nl,k,ctot(3),one,u2(1,ktemp),ldu2,q(ktemp, &
                        1),ldq,zero,u(1,1),ldu)
           else
              call la_qlacpy('F',nl,k,u2,ldu2,u,ldu)
           end if
           call la_qcopy(k,q(1,1),ldq,u(nlp1,1),ldu)
           ktemp = 2 + ctot(1)
           ctemp = ctot(2) + ctot(3)
           call la_qgemm('N','N',nr,k,ctemp,one,u2(nlp2,ktemp),ldu2,q(ktemp,1), &
                     ldq,zero,u(nlp2,1),ldu)
           ! generate the right singular vectors.
           100 continue
           do i = 1,k
              temp = la_qnrm2(k,vt(1,i),1)
              q(i,1) = vt(1,i)/temp
              do j = 2,k
                 jc = idxc(j)
                 q(i,j) = vt(jc,i)/temp
              end do
           end do
           ! update the right singular vector matrix.
           if (k == 2) then
              call la_qgemm('N','N',k,m,k,one,q,ldq,vt2,ldvt2,zero,vt,ldvt)

              return
           end if
           ktemp = 1 + ctot(1)
           call la_qgemm('N','N',k,nlp1,ktemp,one,q(1,1),ldq,vt2(1,1),ldvt2, &
                     zero,vt(1,1),ldvt)
           ktemp = 2 + ctot(1) + ctot(2)
           if (ktemp <= ldvt2) call la_qgemm('N','N',k,nlp1,ctot(3),one,q(1,ktemp), &
                     ldq,vt2(ktemp,1),ldvt2,one,vt(1,1),ldvt)
           ktemp = ctot(1) + 1
           nrp1 = nr + sqre
           if (ktemp > 1) then
              do i = 1,k
                 q(i,ktemp) = q(i,1)
              end do
              do i = nlp2,m
                 vt2(ktemp,i) = vt2(1,i)
              end do
           end if
           ctemp = 1 + ctot(2) + ctot(3)
           call la_qgemm('N','N',k,nrp1,ctemp,one,q(1,ktemp),ldq,vt2(ktemp,nlp2) &
                     ,ldvt2,zero,vt(1,nlp2),ldvt)
           return
     end subroutine la_qlasd3
#endif

     !> SLASD6: computes the SVD of an updated upper bidiagonal matrix B
     !> obtained by merging two smaller ones by appending a row. This
     !> routine is used only for the problem which requires all singular
     !> values and optionally singular vector matrices in factored form.
     !> B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
     !> A related subroutine, SLASD1, handles the case in which all singular
     !> values and singular vectors of the bidiagonal matrix are desired.
     !> SLASD6 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The singular values of B can be computed using D1, D2, the first
     !> components of all the right singular vectors of the lower block, and
     !> the last components of all the right singular vectors of the upper
     !> block. These components are stored and updated in VF and VL,
     !> respectively, in SLASD6. Hence U and VT are not explicitly
     !> referenced.
     !> The singular values are stored in D. The algorithm consists of two
     !> stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or if there is a zero
     !> in the Z vector. For each such occurrence the dimension of the
     !> secular equation problem is reduced by one. This stage is
     !> performed by the routine SLASD7.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the roots of the
     !> secular equation via the routine SLASD4 (as called by SLASD8).
     !> This routine also updates VF and VL and computes the distances
     !> between the updated singular values and the old singular
     !> values.
     !> SLASD6 is called from SLASDA.

     pure subroutine la_slasd6(icompq,nl,nr,sqre,d,vf,vl,alpha,beta,idxq,perm, &
     givptr,givcol,ldgcol,givnum,ldgnum,poles,difl,difr,z,k,c,s,work,iwork,info)
        use la_constants_sp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(sp),intent(inout) :: alpha,beta
           real(sp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),iwork(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(sp),intent(inout) :: d(*),vf(*),vl(*)
           real(sp),intent(out) :: difl(*),difr(*),givnum(ldgnum,*),poles(ldgnum,*),work(*), &
                     z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idx,idxc,idxp,isigma,ivfw,ivlw,iw,m,n,n1,n2
           real(sp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -14
           else if (ldgnum < n) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('SLASD6',-info)
              return
           end if
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_slasd7 and la_slasd8.
           isigma = 1
           iw = isigma + n
           ivfw = iw + m
           ivlw = ivfw + m
           idx = 1
           idxc = idx + n
           idxp = idxc + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_slascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! sort and deflate singular values.
           call la_slasd7(icompq,nl,nr,sqre,k,d,z,work(iw),vf,work(ivfw),vl, &
           work(ivlw),alpha,beta,work(isigma),iwork(idx),iwork(idxp),idxq,perm, &
                     givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
           ! solve secular equation, compute difl, difr, and update vf, vl.
           call la_slasd8(icompq,k,d,z,vf,vl,difl,difr,ldgnum,work(isigma),work( &
                     iw),info)
           ! report the possible convergence failure.
           if (info /= 0) then
              return
           end if
           ! save the poles if icompq = 1.
           if (icompq == 1) then
              call la_scopy(k,d,1,poles(1,1),1)
              call la_scopy(k,work(isigma),1,poles(1,2),1)
           end if
           ! unscale.
           call la_slascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_slamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_slasd6
     !> DLASD6: computes the SVD of an updated upper bidiagonal matrix B
     !> obtained by merging two smaller ones by appending a row. This
     !> routine is used only for the problem which requires all singular
     !> values and optionally singular vector matrices in factored form.
     !> B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
     !> A related subroutine, DLASD1, handles the case in which all singular
     !> values and singular vectors of the bidiagonal matrix are desired.
     !> DLASD6 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The singular values of B can be computed using D1, D2, the first
     !> components of all the right singular vectors of the lower block, and
     !> the last components of all the right singular vectors of the upper
     !> block. These components are stored and updated in VF and VL,
     !> respectively, in DLASD6. Hence U and VT are not explicitly
     !> referenced.
     !> The singular values are stored in D. The algorithm consists of two
     !> stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or if there is a zero
     !> in the Z vector. For each such occurrence the dimension of the
     !> secular equation problem is reduced by one. This stage is
     !> performed by the routine DLASD7.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the roots of the
     !> secular equation via the routine DLASD4 (as called by DLASD8).
     !> This routine also updates VF and VL and computes the distances
     !> between the updated singular values and the old singular
     !> values.
     !> DLASD6 is called from DLASDA.

     pure subroutine la_dlasd6(icompq,nl,nr,sqre,d,vf,vl,alpha,beta,idxq,perm, &
     givptr,givcol,ldgcol,givnum,ldgnum,poles,difl,difr,z,k,c,s,work,iwork,info)
        use la_constants_dp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(dp),intent(inout) :: alpha,beta
           real(dp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),iwork(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(dp),intent(inout) :: d(*),vf(*),vl(*)
           real(dp),intent(out) :: difl(*),difr(*),givnum(ldgnum,*),poles(ldgnum,*),work(*), &
                     z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idx,idxc,idxp,isigma,ivfw,ivlw,iw,m,n,n1,n2
           real(dp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -14
           else if (ldgnum < n) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('DLASD6',-info)
              return
           end if
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_dlasd7 and la_dlasd8.
           isigma = 1
           iw = isigma + n
           ivfw = iw + m
           ivlw = ivfw + m
           idx = 1
           idxc = idx + n
           idxp = idxc + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_dlascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! sort and deflate singular values.
           call la_dlasd7(icompq,nl,nr,sqre,k,d,z,work(iw),vf,work(ivfw),vl, &
           work(ivlw),alpha,beta,work(isigma),iwork(idx),iwork(idxp),idxq,perm, &
                     givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
           ! solve secular equation, compute difl, difr, and update vf, vl.
           call la_dlasd8(icompq,k,d,z,vf,vl,difl,difr,ldgnum,work(isigma),work( &
                     iw),info)
           ! report the possible convergence failure.
           if (info /= 0) then
              return
           end if
           ! save the poles if icompq = 1.
           if (icompq == 1) then
              call la_dcopy(k,d,1,poles(1,1),1)
              call la_dcopy(k,work(isigma),1,poles(1,2),1)
           end if
           ! unscale.
           call la_dlascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_dlamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_dlasd6
#ifdef LA_WITH_XDP
     !> XLASD6: computes the SVD of an updated upper bidiagonal matrix B
     !> obtained by merging two smaller ones by appending a row. This
     !> routine is used only for the problem which requires all singular
     !> values and optionally singular vector matrices in factored form.
     !> B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
     !> A related subroutine, XLASD1, handles the case in which all singular
     !> values and singular vectors of the bidiagonal matrix are desired.
     !> XLASD6 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The singular values of B can be computed using D1, D2, the first
     !> components of all the right singular vectors of the lower block, and
     !> the last components of all the right singular vectors of the upper
     !> block. These components are stored and updated in VF and VL,
     !> respectively, in XLASD6. Hence U and VT are not explicitly
     !> referenced.
     !> The singular values are stored in D. The algorithm consists of two
     !> stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or if there is a zero
     !> in the Z vector. For each such occurrence the dimension of the
     !> secular equation problem is reduced by one. This stage is
     !> performed by the routine XLASD7.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the roots of the
     !> secular equation via the routine XLASD4 (as called by XLASD8).
     !> This routine also updates VF and VL and computes the distances
     !> between the updated singular values and the old singular
     !> values.
     !> XLASD6 is called from XLASDA.

     pure subroutine la_xlasd6(icompq,nl,nr,sqre,d,vf,vl,alpha,beta,idxq,perm, &
     givptr,givcol,ldgcol,givnum,ldgnum,poles,difl,difr,z,k,c,s,work,iwork,info)
        use la_constants_xdp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(xdp),intent(inout) :: alpha,beta
           real(xdp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),iwork(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(xdp),intent(inout) :: d(*),vf(*),vl(*)
           real(xdp),intent(out) :: difl(*),difr(*),givnum(ldgnum,*),poles(ldgnum,*),work(*), &
                     z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idx,idxc,idxp,isigma,ivfw,ivlw,iw,m,n,n1,n2
           real(xdp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -14
           else if (ldgnum < n) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('XLASD6',-info)
              return
           end if
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_xlasd7 and la_xlasd8.
           isigma = 1
           iw = isigma + n
           ivfw = iw + m
           ivlw = ivfw + m
           idx = 1
           idxc = idx + n
           idxp = idxc + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_xlascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! sort and deflate singular values.
           call la_xlasd7(icompq,nl,nr,sqre,k,d,z,work(iw),vf,work(ivfw),vl, &
           work(ivlw),alpha,beta,work(isigma),iwork(idx),iwork(idxp),idxq,perm, &
                     givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
           ! solve secular equation, compute difl, difr, and update vf, vl.
           call la_xlasd8(icompq,k,d,z,vf,vl,difl,difr,ldgnum,work(isigma),work( &
                     iw),info)
           ! report the possible convergence failure.
           if (info /= 0) then
              return
           end if
           ! save the poles if icompq = 1.
           if (icompq == 1) then
              call la_xcopy(k,d,1,poles(1,1),1)
              call la_xcopy(k,work(isigma),1,poles(1,2),1)
           end if
           ! unscale.
           call la_xlascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_xlamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_xlasd6
#endif
#ifdef LA_WITH_QP
     !> QLASD6: computes the SVD of an updated upper bidiagonal matrix B
     !> obtained by merging two smaller ones by appending a row. This
     !> routine is used only for the problem which requires all singular
     !> values and optionally singular vector matrices in factored form.
     !> B is an N-by-M matrix with N = NL + NR + 1 and M = N + SQRE.
     !> A related subroutine, QLASD1, handles the case in which all singular
     !> values and singular vectors of the bidiagonal matrix are desired.
     !> QLASD6 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The singular values of B can be computed using D1, D2, the first
     !> components of all the right singular vectors of the lower block, and
     !> the last components of all the right singular vectors of the upper
     !> block. These components are stored and updated in VF and VL,
     !> respectively, in QLASD6. Hence U and VT are not explicitly
     !> referenced.
     !> The singular values are stored in D. The algorithm consists of two
     !> stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or if there is a zero
     !> in the Z vector. For each such occurrence the dimension of the
     !> secular equation problem is reduced by one. This stage is
     !> performed by the routine QLASD7.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the roots of the
     !> secular equation via the routine QLASD4 (as called by QLASD8).
     !> This routine also updates VF and VL and computes the distances
     !> between the updated singular values and the old singular
     !> values.
     !> QLASD6 is called from QLASDA.

     pure subroutine la_qlasd6(icompq,nl,nr,sqre,d,vf,vl,alpha,beta,idxq,perm, &
     givptr,givcol,ldgcol,givnum,ldgnum,poles,difl,difr,z,k,c,s,work,iwork,info)
        use la_constants_qp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: givptr,info,k
           integer(ilp),intent(in) :: icompq,ldgcol,ldgnum,nl,nr,sqre
           real(qp),intent(inout) :: alpha,beta
           real(qp),intent(out) :: c,s
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),iwork(*),perm(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(qp),intent(inout) :: d(*),vf(*),vl(*)
           real(qp),intent(out) :: difl(*),difr(*),givnum(ldgnum,*),poles(ldgnum,*),work(*), &
                     z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,idx,idxc,idxp,isigma,ivfw,ivlw,iw,m,n,n1,n2
           real(qp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           n = nl + nr + 1
           m = n + sqre
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (nl < 1) then
              info = -2
           else if (nr < 1) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldgcol < n) then
              info = -14
           else if (ldgnum < n) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('QLASD6',-info)
              return
           end if
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_qlasd7 and la_qlasd8.
           isigma = 1
           iw = isigma + n
           ivfw = iw + m
           ivlw = ivfw + m
           idx = 1
           idxc = idx + n
           idxp = idxc + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_qlascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! sort and deflate singular values.
           call la_qlasd7(icompq,nl,nr,sqre,k,d,z,work(iw),vf,work(ivfw),vl, &
           work(ivlw),alpha,beta,work(isigma),iwork(idx),iwork(idxp),idxq,perm, &
                     givptr,givcol,ldgcol,givnum,ldgnum,c,s,info)
           ! solve secular equation, compute difl, difr, and update vf, vl.
           call la_qlasd8(icompq,k,d,z,vf,vl,difl,difr,ldgnum,work(isigma),work( &
                     iw),info)
           ! report the possible convergence failure.
           if (info /= 0) then
              return
           end if
           ! save the poles if icompq = 1.
           if (icompq == 1) then
              call la_qcopy(k,d,1,poles(1,1),1)
              call la_qcopy(k,work(isigma),1,poles(1,2),1)
           end if
           ! unscale.
           call la_qlascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_qlamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_qlasd6
#endif

     !> SLASD2: merges the two sets of singular values together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> singular values are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.
     !> SLASD2 is called from SLASD1.

     pure subroutine la_slasd2(nl,nr,sqre,k,d,z,alpha,beta,u,ldu,vt,ldvt,dsigma, &
               u2,ldu2,vt2,ldvt2,idxp,idx,idxc,idxq,coltyp,info)
        use la_constants_sp,only:zero,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           real(sp),intent(in) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),idx(*),idxc(*),idxp(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(sp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(sp),intent(out) :: dsigma(*),u2(ldu2,*),vt2(ldvt2,*),z(*)
        ! =====================================================================

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(sp) :: c,eps,hlftol,s,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           if (ldu < n) then
              info = -10
           else if (ldvt < m) then
              info = -12
           else if (ldu2 < n) then
              info = -15
           else if (ldvt2 < m) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('SLASD2',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           ! generate the first part of the vector z; and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vt(nlp1,nlp1)
           z(1) = z1
           do i = nl,1,-1
              z(i + 1) = alpha*vt(i,nlp1)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vt(i,nlp2)
           end do
           ! initialize some reference arrays.
           do i = 2,nlp1
              coltyp(i) = 1
           end do
           do i = nlp2,n
              coltyp(i) = 2
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and the first column of u2
           ! are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              u2(i,1) = z(idxq(i))
              idxc(i) = coltyp(idxq(i))
           end do
           call la_slamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = u2(idxi,1)
              coltyp(i) = idxc(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_slamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 coltyp(j) = 4
                 if (j == n) go to 120
              else
                 jprev = j
                 go to 90
              end if
           end do
           90 continue
           j = jprev
           100 continue
           j = j + 1
           if (j > n) go to 110
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
              coltyp(j) = 4
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_slapy2(c,s)
                 c = c/tau
                 s = -s/tau
                 z(j) = tau
                 z(jprev) = zero
                 ! apply back the givens rotation to the left and right
                 ! singular vector matrices.
                 idxjp = idxq(idx(jprev) + 1)
                 idxj = idxq(idx(j) + 1)
                 if (idxjp <= nlp1) then
                    idxjp = idxjp - 1
                 end if
                 if (idxj <= nlp1) then
                    idxj = idxj - 1
                 end if
                 call la_srot(n,u(1,idxjp),1,u(1,idxj),1,c,s)
                 call la_srot(m,vt(idxjp,1),ldvt,vt(idxj,1),ldvt,c,s)
                 if (coltyp(j) /= coltyp(jprev)) then
                    coltyp(j) = 3
                 end if
                 coltyp(jprev) = 4
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 u2(k,1) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 100
           110 continue
           ! record the last singular value.
           k = k + 1
           u2(k,1) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           120 continue
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four groups of uniform structure (although one or more of these
           ! groups may be empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 2,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 2
           psm(2) = 2 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           ! fill out the idxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's, starting from the
           ! second column. this applies similarly to the rows of vt.
           do j = 2,n
              jp = idxp(j)
              ct = coltyp(jp)
              idxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the singular values and corresponding singular vectors into
           ! dsigma, u2, and vt2 respectively.  the singular values/vectors
           ! which were not deflated go into the first k slots of dsigma, u2,
           ! and vt2 respectively, while those which were deflated go into the
           ! last n - k slots, except that the first column/row will be treated
           ! separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              idxj = idxq(idx(idxp(idxc(j))) + 1)
              if (idxj <= nlp1) then
                 idxj = idxj - 1
              end if
              call la_scopy(n,u(1,idxj),1,u2(1,j),1)
              call la_scopy(m,vt(idxj,1),ldvt,vt2(j,1),ldvt2)
           end do
           ! determine dsigma(1), dsigma(2) and z(1)
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_slapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = z(m)/z(1)
              end if
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! move the rest of the updating row to z.
           call la_scopy(k - 1,u2(2,1),1,z(2),1)
           ! determine the first column of u2, the first row of vt2 and the
           ! last row of vt.
           call la_slaset('A',n,1,zero,zero,u2,ldu2)
           u2(nlp1,1) = one
           if (m > n) then
              do i = 1,nlp1
                 vt(m,i) = -s*vt(nlp1,i)
                 vt2(1,i) = c*vt(nlp1,i)
              end do
              do i = nlp2,m
                 vt2(1,i) = s*vt(m,i)
                 vt(m,i) = c*vt(m,i)
              end do
           else
              call la_scopy(m,vt(nlp1,1),ldvt,vt2(1,1),ldvt2)
           end if
           if (m > n) then
              call la_scopy(m,vt(m,1),ldvt,vt2(m,1),ldvt2)
           end if
           ! the deflated singular values and their corresponding vectors go
           ! into the back of d, u, and v respectively.
           if (n > k) then
              call la_scopy(n - k,dsigma(k + 1),1,d(k + 1),1)
              call la_slacpy('A',n,n - k,u2(1,k + 1),ldu2,u(1,k + 1),ldu)
              call la_slacpy('A',n - k,m,vt2(k + 1,1),ldvt2,vt(k + 1,1),ldvt)
           end if
           ! copy ctot into coltyp for referencing in la_slasd3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           return
     end subroutine la_slasd2
     !> DLASD2: merges the two sets of singular values together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> singular values are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.
     !> DLASD2 is called from DLASD1.

     pure subroutine la_dlasd2(nl,nr,sqre,k,d,z,alpha,beta,u,ldu,vt,ldvt,dsigma, &
               u2,ldu2,vt2,ldvt2,idxp,idx,idxc,idxq,coltyp,info)
        use la_constants_dp,only:zero,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           real(dp),intent(in) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),idx(*),idxc(*),idxp(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(dp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(dp),intent(out) :: dsigma(*),u2(ldu2,*),vt2(ldvt2,*),z(*)
        ! =====================================================================

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(dp) :: c,eps,hlftol,s,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           if (ldu < n) then
              info = -10
           else if (ldvt < m) then
              info = -12
           else if (ldu2 < n) then
              info = -15
           else if (ldvt2 < m) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('DLASD2',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           ! generate the first part of the vector z; and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vt(nlp1,nlp1)
           z(1) = z1
           do i = nl,1,-1
              z(i + 1) = alpha*vt(i,nlp1)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vt(i,nlp2)
           end do
           ! initialize some reference arrays.
           do i = 2,nlp1
              coltyp(i) = 1
           end do
           do i = nlp2,n
              coltyp(i) = 2
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and the first column of u2
           ! are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              u2(i,1) = z(idxq(i))
              idxc(i) = coltyp(idxq(i))
           end do
           call la_dlamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = u2(idxi,1)
              coltyp(i) = idxc(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_dlamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 coltyp(j) = 4
                 if (j == n) go to 120
              else
                 jprev = j
                 go to 90
              end if
           end do
           90 continue
           j = jprev
           100 continue
           j = j + 1
           if (j > n) go to 110
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
              coltyp(j) = 4
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_dlapy2(c,s)
                 c = c/tau
                 s = -s/tau
                 z(j) = tau
                 z(jprev) = zero
                 ! apply back the givens rotation to the left and right
                 ! singular vector matrices.
                 idxjp = idxq(idx(jprev) + 1)
                 idxj = idxq(idx(j) + 1)
                 if (idxjp <= nlp1) then
                    idxjp = idxjp - 1
                 end if
                 if (idxj <= nlp1) then
                    idxj = idxj - 1
                 end if
                 call la_drot(n,u(1,idxjp),1,u(1,idxj),1,c,s)
                 call la_drot(m,vt(idxjp,1),ldvt,vt(idxj,1),ldvt,c,s)
                 if (coltyp(j) /= coltyp(jprev)) then
                    coltyp(j) = 3
                 end if
                 coltyp(jprev) = 4
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 u2(k,1) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 100
           110 continue
           ! record the last singular value.
           k = k + 1
           u2(k,1) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           120 continue
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four groups of uniform structure (although one or more of these
           ! groups may be empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 2,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 2
           psm(2) = 2 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           ! fill out the idxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's, starting from the
           ! second column. this applies similarly to the rows of vt.
           do j = 2,n
              jp = idxp(j)
              ct = coltyp(jp)
              idxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the singular values and corresponding singular vectors into
           ! dsigma, u2, and vt2 respectively.  the singular values/vectors
           ! which were not deflated go into the first k slots of dsigma, u2,
           ! and vt2 respectively, while those which were deflated go into the
           ! last n - k slots, except that the first column/row will be treated
           ! separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              idxj = idxq(idx(idxp(idxc(j))) + 1)
              if (idxj <= nlp1) then
                 idxj = idxj - 1
              end if
              call la_dcopy(n,u(1,idxj),1,u2(1,j),1)
              call la_dcopy(m,vt(idxj,1),ldvt,vt2(j,1),ldvt2)
           end do
           ! determine dsigma(1), dsigma(2) and z(1)
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_dlapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = z(m)/z(1)
              end if
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! move the rest of the updating row to z.
           call la_dcopy(k - 1,u2(2,1),1,z(2),1)
           ! determine the first column of u2, the first row of vt2 and the
           ! last row of vt.
           call la_dlaset('A',n,1,zero,zero,u2,ldu2)
           u2(nlp1,1) = one
           if (m > n) then
              do i = 1,nlp1
                 vt(m,i) = -s*vt(nlp1,i)
                 vt2(1,i) = c*vt(nlp1,i)
              end do
              do i = nlp2,m
                 vt2(1,i) = s*vt(m,i)
                 vt(m,i) = c*vt(m,i)
              end do
           else
              call la_dcopy(m,vt(nlp1,1),ldvt,vt2(1,1),ldvt2)
           end if
           if (m > n) then
              call la_dcopy(m,vt(m,1),ldvt,vt2(m,1),ldvt2)
           end if
           ! the deflated singular values and their corresponding vectors go
           ! into the back of d, u, and v respectively.
           if (n > k) then
              call la_dcopy(n - k,dsigma(k + 1),1,d(k + 1),1)
              call la_dlacpy('A',n,n - k,u2(1,k + 1),ldu2,u(1,k + 1),ldu)
              call la_dlacpy('A',n - k,m,vt2(k + 1,1),ldvt2,vt(k + 1,1),ldvt)
           end if
           ! copy ctot into coltyp for referencing in la_dlasd3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           return
     end subroutine la_dlasd2
#ifdef LA_WITH_XDP
     !> XLASD2: merges the two sets of singular values together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> singular values are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.
     !> XLASD2 is called from XLASD1.

     pure subroutine la_xlasd2(nl,nr,sqre,k,d,z,alpha,beta,u,ldu,vt,ldvt,dsigma, &
               u2,ldu2,vt2,ldvt2,idxp,idx,idxc,idxq,coltyp,info)
        use la_constants_xdp,only:zero,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           real(xdp),intent(in) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),idx(*),idxc(*),idxp(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(xdp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(xdp),intent(out) :: dsigma(*),u2(ldu2,*),vt2(ldvt2,*),z(*)
        ! =====================================================================

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(xdp) :: c,eps,hlftol,s,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           if (ldu < n) then
              info = -10
           else if (ldvt < m) then
              info = -12
           else if (ldu2 < n) then
              info = -15
           else if (ldvt2 < m) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('XLASD2',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           ! generate the first part of the vector z; and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vt(nlp1,nlp1)
           z(1) = z1
           do i = nl,1,-1
              z(i + 1) = alpha*vt(i,nlp1)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vt(i,nlp2)
           end do
           ! initialize some reference arrays.
           do i = 2,nlp1
              coltyp(i) = 1
           end do
           do i = nlp2,n
              coltyp(i) = 2
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and the first column of u2
           ! are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              u2(i,1) = z(idxq(i))
              idxc(i) = coltyp(idxq(i))
           end do
           call la_xlamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = u2(idxi,1)
              coltyp(i) = idxc(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_xlamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 coltyp(j) = 4
                 if (j == n) go to 120
              else
                 jprev = j
                 go to 90
              end if
           end do
           90 continue
           j = jprev
           100 continue
           j = j + 1
           if (j > n) go to 110
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
              coltyp(j) = 4
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_xlapy2(c,s)
                 c = c/tau
                 s = -s/tau
                 z(j) = tau
                 z(jprev) = zero
                 ! apply back the givens rotation to the left and right
                 ! singular vector matrices.
                 idxjp = idxq(idx(jprev) + 1)
                 idxj = idxq(idx(j) + 1)
                 if (idxjp <= nlp1) then
                    idxjp = idxjp - 1
                 end if
                 if (idxj <= nlp1) then
                    idxj = idxj - 1
                 end if
                 call la_xrot(n,u(1,idxjp),1,u(1,idxj),1,c,s)
                 call la_xrot(m,vt(idxjp,1),ldvt,vt(idxj,1),ldvt,c,s)
                 if (coltyp(j) /= coltyp(jprev)) then
                    coltyp(j) = 3
                 end if
                 coltyp(jprev) = 4
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 u2(k,1) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 100
           110 continue
           ! record the last singular value.
           k = k + 1
           u2(k,1) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           120 continue
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four groups of uniform structure (although one or more of these
           ! groups may be empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 2,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 2
           psm(2) = 2 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           ! fill out the idxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's, starting from the
           ! second column. this applies similarly to the rows of vt.
           do j = 2,n
              jp = idxp(j)
              ct = coltyp(jp)
              idxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the singular values and corresponding singular vectors into
           ! dsigma, u2, and vt2 respectively.  the singular values/vectors
           ! which were not deflated go into the first k slots of dsigma, u2,
           ! and vt2 respectively, while those which were deflated go into the
           ! last n - k slots, except that the first column/row will be treated
           ! separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              idxj = idxq(idx(idxp(idxc(j))) + 1)
              if (idxj <= nlp1) then
                 idxj = idxj - 1
              end if
              call la_xcopy(n,u(1,idxj),1,u2(1,j),1)
              call la_xcopy(m,vt(idxj,1),ldvt,vt2(j,1),ldvt2)
           end do
           ! determine dsigma(1), dsigma(2) and z(1)
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_xlapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = z(m)/z(1)
              end if
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! move the rest of the updating row to z.
           call la_xcopy(k - 1,u2(2,1),1,z(2),1)
           ! determine the first column of u2, the first row of vt2 and the
           ! last row of vt.
           call la_xlaset('A',n,1,zero,zero,u2,ldu2)
           u2(nlp1,1) = one
           if (m > n) then
              do i = 1,nlp1
                 vt(m,i) = -s*vt(nlp1,i)
                 vt2(1,i) = c*vt(nlp1,i)
              end do
              do i = nlp2,m
                 vt2(1,i) = s*vt(m,i)
                 vt(m,i) = c*vt(m,i)
              end do
           else
              call la_xcopy(m,vt(nlp1,1),ldvt,vt2(1,1),ldvt2)
           end if
           if (m > n) then
              call la_xcopy(m,vt(m,1),ldvt,vt2(m,1),ldvt2)
           end if
           ! the deflated singular values and their corresponding vectors go
           ! into the back of d, u, and v respectively.
           if (n > k) then
              call la_xcopy(n - k,dsigma(k + 1),1,d(k + 1),1)
              call la_xlacpy('A',n,n - k,u2(1,k + 1),ldu2,u(1,k + 1),ldu)
              call la_xlacpy('A',n - k,m,vt2(k + 1,1),ldvt2,vt(k + 1,1),ldvt)
           end if
           ! copy ctot into coltyp for referencing in la_xlasd3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           return
     end subroutine la_xlasd2
#endif
#ifdef LA_WITH_QP
     !> QLASD2: merges the two sets of singular values together into a single
     !> sorted set.  Then it tries to deflate the size of the problem.
     !> There are two ways in which deflation can occur:  when two or more
     !> singular values are close together or if there is a tiny entry in the
     !> Z vector.  For each such occurrence the order of the related secular
     !> equation problem is reduced by one.
     !> QLASD2 is called from QLASD1.

     pure subroutine la_qlasd2(nl,nr,sqre,k,d,z,alpha,beta,u,ldu,vt,ldvt,dsigma, &
               u2,ldu2,vt2,ldvt2,idxp,idx,idxc,idxq,coltyp,info)
        use la_constants_qp,only:zero,one,two,eight
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,k
           integer(ilp),intent(in) :: ldu,ldu2,ldvt,ldvt2,nl,nr,sqre
           real(qp),intent(in) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(out) :: coltyp(*),idx(*),idxc(*),idxp(*)
           integer(ilp),intent(inout) :: idxq(*)
           real(qp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(qp),intent(out) :: dsigma(*),u2(ldu2,*),vt2(ldvt2,*),z(*)
        ! =====================================================================

           ! Local Arrays
           integer(ilp) :: ctot(4),psm(4)
           ! Local Scalars
           integer(ilp) :: ct,i,idxi,idxj,idxjp,j,jp,jprev,k2,m,n,nlp1,nlp2
           real(qp) :: c,eps,hlftol,s,tau,tol,z1
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre /= 1) .and. (sqre /= 0)) then
              info = -3
           end if
           n = nl + nr + 1
           m = n + sqre
           if (ldu < n) then
              info = -10
           else if (ldvt < m) then
              info = -12
           else if (ldu2 < n) then
              info = -15
           else if (ldvt2 < m) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('QLASD2',-info)
              return
           end if
           nlp1 = nl + 1
           nlp2 = nl + 2
           ! generate the first part of the vector z; and move the singular
           ! values in the first part of d one position backward.
           z1 = alpha*vt(nlp1,nlp1)
           z(1) = z1
           do i = nl,1,-1
              z(i + 1) = alpha*vt(i,nlp1)
              d(i + 1) = d(i)
              idxq(i + 1) = idxq(i) + 1
           end do
           ! generate the second part of the vector z.
           do i = nlp2,m
              z(i) = beta*vt(i,nlp2)
           end do
           ! initialize some reference arrays.
           do i = 2,nlp1
              coltyp(i) = 1
           end do
           do i = nlp2,n
              coltyp(i) = 2
           end do
           ! sort the singular values into increasing order
           do i = nlp2,n
              idxq(i) = idxq(i) + nlp1
           end do
           ! dsigma, idxc, idxc, and the first column of u2
           ! are used as storage space.
           do i = 2,n
              dsigma(i) = d(idxq(i))
              u2(i,1) = z(idxq(i))
              idxc(i) = coltyp(idxq(i))
           end do
           call la_qlamrg(nl,nr,dsigma(2),1,1,idx(2))
           do i = 2,n
              idxi = 1 + idx(i)
              d(i) = dsigma(idxi)
              z(i) = u2(idxi,1)
              coltyp(i) = idxc(idxi)
           end do
           ! calculate the allowable deflation tolerance
           eps = la_qlamch('EPSILON')
           tol = max(abs(alpha),abs(beta))
           tol = eight*eps*max(abs(d(n)),tol)
           ! there are 2 kinds of deflation -- first a value in the z-vector
           ! is small, second two (or more) singular values are very close
           ! together (their difference is small).
           ! if the value in the z-vector is small, we simply permute the
           ! array so that the corresponding singular value is moved to the
           ! end.
           ! if two values in the d-vector are close, we perform a two-sided
           ! rotation designed to make one of the corresponding z-vector
           ! entries zero, and then permute the array so that the deflated
           ! singular value is moved to the end.
           ! if there are multiple singular values then the problem deflates.
           ! here the number of equal singular values are found.  as each equal
           ! singular value is found, an elementary reflector is computed to
           ! rotate the corresponding singular subspace so that the
           ! corresponding components of z are zero in this new basis.
           k = 1
           k2 = n + 1
           do j = 2,n
              if (abs(z(j)) <= tol) then
                 ! deflate due to small z component.
                 k2 = k2 - 1
                 idxp(k2) = j
                 coltyp(j) = 4
                 if (j == n) go to 120
              else
                 jprev = j
                 go to 90
              end if
           end do
           90 continue
           j = jprev
           100 continue
           j = j + 1
           if (j > n) go to 110
           if (abs(z(j)) <= tol) then
              ! deflate due to small z component.
              k2 = k2 - 1
              idxp(k2) = j
              coltyp(j) = 4
           else
              ! check if singular values are close enough to allow deflation.
              if (abs(d(j) - d(jprev)) <= tol) then
                 ! deflation is possible.
                 s = z(jprev)
                 c = z(j)
                 ! find sqrt(a**2+b**2) without overflow or
                 ! destructive underflow.
                 tau = la_qlapy2(c,s)
                 c = c/tau
                 s = -s/tau
                 z(j) = tau
                 z(jprev) = zero
                 ! apply back the givens rotation to the left and right
                 ! singular vector matrices.
                 idxjp = idxq(idx(jprev) + 1)
                 idxj = idxq(idx(j) + 1)
                 if (idxjp <= nlp1) then
                    idxjp = idxjp - 1
                 end if
                 if (idxj <= nlp1) then
                    idxj = idxj - 1
                 end if
                 call la_qrot(n,u(1,idxjp),1,u(1,idxj),1,c,s)
                 call la_qrot(m,vt(idxjp,1),ldvt,vt(idxj,1),ldvt,c,s)
                 if (coltyp(j) /= coltyp(jprev)) then
                    coltyp(j) = 3
                 end if
                 coltyp(jprev) = 4
                 k2 = k2 - 1
                 idxp(k2) = jprev
                 jprev = j
              else
                 k = k + 1
                 u2(k,1) = z(jprev)
                 dsigma(k) = d(jprev)
                 idxp(k) = jprev
                 jprev = j
              end if
           end if
           go to 100
           110 continue
           ! record the last singular value.
           k = k + 1
           u2(k,1) = z(jprev)
           dsigma(k) = d(jprev)
           idxp(k) = jprev
           120 continue
           ! count up the total number of the various types of columns, then
           ! form a permutation which positions the four column types into
           ! four groups of uniform structure (although one or more of these
           ! groups may be empty).
           do j = 1,4
              ctot(j) = 0
           end do
           do j = 2,n
              ct = coltyp(j)
              ctot(ct) = ctot(ct) + 1
           end do
           ! psm(*) = position in submatrix (of types 1 through 4)
           psm(1) = 2
           psm(2) = 2 + ctot(1)
           psm(3) = psm(2) + ctot(2)
           psm(4) = psm(3) + ctot(3)
           ! fill out the idxc array so that the permutation which it induces
           ! will place all type-1 columns first, all type-2 columns next,
           ! then all type-3's, and finally all type-4's, starting from the
           ! second column. this applies similarly to the rows of vt.
           do j = 2,n
              jp = idxp(j)
              ct = coltyp(jp)
              idxc(psm(ct)) = j
              psm(ct) = psm(ct) + 1
           end do
           ! sort the singular values and corresponding singular vectors into
           ! dsigma, u2, and vt2 respectively.  the singular values/vectors
           ! which were not deflated go into the first k slots of dsigma, u2,
           ! and vt2 respectively, while those which were deflated go into the
           ! last n - k slots, except that the first column/row will be treated
           ! separately.
           do j = 2,n
              jp = idxp(j)
              dsigma(j) = d(jp)
              idxj = idxq(idx(idxp(idxc(j))) + 1)
              if (idxj <= nlp1) then
                 idxj = idxj - 1
              end if
              call la_qcopy(n,u(1,idxj),1,u2(1,j),1)
              call la_qcopy(m,vt(idxj,1),ldvt,vt2(j,1),ldvt2)
           end do
           ! determine dsigma(1), dsigma(2) and z(1)
           dsigma(1) = zero
           hlftol = tol/two
           if (abs(dsigma(2)) <= hlftol) dsigma(2) = hlftol
           if (m > n) then
              z(1) = la_qlapy2(z1,z(m))
              if (z(1) <= tol) then
                 c = one
                 s = zero
                 z(1) = tol
              else
                 c = z1/z(1)
                 s = z(m)/z(1)
              end if
           else
              if (abs(z1) <= tol) then
                 z(1) = tol
              else
                 z(1) = z1
              end if
           end if
           ! move the rest of the updating row to z.
           call la_qcopy(k - 1,u2(2,1),1,z(2),1)
           ! determine the first column of u2, the first row of vt2 and the
           ! last row of vt.
           call la_qlaset('A',n,1,zero,zero,u2,ldu2)
           u2(nlp1,1) = one
           if (m > n) then
              do i = 1,nlp1
                 vt(m,i) = -s*vt(nlp1,i)
                 vt2(1,i) = c*vt(nlp1,i)
              end do
              do i = nlp2,m
                 vt2(1,i) = s*vt(m,i)
                 vt(m,i) = c*vt(m,i)
              end do
           else
              call la_qcopy(m,vt(nlp1,1),ldvt,vt2(1,1),ldvt2)
           end if
           if (m > n) then
              call la_qcopy(m,vt(m,1),ldvt,vt2(m,1),ldvt2)
           end if
           ! the deflated singular values and their corresponding vectors go
           ! into the back of d, u, and v respectively.
           if (n > k) then
              call la_qcopy(n - k,dsigma(k + 1),1,d(k + 1),1)
              call la_qlacpy('A',n,n - k,u2(1,k + 1),ldu2,u(1,k + 1),ldu)
              call la_qlacpy('A',n - k,m,vt2(k + 1,1),ldvt2,vt(k + 1,1),ldvt)
           end if
           ! copy ctot into coltyp for referencing in la_qlasd3.
           do j = 1,4
              coltyp(j) = ctot(j)
           end do
           return
     end subroutine la_qlasd2
#endif

     !> SLASD1: computes the SVD of an upper bidiagonal N-by-M matrix B,
     !> where N = NL + NR + 1 and M = N + SQRE. SLASD1 is called from SLASD0.
     !> A related subroutine SLASD7 handles the case in which the singular
     !> values (and the singular vectors in factored form) are desired.
     !> SLASD1 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The left singular vectors of the original matrix are stored in U, and
     !> the transpose of the right singular vectors are stored in VT, and the
     !> singular values are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or when there are zeros in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine SLASD2.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the square roots of the
     !> roots of the secular equation via the routine SLASD4 (as called
     !> by SLASD3). This routine also calculates the singular vectors of
     !> the current problem.
     !> The final stage consists of computing the updated singular vectors
     !> directly using the updated singular values.  The singular vectors
     !> for the current problem are multiplied with the singular vectors
     !> from the overall problem.

     pure subroutine la_slasd1(nl,nr,sqre,d,alpha,beta,u,ldu,vt,ldvt,idxq,iwork, &
               work,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,nl,nr,sqre
           real(sp),intent(inout) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(inout) :: idxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,i,idx,idxc,idxp,iq,isigma,iu2,ivt2,iz,k,ldq,ldu2, &
                     ldvt2,m,n,n1,n2
           real(sp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('SLASD1',-info)
              return
           end if
           n = nl + nr + 1
           m = n + sqre
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_slasd2 and la_slasd3.
           ldu2 = n
           ldvt2 = m
           iz = 1
           isigma = iz + m
           iu2 = isigma + n
           ivt2 = iu2 + ldu2*n
           iq = ivt2 + ldvt2*m
           idx = 1
           idxc = idx + n
           coltyp = idxc + n
           idxp = coltyp + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_slascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! deflate singular values.
           call la_slasd2(nl,nr,sqre,k,d,work(iz),alpha,beta,u,ldu,vt,ldvt,work( &
            isigma),work(iu2),ldu2,work(ivt2),ldvt2,iwork(idxp),iwork(idx),iwork( &
                      idxc),idxq,iwork(coltyp),info)
           ! solve secular equation and update singular vectors.
           ldq = k
           call la_slasd3(nl,nr,sqre,k,d,work(iq),ldq,work(isigma),u,ldu,work( &
           iu2),ldu2,vt,ldvt,work(ivt2),ldvt2,iwork(idxc),iwork(coltyp),work(iz), &
                     info)
           ! report the possible convergence failure.
           if (info /= 0) then
              return
           end if
           ! unscale.
           call la_slascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_slamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_slasd1
     !> DLASD1: computes the SVD of an upper bidiagonal N-by-M matrix B,
     !> where N = NL + NR + 1 and M = N + SQRE. DLASD1 is called from DLASD0.
     !> A related subroutine DLASD7 handles the case in which the singular
     !> values (and the singular vectors in factored form) are desired.
     !> DLASD1 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The left singular vectors of the original matrix are stored in U, and
     !> the transpose of the right singular vectors are stored in VT, and the
     !> singular values are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or when there are zeros in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine DLASD2.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the square roots of the
     !> roots of the secular equation via the routine DLASD4 (as called
     !> by DLASD3). This routine also calculates the singular vectors of
     !> the current problem.
     !> The final stage consists of computing the updated singular vectors
     !> directly using the updated singular values.  The singular vectors
     !> for the current problem are multiplied with the singular vectors
     !> from the overall problem.

     pure subroutine la_dlasd1(nl,nr,sqre,d,alpha,beta,u,ldu,vt,ldvt,idxq,iwork, &
               work,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,nl,nr,sqre
           real(dp),intent(inout) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(inout) :: idxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,i,idx,idxc,idxp,iq,isigma,iu2,ivt2,iz,k,ldq,ldu2, &
                     ldvt2,m,n,n1,n2
           real(dp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('DLASD1',-info)
              return
           end if
           n = nl + nr + 1
           m = n + sqre
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_dlasd2 and la_dlasd3.
           ldu2 = n
           ldvt2 = m
           iz = 1
           isigma = iz + m
           iu2 = isigma + n
           ivt2 = iu2 + ldu2*n
           iq = ivt2 + ldvt2*m
           idx = 1
           idxc = idx + n
           coltyp = idxc + n
           idxp = coltyp + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_dlascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! deflate singular values.
           call la_dlasd2(nl,nr,sqre,k,d,work(iz),alpha,beta,u,ldu,vt,ldvt,work( &
            isigma),work(iu2),ldu2,work(ivt2),ldvt2,iwork(idxp),iwork(idx),iwork( &
                      idxc),idxq,iwork(coltyp),info)
           ! solve secular equation and update singular vectors.
           ldq = k
           call la_dlasd3(nl,nr,sqre,k,d,work(iq),ldq,work(isigma),u,ldu,work( &
           iu2),ldu2,vt,ldvt,work(ivt2),ldvt2,iwork(idxc),iwork(coltyp),work(iz), &
                     info)
           ! report the convergence failure.
           if (info /= 0) then
              return
           end if
           ! unscale.
           call la_dlascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_dlamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_dlasd1
#ifdef LA_WITH_XDP
     !> XLASD1: computes the SVD of an upper bidiagonal N-by-M matrix B,
     !> where N = NL + NR + 1 and M = N + SQRE. XLASD1 is called from XLASD0.
     !> A related subroutine XLASD7 handles the case in which the singular
     !> values (and the singular vectors in factored form) are desired.
     !> XLASD1 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The left singular vectors of the original matrix are stored in U, and
     !> the transpose of the right singular vectors are stored in VT, and the
     !> singular values are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or when there are zeros in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine XLASD2.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the square roots of the
     !> roots of the secular equation via the routine XLASD4 (as called
     !> by XLASD3). This routine also calculates the singular vectors of
     !> the current problem.
     !> The final stage consists of computing the updated singular vectors
     !> directly using the updated singular values.  The singular vectors
     !> for the current problem are multiplied with the singular vectors
     !> from the overall problem.

     pure subroutine la_xlasd1(nl,nr,sqre,d,alpha,beta,u,ldu,vt,ldvt,idxq,iwork, &
               work,info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,nl,nr,sqre
           real(xdp),intent(inout) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(inout) :: idxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,i,idx,idxc,idxp,iq,isigma,iu2,ivt2,iz,k,ldq,ldu2, &
                     ldvt2,m,n,n1,n2
           real(xdp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('XLASD1',-info)
              return
           end if
           n = nl + nr + 1
           m = n + sqre
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_xlasd2 and la_xlasd3.
           ldu2 = n
           ldvt2 = m
           iz = 1
           isigma = iz + m
           iu2 = isigma + n
           ivt2 = iu2 + ldu2*n
           iq = ivt2 + ldvt2*m
           idx = 1
           idxc = idx + n
           coltyp = idxc + n
           idxp = coltyp + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_xlascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! deflate singular values.
           call la_xlasd2(nl,nr,sqre,k,d,work(iz),alpha,beta,u,ldu,vt,ldvt,work( &
            isigma),work(iu2),ldu2,work(ivt2),ldvt2,iwork(idxp),iwork(idx),iwork( &
                      idxc),idxq,iwork(coltyp),info)
           ! solve secular equation and update singular vectors.
           ldq = k
           call la_xlasd3(nl,nr,sqre,k,d,work(iq),ldq,work(isigma),u,ldu,work( &
           iu2),ldu2,vt,ldvt,work(ivt2),ldvt2,iwork(idxc),iwork(coltyp),work(iz), &
                     info)
           ! report the convergence failure.
           if (info /= 0) then
              return
           end if
           ! unscale.
           call la_xlascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_xlamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_xlasd1
#endif
#ifdef LA_WITH_QP
     !> QLASD1: computes the SVD of an upper bidiagonal N-by-M matrix B,
     !> where N = NL + NR + 1 and M = N + SQRE. QLASD1 is called from QLASD0.
     !> A related subroutine QLASD7 handles the case in which the singular
     !> values (and the singular vectors in factored form) are desired.
     !> QLASD1 computes the SVD as follows:
     !> ( D1(in)    0    0       0 )
     !> B = U(in) * (   Z1**T   a   Z2**T    b ) * VT(in)
     !> (   0       0   D2(in)   0 )
     !> = U(out) * ( D(out) 0) * VT(out)
     !> where Z**T = (Z1**T a Z2**T b) = u**T VT**T, and u is a vector of dimension M
     !> with ALPHA and BETA in the NL+1 and NL+2 th entries and zeros
     !> elsewhere; and the entry b is empty if SQRE = 0.
     !> The left singular vectors of the original matrix are stored in U, and
     !> the transpose of the right singular vectors are stored in VT, and the
     !> singular values are in D.  The algorithm consists of three stages:
     !> The first stage consists of deflating the size of the problem
     !> when there are multiple singular values or when there are zeros in
     !> the Z vector.  For each such occurrence the dimension of the
     !> secular equation problem is reduced by one.  This stage is
     !> performed by the routine QLASD2.
     !> The second stage consists of calculating the updated
     !> singular values. This is done by finding the square roots of the
     !> roots of the secular equation via the routine QLASD4 (as called
     !> by QLASD3). This routine also calculates the singular vectors of
     !> the current problem.
     !> The final stage consists of computing the updated singular vectors
     !> directly using the updated singular values.  The singular vectors
     !> for the current problem are multiplied with the singular vectors
     !> from the overall problem.

     pure subroutine la_qlasd1(nl,nr,sqre,d,alpha,beta,u,ldu,vt,ldvt,idxq,iwork, &
               work,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,nl,nr,sqre
           real(qp),intent(inout) :: alpha,beta
           ! Array Arguments
           integer(ilp),intent(inout) :: idxq(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),u(ldu,*),vt(ldvt,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: coltyp,i,idx,idxc,idxp,iq,isigma,iu2,ivt2,iz,k,ldq,ldu2, &
                     ldvt2,m,n,n1,n2
           real(qp) :: orgnrm
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (nl < 1) then
              info = -1
           else if (nr < 1) then
              info = -2
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('QLASD1',-info)
              return
           end if
           n = nl + nr + 1
           m = n + sqre
           ! the following values are for bookkeeping purposes only.  they are
           ! integer pointers which indicate the portion of the workspace
           ! used by a particular array in la_qlasd2 and la_qlasd3.
           ldu2 = n
           ldvt2 = m
           iz = 1
           isigma = iz + m
           iu2 = isigma + n
           ivt2 = iu2 + ldu2*n
           iq = ivt2 + ldvt2*m
           idx = 1
           idxc = idx + n
           coltyp = idxc + n
           idxp = coltyp + n
           ! scale.
           orgnrm = max(abs(alpha),abs(beta))
           d(nl + 1) = zero
           do i = 1,n
              if (abs(d(i)) > orgnrm) then
                 orgnrm = abs(d(i))
              end if
           end do
           call la_qlascl('G',0,0,orgnrm,one,n,1,d,n,info)
           alpha = alpha/orgnrm
           beta = beta/orgnrm
           ! deflate singular values.
           call la_qlasd2(nl,nr,sqre,k,d,work(iz),alpha,beta,u,ldu,vt,ldvt,work( &
            isigma),work(iu2),ldu2,work(ivt2),ldvt2,iwork(idxp),iwork(idx),iwork( &
                      idxc),idxq,iwork(coltyp),info)
           ! solve secular equation and update singular vectors.
           ldq = k
           call la_qlasd3(nl,nr,sqre,k,d,work(iq),ldq,work(isigma),u,ldu,work( &
           iu2),ldu2,vt,ldvt,work(ivt2),ldvt2,iwork(idxc),iwork(coltyp),work(iz), &
                     info)
           ! report the convergence failure.
           if (info /= 0) then
              return
           end if
           ! unscale.
           call la_qlascl('G',0,0,one,orgnrm,n,1,d,n,info)
           ! prepare the idxq sorting permutation.
           n1 = k
           n2 = n - k
           call la_qlamrg(n1,n2,d,1,-1,idxq)
           return
     end subroutine la_qlasd1
#endif

     !> SBDSDC: computes the singular value decomposition (SVD) of a real
     !> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
     !> using a divide and conquer method, where S is a diagonal matrix
     !> with non-negative diagonal elements (the singular values of B), and
     !> U and VT are orthogonal matrices of left and right singular vectors,
     !> respectively. SBDSDC can be used to compute all singular values,
     !> and optionally, singular vectors or singular vectors in compact form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See SLASD3 for details.
     !> The code currently calls SLASDQ if singular values only are desired.
     !> However, it can be slightly modified to compute singular values
     !> using the divide and conquer method.

     pure subroutine la_sbdsdc(uplo,compq,n,d,e,u,ldu,vt,ldvt,q,iq,work,iwork, &
               info)
        use la_constants_sp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n
           ! Array Arguments
           integer(ilp),intent(out) :: iq(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: q(*),u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
        ! changed dimension statement in comment describing e from (n) to
        ! (n-1).  sven, 17 feb 05.
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: difl,difr,givcol,givnum,givptr,i,ic,icompq,ierr,ii,is,iu, &
           iuplo,ivt,j,k,kk,mlvl,nm1,nsize,perm,poles,qstart,smlsiz,smlszp,sqre, &
                     start,wstart,z
           real(sp) :: cs,eps,orgnrm,p,r,sn
           ! Intrinsic Functions
           intrinsic :: real,abs,int,log,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (la_lsame(compq,'N')) then
              icompq = 0
           else if (la_lsame(compq,'P')) then
              icompq = 1
           else if (la_lsame(compq,'I')) then
              icompq = 2
           else
              icompq = -1
           end if
           if (iuplo == 0) then
              info = -1
           else if (icompq < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((ldu < 1) .or. ((icompq == 2) .and. (ldu < n))) then
              info = -7
           else if ((ldvt < 1) .or. ((icompq == 2) .and. (ldvt < n))) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SBDSDC',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'SBDSDC',' ',0,0,0,0)
           if (n == 1) then
              if (icompq == 1) then
                 q(1) = sign(one,d(1))
                 q(1 + smlsiz*n) = one
              else if (icompq == 2) then
                 u(1,1) = sign(one,d(1))
                 vt(1,1) = one
              end if
              d(1) = abs(d(1))
              return
           end if
           nm1 = n - 1
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           wstart = 1
           qstart = 3
           if (icompq == 1) then
              call la_scopy(n,d,1,q(1),1)
              call la_scopy(n - 1,e,1,q(n + 1),1)
           end if
           if (iuplo == 2) then
              qstart = 5
              if (icompq == 2) wstart = 2*n - 1
              do i = 1,n - 1
                 call la_slartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (icompq == 1) then
                    q(i + 2*n) = cs
                    q(i + 3*n) = sn
                 else if (icompq == 2) then
                    work(i) = cs
                    work(nm1 + i) = -sn
                 end if
              end do
           end if
           ! if icompq = 0, use la_slasdq to compute the singular values.
           if (icompq == 0) then
              ! ignore wstart, instead using work( 1 ), since the two vectors
              ! for cs and -sn above are added only if icompq == 2,
              ! and adding them exceeds documented work size of 4*n.
              call la_slasdq('U',0,n,0,0,0,d,e,vt,ldvt,u,ldu,u,ldu,work(1), &
                        info)
              go to 40
           end if
           ! if n is smaller than the minimum divide size smlsiz, then solve
           ! the problem with another solver.
           if (n <= smlsiz) then
              if (icompq == 2) then
                 call la_slaset('A',n,n,zero,one,u,ldu)
                 call la_slaset('A',n,n,zero,one,vt,ldvt)
                 call la_slasdq('U',0,n,n,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work( &
                           wstart),info)
              else if (icompq == 1) then
                 iu = 1
                 ivt = iu + n
                 call la_slaset('A',n,n,zero,one,q(iu + (qstart - 1)*n),n)
                 call la_slaset('A',n,n,zero,one,q(ivt + (qstart - 1)*n),n)
                 call la_slasdq('U',0,n,n,n,0,d,e,q(ivt + (qstart - 1)*n),n,q(iu + ( &
                           qstart - 1)*n),n,q(iu + (qstart - 1)*n),n,work(wstart),info)
              end if
              go to 40
           end if
           if (icompq == 2) then
              call la_slaset('A',n,n,zero,one,u,ldu)
              call la_slaset('A',n,n,zero,one,vt,ldvt)
           end if
           ! scale.
           orgnrm = la_slanst('M',n,d,e)
           if (orgnrm == zero) return
           call la_slascl('G',0,0,orgnrm,one,n,1,d,n,ierr)
           call la_slascl('G',0,0,orgnrm,one,nm1,1,e,nm1,ierr)
           eps = la_slamch('EPSILON')
           mlvl = int(log(real(n,KIND=sp)/real(smlsiz + 1,KIND=sp))/log(two),KIND=ilp) + &
                     1
           smlszp = smlsiz + 1
           if (icompq == 1) then
              iu = 1
              ivt = 1 + smlsiz
              difl = ivt + smlszp
              difr = difl + mlvl
              z = difr + mlvl*2
              ic = z + mlvl
              is = ic + 1
              poles = is + 1
              givnum = poles + 2*mlvl
              k = 1
              givptr = 2
              perm = 3
              givcol = perm + mlvl
           end if
           do i = 1,n
              if (abs(d(i)) < eps) then
                 d(i) = sign(eps,d(i))
              end if
           end do
           start = 1
           sqre = 0
           loop_30: do i = 1,nm1
              if ((abs(e(i)) < eps) .or. (i == nm1)) then
                 ! subproblem found. first determine its size and then
                 ! apply divide and conquer on it.
                 if (i < nm1) then
                    ! a subproblem with e(i) small for i < nm1.
                    nsize = i - start + 1
                 else if (abs(e(i)) >= eps) then
                    ! a subproblem with e(nm1) not too small but i = nm1.
                    nsize = n - start + 1
                 else
                    ! a subproblem with e(nm1) small. this implies an
                    ! 1-by-1 subproblem at d(n). solve this 1-by-1 problem
                    ! first.
                    nsize = i - start + 1
                    if (icompq == 2) then
                       u(n,n) = sign(one,d(n))
                       vt(n,n) = one
                    else if (icompq == 1) then
                       q(n + (qstart - 1)*n) = sign(one,d(n))
                       q(n + (smlsiz + qstart - 1)*n) = one
                    end if
                    d(n) = abs(d(n))
                 end if
                 if (icompq == 2) then
                    call la_slasd0(nsize,sqre,d(start),e(start),u(start,start), &
                              ldu,vt(start,start),ldvt,smlsiz,iwork,work(wstart),info)
                 else
                    call la_slasda(icompq,smlsiz,nsize,sqre,d(start),e(start),q( &
                    start + (iu + qstart - 2)*n),n,q(start + (ivt + qstart - 2)*n),iq(start + k*n),q( &
                     start + (difl + qstart - 2)*n),q(start + (difr + qstart - 2)*n),q(start + (z + &
                     qstart - 2)*n),q(start + (poles + qstart - 2)*n),iq(start + givptr*n),iq( &
                     start + givcol*n),n,iq(start + perm*n),q(start + (givnum + qstart - 2)*n),q( &
                     start + (ic + qstart - 2)*n),q(start + (is + qstart - 2)*n),work(wstart),iwork, &
                                info)
                 end if
                 if (info /= 0) then
                    return
                 end if
                 start = i + 1
              end if
           end do loop_30
           ! unscale
           call la_slascl('G',0,0,one,orgnrm,n,1,d,n,ierr)
           40 continue
           ! use selection sort to minimize swaps of singular vectors
           do ii = 2,n
              i = ii - 1
              kk = i
              p = d(i)
              do j = ii,n
                 if (d(j) > p) then
                    kk = j
                    p = d(j)
                 end if
              end do
              if (kk /= i) then
                 d(kk) = d(i)
                 d(i) = p
                 if (icompq == 1) then
                    iq(i) = kk
                 else if (icompq == 2) then
                    call la_sswap(n,u(1,i),1,u(1,kk),1)
                    call la_sswap(n,vt(i,1),ldvt,vt(kk,1),ldvt)
                 end if
              else if (icompq == 1) then
                 iq(i) = i
              end if
           end do
           ! if icompq = 1, use iq(n,1) as the indicator for uplo
           if (icompq == 1) then
              if (iuplo == 1) then
                 iq(n) = 1
              else
                 iq(n) = 0
              end if
           end if
           ! if b is lower bidiagonal, update u by those givens rotations
           ! which rotated b to be upper bidiagonal
           if ((iuplo == 2) .and. (icompq == 2)) call la_slasr('L','V','B',n,n,work(1) &
                     ,work(n),u,ldu)
           return
     end subroutine la_sbdsdc
     !> DBDSDC: computes the singular value decomposition (SVD) of a real
     !> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
     !> using a divide and conquer method, where S is a diagonal matrix
     !> with non-negative diagonal elements (the singular values of B), and
     !> U and VT are orthogonal matrices of left and right singular vectors,
     !> respectively. DBDSDC can be used to compute all singular values,
     !> and optionally, singular vectors or singular vectors in compact form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See DLASD3 for details.
     !> The code currently calls DLASDQ if singular values only are desired.
     !> However, it can be slightly modified to compute singular values
     !> using the divide and conquer method.

     pure subroutine la_dbdsdc(uplo,compq,n,d,e,u,ldu,vt,ldvt,q,iq,work,iwork, &
               info)
        use la_constants_dp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n
           ! Array Arguments
           integer(ilp),intent(out) :: iq(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: q(*),u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
        ! changed dimension statement in comment describing e from (n) to
        ! (n-1).  sven, 17 feb 05.
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: difl,difr,givcol,givnum,givptr,i,ic,icompq,ierr,ii,is,iu, &
           iuplo,ivt,j,k,kk,mlvl,nm1,nsize,perm,poles,qstart,smlsiz,smlszp,sqre, &
                     start,wstart,z
           real(dp) :: cs,eps,orgnrm,p,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (la_lsame(compq,'N')) then
              icompq = 0
           else if (la_lsame(compq,'P')) then
              icompq = 1
           else if (la_lsame(compq,'I')) then
              icompq = 2
           else
              icompq = -1
           end if
           if (iuplo == 0) then
              info = -1
           else if (icompq < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((ldu < 1) .or. ((icompq == 2) .and. (ldu < n))) then
              info = -7
           else if ((ldvt < 1) .or. ((icompq == 2) .and. (ldvt < n))) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DBDSDC',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'DBDSDC',' ',0,0,0,0)
           if (n == 1) then
              if (icompq == 1) then
                 q(1) = sign(one,d(1))
                 q(1 + smlsiz*n) = one
              else if (icompq == 2) then
                 u(1,1) = sign(one,d(1))
                 vt(1,1) = one
              end if
              d(1) = abs(d(1))
              return
           end if
           nm1 = n - 1
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           wstart = 1
           qstart = 3
           if (icompq == 1) then
              call la_dcopy(n,d,1,q(1),1)
              call la_dcopy(n - 1,e,1,q(n + 1),1)
           end if
           if (iuplo == 2) then
              qstart = 5
              if (icompq == 2) wstart = 2*n - 1
              do i = 1,n - 1
                 call la_dlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (icompq == 1) then
                    q(i + 2*n) = cs
                    q(i + 3*n) = sn
                 else if (icompq == 2) then
                    work(i) = cs
                    work(nm1 + i) = -sn
                 end if
              end do
           end if
           ! if icompq = 0, use la_dlasdq to compute the singular values.
           if (icompq == 0) then
              ! ignore wstart, instead using work( 1 ), since the two vectors
              ! for cs and -sn above are added only if icompq == 2,
              ! and adding them exceeds documented work size of 4*n.
              call la_dlasdq('U',0,n,0,0,0,d,e,vt,ldvt,u,ldu,u,ldu,work(1), &
                        info)
              go to 40
           end if
           ! if n is smaller than the minimum divide size smlsiz, then solve
           ! the problem with another solver.
           if (n <= smlsiz) then
              if (icompq == 2) then
                 call la_dlaset('A',n,n,zero,one,u,ldu)
                 call la_dlaset('A',n,n,zero,one,vt,ldvt)
                 call la_dlasdq('U',0,n,n,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work( &
                           wstart),info)
              else if (icompq == 1) then
                 iu = 1
                 ivt = iu + n
                 call la_dlaset('A',n,n,zero,one,q(iu + (qstart - 1)*n),n)
                 call la_dlaset('A',n,n,zero,one,q(ivt + (qstart - 1)*n),n)
                 call la_dlasdq('U',0,n,n,n,0,d,e,q(ivt + (qstart - 1)*n),n,q(iu + ( &
                           qstart - 1)*n),n,q(iu + (qstart - 1)*n),n,work(wstart),info)
              end if
              go to 40
           end if
           if (icompq == 2) then
              call la_dlaset('A',n,n,zero,one,u,ldu)
              call la_dlaset('A',n,n,zero,one,vt,ldvt)
           end if
           ! scale.
           orgnrm = la_dlanst('M',n,d,e)
           if (orgnrm == zero) return
           call la_dlascl('G',0,0,orgnrm,one,n,1,d,n,ierr)
           call la_dlascl('G',0,0,orgnrm,one,nm1,1,e,nm1,ierr)
           eps = (0.9e+0_dp)*la_dlamch('EPSILON')
           mlvl = int(log(real(n,KIND=dp)/real(smlsiz + 1,KIND=dp))/log(two),KIND=ilp) + &
                     1
           smlszp = smlsiz + 1
           if (icompq == 1) then
              iu = 1
              ivt = 1 + smlsiz
              difl = ivt + smlszp
              difr = difl + mlvl
              z = difr + mlvl*2
              ic = z + mlvl
              is = ic + 1
              poles = is + 1
              givnum = poles + 2*mlvl
              k = 1
              givptr = 2
              perm = 3
              givcol = perm + mlvl
           end if
           do i = 1,n
              if (abs(d(i)) < eps) then
                 d(i) = sign(eps,d(i))
              end if
           end do
           start = 1
           sqre = 0
           loop_30: do i = 1,nm1
              if ((abs(e(i)) < eps) .or. (i == nm1)) then
                 ! subproblem found. first determine its size and then
                 ! apply divide and conquer on it.
                 if (i < nm1) then
                    ! a subproblem with e(i) small for i < nm1.
                    nsize = i - start + 1
                 else if (abs(e(i)) >= eps) then
                    ! a subproblem with e(nm1) not too small but i = nm1.
                    nsize = n - start + 1
                 else
                    ! a subproblem with e(nm1) small. this implies an
                    ! 1-by-1 subproblem at d(n). solve this 1-by-1 problem
                    ! first.
                    nsize = i - start + 1
                    if (icompq == 2) then
                       u(n,n) = sign(one,d(n))
                       vt(n,n) = one
                    else if (icompq == 1) then
                       q(n + (qstart - 1)*n) = sign(one,d(n))
                       q(n + (smlsiz + qstart - 1)*n) = one
                    end if
                    d(n) = abs(d(n))
                 end if
                 if (icompq == 2) then
                    call la_dlasd0(nsize,sqre,d(start),e(start),u(start,start), &
                              ldu,vt(start,start),ldvt,smlsiz,iwork,work(wstart),info)
                 else
                    call la_dlasda(icompq,smlsiz,nsize,sqre,d(start),e(start),q( &
                    start + (iu + qstart - 2)*n),n,q(start + (ivt + qstart - 2)*n),iq(start + k*n),q( &
                     start + (difl + qstart - 2)*n),q(start + (difr + qstart - 2)*n),q(start + (z + &
                     qstart - 2)*n),q(start + (poles + qstart - 2)*n),iq(start + givptr*n),iq( &
                     start + givcol*n),n,iq(start + perm*n),q(start + (givnum + qstart - 2)*n),q( &
                     start + (ic + qstart - 2)*n),q(start + (is + qstart - 2)*n),work(wstart),iwork, &
                                info)
                 end if
                 if (info /= 0) then
                    return
                 end if
                 start = i + 1
              end if
           end do loop_30
           ! unscale
           call la_dlascl('G',0,0,one,orgnrm,n,1,d,n,ierr)
           40 continue
           ! use selection sort to minimize swaps of singular vectors
           do ii = 2,n
              i = ii - 1
              kk = i
              p = d(i)
              do j = ii,n
                 if (d(j) > p) then
                    kk = j
                    p = d(j)
                 end if
              end do
              if (kk /= i) then
                 d(kk) = d(i)
                 d(i) = p
                 if (icompq == 1) then
                    iq(i) = kk
                 else if (icompq == 2) then
                    call la_dswap(n,u(1,i),1,u(1,kk),1)
                    call la_dswap(n,vt(i,1),ldvt,vt(kk,1),ldvt)
                 end if
              else if (icompq == 1) then
                 iq(i) = i
              end if
           end do
           ! if icompq = 1, use iq(n,1) as the indicator for uplo
           if (icompq == 1) then
              if (iuplo == 1) then
                 iq(n) = 1
              else
                 iq(n) = 0
              end if
           end if
           ! if b is lower bidiagonal, update u by those givens rotations
           ! which rotated b to be upper bidiagonal
           if ((iuplo == 2) .and. (icompq == 2)) call la_dlasr('L','V','B',n,n,work(1) &
                     ,work(n),u,ldu)
           return
     end subroutine la_dbdsdc
#ifdef LA_WITH_XDP
     !> XBDSDC: computes the singular value decomposition (SVD) of a real
     !> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
     !> using a divide and conquer method, where S is a diagonal matrix
     !> with non-negative diagonal elements (the singular values of B), and
     !> U and VT are orthogonal matrices of left and right singular vectors,
     !> respectively. XBDSDC can be used to compute all singular values,
     !> and optionally, singular vectors or singular vectors in compact form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See XLASD3 for details.
     !> The code currently calls XLASDQ if singular values only are desired.
     !> However, it can be slightly modified to compute singular values
     !> using the divide and conquer method.

     pure subroutine la_xbdsdc(uplo,compq,n,d,e,u,ldu,vt,ldvt,q,iq,work,iwork, &
               info)
        use la_constants_xdp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n
           ! Array Arguments
           integer(ilp),intent(out) :: iq(*),iwork(*)
           real(xdp),intent(inout) :: d(*),e(*)
           real(xdp),intent(out) :: q(*),u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
        ! changed dimension statement in comment describing e from (n) to
        ! (n-1).  sven, 17 feb 05.
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: difl,difr,givcol,givnum,givptr,i,ic,icompq,ierr,ii,is,iu, &
           iuplo,ivt,j,k,kk,mlvl,nm1,nsize,perm,poles,qstart,smlsiz,smlszp,sqre, &
                     start,wstart,z
           real(xdp) :: cs,eps,orgnrm,p,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (la_lsame(compq,'N')) then
              icompq = 0
           else if (la_lsame(compq,'P')) then
              icompq = 1
           else if (la_lsame(compq,'I')) then
              icompq = 2
           else
              icompq = -1
           end if
           if (iuplo == 0) then
              info = -1
           else if (icompq < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((ldu < 1) .or. ((icompq == 2) .and. (ldu < n))) then
              info = -7
           else if ((ldvt < 1) .or. ((icompq == 2) .and. (ldvt < n))) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XBDSDC',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'XBDSDC',' ',0,0,0,0)
           if (n == 1) then
              if (icompq == 1) then
                 q(1) = sign(one,d(1))
                 q(1 + smlsiz*n) = one
              else if (icompq == 2) then
                 u(1,1) = sign(one,d(1))
                 vt(1,1) = one
              end if
              d(1) = abs(d(1))
              return
           end if
           nm1 = n - 1
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           wstart = 1
           qstart = 3
           if (icompq == 1) then
              call la_xcopy(n,d,1,q(1),1)
              call la_xcopy(n - 1,e,1,q(n + 1),1)
           end if
           if (iuplo == 2) then
              qstart = 5
              if (icompq == 2) wstart = 2*n - 1
              do i = 1,n - 1
                 call la_xlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (icompq == 1) then
                    q(i + 2*n) = cs
                    q(i + 3*n) = sn
                 else if (icompq == 2) then
                    work(i) = cs
                    work(nm1 + i) = -sn
                 end if
              end do
           end if
           ! if icompq = 0, use la_xlasdq to compute the singular values.
           if (icompq == 0) then
              ! ignore wstart, instead using work( 1 ), since the two vectors
              ! for cs and -sn above are added only if icompq == 2,
              ! and adding them exceeds documented work size of 4*n.
              call la_xlasdq('U',0,n,0,0,0,d,e,vt,ldvt,u,ldu,u,ldu,work(1), &
                        info)
              go to 40
           end if
           ! if n is smaller than the minimum divide size smlsiz, then solve
           ! the problem with another solver.
           if (n <= smlsiz) then
              if (icompq == 2) then
                 call la_xlaset('A',n,n,zero,one,u,ldu)
                 call la_xlaset('A',n,n,zero,one,vt,ldvt)
                 call la_xlasdq('U',0,n,n,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work( &
                           wstart),info)
              else if (icompq == 1) then
                 iu = 1
                 ivt = iu + n
                 call la_xlaset('A',n,n,zero,one,q(iu + (qstart - 1)*n),n)
                 call la_xlaset('A',n,n,zero,one,q(ivt + (qstart - 1)*n),n)
                 call la_xlasdq('U',0,n,n,n,0,d,e,q(ivt + (qstart - 1)*n),n,q(iu + ( &
                           qstart - 1)*n),n,q(iu + (qstart - 1)*n),n,work(wstart),info)
              end if
              go to 40
           end if
           if (icompq == 2) then
              call la_xlaset('A',n,n,zero,one,u,ldu)
              call la_xlaset('A',n,n,zero,one,vt,ldvt)
           end if
           ! scale.
           orgnrm = la_xlanst('M',n,d,e)
           if (orgnrm == zero) return
           call la_xlascl('G',0,0,orgnrm,one,n,1,d,n,ierr)
           call la_xlascl('G',0,0,orgnrm,one,nm1,1,e,nm1,ierr)
           eps = (0.9e+0_xdp)*la_xlamch('EPSILON')
           mlvl = int(log(real(n,KIND=xdp)/real(smlsiz + 1,KIND=xdp))/log(two),KIND=ilp) + &
                     1
           smlszp = smlsiz + 1
           if (icompq == 1) then
              iu = 1
              ivt = 1 + smlsiz
              difl = ivt + smlszp
              difr = difl + mlvl
              z = difr + mlvl*2
              ic = z + mlvl
              is = ic + 1
              poles = is + 1
              givnum = poles + 2*mlvl
              k = 1
              givptr = 2
              perm = 3
              givcol = perm + mlvl
           end if
           do i = 1,n
              if (abs(d(i)) < eps) then
                 d(i) = sign(eps,d(i))
              end if
           end do
           start = 1
           sqre = 0
           loop_30: do i = 1,nm1
              if ((abs(e(i)) < eps) .or. (i == nm1)) then
                 ! subproblem found. first determine its size and then
                 ! apply divide and conquer on it.
                 if (i < nm1) then
                    ! a subproblem with e(i) small for i < nm1.
                    nsize = i - start + 1
                 else if (abs(e(i)) >= eps) then
                    ! a subproblem with e(nm1) not too small but i = nm1.
                    nsize = n - start + 1
                 else
                    ! a subproblem with e(nm1) small. this implies an
                    ! 1-by-1 subproblem at d(n). solve this 1-by-1 problem
                    ! first.
                    nsize = i - start + 1
                    if (icompq == 2) then
                       u(n,n) = sign(one,d(n))
                       vt(n,n) = one
                    else if (icompq == 1) then
                       q(n + (qstart - 1)*n) = sign(one,d(n))
                       q(n + (smlsiz + qstart - 1)*n) = one
                    end if
                    d(n) = abs(d(n))
                 end if
                 if (icompq == 2) then
                    call la_xlasd0(nsize,sqre,d(start),e(start),u(start,start), &
                              ldu,vt(start,start),ldvt,smlsiz,iwork,work(wstart),info)
                 else
                    call la_xlasda(icompq,smlsiz,nsize,sqre,d(start),e(start),q( &
                    start + (iu + qstart - 2)*n),n,q(start + (ivt + qstart - 2)*n),iq(start + k*n),q( &
                     start + (difl + qstart - 2)*n),q(start + (difr + qstart - 2)*n),q(start + (z + &
                     qstart - 2)*n),q(start + (poles + qstart - 2)*n),iq(start + givptr*n),iq( &
                     start + givcol*n),n,iq(start + perm*n),q(start + (givnum + qstart - 2)*n),q( &
                     start + (ic + qstart - 2)*n),q(start + (is + qstart - 2)*n),work(wstart),iwork, &
                                info)
                 end if
                 if (info /= 0) then
                    return
                 end if
                 start = i + 1
              end if
           end do loop_30
           ! unscale
           call la_xlascl('G',0,0,one,orgnrm,n,1,d,n,ierr)
           40 continue
           ! use selection sort to minimize swaps of singular vectors
           do ii = 2,n
              i = ii - 1
              kk = i
              p = d(i)
              do j = ii,n
                 if (d(j) > p) then
                    kk = j
                    p = d(j)
                 end if
              end do
              if (kk /= i) then
                 d(kk) = d(i)
                 d(i) = p
                 if (icompq == 1) then
                    iq(i) = kk
                 else if (icompq == 2) then
                    call la_xswap(n,u(1,i),1,u(1,kk),1)
                    call la_xswap(n,vt(i,1),ldvt,vt(kk,1),ldvt)
                 end if
              else if (icompq == 1) then
                 iq(i) = i
              end if
           end do
           ! if icompq = 1, use iq(n,1) as the indicator for uplo
           if (icompq == 1) then
              if (iuplo == 1) then
                 iq(n) = 1
              else
                 iq(n) = 0
              end if
           end if
           ! if b is lower bidiagonal, update u by those givens rotations
           ! which rotated b to be upper bidiagonal
           if ((iuplo == 2) .and. (icompq == 2)) call la_xlasr('L','V','B',n,n,work(1) &
                     ,work(n),u,ldu)
           return
     end subroutine la_xbdsdc
#endif
#ifdef LA_WITH_QP
     !> QBDSDC: computes the singular value decomposition (SVD) of a real
     !> N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
     !> using a divide and conquer method, where S is a diagonal matrix
     !> with non-negative diagonal elements (the singular values of B), and
     !> U and VT are orthogonal matrices of left and right singular vectors,
     !> respectively. QBDSDC can be used to compute all singular values,
     !> and optionally, singular vectors or singular vectors in compact form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See QLASD3 for details.
     !> The code currently calls QLASDQ if singular values only are desired.
     !> However, it can be slightly modified to compute singular values
     !> using the divide and conquer method.

     pure subroutine la_qbdsdc(uplo,compq,n,d,e,u,ldu,vt,ldvt,q,iq,work,iwork, &
               info)
        use la_constants_qp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compq,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n
           ! Array Arguments
           integer(ilp),intent(out) :: iq(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: q(*),u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
        ! changed dimension statement in comment describing e from (n) to
        ! (n-1).  sven, 17 feb 05.
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: difl,difr,givcol,givnum,givptr,i,ic,icompq,ierr,ii,is,iu, &
           iuplo,ivt,j,k,kk,mlvl,nm1,nsize,perm,poles,qstart,smlsiz,smlszp,sqre, &
                     start,wstart,z
           real(qp) :: cs,eps,orgnrm,p,r,sn
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,sign
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (la_lsame(compq,'N')) then
              icompq = 0
           else if (la_lsame(compq,'P')) then
              icompq = 1
           else if (la_lsame(compq,'I')) then
              icompq = 2
           else
              icompq = -1
           end if
           if (iuplo == 0) then
              info = -1
           else if (icompq < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((ldu < 1) .or. ((icompq == 2) .and. (ldu < n))) then
              info = -7
           else if ((ldvt < 1) .or. ((icompq == 2) .and. (ldvt < n))) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QBDSDC',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           smlsiz = la_ilaenv(9,'QBDSDC',' ',0,0,0,0)
           if (n == 1) then
              if (icompq == 1) then
                 q(1) = sign(one,d(1))
                 q(1 + smlsiz*n) = one
              else if (icompq == 2) then
                 u(1,1) = sign(one,d(1))
                 vt(1,1) = one
              end if
              d(1) = abs(d(1))
              return
           end if
           nm1 = n - 1
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           wstart = 1
           qstart = 3
           if (icompq == 1) then
              call la_qcopy(n,d,1,q(1),1)
              call la_qcopy(n - 1,e,1,q(n + 1),1)
           end if
           if (iuplo == 2) then
              qstart = 5
              if (icompq == 2) wstart = 2*n - 1
              do i = 1,n - 1
                 call la_qlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (icompq == 1) then
                    q(i + 2*n) = cs
                    q(i + 3*n) = sn
                 else if (icompq == 2) then
                    work(i) = cs
                    work(nm1 + i) = -sn
                 end if
              end do
           end if
           ! if icompq = 0, use la_qlasdq to compute the singular values.
           if (icompq == 0) then
              ! ignore wstart, instead using work( 1 ), since the two vectors
              ! for cs and -sn above are added only if icompq == 2,
              ! and adding them exceeds documented work size of 4*n.
              call la_qlasdq('U',0,n,0,0,0,d,e,vt,ldvt,u,ldu,u,ldu,work(1), &
                        info)
              go to 40
           end if
           ! if n is smaller than the minimum divide size smlsiz, then solve
           ! the problem with another solver.
           if (n <= smlsiz) then
              if (icompq == 2) then
                 call la_qlaset('A',n,n,zero,one,u,ldu)
                 call la_qlaset('A',n,n,zero,one,vt,ldvt)
                 call la_qlasdq('U',0,n,n,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work( &
                           wstart),info)
              else if (icompq == 1) then
                 iu = 1
                 ivt = iu + n
                 call la_qlaset('A',n,n,zero,one,q(iu + (qstart - 1)*n),n)
                 call la_qlaset('A',n,n,zero,one,q(ivt + (qstart - 1)*n),n)
                 call la_qlasdq('U',0,n,n,n,0,d,e,q(ivt + (qstart - 1)*n),n,q(iu + ( &
                           qstart - 1)*n),n,q(iu + (qstart - 1)*n),n,work(wstart),info)
              end if
              go to 40
           end if
           if (icompq == 2) then
              call la_qlaset('A',n,n,zero,one,u,ldu)
              call la_qlaset('A',n,n,zero,one,vt,ldvt)
           end if
           ! scale.
           orgnrm = la_qlanst('M',n,d,e)
           if (orgnrm == zero) return
           call la_qlascl('G',0,0,orgnrm,one,n,1,d,n,ierr)
           call la_qlascl('G',0,0,orgnrm,one,nm1,1,e,nm1,ierr)
           eps = (0.9e+0_qp)*la_qlamch('EPSILON')
           mlvl = int(log(real(n,KIND=qp)/real(smlsiz + 1,KIND=qp))/log(two),KIND=ilp) + &
                     1
           smlszp = smlsiz + 1
           if (icompq == 1) then
              iu = 1
              ivt = 1 + smlsiz
              difl = ivt + smlszp
              difr = difl + mlvl
              z = difr + mlvl*2
              ic = z + mlvl
              is = ic + 1
              poles = is + 1
              givnum = poles + 2*mlvl
              k = 1
              givptr = 2
              perm = 3
              givcol = perm + mlvl
           end if
           do i = 1,n
              if (abs(d(i)) < eps) then
                 d(i) = sign(eps,d(i))
              end if
           end do
           start = 1
           sqre = 0
           loop_30: do i = 1,nm1
              if ((abs(e(i)) < eps) .or. (i == nm1)) then
                 ! subproblem found. first determine its size and then
                 ! apply divide and conquer on it.
                 if (i < nm1) then
                    ! a subproblem with e(i) small for i < nm1.
                    nsize = i - start + 1
                 else if (abs(e(i)) >= eps) then
                    ! a subproblem with e(nm1) not too small but i = nm1.
                    nsize = n - start + 1
                 else
                    ! a subproblem with e(nm1) small. this implies an
                    ! 1-by-1 subproblem at d(n). solve this 1-by-1 problem
                    ! first.
                    nsize = i - start + 1
                    if (icompq == 2) then
                       u(n,n) = sign(one,d(n))
                       vt(n,n) = one
                    else if (icompq == 1) then
                       q(n + (qstart - 1)*n) = sign(one,d(n))
                       q(n + (smlsiz + qstart - 1)*n) = one
                    end if
                    d(n) = abs(d(n))
                 end if
                 if (icompq == 2) then
                    call la_qlasd0(nsize,sqre,d(start),e(start),u(start,start), &
                              ldu,vt(start,start),ldvt,smlsiz,iwork,work(wstart),info)
                 else
                    call la_qlasda(icompq,smlsiz,nsize,sqre,d(start),e(start),q( &
                    start + (iu + qstart - 2)*n),n,q(start + (ivt + qstart - 2)*n),iq(start + k*n),q( &
                     start + (difl + qstart - 2)*n),q(start + (difr + qstart - 2)*n),q(start + (z + &
                     qstart - 2)*n),q(start + (poles + qstart - 2)*n),iq(start + givptr*n),iq( &
                     start + givcol*n),n,iq(start + perm*n),q(start + (givnum + qstart - 2)*n),q( &
                     start + (ic + qstart - 2)*n),q(start + (is + qstart - 2)*n),work(wstart),iwork, &
                                info)
                 end if
                 if (info /= 0) then
                    return
                 end if
                 start = i + 1
              end if
           end do loop_30
           ! unscale
           call la_qlascl('G',0,0,one,orgnrm,n,1,d,n,ierr)
           40 continue
           ! use selection sort to minimize swaps of singular vectors
           do ii = 2,n
              i = ii - 1
              kk = i
              p = d(i)
              do j = ii,n
                 if (d(j) > p) then
                    kk = j
                    p = d(j)
                 end if
              end do
              if (kk /= i) then
                 d(kk) = d(i)
                 d(i) = p
                 if (icompq == 1) then
                    iq(i) = kk
                 else if (icompq == 2) then
                    call la_qswap(n,u(1,i),1,u(1,kk),1)
                    call la_qswap(n,vt(i,1),ldvt,vt(kk,1),ldvt)
                 end if
              else if (icompq == 1) then
                 iq(i) = i
              end if
           end do
           ! if icompq = 1, use iq(n,1) as the indicator for uplo
           if (icompq == 1) then
              if (iuplo == 1) then
                 iq(n) = 1
              else
                 iq(n) = 0
              end if
           end if
           ! if b is lower bidiagonal, update u by those givens rotations
           ! which rotated b to be upper bidiagonal
           if ((iuplo == 2) .and. (icompq == 2)) call la_qlasr('L','V','B',n,n,work(1) &
                     ,work(n),u,ldu)
           return
     end subroutine la_qbdsdc
#endif

     !> Using a divide and conquer approach, SLASD0: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M
     !> matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
     !> The algorithm computes orthogonal matrices U and VT such that
     !> B = U * S * VT. The singular values S are overwritten on D.
     !> A related subroutine, SLASDA, computes only the singular values,
     !> and optionally, the singular vectors in compact form.

     pure subroutine la_slasd0(n,sqre,d,e,u,ldu,vt,ldvt,smlsiz,iwork,work,info)

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n,smlsiz,sqre
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqc,im1,inode,itemp,iwk,j,lf,ll,lvl,m,ncc, &
                      nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,sqrei
           real(sp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           end if
           m = n + sqre
           if (ldu < n) then
              info = -6
           else if (ldvt < m) then
              info = -8
           else if (smlsiz < 3) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SLASD0',-info)
              return
           end if
           ! if the input matrix is too small, call la_slasdq to find the svd.
           if (n <= smlsiz) then
              call la_slasdq('U',sqre,n,m,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work, &
                        info)
              return
           end if
           ! set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           call la_slasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_slasdq.
           ndb1 = (nd + 1)/2
           ncc = 0
           loop_30: do i = ndb1,nd
           ! ic : center row of each node
           ! nl : number of rows of left  subproblem
           ! nr : number of rows of right subproblem
           ! nlf: starting row of the left   subproblem
           ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nrp1 = nr + 1
              nlf = ic - nl
              nrf = ic + 1
              sqrei = 1
              call la_slasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,nlf) &
                        ,ldvt,u(nlf,nlf),ldu,u(nlf,nlf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + nlf - 2
              do j = 1,nl
                 iwork(itemp + j) = j
              end do
              if (i == nd) then
                 sqrei = sqre
              else
                 sqrei = 1
              end if
              nrp1 = nr + sqrei
              call la_slasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,nrf) &
                        ,ldvt,u(nrf,nrf),ldu,u(nrf,nrf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + ic
              do j = 1,nr
                 iwork(itemp + j - 1) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           loop_50: do lvl = nlvl,1,-1
              ! find the first node lf and last node ll on the
              ! current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 if ((sqre == 0) .and. (i == ll)) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 idxqc = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 call la_slasd1(nl,nr,sqrei,d(nlf),alpha,beta,u(nlf,nlf),ldu,vt( &
                           nlf,nlf),ldvt,iwork(idxqc),iwork(iwk),work,info)
              ! report the possible convergence failure.
                 if (info /= 0) then
                    return
                 end if
              end do
           end do loop_50
           return
     end subroutine la_slasd0
     !> Using a divide and conquer approach, DLASD0: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M
     !> matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
     !> The algorithm computes orthogonal matrices U and VT such that
     !> B = U * S * VT. The singular values S are overwritten on D.
     !> A related subroutine, DLASDA, computes only the singular values,
     !> and optionally, the singular vectors in compact form.

     pure subroutine la_dlasd0(n,sqre,d,e,u,ldu,vt,ldvt,smlsiz,iwork,work,info)

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n,smlsiz,sqre
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqc,im1,inode,itemp,iwk,j,lf,ll,lvl,m,ncc, &
                      nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,sqrei
           real(dp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           end if
           m = n + sqre
           if (ldu < n) then
              info = -6
           else if (ldvt < m) then
              info = -8
           else if (smlsiz < 3) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DLASD0',-info)
              return
           end if
           ! if the input matrix is too small, call la_dlasdq to find the svd.
           if (n <= smlsiz) then
              call la_dlasdq('U',sqre,n,m,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work, &
                        info)
              return
           end if
           ! set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           call la_dlasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_dlasdq.
           ndb1 = (nd + 1)/2
           ncc = 0
           loop_30: do i = ndb1,nd
           ! ic : center row of each node
           ! nl : number of rows of left  subproblem
           ! nr : number of rows of right subproblem
           ! nlf: starting row of the left   subproblem
           ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nrp1 = nr + 1
              nlf = ic - nl
              nrf = ic + 1
              sqrei = 1
              call la_dlasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,nlf) &
                        ,ldvt,u(nlf,nlf),ldu,u(nlf,nlf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + nlf - 2
              do j = 1,nl
                 iwork(itemp + j) = j
              end do
              if (i == nd) then
                 sqrei = sqre
              else
                 sqrei = 1
              end if
              nrp1 = nr + sqrei
              call la_dlasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,nrf) &
                        ,ldvt,u(nrf,nrf),ldu,u(nrf,nrf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + ic
              do j = 1,nr
                 iwork(itemp + j - 1) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           loop_50: do lvl = nlvl,1,-1
              ! find the first node lf and last node ll on the
              ! current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 if ((sqre == 0) .and. (i == ll)) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 idxqc = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 call la_dlasd1(nl,nr,sqrei,d(nlf),alpha,beta,u(nlf,nlf),ldu,vt( &
                           nlf,nlf),ldvt,iwork(idxqc),iwork(iwk),work,info)
              ! report the possible convergence failure.
                 if (info /= 0) then
                    return
                 end if
              end do
           end do loop_50
           return
     end subroutine la_dlasd0
#ifdef LA_WITH_XDP
     !> Using a divide and conquer approach, XLASD0: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M
     !> matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
     !> The algorithm computes orthogonal matrices U and VT such that
     !> B = U * S * VT. The singular values S are overwritten on D.
     !> A related subroutine, XLASDA, computes only the singular values,
     !> and optionally, the singular vectors in compact form.

     pure subroutine la_xlasd0(n,sqre,d,e,u,ldu,vt,ldvt,smlsiz,iwork,work,info)

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n,smlsiz,sqre
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(inout) :: d(*),e(*)
           real(xdp),intent(out) :: u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqc,im1,inode,itemp,iwk,j,lf,ll,lvl,m,ncc, &
                      nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,sqrei
           real(xdp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           end if
           m = n + sqre
           if (ldu < n) then
              info = -6
           else if (ldvt < m) then
              info = -8
           else if (smlsiz < 3) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XLASD0',-info)
              return
           end if
           ! if the input matrix is too small, call la_xlasdq to find the svd.
           if (n <= smlsiz) then
              call la_xlasdq('U',sqre,n,m,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work, &
                        info)
              return
           end if
           ! set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           call la_xlasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_xlasdq.
           ndb1 = (nd + 1)/2
           ncc = 0
           loop_30: do i = ndb1,nd
           ! ic : center row of each node
           ! nl : number of rows of left  subproblem
           ! nr : number of rows of right subproblem
           ! nlf: starting row of the left   subproblem
           ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nrp1 = nr + 1
              nlf = ic - nl
              nrf = ic + 1
              sqrei = 1
              call la_xlasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,nlf) &
                        ,ldvt,u(nlf,nlf),ldu,u(nlf,nlf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + nlf - 2
              do j = 1,nl
                 iwork(itemp + j) = j
              end do
              if (i == nd) then
                 sqrei = sqre
              else
                 sqrei = 1
              end if
              nrp1 = nr + sqrei
              call la_xlasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,nrf) &
                        ,ldvt,u(nrf,nrf),ldu,u(nrf,nrf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + ic
              do j = 1,nr
                 iwork(itemp + j - 1) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           loop_50: do lvl = nlvl,1,-1
              ! find the first node lf and last node ll on the
              ! current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 if ((sqre == 0) .and. (i == ll)) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 idxqc = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 call la_xlasd1(nl,nr,sqrei,d(nlf),alpha,beta,u(nlf,nlf),ldu,vt( &
                           nlf,nlf),ldvt,iwork(idxqc),iwork(iwk),work,info)
              ! report the possible convergence failure.
                 if (info /= 0) then
                    return
                 end if
              end do
           end do loop_50
           return
     end subroutine la_xlasd0
#endif
#ifdef LA_WITH_QP
     !> Using a divide and conquer approach, QLASD0: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M
     !> matrix B with diagonal D and offdiagonal E, where M = N + SQRE.
     !> The algorithm computes orthogonal matrices U and VT such that
     !> B = U * S * VT. The singular values S are overwritten on D.
     !> A related subroutine, QLASDA, computes only the singular values,
     !> and optionally, the singular vectors in compact form.

     pure subroutine la_qlasd0(n,sqre,d,e,u,ldu,vt,ldvt,smlsiz,iwork,work,info)

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldu,ldvt,n,smlsiz,sqre
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: u(ldu,*),vt(ldvt,*),work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqc,im1,inode,itemp,iwk,j,lf,ll,lvl,m,ncc, &
                      nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,sqrei
           real(qp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           end if
           m = n + sqre
           if (ldu < n) then
              info = -6
           else if (ldvt < m) then
              info = -8
           else if (smlsiz < 3) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QLASD0',-info)
              return
           end if
           ! if the input matrix is too small, call la_qlasdq to find the svd.
           if (n <= smlsiz) then
              call la_qlasdq('U',sqre,n,m,n,0,d,e,vt,ldvt,u,ldu,u,ldu,work, &
                        info)
              return
           end if
           ! set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           call la_qlasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_qlasdq.
           ndb1 = (nd + 1)/2
           ncc = 0
           loop_30: do i = ndb1,nd
           ! ic : center row of each node
           ! nl : number of rows of left  subproblem
           ! nr : number of rows of right subproblem
           ! nlf: starting row of the left   subproblem
           ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nrp1 = nr + 1
              nlf = ic - nl
              nrf = ic + 1
              sqrei = 1
              call la_qlasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,nlf) &
                        ,ldvt,u(nlf,nlf),ldu,u(nlf,nlf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + nlf - 2
              do j = 1,nl
                 iwork(itemp + j) = j
              end do
              if (i == nd) then
                 sqrei = sqre
              else
                 sqrei = 1
              end if
              nrp1 = nr + sqrei
              call la_qlasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,nrf) &
                        ,ldvt,u(nrf,nrf),ldu,u(nrf,nrf),ldu,work,info)
              if (info /= 0) then
                 return
              end if
              itemp = idxq + ic
              do j = 1,nr
                 iwork(itemp + j - 1) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           loop_50: do lvl = nlvl,1,-1
              ! find the first node lf and last node ll on the
              ! current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 if ((sqre == 0) .and. (i == ll)) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 idxqc = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 call la_qlasd1(nl,nr,sqrei,d(nlf),alpha,beta,u(nlf,nlf),ldu,vt( &
                           nlf,nlf),ldvt,iwork(idxqc),iwork(iwk),work,info)
              ! report the possible convergence failure.
                 if (info /= 0) then
                    return
                 end if
              end do
           end do loop_50
           return
     end subroutine la_qlasd0
#endif

     !> Using a divide and conquer approach, SLASDA: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
     !> B with diagonal D and offdiagonal E, where M = N + SQRE. The
     !> algorithm computes the singular values in the SVD B = U * S * VT.
     !> The orthogonal matrices U and VT are optionally computed in
     !> compact form.
     !> A related subroutine, SLASD0, computes the singular values and
     !> the singular vectors in explicit form.

     pure subroutine la_slasda(icompq,smlsiz,n,sqre,d,e,u,ldu,vt,k,difl,difr,z, &
               poles,givptr,givcol,ldgcol,perm,givnum,c,s,work,iwork,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldgcol,ldu,n,smlsiz,sqre
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),givptr(*),iwork(*),k(*),perm(ldgcol, &
                     *)
           real(sp),intent(out) :: c(*),difl(ldu,*),difr(ldu,*),givnum(ldu,*),poles(ldu,*), &
                     s(*),u(ldu,*),vt(ldu,*),work(*),z(ldu,*)
           real(sp),intent(inout) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqi,im1,inode,itemp,iwk,j,lf,ll,lvl,lvl2, &
           m,ncc,nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,nru,nwork1, &
                     nwork2,smlszp,sqrei,vf,vfi,vl,vli
           real(sp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (smlsiz < 3) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldu < (n + sqre)) then
              info = -8
           else if (ldgcol < n) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('SLASDA',-info)
              return
           end if
           m = n + sqre
           ! if the input matrix is too small, call la_slasdq to find the svd.
           if (n <= smlsiz) then
              if (icompq == 0) then
                 call la_slasdq('U',sqre,n,0,0,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              else
                 call la_slasdq('U',sqre,n,m,n,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              end if
              return
           end if
           ! book-keeping and  set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           ncc = 0
           nru = 0
           smlszp = smlsiz + 1
           vf = 1
           vl = vf + m
           nwork1 = vl + m
           nwork2 = nwork1 + smlszp*smlszp
           call la_slasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_slasdq.
           ndb1 = (nd + 1)/2
           loop_30: do i = ndb1,nd
              ! ic : center row of each node
              ! nl : number of rows of left  subproblem
              ! nr : number of rows of right subproblem
              ! nlf: starting row of the left   subproblem
              ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nlf = ic - nl
              nrf = ic + 1
              idxqi = idxq + nlf - 2
              vfi = vf + nlf - 1
              vli = vl + nlf - 1
              sqrei = 1
              if (icompq == 0) then
                 call la_slaset('A',nlp1,nlp1,zero,one,work(nwork1),smlszp)
                 call la_slasdq('U',sqrei,nl,nlp1,nru,ncc,d(nlf),e(nlf),work( &
                 nwork1),smlszp,work(nwork2),nl,work(nwork2),nl,work(nwork2),info)

                 itemp = nwork1 + nl*smlszp
                 call la_scopy(nlp1,work(nwork1),1,work(vfi),1)
                 call la_scopy(nlp1,work(itemp),1,work(vli),1)
              else
                 call la_slaset('A',nl,nl,zero,one,u(nlf,1),ldu)
                 call la_slaset('A',nlp1,nlp1,zero,one,vt(nlf,1),ldu)
                 call la_slasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,1 &
                           ),ldu,u(nlf,1),ldu,u(nlf,1),ldu,work(nwork1),info)
                 call la_scopy(nlp1,vt(nlf,1),1,work(vfi),1)
                 call la_scopy(nlp1,vt(nlf,nlp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nl
                 iwork(idxqi + j) = j
              end do
              if ((i == nd) .and. (sqre == 0)) then
                 sqrei = 0
              else
                 sqrei = 1
              end if
              idxqi = idxqi + nlp1
              vfi = vfi + nlp1
              vli = vli + nlp1
              nrp1 = nr + sqrei
              if (icompq == 0) then
                 call la_slaset('A',nrp1,nrp1,zero,one,work(nwork1),smlszp)
                 call la_slasdq('U',sqrei,nr,nrp1,nru,ncc,d(nrf),e(nrf),work( &
                 nwork1),smlszp,work(nwork2),nr,work(nwork2),nr,work(nwork2),info)

                 itemp = nwork1 + (nrp1 - 1)*smlszp
                 call la_scopy(nrp1,work(nwork1),1,work(vfi),1)
                 call la_scopy(nrp1,work(itemp),1,work(vli),1)
              else
                 call la_slaset('A',nr,nr,zero,one,u(nrf,1),ldu)
                 call la_slaset('A',nrp1,nrp1,zero,one,vt(nrf,1),ldu)
                 call la_slasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,1 &
                           ),ldu,u(nrf,1),ldu,u(nrf,1),ldu,work(nwork1),info)
                 call la_scopy(nrp1,vt(nrf,1),1,work(vfi),1)
                 call la_scopy(nrp1,vt(nrf,nrp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nr
                 iwork(idxqi + j) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           j = 2**nlvl
           loop_50: do lvl = nlvl,1,-1
              lvl2 = lvl*2 - 1
              ! find the first node lf and last node ll on
              ! the current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              loop_40: do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 nrf = ic + 1
                 if (i == ll) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 vfi = vf + nlf - 1
                 vli = vl + nlf - 1
                 idxqi = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 if (icompq == 0) then
                    call la_slasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm,givptr(1),givcol,ldgcol,givnum,ldu, &
                    poles,difl,difr,z,k(1),c(1),s(1),work(nwork1),iwork(iwk), &
                              info)
                 else
                    j = j - 1
                    call la_slasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm(nlf,lvl),givptr(j),givcol(nlf,lvl2), &
                     ldgcol,givnum(nlf,lvl2),ldu,poles(nlf,lvl2),difl(nlf,lvl),difr( &
                     nlf,lvl2),z(nlf,lvl),k(j),c(j),s(j),work(nwork1),iwork(iwk &
                               ),info)
                 end if
                 if (info /= 0) then
                    return
                 end if
              end do loop_40
           end do loop_50
           return
     end subroutine la_slasda
     !> Using a divide and conquer approach, DLASDA: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
     !> B with diagonal D and offdiagonal E, where M = N + SQRE. The
     !> algorithm computes the singular values in the SVD B = U * S * VT.
     !> The orthogonal matrices U and VT are optionally computed in
     !> compact form.
     !> A related subroutine, DLASD0, computes the singular values and
     !> the singular vectors in explicit form.

     pure subroutine la_dlasda(icompq,smlsiz,n,sqre,d,e,u,ldu,vt,k,difl,difr,z, &
               poles,givptr,givcol,ldgcol,perm,givnum,c,s,work,iwork,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldgcol,ldu,n,smlsiz,sqre
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),givptr(*),iwork(*),k(*),perm(ldgcol, &
                     *)
           real(dp),intent(out) :: c(*),difl(ldu,*),difr(ldu,*),givnum(ldu,*),poles(ldu,*), &
                     s(*),u(ldu,*),vt(ldu,*),work(*),z(ldu,*)
           real(dp),intent(inout) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqi,im1,inode,itemp,iwk,j,lf,ll,lvl,lvl2, &
           m,ncc,nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,nru,nwork1, &
                     nwork2,smlszp,sqrei,vf,vfi,vl,vli
           real(dp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (smlsiz < 3) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldu < (n + sqre)) then
              info = -8
           else if (ldgcol < n) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('DLASDA',-info)
              return
           end if
           m = n + sqre
           ! if the input matrix is too small, call la_dlasdq to find the svd.
           if (n <= smlsiz) then
              if (icompq == 0) then
                 call la_dlasdq('U',sqre,n,0,0,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              else
                 call la_dlasdq('U',sqre,n,m,n,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              end if
              return
           end if
           ! book-keeping and  set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           ncc = 0
           nru = 0
           smlszp = smlsiz + 1
           vf = 1
           vl = vf + m
           nwork1 = vl + m
           nwork2 = nwork1 + smlszp*smlszp
           call la_dlasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_dlasdq.
           ndb1 = (nd + 1)/2
           loop_30: do i = ndb1,nd
              ! ic : center row of each node
              ! nl : number of rows of left  subproblem
              ! nr : number of rows of right subproblem
              ! nlf: starting row of the left   subproblem
              ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nlf = ic - nl
              nrf = ic + 1
              idxqi = idxq + nlf - 2
              vfi = vf + nlf - 1
              vli = vl + nlf - 1
              sqrei = 1
              if (icompq == 0) then
                 call la_dlaset('A',nlp1,nlp1,zero,one,work(nwork1),smlszp)
                 call la_dlasdq('U',sqrei,nl,nlp1,nru,ncc,d(nlf),e(nlf),work( &
                 nwork1),smlszp,work(nwork2),nl,work(nwork2),nl,work(nwork2),info)

                 itemp = nwork1 + nl*smlszp
                 call la_dcopy(nlp1,work(nwork1),1,work(vfi),1)
                 call la_dcopy(nlp1,work(itemp),1,work(vli),1)
              else
                 call la_dlaset('A',nl,nl,zero,one,u(nlf,1),ldu)
                 call la_dlaset('A',nlp1,nlp1,zero,one,vt(nlf,1),ldu)
                 call la_dlasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,1 &
                           ),ldu,u(nlf,1),ldu,u(nlf,1),ldu,work(nwork1),info)
                 call la_dcopy(nlp1,vt(nlf,1),1,work(vfi),1)
                 call la_dcopy(nlp1,vt(nlf,nlp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nl
                 iwork(idxqi + j) = j
              end do
              if ((i == nd) .and. (sqre == 0)) then
                 sqrei = 0
              else
                 sqrei = 1
              end if
              idxqi = idxqi + nlp1
              vfi = vfi + nlp1
              vli = vli + nlp1
              nrp1 = nr + sqrei
              if (icompq == 0) then
                 call la_dlaset('A',nrp1,nrp1,zero,one,work(nwork1),smlszp)
                 call la_dlasdq('U',sqrei,nr,nrp1,nru,ncc,d(nrf),e(nrf),work( &
                 nwork1),smlszp,work(nwork2),nr,work(nwork2),nr,work(nwork2),info)

                 itemp = nwork1 + (nrp1 - 1)*smlszp
                 call la_dcopy(nrp1,work(nwork1),1,work(vfi),1)
                 call la_dcopy(nrp1,work(itemp),1,work(vli),1)
              else
                 call la_dlaset('A',nr,nr,zero,one,u(nrf,1),ldu)
                 call la_dlaset('A',nrp1,nrp1,zero,one,vt(nrf,1),ldu)
                 call la_dlasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,1 &
                           ),ldu,u(nrf,1),ldu,u(nrf,1),ldu,work(nwork1),info)
                 call la_dcopy(nrp1,vt(nrf,1),1,work(vfi),1)
                 call la_dcopy(nrp1,vt(nrf,nrp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nr
                 iwork(idxqi + j) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           j = 2**nlvl
           loop_50: do lvl = nlvl,1,-1
              lvl2 = lvl*2 - 1
              ! find the first node lf and last node ll on
              ! the current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              loop_40: do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 nrf = ic + 1
                 if (i == ll) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 vfi = vf + nlf - 1
                 vli = vl + nlf - 1
                 idxqi = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 if (icompq == 0) then
                    call la_dlasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm,givptr(1),givcol,ldgcol,givnum,ldu, &
                    poles,difl,difr,z,k(1),c(1),s(1),work(nwork1),iwork(iwk), &
                              info)
                 else
                    j = j - 1
                    call la_dlasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm(nlf,lvl),givptr(j),givcol(nlf,lvl2), &
                     ldgcol,givnum(nlf,lvl2),ldu,poles(nlf,lvl2),difl(nlf,lvl),difr( &
                     nlf,lvl2),z(nlf,lvl),k(j),c(j),s(j),work(nwork1),iwork(iwk &
                               ),info)
                 end if
                 if (info /= 0) then
                    return
                 end if
              end do loop_40
           end do loop_50
           return
     end subroutine la_dlasda
#ifdef LA_WITH_XDP
     !> Using a divide and conquer approach, XLASDA: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
     !> B with diagonal D and offdiagonal E, where M = N + SQRE. The
     !> algorithm computes the singular values in the SVD B = U * S * VT.
     !> The orthogonal matrices U and VT are optionally computed in
     !> compact form.
     !> A related subroutine, XLASD0, computes the singular values and
     !> the singular vectors in explicit form.

     pure subroutine la_xlasda(icompq,smlsiz,n,sqre,d,e,u,ldu,vt,k,difl,difr,z, &
               poles,givptr,givcol,ldgcol,perm,givnum,c,s,work,iwork,info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldgcol,ldu,n,smlsiz,sqre
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),givptr(*),iwork(*),k(*),perm(ldgcol, &
                     *)
           real(xdp),intent(out) :: c(*),difl(ldu,*),difr(ldu,*),givnum(ldu,*),poles(ldu,*), &
                     s(*),u(ldu,*),vt(ldu,*),work(*),z(ldu,*)
           real(xdp),intent(inout) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqi,im1,inode,itemp,iwk,j,lf,ll,lvl,lvl2, &
           m,ncc,nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,nru,nwork1, &
                     nwork2,smlszp,sqrei,vf,vfi,vl,vli
           real(xdp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (smlsiz < 3) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldu < (n + sqre)) then
              info = -8
           else if (ldgcol < n) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('XLASDA',-info)
              return
           end if
           m = n + sqre
           ! if the input matrix is too small, call la_xlasdq to find the svd.
           if (n <= smlsiz) then
              if (icompq == 0) then
                 call la_xlasdq('U',sqre,n,0,0,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              else
                 call la_xlasdq('U',sqre,n,m,n,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              end if
              return
           end if
           ! book-keeping and  set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           ncc = 0
           nru = 0
           smlszp = smlsiz + 1
           vf = 1
           vl = vf + m
           nwork1 = vl + m
           nwork2 = nwork1 + smlszp*smlszp
           call la_xlasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_xlasdq.
           ndb1 = (nd + 1)/2
           loop_30: do i = ndb1,nd
              ! ic : center row of each node
              ! nl : number of rows of left  subproblem
              ! nr : number of rows of right subproblem
              ! nlf: starting row of the left   subproblem
              ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nlf = ic - nl
              nrf = ic + 1
              idxqi = idxq + nlf - 2
              vfi = vf + nlf - 1
              vli = vl + nlf - 1
              sqrei = 1
              if (icompq == 0) then
                 call la_xlaset('A',nlp1,nlp1,zero,one,work(nwork1),smlszp)
                 call la_xlasdq('U',sqrei,nl,nlp1,nru,ncc,d(nlf),e(nlf),work( &
                 nwork1),smlszp,work(nwork2),nl,work(nwork2),nl,work(nwork2),info)

                 itemp = nwork1 + nl*smlszp
                 call la_xcopy(nlp1,work(nwork1),1,work(vfi),1)
                 call la_xcopy(nlp1,work(itemp),1,work(vli),1)
              else
                 call la_xlaset('A',nl,nl,zero,one,u(nlf,1),ldu)
                 call la_xlaset('A',nlp1,nlp1,zero,one,vt(nlf,1),ldu)
                 call la_xlasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,1 &
                           ),ldu,u(nlf,1),ldu,u(nlf,1),ldu,work(nwork1),info)
                 call la_xcopy(nlp1,vt(nlf,1),1,work(vfi),1)
                 call la_xcopy(nlp1,vt(nlf,nlp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nl
                 iwork(idxqi + j) = j
              end do
              if ((i == nd) .and. (sqre == 0)) then
                 sqrei = 0
              else
                 sqrei = 1
              end if
              idxqi = idxqi + nlp1
              vfi = vfi + nlp1
              vli = vli + nlp1
              nrp1 = nr + sqrei
              if (icompq == 0) then
                 call la_xlaset('A',nrp1,nrp1,zero,one,work(nwork1),smlszp)
                 call la_xlasdq('U',sqrei,nr,nrp1,nru,ncc,d(nrf),e(nrf),work( &
                 nwork1),smlszp,work(nwork2),nr,work(nwork2),nr,work(nwork2),info)

                 itemp = nwork1 + (nrp1 - 1)*smlszp
                 call la_xcopy(nrp1,work(nwork1),1,work(vfi),1)
                 call la_xcopy(nrp1,work(itemp),1,work(vli),1)
              else
                 call la_xlaset('A',nr,nr,zero,one,u(nrf,1),ldu)
                 call la_xlaset('A',nrp1,nrp1,zero,one,vt(nrf,1),ldu)
                 call la_xlasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,1 &
                           ),ldu,u(nrf,1),ldu,u(nrf,1),ldu,work(nwork1),info)
                 call la_xcopy(nrp1,vt(nrf,1),1,work(vfi),1)
                 call la_xcopy(nrp1,vt(nrf,nrp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nr
                 iwork(idxqi + j) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           j = 2**nlvl
           loop_50: do lvl = nlvl,1,-1
              lvl2 = lvl*2 - 1
              ! find the first node lf and last node ll on
              ! the current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              loop_40: do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 nrf = ic + 1
                 if (i == ll) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 vfi = vf + nlf - 1
                 vli = vl + nlf - 1
                 idxqi = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 if (icompq == 0) then
                    call la_xlasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm,givptr(1),givcol,ldgcol,givnum,ldu, &
                    poles,difl,difr,z,k(1),c(1),s(1),work(nwork1),iwork(iwk), &
                              info)
                 else
                    j = j - 1
                    call la_xlasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm(nlf,lvl),givptr(j),givcol(nlf,lvl2), &
                     ldgcol,givnum(nlf,lvl2),ldu,poles(nlf,lvl2),difl(nlf,lvl),difr( &
                     nlf,lvl2),z(nlf,lvl),k(j),c(j),s(j),work(nwork1),iwork(iwk &
                               ),info)
                 end if
                 if (info /= 0) then
                    return
                 end if
              end do loop_40
           end do loop_50
           return
     end subroutine la_xlasda
#endif
#ifdef LA_WITH_QP
     !> Using a divide and conquer approach, QLASDA: computes the singular
     !> value decomposition (SVD) of a real upper bidiagonal N-by-M matrix
     !> B with diagonal D and offdiagonal E, where M = N + SQRE. The
     !> algorithm computes the singular values in the SVD B = U * S * VT.
     !> The orthogonal matrices U and VT are optionally computed in
     !> compact form.
     !> A related subroutine, QLASD0, computes the singular values and
     !> the singular vectors in explicit form.

     pure subroutine la_qlasda(icompq,smlsiz,n,sqre,d,e,u,ldu,vt,k,difl,difr,z, &
               poles,givptr,givcol,ldgcol,perm,givnum,c,s,work,iwork,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: icompq,ldgcol,ldu,n,smlsiz,sqre
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: givcol(ldgcol,*),givptr(*),iwork(*),k(*),perm(ldgcol, &
                     *)
           real(qp),intent(out) :: c(*),difl(ldu,*),difr(ldu,*),givnum(ldu,*),poles(ldu,*), &
                     s(*),u(ldu,*),vt(ldu,*),work(*),z(ldu,*)
           real(qp),intent(inout) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,i1,ic,idxq,idxqi,im1,inode,itemp,iwk,j,lf,ll,lvl,lvl2, &
           m,ncc,nd,ndb1,ndiml,ndimr,nl,nlf,nlp1,nlvl,nr,nrf,nrp1,nru,nwork1, &
                     nwork2,smlszp,sqrei,vf,vfi,vl,vli
           real(qp) :: alpha,beta
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if ((icompq < 0) .or. (icompq > 1)) then
              info = -1
           else if (smlsiz < 3) then
              info = -2
           else if (n < 0) then
              info = -3
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -4
           else if (ldu < (n + sqre)) then
              info = -8
           else if (ldgcol < n) then
              info = -17
           end if
           if (info /= 0) then
              call la_xerbla('QLASDA',-info)
              return
           end if
           m = n + sqre
           ! if the input matrix is too small, call la_qlasdq to find the svd.
           if (n <= smlsiz) then
              if (icompq == 0) then
                 call la_qlasdq('U',sqre,n,0,0,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              else
                 call la_qlasdq('U',sqre,n,m,n,0,d,e,vt,ldu,u,ldu,u,ldu,work, &
                           info)
              end if
              return
           end if
           ! book-keeping and  set up the computation tree.
           inode = 1
           ndiml = inode + n
           ndimr = ndiml + n
           idxq = ndimr + n
           iwk = idxq + n
           ncc = 0
           nru = 0
           smlszp = smlsiz + 1
           vf = 1
           vl = vf + m
           nwork1 = vl + m
           nwork2 = nwork1 + smlszp*smlszp
           call la_qlasdt(n,nlvl,nd,iwork(inode),iwork(ndiml),iwork(ndimr),smlsiz &
                     )
           ! for the nodes on bottom level of the tree, solve
           ! their subproblems by la_qlasdq.
           ndb1 = (nd + 1)/2
           loop_30: do i = ndb1,nd
              ! ic : center row of each node
              ! nl : number of rows of left  subproblem
              ! nr : number of rows of right subproblem
              ! nlf: starting row of the left   subproblem
              ! nrf: starting row of the right  subproblem
              i1 = i - 1
              ic = iwork(inode + i1)
              nl = iwork(ndiml + i1)
              nlp1 = nl + 1
              nr = iwork(ndimr + i1)
              nlf = ic - nl
              nrf = ic + 1
              idxqi = idxq + nlf - 2
              vfi = vf + nlf - 1
              vli = vl + nlf - 1
              sqrei = 1
              if (icompq == 0) then
                 call la_qlaset('A',nlp1,nlp1,zero,one,work(nwork1),smlszp)
                 call la_qlasdq('U',sqrei,nl,nlp1,nru,ncc,d(nlf),e(nlf),work( &
                 nwork1),smlszp,work(nwork2),nl,work(nwork2),nl,work(nwork2),info)

                 itemp = nwork1 + nl*smlszp
                 call la_qcopy(nlp1,work(nwork1),1,work(vfi),1)
                 call la_qcopy(nlp1,work(itemp),1,work(vli),1)
              else
                 call la_qlaset('A',nl,nl,zero,one,u(nlf,1),ldu)
                 call la_qlaset('A',nlp1,nlp1,zero,one,vt(nlf,1),ldu)
                 call la_qlasdq('U',sqrei,nl,nlp1,nl,ncc,d(nlf),e(nlf),vt(nlf,1 &
                           ),ldu,u(nlf,1),ldu,u(nlf,1),ldu,work(nwork1),info)
                 call la_qcopy(nlp1,vt(nlf,1),1,work(vfi),1)
                 call la_qcopy(nlp1,vt(nlf,nlp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nl
                 iwork(idxqi + j) = j
              end do
              if ((i == nd) .and. (sqre == 0)) then
                 sqrei = 0
              else
                 sqrei = 1
              end if
              idxqi = idxqi + nlp1
              vfi = vfi + nlp1
              vli = vli + nlp1
              nrp1 = nr + sqrei
              if (icompq == 0) then
                 call la_qlaset('A',nrp1,nrp1,zero,one,work(nwork1),smlszp)
                 call la_qlasdq('U',sqrei,nr,nrp1,nru,ncc,d(nrf),e(nrf),work( &
                 nwork1),smlszp,work(nwork2),nr,work(nwork2),nr,work(nwork2),info)

                 itemp = nwork1 + (nrp1 - 1)*smlszp
                 call la_qcopy(nrp1,work(nwork1),1,work(vfi),1)
                 call la_qcopy(nrp1,work(itemp),1,work(vli),1)
              else
                 call la_qlaset('A',nr,nr,zero,one,u(nrf,1),ldu)
                 call la_qlaset('A',nrp1,nrp1,zero,one,vt(nrf,1),ldu)
                 call la_qlasdq('U',sqrei,nr,nrp1,nr,ncc,d(nrf),e(nrf),vt(nrf,1 &
                           ),ldu,u(nrf,1),ldu,u(nrf,1),ldu,work(nwork1),info)
                 call la_qcopy(nrp1,vt(nrf,1),1,work(vfi),1)
                 call la_qcopy(nrp1,vt(nrf,nrp1),1,work(vli),1)
              end if
              if (info /= 0) then
                 return
              end if
              do j = 1,nr
                 iwork(idxqi + j) = j
              end do
           end do loop_30
           ! now conquer each subproblem bottom-up.
           j = 2**nlvl
           loop_50: do lvl = nlvl,1,-1
              lvl2 = lvl*2 - 1
              ! find the first node lf and last node ll on
              ! the current level lvl.
              if (lvl == 1) then
                 lf = 1
                 ll = 1
              else
                 lf = 2**(lvl - 1)
                 ll = 2*lf - 1
              end if
              loop_40: do i = lf,ll
                 im1 = i - 1
                 ic = iwork(inode + im1)
                 nl = iwork(ndiml + im1)
                 nr = iwork(ndimr + im1)
                 nlf = ic - nl
                 nrf = ic + 1
                 if (i == ll) then
                    sqrei = sqre
                 else
                    sqrei = 1
                 end if
                 vfi = vf + nlf - 1
                 vli = vl + nlf - 1
                 idxqi = idxq + nlf - 1
                 alpha = d(ic)
                 beta = e(ic)
                 if (icompq == 0) then
                    call la_qlasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm,givptr(1),givcol,ldgcol,givnum,ldu, &
                    poles,difl,difr,z,k(1),c(1),s(1),work(nwork1),iwork(iwk), &
                              info)
                 else
                    j = j - 1
                    call la_qlasd6(icompq,nl,nr,sqrei,d(nlf),work(vfi),work(vli), &
                    alpha,beta,iwork(idxqi),perm(nlf,lvl),givptr(j),givcol(nlf,lvl2), &
                     ldgcol,givnum(nlf,lvl2),ldu,poles(nlf,lvl2),difl(nlf,lvl),difr( &
                     nlf,lvl2),z(nlf,lvl),k(j),c(j),s(j),work(nwork1),iwork(iwk &
                               ),info)
                 end if
                 if (info /= 0) then
                    return
                 end if
              end do loop_40
           end do loop_50
           return
     end subroutine la_qlasda
#endif

     !> SLASDQ: computes the singular value decomposition (SVD) of a real
     !> (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
     !> E, accumulating the transformations if desired. Letting B denote
     !> the input bidiagonal matrix, the algorithm computes orthogonal
     !> matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose
     !> of P). The singular values S are overwritten on D.
     !> The input matrix U  is changed to U  * Q  if desired.
     !> The input matrix VT is changed to P**T * VT if desired.
     !> The input matrix C  is changed to Q**T * C  if desired.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3, for a detailed description of the algorithm.

     pure subroutine la_slasdq(uplo,sqre,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc, &
               work,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru,sqre
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: rotate
           integer(ilp) :: i,isub,iuplo,j,np1,sqre1
           real(sp) :: cs,r,smin,sn
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (iuplo == 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncvt < 0) then
              info = -4
           else if (nru < 0) then
              info = -5
           else if (ncc < 0) then
              info = -6
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -10
           else if (ldu < max(1,nru)) then
              info = -12
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('SLASDQ',-info)
              return
           end if
           if (n == 0) return
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           np1 = n + 1
           sqre1 = sqre
           ! if matrix non-square upper bidiagonal, rotate to be lower
           ! bidiagonal.  the rotations are on the right.
           if ((iuplo == 1) .and. (sqre1 == 1)) then
              do i = 1,n - 1
                 call la_slartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              call la_slartg(d(n),e(n),cs,sn,r)
              d(n) = r
              e(n) = zero
              if (rotate) then
                 work(n) = cs
                 work(n + n) = sn
              end if
              iuplo = 2
              sqre1 = 0
              ! update singular vectors if desired.
              if (ncvt > 0) call la_slasr('L','V','F',np1,ncvt,work(1),work(np1),vt, &
                        ldvt)
           end if
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left.
           if (iuplo == 2) then
              do i = 1,n - 1
                 call la_slartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              ! if matrix (n+1)-by-n lower bidiagonal, one additional
              ! rotation is needed.
              if (sqre1 == 1) then
                 call la_slartg(d(n),e(n),cs,sn,r)
                 d(n) = r
                 if (rotate) then
                    work(n) = cs
                    work(n + n) = sn
                 end if
              end if
              ! update singular vectors if desired.
              if (nru > 0) then
                 if (sqre1 == 0) then
                    call la_slasr('R','V','F',nru,n,work(1),work(np1),u,ldu)

                 else
                    call la_slasr('R','V','F',nru,np1,work(1),work(np1),u,ldu)

                 end if
              end if
              if (ncc > 0) then
                 if (sqre1 == 0) then
                    call la_slasr('L','V','F',n,ncc,work(1),work(np1),c,ldc)

                 else
                    call la_slasr('L','V','F',np1,ncc,work(1),work(np1),c,ldc)

                 end if
              end if
           end if
           ! call la_sbdsqr to compute the svd of the reduced real
           ! n-by-n upper bidiagonal matrix.
           call la_sbdsqr('U',n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)

           ! sort the singular values into ascending order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n
              ! scan for smallest d(i).
              isub = i
              smin = d(i)
              do j = i + 1,n
                 if (d(j) < smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= i) then
                 ! swap singular values and vectors.
                 d(isub) = d(i)
                 d(i) = smin
                 if (ncvt > 0) call la_sswap(ncvt,vt(isub,1),ldvt,vt(i,1),ldvt)

                 if (nru > 0) call la_sswap(nru,u(1,isub),1,u(1,i),1)
                 if (ncc > 0) call la_sswap(ncc,c(isub,1),ldc,c(i,1),ldc)
              end if
           end do
           return
     end subroutine la_slasdq
     !> DLASDQ: computes the singular value decomposition (SVD) of a real
     !> (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
     !> E, accumulating the transformations if desired. Letting B denote
     !> the input bidiagonal matrix, the algorithm computes orthogonal
     !> matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose
     !> of P). The singular values S are overwritten on D.
     !> The input matrix U  is changed to U  * Q  if desired.
     !> The input matrix VT is changed to P**T * VT if desired.
     !> The input matrix C  is changed to Q**T * C  if desired.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3, for a detailed description of the algorithm.

     pure subroutine la_dlasdq(uplo,sqre,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc, &
               work,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru,sqre
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: rotate
           integer(ilp) :: i,isub,iuplo,j,np1,sqre1
           real(dp) :: cs,r,smin,sn
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (iuplo == 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncvt < 0) then
              info = -4
           else if (nru < 0) then
              info = -5
           else if (ncc < 0) then
              info = -6
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -10
           else if (ldu < max(1,nru)) then
              info = -12
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('DLASDQ',-info)
              return
           end if
           if (n == 0) return
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           np1 = n + 1
           sqre1 = sqre
           ! if matrix non-square upper bidiagonal, rotate to be lower
           ! bidiagonal.  the rotations are on the right.
           if ((iuplo == 1) .and. (sqre1 == 1)) then
              do i = 1,n - 1
                 call la_dlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              call la_dlartg(d(n),e(n),cs,sn,r)
              d(n) = r
              e(n) = zero
              if (rotate) then
                 work(n) = cs
                 work(n + n) = sn
              end if
              iuplo = 2
              sqre1 = 0
              ! update singular vectors if desired.
              if (ncvt > 0) call la_dlasr('L','V','F',np1,ncvt,work(1),work(np1),vt, &
                        ldvt)
           end if
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left.
           if (iuplo == 2) then
              do i = 1,n - 1
                 call la_dlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              ! if matrix (n+1)-by-n lower bidiagonal, one additional
              ! rotation is needed.
              if (sqre1 == 1) then
                 call la_dlartg(d(n),e(n),cs,sn,r)
                 d(n) = r
                 if (rotate) then
                    work(n) = cs
                    work(n + n) = sn
                 end if
              end if
              ! update singular vectors if desired.
              if (nru > 0) then
                 if (sqre1 == 0) then
                    call la_dlasr('R','V','F',nru,n,work(1),work(np1),u,ldu)

                 else
                    call la_dlasr('R','V','F',nru,np1,work(1),work(np1),u,ldu)

                 end if
              end if
              if (ncc > 0) then
                 if (sqre1 == 0) then
                    call la_dlasr('L','V','F',n,ncc,work(1),work(np1),c,ldc)

                 else
                    call la_dlasr('L','V','F',np1,ncc,work(1),work(np1),c,ldc)

                 end if
              end if
           end if
           ! call la_dbdsqr to compute the svd of the reduced real
           ! n-by-n upper bidiagonal matrix.
           call la_dbdsqr('U',n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)

           ! sort the singular values into ascending order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n
              ! scan for smallest d(i).
              isub = i
              smin = d(i)
              do j = i + 1,n
                 if (d(j) < smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= i) then
                 ! swap singular values and vectors.
                 d(isub) = d(i)
                 d(i) = smin
                 if (ncvt > 0) call la_dswap(ncvt,vt(isub,1),ldvt,vt(i,1),ldvt)

                 if (nru > 0) call la_dswap(nru,u(1,isub),1,u(1,i),1)
                 if (ncc > 0) call la_dswap(ncc,c(isub,1),ldc,c(i,1),ldc)
              end if
           end do
           return
     end subroutine la_dlasdq
#ifdef LA_WITH_XDP
     !> XLASDQ: computes the singular value decomposition (SVD) of a real
     !> (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
     !> E, accumulating the transformations if desired. Letting B denote
     !> the input bidiagonal matrix, the algorithm computes orthogonal
     !> matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose
     !> of P). The singular values S are overwritten on D.
     !> The input matrix U  is changed to U  * Q  if desired.
     !> The input matrix VT is changed to P**T * VT if desired.
     !> The input matrix C  is changed to Q**T * C  if desired.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3, for a detailed description of the algorithm.

     pure subroutine la_xlasdq(uplo,sqre,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc, &
               work,info)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru,sqre
           ! Array Arguments
           real(xdp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: rotate
           integer(ilp) :: i,isub,iuplo,j,np1,sqre1
           real(xdp) :: cs,r,smin,sn
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (iuplo == 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncvt < 0) then
              info = -4
           else if (nru < 0) then
              info = -5
           else if (ncc < 0) then
              info = -6
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -10
           else if (ldu < max(1,nru)) then
              info = -12
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('XLASDQ',-info)
              return
           end if
           if (n == 0) return
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           np1 = n + 1
           sqre1 = sqre
           ! if matrix non-square upper bidiagonal, rotate to be lower
           ! bidiagonal.  the rotations are on the right.
           if ((iuplo == 1) .and. (sqre1 == 1)) then
              do i = 1,n - 1
                 call la_xlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              call la_xlartg(d(n),e(n),cs,sn,r)
              d(n) = r
              e(n) = zero
              if (rotate) then
                 work(n) = cs
                 work(n + n) = sn
              end if
              iuplo = 2
              sqre1 = 0
              ! update singular vectors if desired.
              if (ncvt > 0) call la_xlasr('L','V','F',np1,ncvt,work(1),work(np1),vt, &
                        ldvt)
           end if
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left.
           if (iuplo == 2) then
              do i = 1,n - 1
                 call la_xlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              ! if matrix (n+1)-by-n lower bidiagonal, one additional
              ! rotation is needed.
              if (sqre1 == 1) then
                 call la_xlartg(d(n),e(n),cs,sn,r)
                 d(n) = r
                 if (rotate) then
                    work(n) = cs
                    work(n + n) = sn
                 end if
              end if
              ! update singular vectors if desired.
              if (nru > 0) then
                 if (sqre1 == 0) then
                    call la_xlasr('R','V','F',nru,n,work(1),work(np1),u,ldu)

                 else
                    call la_xlasr('R','V','F',nru,np1,work(1),work(np1),u,ldu)

                 end if
              end if
              if (ncc > 0) then
                 if (sqre1 == 0) then
                    call la_xlasr('L','V','F',n,ncc,work(1),work(np1),c,ldc)

                 else
                    call la_xlasr('L','V','F',np1,ncc,work(1),work(np1),c,ldc)

                 end if
              end if
           end if
           ! call la_xbdsqr to compute the svd of the reduced real
           ! n-by-n upper bidiagonal matrix.
           call la_xbdsqr('U',n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)

           ! sort the singular values into ascending order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n
              ! scan for smallest d(i).
              isub = i
              smin = d(i)
              do j = i + 1,n
                 if (d(j) < smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= i) then
                 ! swap singular values and vectors.
                 d(isub) = d(i)
                 d(i) = smin
                 if (ncvt > 0) call la_xswap(ncvt,vt(isub,1),ldvt,vt(i,1),ldvt)

                 if (nru > 0) call la_xswap(nru,u(1,isub),1,u(1,i),1)
                 if (ncc > 0) call la_xswap(ncc,c(isub,1),ldc,c(i,1),ldc)
              end if
           end do
           return
     end subroutine la_xlasdq
#endif
#ifdef LA_WITH_QP
     !> QLASDQ: computes the singular value decomposition (SVD) of a real
     !> (upper or lower) bidiagonal matrix with diagonal D and offdiagonal
     !> E, accumulating the transformations if desired. Letting B denote
     !> the input bidiagonal matrix, the algorithm computes orthogonal
     !> matrices Q and P such that B = Q * S * P**T (P**T denotes the transpose
     !> of P). The singular values S are overwritten on D.
     !> The input matrix U  is changed to U  * Q  if desired.
     !> The input matrix VT is changed to P**T * VT if desired.
     !> The input matrix C  is changed to Q**T * C  if desired.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3, for a detailed description of the algorithm.

     pure subroutine la_qlasdq(uplo,sqre,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc, &
               work,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru,sqre
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: rotate
           integer(ilp) :: i,isub,iuplo,j,np1,sqre1
           real(qp) :: cs,r,smin,sn
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           iuplo = 0
           if (la_lsame(uplo,'U')) iuplo = 1
           if (la_lsame(uplo,'L')) iuplo = 2
           if (iuplo == 0) then
              info = -1
           else if ((sqre < 0) .or. (sqre > 1)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncvt < 0) then
              info = -4
           else if (nru < 0) then
              info = -5
           else if (ncc < 0) then
              info = -6
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -10
           else if (ldu < max(1,nru)) then
              info = -12
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -14
           end if
           if (info /= 0) then
              call la_xerbla('QLASDQ',-info)
              return
           end if
           if (n == 0) return
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           np1 = n + 1
           sqre1 = sqre
           ! if matrix non-square upper bidiagonal, rotate to be lower
           ! bidiagonal.  the rotations are on the right.
           if ((iuplo == 1) .and. (sqre1 == 1)) then
              do i = 1,n - 1
                 call la_qlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              call la_qlartg(d(n),e(n),cs,sn,r)
              d(n) = r
              e(n) = zero
              if (rotate) then
                 work(n) = cs
                 work(n + n) = sn
              end if
              iuplo = 2
              sqre1 = 0
              ! update singular vectors if desired.
              if (ncvt > 0) call la_qlasr('L','V','F',np1,ncvt,work(1),work(np1),vt, &
                        ldvt)
           end if
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left.
           if (iuplo == 2) then
              do i = 1,n - 1
                 call la_qlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 if (rotate) then
                    work(i) = cs
                    work(n + i) = sn
                 end if
              end do
              ! if matrix (n+1)-by-n lower bidiagonal, one additional
              ! rotation is needed.
              if (sqre1 == 1) then
                 call la_qlartg(d(n),e(n),cs,sn,r)
                 d(n) = r
                 if (rotate) then
                    work(n) = cs
                    work(n + n) = sn
                 end if
              end if
              ! update singular vectors if desired.
              if (nru > 0) then
                 if (sqre1 == 0) then
                    call la_qlasr('R','V','F',nru,n,work(1),work(np1),u,ldu)

                 else
                    call la_qlasr('R','V','F',nru,np1,work(1),work(np1),u,ldu)

                 end if
              end if
              if (ncc > 0) then
                 if (sqre1 == 0) then
                    call la_qlasr('L','V','F',n,ncc,work(1),work(np1),c,ldc)

                 else
                    call la_qlasr('L','V','F',np1,ncc,work(1),work(np1),c,ldc)

                 end if
              end if
           end if
           ! call la_qbdsqr to compute the svd of the reduced real
           ! n-by-n upper bidiagonal matrix.
           call la_qbdsqr('U',n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work,info)

           ! sort the singular values into ascending order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n
              ! scan for smallest d(i).
              isub = i
              smin = d(i)
              do j = i + 1,n
                 if (d(j) < smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= i) then
                 ! swap singular values and vectors.
                 d(isub) = d(i)
                 d(i) = smin
                 if (ncvt > 0) call la_qswap(ncvt,vt(isub,1),ldvt,vt(i,1),ldvt)

                 if (nru > 0) call la_qswap(nru,u(1,isub),1,u(1,i),1)
                 if (ncc > 0) call la_qswap(ncc,c(isub,1),ldc,c(i,1),ldc)
              end if
           end do
           return
     end subroutine la_qlasdq
#endif

end module la_lapack_eigv_svd_bidiag_dc
