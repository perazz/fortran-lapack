!> Bidiagonal singular values: implicit QR sweep and the dqds algorithm
module la_lapack_svd_bidiag_qr
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_l2
     use la_lapack_blas_like_scalar
     use la_lapack_givens_jacobi_rot
     use la_lapack_svd_comp2
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slasq4
     public :: la_slasq5
     public :: la_slasq6
     public :: la_slasq3
     public :: la_sbdsqr
     public :: la_slasq1
     public :: la_slasq2
     public :: la_dlasq4
     public :: la_dlasq5
     public :: la_dlasq6
     public :: la_dlasq3
     public :: la_dbdsqr
     public :: la_dlasq1
     public :: la_dlasq2
#ifdef LA_WITH_XDP
     public :: la_xlasq4
     public :: la_xlasq5
     public :: la_xlasq6
     public :: la_xlasq3
     public :: la_xbdsqr
     public :: la_xlasq1
     public :: la_xlasq2
#endif
#ifdef LA_WITH_QP
     public :: la_qlasq4
     public :: la_qlasq5
     public :: la_qlasq6
     public :: la_qlasq3
     public :: la_qbdsqr
     public :: la_qlasq1
     public :: la_qlasq2
#endif
     public :: la_cbdsqr
     public :: la_zbdsqr
#ifdef LA_WITH_XDP
     public :: la_ybdsqr
#endif
#ifdef LA_WITH_QP
     public :: la_wbdsqr
#endif

     contains

     !> SLASQ4: computes an approximation TAU to the smallest eigenvalue
     !> using values of d from the previous transform.

     pure subroutine la_slasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau, &
               ttype,g)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,n0in,pp
           integer(ilp),intent(out) :: ttype
           real(sp),intent(in) :: dmin,dmin1,dmin2,dn,dn1,dn2
           real(sp),intent(inout) :: g
           real(sp),intent(out) :: tau
           ! Array Arguments
           real(sp),intent(in) :: z(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: cnst1 = 0.5630_sp
           real(sp),parameter :: cnst2 = 1.010_sp
           real(sp),parameter :: cnst3 = 1.050_sp
           real(sp),parameter :: qurtr = 0.250_sp
           real(sp),parameter :: third = 0.3330_sp
           real(sp),parameter :: hundrd = 100.0_sp

           ! Local Scalars
           integer(ilp) :: i4,nn,np
           real(sp) :: a2,b1,b2,gam,gap1,gap2,s
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! a negative dmin forces the shift to take that absolute value
           ! ttype records the type of shift.
           if (dmin <= zero) then
              tau = -dmin
              ttype = -1
              return
           end if
           nn = 4*n0 + pp
           if (n0in == n0) then
              ! no eigenvalues deflated.
              if (dmin == dn .or. dmin == dn1) then
                 b1 = sqrt(z(nn - 3))*sqrt(z(nn - 5))
                 b2 = sqrt(z(nn - 7))*sqrt(z(nn - 9))
                 a2 = z(nn - 7) + z(nn - 5)
                 ! cases 2 and 3.
                 if (dmin == dn .and. dmin1 == dn1) then
                    gap2 = dmin2 - a2 - dmin2*qurtr
                    if (gap2 > zero .and. gap2 > b2) then
                       gap1 = a2 - dn - (b2/gap2)*b2
                    else
                       gap1 = a2 - dn - (b1 + b2)
                    end if
                    if (gap1 > zero .and. gap1 > b1) then
                       s = max(dn - (b1/gap1)*b1,half*dmin)
                       ttype = -2
                    else
                       s = zero
                       if (dn > b1) s = dn - b1
                       if (a2 > (b1 + b2)) s = min(s,a2 - (b1 + b2))
                       s = max(s,third*dmin)
                       ttype = -3
                    end if
                 else
                    ! case 4.
                    ttype = -4
                    s = qurtr*dmin
                    if (dmin == dn) then
                       gam = dn
                       a2 = zero
                       if (z(nn - 5) > z(nn - 7)) return
                       b2 = z(nn - 5)/z(nn - 7)
                       np = nn - 9
                    else
                       np = nn - 2*pp
                       gam = dn1
                       if (z(np - 4) > z(np - 2)) return
                       a2 = z(np - 4)/z(np - 2)
                       if (z(nn - 9) > z(nn - 11)) return
                       b2 = z(nn - 9)/z(nn - 11)
                       np = nn - 13
                    end if
                    ! approximate contribution to norm squared from i < nn-1.
                    a2 = a2 + b2
                    do i4 = np,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 20
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 20
                    end do
                    20 continue
                    a2 = cnst3*a2
                    ! rayleigh quotient residual bound.
                    if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
                 end if
              else if (dmin == dn2) then
                 ! case 5.
                 ttype = -5
                 s = qurtr*dmin
                 ! compute contribution to norm squared from i > nn-2.
                 np = nn - 2*pp
                 b1 = z(np - 2)
                 b2 = z(np - 6)
                 gam = dn2
                 if (z(np - 8) > b2 .or. z(np - 4) > b1) return
                 a2 = (z(np - 8)/b2)*(one + z(np - 4)/b1)
                 ! approximate contribution to norm squared from i < nn-2.
                 if (n0 - i0 > 2) then
                    b2 = z(nn - 13)/z(nn - 15)
                    a2 = a2 + b2
                    do i4 = nn - 17,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 40
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 40
                    end do
                    40 continue
                    a2 = cnst3*a2
                 end if
                 if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
              else
                 ! case 6, no information to guide us.
                 if (ttype == -6) then
                    g = g + third*(one - g)
                 else if (ttype == -18) then
                    g = qurtr*third
                 else
                    g = qurtr
                 end if
                 s = g*dmin
                 ttype = -6
              end if
           else if (n0in == (n0 + 1)) then
              ! one eigenvalue just deflated. use dmin1, dn1 for dmin and dn.
              if (dmin1 == dn1 .and. dmin2 == dn2) then
                 ! cases 7 and 8.
                 ttype = -7
                 s = third*dmin1
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 60
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    a2 = b1
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*max(b1,a2) < b2) go to 60
                 end do
                 60 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin1/(one + b2**2)
                 gap2 = half*dmin2 - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                    ttype = -8
                 end if
              else
                 ! case 9.
                 s = qurtr*dmin1
                 if (dmin1 == dn1) s = half*dmin1
                 ttype = -9
              end if
           else if (n0in == (n0 + 2)) then
              ! two eigenvalues deflated. use dmin2, dn2 for dmin and dn.
              ! cases 10 and 11.
              if (dmin2 == dn2 .and. two*z(nn - 5) < z(nn - 7)) then
                 ttype = -10
                 s = third*dmin2
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 80
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*b1 < b2) go to 80
                 end do
                 80 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin2/(one + b2**2)
                 gap2 = z(nn - 7) + z(nn - 9) - sqrt(z(nn - 11))*sqrt(z(nn - 9)) - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                 end if
              else
                 s = qurtr*dmin2
                 ttype = -11
              end if
           else if (n0in > (n0 + 2)) then
              ! case 12, more than two eigenvalues deflated. no information.
              s = zero
              ttype = -12
           end if
           tau = s
           return
     end subroutine la_slasq4
     !> DLASQ4: computes an approximation TAU to the smallest eigenvalue
     !> using values of d from the previous transform.

     pure subroutine la_dlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau, &
               ttype,g)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,n0in,pp
           integer(ilp),intent(out) :: ttype
           real(dp),intent(in) :: dmin,dmin1,dmin2,dn,dn1,dn2
           real(dp),intent(inout) :: g
           real(dp),intent(out) :: tau
           ! Array Arguments
           real(dp),intent(in) :: z(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: cnst1 = 0.5630_dp
           real(dp),parameter :: cnst2 = 1.010_dp
           real(dp),parameter :: cnst3 = 1.050_dp
           real(dp),parameter :: qurtr = 0.250_dp
           real(dp),parameter :: third = 0.3330_dp
           real(dp),parameter :: hundrd = 100.0_dp

           ! Local Scalars
           integer(ilp) :: i4,nn,np
           real(dp) :: a2,b1,b2,gam,gap1,gap2,s
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! a negative dmin forces the shift to take that absolute value
           ! ttype records the type of shift.
           if (dmin <= zero) then
              tau = -dmin
              ttype = -1
              return
           end if
           nn = 4*n0 + pp
           if (n0in == n0) then
              ! no eigenvalues deflated.
              if (dmin == dn .or. dmin == dn1) then
                 b1 = sqrt(z(nn - 3))*sqrt(z(nn - 5))
                 b2 = sqrt(z(nn - 7))*sqrt(z(nn - 9))
                 a2 = z(nn - 7) + z(nn - 5)
                 ! cases 2 and 3.
                 if (dmin == dn .and. dmin1 == dn1) then
                    gap2 = dmin2 - a2 - dmin2*qurtr
                    if (gap2 > zero .and. gap2 > b2) then
                       gap1 = a2 - dn - (b2/gap2)*b2
                    else
                       gap1 = a2 - dn - (b1 + b2)
                    end if
                    if (gap1 > zero .and. gap1 > b1) then
                       s = max(dn - (b1/gap1)*b1,half*dmin)
                       ttype = -2
                    else
                       s = zero
                       if (dn > b1) s = dn - b1
                       if (a2 > (b1 + b2)) s = min(s,a2 - (b1 + b2))
                       s = max(s,third*dmin)
                       ttype = -3
                    end if
                 else
                    ! case 4.
                    ttype = -4
                    s = qurtr*dmin
                    if (dmin == dn) then
                       gam = dn
                       a2 = zero
                       if (z(nn - 5) > z(nn - 7)) return
                       b2 = z(nn - 5)/z(nn - 7)
                       np = nn - 9
                    else
                       np = nn - 2*pp
                       gam = dn1
                       if (z(np - 4) > z(np - 2)) return
                       a2 = z(np - 4)/z(np - 2)
                       if (z(nn - 9) > z(nn - 11)) return
                       b2 = z(nn - 9)/z(nn - 11)
                       np = nn - 13
                    end if
                    ! approximate contribution to norm squared from i < nn-1.
                    a2 = a2 + b2
                    do i4 = np,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 20
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 20
                    end do
                    20 continue
                    a2 = cnst3*a2
                    ! rayleigh quotient residual bound.
                    if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
                 end if
              else if (dmin == dn2) then
                 ! case 5.
                 ttype = -5
                 s = qurtr*dmin
                 ! compute contribution to norm squared from i > nn-2.
                 np = nn - 2*pp
                 b1 = z(np - 2)
                 b2 = z(np - 6)
                 gam = dn2
                 if (z(np - 8) > b2 .or. z(np - 4) > b1) return
                 a2 = (z(np - 8)/b2)*(one + z(np - 4)/b1)
                 ! approximate contribution to norm squared from i < nn-2.
                 if (n0 - i0 > 2) then
                    b2 = z(nn - 13)/z(nn - 15)
                    a2 = a2 + b2
                    do i4 = nn - 17,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 40
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 40
                    end do
                    40 continue
                    a2 = cnst3*a2
                 end if
                 if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
              else
                 ! case 6, no information to guide us.
                 if (ttype == -6) then
                    g = g + third*(one - g)
                 else if (ttype == -18) then
                    g = qurtr*third
                 else
                    g = qurtr
                 end if
                 s = g*dmin
                 ttype = -6
              end if
           else if (n0in == (n0 + 1)) then
              ! one eigenvalue just deflated. use dmin1, dn1 for dmin and dn.
              if (dmin1 == dn1 .and. dmin2 == dn2) then
                 ! cases 7 and 8.
                 ttype = -7
                 s = third*dmin1
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 60
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    a2 = b1
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*max(b1,a2) < b2) go to 60
                 end do
                 60 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin1/(one + b2**2)
                 gap2 = half*dmin2 - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                    ttype = -8
                 end if
              else
                 ! case 9.
                 s = qurtr*dmin1
                 if (dmin1 == dn1) s = half*dmin1
                 ttype = -9
              end if
           else if (n0in == (n0 + 2)) then
              ! two eigenvalues deflated. use dmin2, dn2 for dmin and dn.
              ! cases 10 and 11.
              if (dmin2 == dn2 .and. two*z(nn - 5) < z(nn - 7)) then
                 ttype = -10
                 s = third*dmin2
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 80
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*b1 < b2) go to 80
                 end do
                 80 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin2/(one + b2**2)
                 gap2 = z(nn - 7) + z(nn - 9) - sqrt(z(nn - 11))*sqrt(z(nn - 9)) - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                 end if
              else
                 s = qurtr*dmin2
                 ttype = -11
              end if
           else if (n0in > (n0 + 2)) then
              ! case 12, more than two eigenvalues deflated. no information.
              s = zero
              ttype = -12
           end if
           tau = s
           return
     end subroutine la_dlasq4
#ifdef LA_WITH_XDP
     !> XLASQ4: computes an approximation TAU to the smallest eigenvalue
     !> using values of d from the previous transform.

     pure subroutine la_xlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau, &
               ttype,g)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,n0in,pp
           integer(ilp),intent(out) :: ttype
           real(xdp),intent(in) :: dmin,dmin1,dmin2,dn,dn1,dn2
           real(xdp),intent(inout) :: g
           real(xdp),intent(out) :: tau
           ! Array Arguments
           real(xdp),intent(in) :: z(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: cnst1 = 0.5630_xdp
           real(xdp),parameter :: cnst2 = 1.010_xdp
           real(xdp),parameter :: cnst3 = 1.050_xdp
           real(xdp),parameter :: qurtr = 0.250_xdp
           real(xdp),parameter :: third = 0.3330_xdp
           real(xdp),parameter :: hundrd = 100.0_xdp

           ! Local Scalars
           integer(ilp) :: i4,nn,np
           real(xdp) :: a2,b1,b2,gam,gap1,gap2,s
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! a negative dmin forces the shift to take that absolute value
           ! ttype records the type of shift.
           if (dmin <= zero) then
              tau = -dmin
              ttype = -1
              return
           end if
           nn = 4*n0 + pp
           if (n0in == n0) then
              ! no eigenvalues deflated.
              if (dmin == dn .or. dmin == dn1) then
                 b1 = sqrt(z(nn - 3))*sqrt(z(nn - 5))
                 b2 = sqrt(z(nn - 7))*sqrt(z(nn - 9))
                 a2 = z(nn - 7) + z(nn - 5)
                 ! cases 2 and 3.
                 if (dmin == dn .and. dmin1 == dn1) then
                    gap2 = dmin2 - a2 - dmin2*qurtr
                    if (gap2 > zero .and. gap2 > b2) then
                       gap1 = a2 - dn - (b2/gap2)*b2
                    else
                       gap1 = a2 - dn - (b1 + b2)
                    end if
                    if (gap1 > zero .and. gap1 > b1) then
                       s = max(dn - (b1/gap1)*b1,half*dmin)
                       ttype = -2
                    else
                       s = zero
                       if (dn > b1) s = dn - b1
                       if (a2 > (b1 + b2)) s = min(s,a2 - (b1 + b2))
                       s = max(s,third*dmin)
                       ttype = -3
                    end if
                 else
                    ! case 4.
                    ttype = -4
                    s = qurtr*dmin
                    if (dmin == dn) then
                       gam = dn
                       a2 = zero
                       if (z(nn - 5) > z(nn - 7)) return
                       b2 = z(nn - 5)/z(nn - 7)
                       np = nn - 9
                    else
                       np = nn - 2*pp
                       gam = dn1
                       if (z(np - 4) > z(np - 2)) return
                       a2 = z(np - 4)/z(np - 2)
                       if (z(nn - 9) > z(nn - 11)) return
                       b2 = z(nn - 9)/z(nn - 11)
                       np = nn - 13
                    end if
                    ! approximate contribution to norm squared from i < nn-1.
                    a2 = a2 + b2
                    do i4 = np,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 20
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 20
                    end do
                    20 continue
                    a2 = cnst3*a2
                    ! rayleigh quotient residual bound.
                    if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
                 end if
              else if (dmin == dn2) then
                 ! case 5.
                 ttype = -5
                 s = qurtr*dmin
                 ! compute contribution to norm squared from i > nn-2.
                 np = nn - 2*pp
                 b1 = z(np - 2)
                 b2 = z(np - 6)
                 gam = dn2
                 if (z(np - 8) > b2 .or. z(np - 4) > b1) return
                 a2 = (z(np - 8)/b2)*(one + z(np - 4)/b1)
                 ! approximate contribution to norm squared from i < nn-2.
                 if (n0 - i0 > 2) then
                    b2 = z(nn - 13)/z(nn - 15)
                    a2 = a2 + b2
                    do i4 = nn - 17,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 40
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 40
                    end do
                    40 continue
                    a2 = cnst3*a2
                 end if
                 if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
              else
                 ! case 6, no information to guide us.
                 if (ttype == -6) then
                    g = g + third*(one - g)
                 else if (ttype == -18) then
                    g = qurtr*third
                 else
                    g = qurtr
                 end if
                 s = g*dmin
                 ttype = -6
              end if
           else if (n0in == (n0 + 1)) then
              ! one eigenvalue just deflated. use dmin1, dn1 for dmin and dn.
              if (dmin1 == dn1 .and. dmin2 == dn2) then
                 ! cases 7 and 8.
                 ttype = -7
                 s = third*dmin1
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 60
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    a2 = b1
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*max(b1,a2) < b2) go to 60
                 end do
                 60 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin1/(one + b2**2)
                 gap2 = half*dmin2 - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                    ttype = -8
                 end if
              else
                 ! case 9.
                 s = qurtr*dmin1
                 if (dmin1 == dn1) s = half*dmin1
                 ttype = -9
              end if
           else if (n0in == (n0 + 2)) then
              ! two eigenvalues deflated. use dmin2, dn2 for dmin and dn.
              ! cases 10 and 11.
              if (dmin2 == dn2 .and. two*z(nn - 5) < z(nn - 7)) then
                 ttype = -10
                 s = third*dmin2
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 80
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*b1 < b2) go to 80
                 end do
                 80 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin2/(one + b2**2)
                 gap2 = z(nn - 7) + z(nn - 9) - sqrt(z(nn - 11))*sqrt(z(nn - 9)) - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                 end if
              else
                 s = qurtr*dmin2
                 ttype = -11
              end if
           else if (n0in > (n0 + 2)) then
              ! case 12, more than two eigenvalues deflated. no information.
              s = zero
              ttype = -12
           end if
           tau = s
           return
     end subroutine la_xlasq4
#endif
#ifdef LA_WITH_QP
     !> QLASQ4: computes an approximation TAU to the smallest eigenvalue
     !> using values of d from the previous transform.

     pure subroutine la_qlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau, &
               ttype,g)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,n0in,pp
           integer(ilp),intent(out) :: ttype
           real(qp),intent(in) :: dmin,dmin1,dmin2,dn,dn1,dn2
           real(qp),intent(inout) :: g
           real(qp),intent(out) :: tau
           ! Array Arguments
           real(qp),intent(in) :: z(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: cnst1 = 0.5630_qp
           real(qp),parameter :: cnst2 = 1.010_qp
           real(qp),parameter :: cnst3 = 1.050_qp
           real(qp),parameter :: qurtr = 0.250_qp
           real(qp),parameter :: third = 0.3330_qp
           real(qp),parameter :: hundrd = 100.0_qp

           ! Local Scalars
           integer(ilp) :: i4,nn,np
           real(qp) :: a2,b1,b2,gam,gap1,gap2,s
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! a negative dmin forces the shift to take that absolute value
           ! ttype records the type of shift.
           if (dmin <= zero) then
              tau = -dmin
              ttype = -1
              return
           end if
           nn = 4*n0 + pp
           if (n0in == n0) then
              ! no eigenvalues deflated.
              if (dmin == dn .or. dmin == dn1) then
                 b1 = sqrt(z(nn - 3))*sqrt(z(nn - 5))
                 b2 = sqrt(z(nn - 7))*sqrt(z(nn - 9))
                 a2 = z(nn - 7) + z(nn - 5)
                 ! cases 2 and 3.
                 if (dmin == dn .and. dmin1 == dn1) then
                    gap2 = dmin2 - a2 - dmin2*qurtr
                    if (gap2 > zero .and. gap2 > b2) then
                       gap1 = a2 - dn - (b2/gap2)*b2
                    else
                       gap1 = a2 - dn - (b1 + b2)
                    end if
                    if (gap1 > zero .and. gap1 > b1) then
                       s = max(dn - (b1/gap1)*b1,half*dmin)
                       ttype = -2
                    else
                       s = zero
                       if (dn > b1) s = dn - b1
                       if (a2 > (b1 + b2)) s = min(s,a2 - (b1 + b2))
                       s = max(s,third*dmin)
                       ttype = -3
                    end if
                 else
                    ! case 4.
                    ttype = -4
                    s = qurtr*dmin
                    if (dmin == dn) then
                       gam = dn
                       a2 = zero
                       if (z(nn - 5) > z(nn - 7)) return
                       b2 = z(nn - 5)/z(nn - 7)
                       np = nn - 9
                    else
                       np = nn - 2*pp
                       gam = dn1
                       if (z(np - 4) > z(np - 2)) return
                       a2 = z(np - 4)/z(np - 2)
                       if (z(nn - 9) > z(nn - 11)) return
                       b2 = z(nn - 9)/z(nn - 11)
                       np = nn - 13
                    end if
                    ! approximate contribution to norm squared from i < nn-1.
                    a2 = a2 + b2
                    do i4 = np,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 20
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 20
                    end do
                    20 continue
                    a2 = cnst3*a2
                    ! rayleigh quotient residual bound.
                    if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
                 end if
              else if (dmin == dn2) then
                 ! case 5.
                 ttype = -5
                 s = qurtr*dmin
                 ! compute contribution to norm squared from i > nn-2.
                 np = nn - 2*pp
                 b1 = z(np - 2)
                 b2 = z(np - 6)
                 gam = dn2
                 if (z(np - 8) > b2 .or. z(np - 4) > b1) return
                 a2 = (z(np - 8)/b2)*(one + z(np - 4)/b1)
                 ! approximate contribution to norm squared from i < nn-2.
                 if (n0 - i0 > 2) then
                    b2 = z(nn - 13)/z(nn - 15)
                    a2 = a2 + b2
                    do i4 = nn - 17,4*i0 - 1 + pp,-4
                       if (b2 == zero) go to 40
                       b1 = b2
                       if (z(i4) > z(i4 - 2)) return
                       b2 = b2*(z(i4)/z(i4 - 2))
                       a2 = a2 + b2
                       if (hundrd*max(b2,b1) < a2 .or. cnst1 < a2) go to 40
                    end do
                    40 continue
                    a2 = cnst3*a2
                 end if
                 if (a2 < cnst1) s = gam*(one - sqrt(a2))/(one + a2)
              else
                 ! case 6, no information to guide us.
                 if (ttype == -6) then
                    g = g + third*(one - g)
                 else if (ttype == -18) then
                    g = qurtr*third
                 else
                    g = qurtr
                 end if
                 s = g*dmin
                 ttype = -6
              end if
           else if (n0in == (n0 + 1)) then
              ! one eigenvalue just deflated. use dmin1, dn1 for dmin and dn.
              if (dmin1 == dn1 .and. dmin2 == dn2) then
                 ! cases 7 and 8.
                 ttype = -7
                 s = third*dmin1
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 60
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    a2 = b1
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*max(b1,a2) < b2) go to 60
                 end do
                 60 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin1/(one + b2**2)
                 gap2 = half*dmin2 - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                    ttype = -8
                 end if
              else
                 ! case 9.
                 s = qurtr*dmin1
                 if (dmin1 == dn1) s = half*dmin1
                 ttype = -9
              end if
           else if (n0in == (n0 + 2)) then
              ! two eigenvalues deflated. use dmin2, dn2 for dmin and dn.
              ! cases 10 and 11.
              if (dmin2 == dn2 .and. two*z(nn - 5) < z(nn - 7)) then
                 ttype = -10
                 s = third*dmin2
                 if (z(nn - 5) > z(nn - 7)) return
                 b1 = z(nn - 5)/z(nn - 7)
                 b2 = b1
                 if (b2 == zero) go to 80
                 do i4 = 4*n0 - 9 + pp,4*i0 - 1 + pp,-4
                    if (z(i4) > z(i4 - 2)) return
                    b1 = b1*(z(i4)/z(i4 - 2))
                    b2 = b2 + b1
                    if (hundrd*b1 < b2) go to 80
                 end do
                 80 continue
                 b2 = sqrt(cnst3*b2)
                 a2 = dmin2/(one + b2**2)
                 gap2 = z(nn - 7) + z(nn - 9) - sqrt(z(nn - 11))*sqrt(z(nn - 9)) - a2
                 if (gap2 > zero .and. gap2 > b2*a2) then
                    s = max(s,a2*(one - cnst2*a2*(b2/gap2)*b2))
                 else
                    s = max(s,a2*(one - cnst2*b2))
                 end if
              else
                 s = qurtr*dmin2
                 ttype = -11
              end if
           else if (n0in > (n0 + 2)) then
              ! case 12, more than two eigenvalues deflated. no information.
              s = zero
              ttype = -12
           end if
           tau = s
           return
     end subroutine la_qlasq4
#endif

     !> SLASQ5: computes one dqds transform in ping-pong form, one
     !> version for IEEE machines another for non IEEE machines.

     pure subroutine la_slasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dnm1,dnm2, &
               ieee,eps)
        use la_constants_sp,only:zero,half
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0,n0,pp
           real(sp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           real(sp),intent(inout) :: tau
           real(sp),intent(in) :: sigma,eps
           ! Array Arguments
           real(sp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(sp) :: d,emin,temp,dthresh
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           dthresh = eps*(sigma + tau)
           if (tau < dthresh*half) tau = zero
           if (tau /= zero) then
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4) - tau
           dmin = d
           dmin1 = -z(j4)
           if (ieee) then
              ! code for ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    temp = z(j4 + 1)/z(j4 - 2)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4) = z(j4 - 1)*temp
                    emin = min(z(j4),emin)
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    temp = z(j4 + 2)/z(j4 - 3)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4 - 1) = z(j4)*temp
                    emin = min(z(j4 - 1),emin)
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              dmin = min(dmin,dn)
           else
              ! code for non ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    if (d < zero) then
                       return
                    else
                       z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                       d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4))
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    if (d < zero) then
                       return
                    else
                       z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                       d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4 - 1))
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              if (dnm2 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              if (dnm1 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dn)
           end if
           else
           ! this is the version that sets d's to zero if they are small enough
              j4 = 4*i0 + pp - 3
              emin = z(j4 + 4)
              d = z(j4) - tau
              dmin = d
              dmin1 = -z(j4)
              if (ieee) then
           ! code for ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       temp = z(j4 + 1)/z(j4 - 2)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4) = z(j4 - 1)*temp
                       emin = min(z(j4),emin)
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       temp = z(j4 + 2)/z(j4 - 3)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4 - 1) = z(j4)*temp
                       emin = min(z(j4 - 1),emin)
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 dmin = min(dmin,dn)
              else
           ! code for non ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       if (d < zero) then
                          return
                       else
                          z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                          d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4))
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       if (d < zero) then
                          return
                       else
                          z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                          d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4 - 1))
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 if (dnm2 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 if (dnm1 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dn)
              end if
           end if
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_slasq5
     !> DLASQ5: computes one dqds transform in ping-pong form, one
     !> version for IEEE machines another for non IEEE machines.

     pure subroutine la_dlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dnm1,dnm2, &
               ieee,eps)
        use la_constants_dp,only:zero,half
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0,n0,pp
           real(dp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           real(dp),intent(inout) :: tau
           real(dp),intent(in) :: sigma,eps
           ! Array Arguments
           real(dp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(dp) :: d,emin,temp,dthresh
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           dthresh = eps*(sigma + tau)
           if (tau < dthresh*half) tau = zero
           if (tau /= zero) then
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4) - tau
           dmin = d
           dmin1 = -z(j4)
           if (ieee) then
              ! code for ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    temp = z(j4 + 1)/z(j4 - 2)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4) = z(j4 - 1)*temp
                    emin = min(z(j4),emin)
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    temp = z(j4 + 2)/z(j4 - 3)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4 - 1) = z(j4)*temp
                    emin = min(z(j4 - 1),emin)
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              dmin = min(dmin,dn)
           else
              ! code for non ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    if (d < zero) then
                       return
                    else
                       z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                       d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4))
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    if (d < zero) then
                       return
                    else
                       z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                       d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4 - 1))
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              if (dnm2 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              if (dnm1 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dn)
           end if
           else
           ! this is the version that sets d's to zero if they are small enough
              j4 = 4*i0 + pp - 3
              emin = z(j4 + 4)
              d = z(j4) - tau
              dmin = d
              dmin1 = -z(j4)
              if (ieee) then
           ! code for ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       temp = z(j4 + 1)/z(j4 - 2)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4) = z(j4 - 1)*temp
                       emin = min(z(j4),emin)
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       temp = z(j4 + 2)/z(j4 - 3)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4 - 1) = z(j4)*temp
                       emin = min(z(j4 - 1),emin)
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 dmin = min(dmin,dn)
              else
           ! code for non ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       if (d < zero) then
                          return
                       else
                          z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                          d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4))
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       if (d < zero) then
                          return
                       else
                          z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                          d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4 - 1))
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 if (dnm2 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 if (dnm1 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dn)
              end if
           end if
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_dlasq5
#ifdef LA_WITH_XDP
     !> XLASQ5: computes one dqds transform in ping-pong form, one
     !> version for IEEE machines another for non IEEE machines.

     pure subroutine la_xlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dnm1,dnm2, &
               ieee,eps)
        use la_constants_xdp,only:zero,half
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0,n0,pp
           real(xdp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           real(xdp),intent(inout) :: tau
           real(xdp),intent(in) :: sigma,eps
           ! Array Arguments
           real(xdp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(xdp) :: d,emin,temp,dthresh
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           dthresh = eps*(sigma + tau)
           if (tau < dthresh*half) tau = zero
           if (tau /= zero) then
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4) - tau
           dmin = d
           dmin1 = -z(j4)
           if (ieee) then
              ! code for ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    temp = z(j4 + 1)/z(j4 - 2)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4) = z(j4 - 1)*temp
                    emin = min(z(j4),emin)
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    temp = z(j4 + 2)/z(j4 - 3)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4 - 1) = z(j4)*temp
                    emin = min(z(j4 - 1),emin)
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              dmin = min(dmin,dn)
           else
              ! code for non ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    if (d < zero) then
                       return
                    else
                       z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                       d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4))
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    if (d < zero) then
                       return
                    else
                       z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                       d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4 - 1))
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              if (dnm2 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              if (dnm1 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dn)
           end if
           else
           ! this is the version that sets d's to zero if they are small enough
              j4 = 4*i0 + pp - 3
              emin = z(j4 + 4)
              d = z(j4) - tau
              dmin = d
              dmin1 = -z(j4)
              if (ieee) then
           ! code for ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       temp = z(j4 + 1)/z(j4 - 2)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4) = z(j4 - 1)*temp
                       emin = min(z(j4),emin)
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       temp = z(j4 + 2)/z(j4 - 3)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4 - 1) = z(j4)*temp
                       emin = min(z(j4 - 1),emin)
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 dmin = min(dmin,dn)
              else
           ! code for non ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       if (d < zero) then
                          return
                       else
                          z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                          d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4))
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       if (d < zero) then
                          return
                       else
                          z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                          d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4 - 1))
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 if (dnm2 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 if (dnm1 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dn)
              end if
           end if
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_xlasq5
#endif
#ifdef LA_WITH_QP
     !> QLASQ5: computes one dqds transform in ping-pong form, one
     !> version for IEEE machines another for non IEEE machines.

     pure subroutine la_qlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dnm1,dnm2, &
               ieee,eps)
        use la_constants_qp,only:zero,half
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0,n0,pp
           real(qp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           real(qp),intent(inout) :: tau
           real(qp),intent(in) :: sigma,eps
           ! Array Arguments
           real(qp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(qp) :: d,emin,temp,dthresh
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           dthresh = eps*(sigma + tau)
           if (tau < dthresh*half) tau = zero
           if (tau /= zero) then
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4) - tau
           dmin = d
           dmin1 = -z(j4)
           if (ieee) then
              ! code for ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    temp = z(j4 + 1)/z(j4 - 2)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4) = z(j4 - 1)*temp
                    emin = min(z(j4),emin)
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    temp = z(j4 + 2)/z(j4 - 3)
                    d = d*temp - tau
                    dmin = min(dmin,d)
                    z(j4 - 1) = z(j4)*temp
                    emin = min(z(j4 - 1),emin)
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              dmin = min(dmin,dn)
           else
              ! code for non ieee arithmetic.
              if (pp == 0) then
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 2) = d + z(j4 - 1)
                    if (d < zero) then
                       return
                    else
                       z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                       d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4))
                 end do
              else
                 do j4 = 4*i0,4*(n0 - 3),4
                    z(j4 - 3) = d + z(j4)
                    if (d < zero) then
                       return
                    else
                       z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                       d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                    end if
                    dmin = min(dmin,d)
                    emin = min(emin,z(j4 - 1))
                 end do
              end if
              ! unroll last two steps.
              dnm2 = d
              dmin2 = dmin
              j4 = 4*(n0 - 2) - pp
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm2 + z(j4p2)
              if (dnm2 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dnm1)
              dmin1 = dmin
              j4 = j4 + 4
              j4p2 = j4 + 2*pp - 1
              z(j4 - 2) = dnm1 + z(j4p2)
              if (dnm1 < zero) then
                 return
              else
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
              end if
              dmin = min(dmin,dn)
           end if
           else
           ! this is the version that sets d's to zero if they are small enough
              j4 = 4*i0 + pp - 3
              emin = z(j4 + 4)
              d = z(j4) - tau
              dmin = d
              dmin1 = -z(j4)
              if (ieee) then
           ! code for ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       temp = z(j4 + 1)/z(j4 - 2)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4) = z(j4 - 1)*temp
                       emin = min(z(j4),emin)
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       temp = z(j4 + 2)/z(j4 - 3)
                       d = d*temp - tau
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       z(j4 - 1) = z(j4)*temp
                       emin = min(z(j4 - 1),emin)
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                 dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 dmin = min(dmin,dn)
              else
           ! code for non ieee arithmetic.
                 if (pp == 0) then
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 2) = d + z(j4 - 1)
                       if (d < zero) then
                          return
                       else
                          z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                          d = z(j4 + 1)*(d/z(j4 - 2)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4))
                    end do
                 else
                    do j4 = 4*i0,4*(n0 - 3),4
                       z(j4 - 3) = d + z(j4)
                       if (d < zero) then
                          return
                       else
                          z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                          d = z(j4 + 2)*(d/z(j4 - 3)) - tau
                       end if
                       if (d < dthresh) d = zero
                       dmin = min(dmin,d)
                       emin = min(emin,z(j4 - 1))
                    end do
                 end if
           ! unroll last two steps.
                 dnm2 = d
                 dmin2 = dmin
                 j4 = 4*(n0 - 2) - pp
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm2 + z(j4p2)
                 if (dnm2 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dnm1)
                 dmin1 = dmin
                 j4 = j4 + 4
                 j4p2 = j4 + 2*pp - 1
                 z(j4 - 2) = dnm1 + z(j4p2)
                 if (dnm1 < zero) then
                    return
                 else
                    z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
                    dn = z(j4p2 + 2)*(dnm1/z(j4 - 2)) - tau
                 end if
                 dmin = min(dmin,dn)
              end if
           end if
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_qlasq5
#endif

     !> SLASQ6: computes one dqd (shift equal to zero) transform in
     !> ping-pong form, with protection against underflow and overflow.

     pure subroutine la_slasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dnm1,dnm2)
        use la_constants_sp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,pp
           real(sp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           ! Array Arguments
           real(sp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(sp) :: d,emin,safmin,temp
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           safmin = la_slamch('SAFE MINIMUM')
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4)
           dmin = d
           if (pp == 0) then
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 2) = d + z(j4 - 1)
                 if (z(j4 - 2) == zero) then
                    z(j4) = zero
                    d = z(j4 + 1)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 1) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4 + 1)) &
                           then
                    temp = z(j4 + 1)/z(j4 - 2)
                    z(j4) = z(j4 - 1)*temp
                    d = d*temp
                 else
                    z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                    d = z(j4 + 1)*(d/z(j4 - 2))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4))
              end do
           else
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 3) = d + z(j4)
                 if (z(j4 - 3) == zero) then
                    z(j4 - 1) = zero
                    d = z(j4 + 2)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 2) < z(j4 - 3) .and. safmin*z(j4 - 3) < z(j4 + 2)) &
                           then
                    temp = z(j4 + 2)/z(j4 - 3)
                    z(j4 - 1) = z(j4)*temp
                    d = d*temp
                 else
                    z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                    d = z(j4 + 2)*(d/z(j4 - 3))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4 - 1))
              end do
           end if
           ! unroll last two steps.
           dnm2 = d
           dmin2 = dmin
           j4 = 4*(n0 - 2) - pp
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm2 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dnm1 = z(j4p2 + 2)
              dmin = dnm1
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dnm1 = dnm2*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2))
           end if
           dmin = min(dmin,dnm1)
           dmin1 = dmin
           j4 = j4 + 4
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm1 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dn = z(j4p2 + 2)
              dmin = dn
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dn = dnm1*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2))
           end if
           dmin = min(dmin,dn)
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_slasq6
     !> DLASQ6: computes one dqd (shift equal to zero) transform in
     !> ping-pong form, with protection against underflow and overflow.

     pure subroutine la_dlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dnm1,dnm2)
        use la_constants_dp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,pp
           real(dp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           ! Array Arguments
           real(dp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(dp) :: d,emin,safmin,temp
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           safmin = la_dlamch('SAFE MINIMUM')
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4)
           dmin = d
           if (pp == 0) then
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 2) = d + z(j4 - 1)
                 if (z(j4 - 2) == zero) then
                    z(j4) = zero
                    d = z(j4 + 1)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 1) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4 + 1)) &
                           then
                    temp = z(j4 + 1)/z(j4 - 2)
                    z(j4) = z(j4 - 1)*temp
                    d = d*temp
                 else
                    z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                    d = z(j4 + 1)*(d/z(j4 - 2))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4))
              end do
           else
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 3) = d + z(j4)
                 if (z(j4 - 3) == zero) then
                    z(j4 - 1) = zero
                    d = z(j4 + 2)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 2) < z(j4 - 3) .and. safmin*z(j4 - 3) < z(j4 + 2)) &
                           then
                    temp = z(j4 + 2)/z(j4 - 3)
                    z(j4 - 1) = z(j4)*temp
                    d = d*temp
                 else
                    z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                    d = z(j4 + 2)*(d/z(j4 - 3))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4 - 1))
              end do
           end if
           ! unroll last two steps.
           dnm2 = d
           dmin2 = dmin
           j4 = 4*(n0 - 2) - pp
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm2 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dnm1 = z(j4p2 + 2)
              dmin = dnm1
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dnm1 = dnm2*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2))
           end if
           dmin = min(dmin,dnm1)
           dmin1 = dmin
           j4 = j4 + 4
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm1 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dn = z(j4p2 + 2)
              dmin = dn
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dn = dnm1*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2))
           end if
           dmin = min(dmin,dn)
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_dlasq6
#ifdef LA_WITH_XDP
     !> XLASQ6: computes one dqd (shift equal to zero) transform in
     !> ping-pong form, with protection against underflow and overflow.

     pure subroutine la_xlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dnm1,dnm2)
        use la_constants_xdp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,pp
           real(xdp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           ! Array Arguments
           real(xdp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(xdp) :: d,emin,safmin,temp
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           safmin = la_xlamch('SAFE MINIMUM')
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4)
           dmin = d
           if (pp == 0) then
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 2) = d + z(j4 - 1)
                 if (z(j4 - 2) == zero) then
                    z(j4) = zero
                    d = z(j4 + 1)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 1) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4 + 1)) &
                           then
                    temp = z(j4 + 1)/z(j4 - 2)
                    z(j4) = z(j4 - 1)*temp
                    d = d*temp
                 else
                    z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                    d = z(j4 + 1)*(d/z(j4 - 2))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4))
              end do
           else
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 3) = d + z(j4)
                 if (z(j4 - 3) == zero) then
                    z(j4 - 1) = zero
                    d = z(j4 + 2)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 2) < z(j4 - 3) .and. safmin*z(j4 - 3) < z(j4 + 2)) &
                           then
                    temp = z(j4 + 2)/z(j4 - 3)
                    z(j4 - 1) = z(j4)*temp
                    d = d*temp
                 else
                    z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                    d = z(j4 + 2)*(d/z(j4 - 3))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4 - 1))
              end do
           end if
           ! unroll last two steps.
           dnm2 = d
           dmin2 = dmin
           j4 = 4*(n0 - 2) - pp
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm2 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dnm1 = z(j4p2 + 2)
              dmin = dnm1
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dnm1 = dnm2*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2))
           end if
           dmin = min(dmin,dnm1)
           dmin1 = dmin
           j4 = j4 + 4
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm1 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dn = z(j4p2 + 2)
              dmin = dn
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dn = dnm1*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2))
           end if
           dmin = min(dmin,dn)
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_xlasq6
#endif
#ifdef LA_WITH_QP
     !> QLASQ6: computes one dqd (shift equal to zero) transform in
     !> ping-pong form, with protection against underflow and overflow.

     pure subroutine la_qlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dnm1,dnm2)
        use la_constants_qp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: i0,n0,pp
           real(qp),intent(out) :: dmin,dmin1,dmin2,dn,dnm1,dnm2
           ! Array Arguments
           real(qp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j4,j4p2
           real(qp) :: d,emin,safmin,temp
           ! Intrinsic Functions
           intrinsic :: min
           ! Executable Statements
           if ((n0 - i0 - 1) <= 0) return
           safmin = la_qlamch('SAFE MINIMUM')
           j4 = 4*i0 + pp - 3
           emin = z(j4 + 4)
           d = z(j4)
           dmin = d
           if (pp == 0) then
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 2) = d + z(j4 - 1)
                 if (z(j4 - 2) == zero) then
                    z(j4) = zero
                    d = z(j4 + 1)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 1) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4 + 1)) &
                           then
                    temp = z(j4 + 1)/z(j4 - 2)
                    z(j4) = z(j4 - 1)*temp
                    d = d*temp
                 else
                    z(j4) = z(j4 + 1)*(z(j4 - 1)/z(j4 - 2))
                    d = z(j4 + 1)*(d/z(j4 - 2))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4))
              end do
           else
              do j4 = 4*i0,4*(n0 - 3),4
                 z(j4 - 3) = d + z(j4)
                 if (z(j4 - 3) == zero) then
                    z(j4 - 1) = zero
                    d = z(j4 + 2)
                    dmin = d
                    emin = zero
                 else if (safmin*z(j4 + 2) < z(j4 - 3) .and. safmin*z(j4 - 3) < z(j4 + 2)) &
                           then
                    temp = z(j4 + 2)/z(j4 - 3)
                    z(j4 - 1) = z(j4)*temp
                    d = d*temp
                 else
                    z(j4 - 1) = z(j4 + 2)*(z(j4)/z(j4 - 3))
                    d = z(j4 + 2)*(d/z(j4 - 3))
                 end if
                 dmin = min(dmin,d)
                 emin = min(emin,z(j4 - 1))
              end do
           end if
           ! unroll last two steps.
           dnm2 = d
           dmin2 = dmin
           j4 = 4*(n0 - 2) - pp
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm2 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dnm1 = z(j4p2 + 2)
              dmin = dnm1
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dnm1 = dnm2*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dnm1 = z(j4p2 + 2)*(dnm2/z(j4 - 2))
           end if
           dmin = min(dmin,dnm1)
           dmin1 = dmin
           j4 = j4 + 4
           j4p2 = j4 + 2*pp - 1
           z(j4 - 2) = dnm1 + z(j4p2)
           if (z(j4 - 2) == zero) then
              z(j4) = zero
              dn = z(j4p2 + 2)
              dmin = dn
              emin = zero
           else if (safmin*z(j4p2 + 2) < z(j4 - 2) .and. safmin*z(j4 - 2) < z(j4p2 + 2)) then
              temp = z(j4p2 + 2)/z(j4 - 2)
              z(j4) = z(j4p2)*temp
              dn = dnm1*temp
           else
              z(j4) = z(j4p2 + 2)*(z(j4p2)/z(j4 - 2))
              dn = z(j4p2 + 2)*(dnm1/z(j4 - 2))
           end if
           dmin = min(dmin,dn)
           z(j4 + 2) = dn
           z(4*n0 - pp) = emin
           return
     end subroutine la_qlasq6
#endif

     !> SLASQ3: checks for deflation, computes a shift (TAU) and calls dqds.
     !> In case of failure it changes shifts, and tries again until output
     !> is positive.

     pure subroutine la_slasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
               ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
        use la_constants_sp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0
           integer(ilp),intent(inout) :: iter,n0,ndiv,nfail,pp
           real(sp),intent(inout) :: desig,dmin1,dmin2,dn,dn1,dn2,g,qmax,tau
           real(sp),intent(out) :: dmin,sigma
           ! Array Arguments
           real(sp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: cbias = 1.50_sp
           real(sp),parameter :: qurtr = 0.250_sp
           real(sp),parameter :: hundrd = 100.0_sp

           ! Local Scalars
           integer(ilp) :: ipn4,j4,n0in,nn
           integer(ilp),intent(inout) :: ttype
           real(sp) :: eps,s,t,temp,tol,tol2
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           n0in = n0
           eps = la_slamch('PRECISION')
           tol = eps*hundrd
           tol2 = tol**2
           ! check for deflation.
           10 continue
           if (n0 < i0) return
           if (n0 == i0) go to 20
           nn = 4*n0 + pp
           if (n0 == (i0 + 1)) go to 40
           ! check whether e(n0-1) is negligible, 1 eigenvalue.
           if (z(nn - 5) > tol2*(sigma + z(nn - 3)) .and. z(nn - 2*pp - 4) > tol2*z(nn - 7)) go to &
                     30
                     20 continue
           z(4*n0 - 3) = z(4*n0 + pp - 3) + sigma
           n0 = n0 - 1
           go to 10
           ! check  whether e(n0-2) is negligible, 2 eigenvalues.
           30 continue
           if (z(nn - 9) > tol2*sigma .and. z(nn - 2*pp - 8) > tol2*z(nn - 11)) go to 50
           40 continue
           if (z(nn - 3) > z(nn - 7)) then
              s = z(nn - 3)
              z(nn - 3) = z(nn - 7)
              z(nn - 7) = s
           end if
           t = half*((z(nn - 7) - z(nn - 3)) + z(nn - 5))
           if (z(nn - 5) > z(nn - 3)*tol2 .and. t /= zero) then
              s = z(nn - 3)*(z(nn - 5)/t)
              if (s <= t) then
                 s = z(nn - 3)*(z(nn - 5)/(t*(one + sqrt(one + s/t))))
              else
                 s = z(nn - 3)*(z(nn - 5)/(t + sqrt(t)*sqrt(t + s)))
              end if
              t = z(nn - 7) + (s + z(nn - 5))
              z(nn - 3) = z(nn - 3)*(z(nn - 7)/t)
              z(nn - 7) = t
           end if
           z(4*n0 - 7) = z(nn - 7) + sigma
           z(4*n0 - 3) = z(nn - 3) + sigma
           n0 = n0 - 2
           go to 10
           50 continue
           if (pp == 2) pp = 0
           ! reverse the qd-array, if warranted.
           if (dmin <= zero .or. n0 < n0in) then
              if (cbias*z(4*i0 + pp - 3) < z(4*n0 + pp - 3)) then
                 ipn4 = 4*(i0 + n0)
                 do j4 = 4*i0,2*(i0 + n0 - 1),4
                    temp = z(j4 - 3)
                    z(j4 - 3) = z(ipn4 - j4 - 3)
                    z(ipn4 - j4 - 3) = temp
                    temp = z(j4 - 2)
                    z(j4 - 2) = z(ipn4 - j4 - 2)
                    z(ipn4 - j4 - 2) = temp
                    temp = z(j4 - 1)
                    z(j4 - 1) = z(ipn4 - j4 - 5)
                    z(ipn4 - j4 - 5) = temp
                    temp = z(j4)
                    z(j4) = z(ipn4 - j4 - 4)
                    z(ipn4 - j4 - 4) = temp
                 end do
                 if (n0 - i0 <= 4) then
                    z(4*n0 + pp - 1) = z(4*i0 + pp - 1)
                    z(4*n0 - pp) = z(4*i0 - pp)
                 end if
                 dmin2 = min(dmin2,z(4*n0 + pp - 1))
                 z(4*n0 + pp - 1) = min(z(4*n0 + pp - 1),z(4*i0 + pp - 1),z(4*i0 + pp + 3))
                 z(4*n0 - pp) = min(z(4*n0 - pp),z(4*i0 - pp),z(4*i0 - pp + 4))
                 qmax = max(qmax,z(4*i0 + pp - 3),z(4*i0 + pp + 1))
                 dmin = -zero
              end if
           end if
           ! choose a shift.
           call la_slasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau,ttype, &
                     g)
           ! call dqds until dmin > 0.
           70 continue
           call la_slasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dn1,dn2,ieee, &
                     eps)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           ! check status.
           if (dmin >= zero .and. dmin1 >= zero) then
              ! success.
              go to 90
           else if (dmin < zero .and. dmin1 > zero .and. z(4*(n0 - 1) - pp) < tol*(sigma + dn1) .and. abs( &
                      dn) < tol*sigma) then
              ! convergence hidden by negative dn.
              z(4*(n0 - 1) - pp + 2) = zero
              dmin = zero
              go to 90
           else if (dmin < zero) then
              ! tau too big. select new tau and try again.
              nfail = nfail + 1
              if (ttype < -22) then
                 ! failed twice. play it safe.
                 tau = zero
              else if (dmin1 > zero) then
                 ! late failure. gives excellent shift.
                 tau = (tau + dmin)*(one - two*eps)
                 ttype = ttype - 11
              else
                 ! early failure. divide by 4.
                 tau = qurtr*tau
                 ttype = ttype - 12
              end if
              go to 70
           else if (la_sisnan(dmin)) then
              ! nan.
              if (tau == zero) then
                 go to 80
              else
                 tau = zero
                 go to 70
              end if
           else
              ! possible underflow. play it safe.
              go to 80
           end if
           ! risk of underflow.
           80 continue
           call la_slasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dn1,dn2)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           tau = zero
           90 continue
           if (tau < sigma) then
              desig = desig + tau
              t = sigma + desig
              desig = desig - (t - sigma)
           else
              t = sigma + tau
              desig = sigma - (t - tau) + desig
           end if
           sigma = t
           return
     end subroutine la_slasq3
     !> DLASQ3: checks for deflation, computes a shift (TAU) and calls dqds.
     !> In case of failure it changes shifts, and tries again until output
     !> is positive.

     pure subroutine la_dlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
               ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
        use la_constants_dp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0
           integer(ilp),intent(inout) :: iter,n0,ndiv,nfail,pp
           real(dp),intent(inout) :: desig,dmin1,dmin2,dn,dn1,dn2,g,qmax,tau
           real(dp),intent(out) :: dmin,sigma
           ! Array Arguments
           real(dp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: cbias = 1.50_dp
           real(dp),parameter :: qurtr = 0.250_dp
           real(dp),parameter :: hundrd = 100.0_dp

           ! Local Scalars
           integer(ilp) :: ipn4,j4,n0in,nn
           integer(ilp),intent(inout) :: ttype
           real(dp) :: eps,s,t,temp,tol,tol2
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           n0in = n0
           eps = la_dlamch('PRECISION')
           tol = eps*hundrd
           tol2 = tol**2
           ! check for deflation.
           10 continue
           if (n0 < i0) return
           if (n0 == i0) go to 20
           nn = 4*n0 + pp
           if (n0 == (i0 + 1)) go to 40
           ! check whether e(n0-1) is negligible, 1 eigenvalue.
           if (z(nn - 5) > tol2*(sigma + z(nn - 3)) .and. z(nn - 2*pp - 4) > tol2*z(nn - 7)) go to &
                     30
                     20 continue
           z(4*n0 - 3) = z(4*n0 + pp - 3) + sigma
           n0 = n0 - 1
           go to 10
           ! check  whether e(n0-2) is negligible, 2 eigenvalues.
           30 continue
           if (z(nn - 9) > tol2*sigma .and. z(nn - 2*pp - 8) > tol2*z(nn - 11)) go to 50
           40 continue
           if (z(nn - 3) > z(nn - 7)) then
              s = z(nn - 3)
              z(nn - 3) = z(nn - 7)
              z(nn - 7) = s
           end if
           t = half*((z(nn - 7) - z(nn - 3)) + z(nn - 5))
           if (z(nn - 5) > z(nn - 3)*tol2 .and. t /= zero) then
              s = z(nn - 3)*(z(nn - 5)/t)
              if (s <= t) then
                 s = z(nn - 3)*(z(nn - 5)/(t*(one + sqrt(one + s/t))))
              else
                 s = z(nn - 3)*(z(nn - 5)/(t + sqrt(t)*sqrt(t + s)))
              end if
              t = z(nn - 7) + (s + z(nn - 5))
              z(nn - 3) = z(nn - 3)*(z(nn - 7)/t)
              z(nn - 7) = t
           end if
           z(4*n0 - 7) = z(nn - 7) + sigma
           z(4*n0 - 3) = z(nn - 3) + sigma
           n0 = n0 - 2
           go to 10
           50 continue
           if (pp == 2) pp = 0
           ! reverse the qd-array, if warranted.
           if (dmin <= zero .or. n0 < n0in) then
              if (cbias*z(4*i0 + pp - 3) < z(4*n0 + pp - 3)) then
                 ipn4 = 4*(i0 + n0)
                 do j4 = 4*i0,2*(i0 + n0 - 1),4
                    temp = z(j4 - 3)
                    z(j4 - 3) = z(ipn4 - j4 - 3)
                    z(ipn4 - j4 - 3) = temp
                    temp = z(j4 - 2)
                    z(j4 - 2) = z(ipn4 - j4 - 2)
                    z(ipn4 - j4 - 2) = temp
                    temp = z(j4 - 1)
                    z(j4 - 1) = z(ipn4 - j4 - 5)
                    z(ipn4 - j4 - 5) = temp
                    temp = z(j4)
                    z(j4) = z(ipn4 - j4 - 4)
                    z(ipn4 - j4 - 4) = temp
                 end do
                 if (n0 - i0 <= 4) then
                    z(4*n0 + pp - 1) = z(4*i0 + pp - 1)
                    z(4*n0 - pp) = z(4*i0 - pp)
                 end if
                 dmin2 = min(dmin2,z(4*n0 + pp - 1))
                 z(4*n0 + pp - 1) = min(z(4*n0 + pp - 1),z(4*i0 + pp - 1),z(4*i0 + pp + 3))
                 z(4*n0 - pp) = min(z(4*n0 - pp),z(4*i0 - pp),z(4*i0 - pp + 4))
                 qmax = max(qmax,z(4*i0 + pp - 3),z(4*i0 + pp + 1))
                 dmin = -zero
              end if
           end if
           ! choose a shift.
           call la_dlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau,ttype, &
                     g)
           ! call dqds until dmin > 0.
           70 continue
           call la_dlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dn1,dn2,ieee, &
                     eps)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           ! check status.
           if (dmin >= zero .and. dmin1 >= zero) then
              ! success.
              go to 90
           else if (dmin < zero .and. dmin1 > zero .and. z(4*(n0 - 1) - pp) < tol*(sigma + dn1) .and. abs( &
                      dn) < tol*sigma) then
              ! convergence hidden by negative dn.
              z(4*(n0 - 1) - pp + 2) = zero
              dmin = zero
              go to 90
           else if (dmin < zero) then
              ! tau too big. select new tau and try again.
              nfail = nfail + 1
              if (ttype < -22) then
                 ! failed twice. play it safe.
                 tau = zero
              else if (dmin1 > zero) then
                 ! late failure. gives excellent shift.
                 tau = (tau + dmin)*(one - two*eps)
                 ttype = ttype - 11
              else
                 ! early failure. divide by 4.
                 tau = qurtr*tau
                 ttype = ttype - 12
              end if
              go to 70
           else if (la_disnan(dmin)) then
              ! nan.
              if (tau == zero) then
                 go to 80
              else
                 tau = zero
                 go to 70
              end if
           else
              ! possible underflow. play it safe.
              go to 80
           end if
           ! risk of underflow.
           80 continue
           call la_dlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dn1,dn2)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           tau = zero
           90 continue
           if (tau < sigma) then
              desig = desig + tau
              t = sigma + desig
              desig = desig - (t - sigma)
           else
              t = sigma + tau
              desig = sigma - (t - tau) + desig
           end if
           sigma = t
           return
     end subroutine la_dlasq3
#ifdef LA_WITH_XDP
     !> XLASQ3: checks for deflation, computes a shift (TAU) and calls dqds.
     !> In case of failure it changes shifts, and tries again until output
     !> is positive.

     pure subroutine la_xlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
               ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
        use la_constants_xdp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0
           integer(ilp),intent(inout) :: iter,n0,ndiv,nfail,pp
           real(xdp),intent(inout) :: desig,dmin1,dmin2,dn,dn1,dn2,g,qmax,tau
           real(xdp),intent(out) :: dmin,sigma
           ! Array Arguments
           real(xdp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: cbias = 1.50_xdp
           real(xdp),parameter :: qurtr = 0.250_xdp
           real(xdp),parameter :: hundrd = 100.0_xdp

           ! Local Scalars
           integer(ilp) :: ipn4,j4,n0in,nn
           integer(ilp),intent(inout) :: ttype
           real(xdp) :: eps,s,t,temp,tol,tol2
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           n0in = n0
           eps = la_xlamch('PRECISION')
           tol = eps*hundrd
           tol2 = tol**2
           ! check for deflation.
           10 continue
           if (n0 < i0) return
           if (n0 == i0) go to 20
           nn = 4*n0 + pp
           if (n0 == (i0 + 1)) go to 40
           ! check whether e(n0-1) is negligible, 1 eigenvalue.
           if (z(nn - 5) > tol2*(sigma + z(nn - 3)) .and. z(nn - 2*pp - 4) > tol2*z(nn - 7)) go to &
                     30
                     20 continue
           z(4*n0 - 3) = z(4*n0 + pp - 3) + sigma
           n0 = n0 - 1
           go to 10
           ! check  whether e(n0-2) is negligible, 2 eigenvalues.
           30 continue
           if (z(nn - 9) > tol2*sigma .and. z(nn - 2*pp - 8) > tol2*z(nn - 11)) go to 50
           40 continue
           if (z(nn - 3) > z(nn - 7)) then
              s = z(nn - 3)
              z(nn - 3) = z(nn - 7)
              z(nn - 7) = s
           end if
           t = half*((z(nn - 7) - z(nn - 3)) + z(nn - 5))
           if (z(nn - 5) > z(nn - 3)*tol2 .and. t /= zero) then
              s = z(nn - 3)*(z(nn - 5)/t)
              if (s <= t) then
                 s = z(nn - 3)*(z(nn - 5)/(t*(one + sqrt(one + s/t))))
              else
                 s = z(nn - 3)*(z(nn - 5)/(t + sqrt(t)*sqrt(t + s)))
              end if
              t = z(nn - 7) + (s + z(nn - 5))
              z(nn - 3) = z(nn - 3)*(z(nn - 7)/t)
              z(nn - 7) = t
           end if
           z(4*n0 - 7) = z(nn - 7) + sigma
           z(4*n0 - 3) = z(nn - 3) + sigma
           n0 = n0 - 2
           go to 10
           50 continue
           if (pp == 2) pp = 0
           ! reverse the qd-array, if warranted.
           if (dmin <= zero .or. n0 < n0in) then
              if (cbias*z(4*i0 + pp - 3) < z(4*n0 + pp - 3)) then
                 ipn4 = 4*(i0 + n0)
                 do j4 = 4*i0,2*(i0 + n0 - 1),4
                    temp = z(j4 - 3)
                    z(j4 - 3) = z(ipn4 - j4 - 3)
                    z(ipn4 - j4 - 3) = temp
                    temp = z(j4 - 2)
                    z(j4 - 2) = z(ipn4 - j4 - 2)
                    z(ipn4 - j4 - 2) = temp
                    temp = z(j4 - 1)
                    z(j4 - 1) = z(ipn4 - j4 - 5)
                    z(ipn4 - j4 - 5) = temp
                    temp = z(j4)
                    z(j4) = z(ipn4 - j4 - 4)
                    z(ipn4 - j4 - 4) = temp
                 end do
                 if (n0 - i0 <= 4) then
                    z(4*n0 + pp - 1) = z(4*i0 + pp - 1)
                    z(4*n0 - pp) = z(4*i0 - pp)
                 end if
                 dmin2 = min(dmin2,z(4*n0 + pp - 1))
                 z(4*n0 + pp - 1) = min(z(4*n0 + pp - 1),z(4*i0 + pp - 1),z(4*i0 + pp + 3))
                 z(4*n0 - pp) = min(z(4*n0 - pp),z(4*i0 - pp),z(4*i0 - pp + 4))
                 qmax = max(qmax,z(4*i0 + pp - 3),z(4*i0 + pp + 1))
                 dmin = -zero
              end if
           end if
           ! choose a shift.
           call la_xlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau,ttype, &
                     g)
           ! call dqds until dmin > 0.
           70 continue
           call la_xlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dn1,dn2,ieee, &
                     eps)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           ! check status.
           if (dmin >= zero .and. dmin1 >= zero) then
              ! success.
              go to 90
           else if (dmin < zero .and. dmin1 > zero .and. z(4*(n0 - 1) - pp) < tol*(sigma + dn1) .and. abs( &
                      dn) < tol*sigma) then
              ! convergence hidden by negative dn.
              z(4*(n0 - 1) - pp + 2) = zero
              dmin = zero
              go to 90
           else if (dmin < zero) then
              ! tau too big. select new tau and try again.
              nfail = nfail + 1
              if (ttype < -22) then
                 ! failed twice. play it safe.
                 tau = zero
              else if (dmin1 > zero) then
                 ! late failure. gives excellent shift.
                 tau = (tau + dmin)*(one - two*eps)
                 ttype = ttype - 11
              else
                 ! early failure. divide by 4.
                 tau = qurtr*tau
                 ttype = ttype - 12
              end if
              go to 70
           else if (la_xisnan(dmin)) then
              ! nan.
              if (tau == zero) then
                 go to 80
              else
                 tau = zero
                 go to 70
              end if
           else
              ! possible underflow. play it safe.
              go to 80
           end if
           ! risk of underflow.
           80 continue
           call la_xlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dn1,dn2)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           tau = zero
           90 continue
           if (tau < sigma) then
              desig = desig + tau
              t = sigma + desig
              desig = desig - (t - sigma)
           else
              t = sigma + tau
              desig = sigma - (t - tau) + desig
           end if
           sigma = t
           return
     end subroutine la_xlasq3
#endif
#ifdef LA_WITH_QP
     !> QLASQ3: checks for deflation, computes a shift (TAU) and calls dqds.
     !> In case of failure it changes shifts, and tries again until output
     !> is positive.

     pure subroutine la_qlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
               ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
        use la_constants_qp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: ieee
           integer(ilp),intent(in) :: i0
           integer(ilp),intent(inout) :: iter,n0,ndiv,nfail,pp
           real(qp),intent(inout) :: desig,dmin1,dmin2,dn,dn1,dn2,g,qmax,tau
           real(qp),intent(out) :: dmin,sigma
           ! Array Arguments
           real(qp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: cbias = 1.50_qp
           real(qp),parameter :: qurtr = 0.250_qp
           real(qp),parameter :: hundrd = 100.0_qp

           ! Local Scalars
           integer(ilp) :: ipn4,j4,n0in,nn
           integer(ilp),intent(inout) :: ttype
           real(qp) :: eps,s,t,temp,tol,tol2
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           n0in = n0
           eps = la_qlamch('PRECISION')
           tol = eps*hundrd
           tol2 = tol**2
           ! check for deflation.
           10 continue
           if (n0 < i0) return
           if (n0 == i0) go to 20
           nn = 4*n0 + pp
           if (n0 == (i0 + 1)) go to 40
           ! check whether e(n0-1) is negligible, 1 eigenvalue.
           if (z(nn - 5) > tol2*(sigma + z(nn - 3)) .and. z(nn - 2*pp - 4) > tol2*z(nn - 7)) go to &
                     30
                     20 continue
           z(4*n0 - 3) = z(4*n0 + pp - 3) + sigma
           n0 = n0 - 1
           go to 10
           ! check  whether e(n0-2) is negligible, 2 eigenvalues.
           30 continue
           if (z(nn - 9) > tol2*sigma .and. z(nn - 2*pp - 8) > tol2*z(nn - 11)) go to 50
           40 continue
           if (z(nn - 3) > z(nn - 7)) then
              s = z(nn - 3)
              z(nn - 3) = z(nn - 7)
              z(nn - 7) = s
           end if
           t = half*((z(nn - 7) - z(nn - 3)) + z(nn - 5))
           if (z(nn - 5) > z(nn - 3)*tol2 .and. t /= zero) then
              s = z(nn - 3)*(z(nn - 5)/t)
              if (s <= t) then
                 s = z(nn - 3)*(z(nn - 5)/(t*(one + sqrt(one + s/t))))
              else
                 s = z(nn - 3)*(z(nn - 5)/(t + sqrt(t)*sqrt(t + s)))
              end if
              t = z(nn - 7) + (s + z(nn - 5))
              z(nn - 3) = z(nn - 3)*(z(nn - 7)/t)
              z(nn - 7) = t
           end if
           z(4*n0 - 7) = z(nn - 7) + sigma
           z(4*n0 - 3) = z(nn - 3) + sigma
           n0 = n0 - 2
           go to 10
           50 continue
           if (pp == 2) pp = 0
           ! reverse the qd-array, if warranted.
           if (dmin <= zero .or. n0 < n0in) then
              if (cbias*z(4*i0 + pp - 3) < z(4*n0 + pp - 3)) then
                 ipn4 = 4*(i0 + n0)
                 do j4 = 4*i0,2*(i0 + n0 - 1),4
                    temp = z(j4 - 3)
                    z(j4 - 3) = z(ipn4 - j4 - 3)
                    z(ipn4 - j4 - 3) = temp
                    temp = z(j4 - 2)
                    z(j4 - 2) = z(ipn4 - j4 - 2)
                    z(ipn4 - j4 - 2) = temp
                    temp = z(j4 - 1)
                    z(j4 - 1) = z(ipn4 - j4 - 5)
                    z(ipn4 - j4 - 5) = temp
                    temp = z(j4)
                    z(j4) = z(ipn4 - j4 - 4)
                    z(ipn4 - j4 - 4) = temp
                 end do
                 if (n0 - i0 <= 4) then
                    z(4*n0 + pp - 1) = z(4*i0 + pp - 1)
                    z(4*n0 - pp) = z(4*i0 - pp)
                 end if
                 dmin2 = min(dmin2,z(4*n0 + pp - 1))
                 z(4*n0 + pp - 1) = min(z(4*n0 + pp - 1),z(4*i0 + pp - 1),z(4*i0 + pp + 3))
                 z(4*n0 - pp) = min(z(4*n0 - pp),z(4*i0 - pp),z(4*i0 - pp + 4))
                 qmax = max(qmax,z(4*i0 + pp - 3),z(4*i0 + pp + 1))
                 dmin = -zero
              end if
           end if
           ! choose a shift.
           call la_qlasq4(i0,n0,z,pp,n0in,dmin,dmin1,dmin2,dn,dn1,dn2,tau,ttype, &
                     g)
           ! call dqds until dmin > 0.
           70 continue
           call la_qlasq5(i0,n0,z,pp,tau,sigma,dmin,dmin1,dmin2,dn,dn1,dn2,ieee, &
                     eps)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           ! check status.
           if (dmin >= zero .and. dmin1 >= zero) then
              ! success.
              go to 90
           else if (dmin < zero .and. dmin1 > zero .and. z(4*(n0 - 1) - pp) < tol*(sigma + dn1) .and. abs( &
                      dn) < tol*sigma) then
              ! convergence hidden by negative dn.
              z(4*(n0 - 1) - pp + 2) = zero
              dmin = zero
              go to 90
           else if (dmin < zero) then
              ! tau too big. select new tau and try again.
              nfail = nfail + 1
              if (ttype < -22) then
                 ! failed twice. play it safe.
                 tau = zero
              else if (dmin1 > zero) then
                 ! late failure. gives excellent shift.
                 tau = (tau + dmin)*(one - two*eps)
                 ttype = ttype - 11
              else
                 ! early failure. divide by 4.
                 tau = qurtr*tau
                 ttype = ttype - 12
              end if
              go to 70
           else if (la_qisnan(dmin)) then
              ! nan.
              if (tau == zero) then
                 go to 80
              else
                 tau = zero
                 go to 70
              end if
           else
              ! possible underflow. play it safe.
              go to 80
           end if
           ! risk of underflow.
           80 continue
           call la_qlasq6(i0,n0,z,pp,dmin,dmin1,dmin2,dn,dn1,dn2)
           ndiv = ndiv + (n0 - i0 + 2)
           iter = iter + 1
           tau = zero
           90 continue
           if (tau < sigma) then
              desig = desig + tau
              t = sigma + desig
              desig = desig - (t - sigma)
           else
              t = sigma + tau
              desig = sigma - (t - tau) + desig
           end if
           sigma = t
           return
     end subroutine la_qlasq3
#endif

     !> SBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**T
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**T*VT instead of
     !> P**T, for given real input matrices U and VT.  When U and VT are the
     !> orthogonal matrices that reduce a general matrix A to bidiagonal
     !> form:  A = U*B*VT, as computed by SGEBRD, then
     !> A = (U*Q) * S * (P**T*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
     !> for a given real input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_sbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work, &
               info)
        use la_constants_sp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(sp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: hndrth = 0.01_sp
           real(sp),parameter :: hndrd = 100.0_sp
           real(sp),parameter :: meigth = -0.125_sp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,iterdivn,j,ll,lll,m,maxitdivn,nm1,nm12, &
                     nm13,oldll,oldm
           real(sp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,max,min,real,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('SBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_slasq1(n,d,e,work,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_slamch('EPSILON')
           unfl = la_slamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_slartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 work(i) = cs
                 work(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_slasr('R','V','F',nru,n,work(1),work(n),u,ldu)

              if (ncc > 0) call la_slasr('L','V','F',n,ncc,work(1),work(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=sp))
              thresh = max(tol*sminoa,maxitr*(n*(n*unfl)))
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*(n*(n*unfl)))
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxitdivn = maxitr*n
           iterdivn = 0
           iter = -1
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter >= n) then
              iter = iter - n
              iterdivn = iterdivn + 1
              if (iterdivn >= maxitdivn) go to 200
           end if
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_slasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_srot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr,sinr &
                        )
              if (nru > 0) call la_srot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)
              if (ncc > 0) call la_srot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_slas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_slas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_slartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_slartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll + 1) = cs
                    work(i - ll + 1 + nm1) = sn
                    work(i - ll + 1 + nm12) = oldcs
                    work(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_slasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_slasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_slasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_slartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_slartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll) = cs
                    work(i - ll + nm1) = -sn
                    work(i - ll + nm12) = oldcs
                    work(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_slasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_slasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_slasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_slartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_slartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    work(i - ll + 1) = cosr
                    work(i - ll + 1 + nm1) = sinr
                    work(i - ll + 1 + nm12) = cosl
                    work(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_slasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_slasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_slasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_slartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_slartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    work(i - ll) = cosr
                    work(i - ll + nm1) = -sinr
                    work(i - ll + nm12) = cosl
                    work(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_slasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_slasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_slasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_sscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_sswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_sswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_sswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_sbdsqr
     !> DBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**T
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**T*VT instead of
     !> P**T, for given real input matrices U and VT.  When U and VT are the
     !> orthogonal matrices that reduce a general matrix A to bidiagonal
     !> form:  A = U*B*VT, as computed by DGEBRD, then
     !> A = (U*Q) * S * (P**T*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
     !> for a given real input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_dbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work, &
               info)
        use la_constants_dp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(dp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: hndrth = 0.01_dp
           real(dp),parameter :: hndrd = 100.0_dp
           real(dp),parameter :: meigth = -0.125_dp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,iterdivn,j,ll,lll,m,maxitdivn,nm1,nm12, &
                     nm13,oldll,oldm
           real(dp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('DBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_dlasq1(n,d,e,work,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_dlamch('EPSILON')
           unfl = la_dlamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_dlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 work(i) = cs
                 work(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_dlasr('R','V','F',nru,n,work(1),work(n),u,ldu)

              if (ncc > 0) call la_dlasr('L','V','F',n,ncc,work(1),work(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=dp))
              thresh = max(tol*sminoa,maxitr*(n*(n*unfl)))
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*(n*(n*unfl)))
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxitdivn = maxitr*n
           iterdivn = 0
           iter = -1
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter >= n) then
              iter = iter - n
              iterdivn = iterdivn + 1
              if (iterdivn >= maxitdivn) go to 200
           end if
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_dlasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_drot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr,sinr &
                        )
              if (nru > 0) call la_drot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)
              if (ncc > 0) call la_drot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_dlas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_dlas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_dlartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_dlartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll + 1) = cs
                    work(i - ll + 1 + nm1) = sn
                    work(i - ll + 1 + nm12) = oldcs
                    work(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_dlasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_dlasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_dlasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_dlartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_dlartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll) = cs
                    work(i - ll + nm1) = -sn
                    work(i - ll + nm12) = oldcs
                    work(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_dlasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_dlasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_dlasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_dlartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_dlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    work(i - ll + 1) = cosr
                    work(i - ll + 1 + nm1) = sinr
                    work(i - ll + 1 + nm12) = cosl
                    work(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_dlasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_dlasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_dlasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_dlartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_dlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    work(i - ll) = cosr
                    work(i - ll + nm1) = -sinr
                    work(i - ll + nm12) = cosl
                    work(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_dlasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_dlasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_dlasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_dscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_dswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_dswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_dswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_dbdsqr
#ifdef LA_WITH_XDP
     !> XBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**T
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**T*VT instead of
     !> P**T, for given real input matrices U and VT.  When U and VT are the
     !> orthogonal matrices that reduce a general matrix A to bidiagonal
     !> form:  A = U*B*VT, as computed by XGEBRD, then
     !> A = (U*Q) * S * (P**T*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
     !> for a given real input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_xbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work, &
               info)
        use la_constants_xdp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(xdp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: hndrth = 0.01_xdp
           real(xdp),parameter :: hndrd = 100.0_xdp
           real(xdp),parameter :: meigth = -0.125_xdp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,iterdivn,j,ll,lll,m,maxitdivn,nm1,nm12, &
                     nm13,oldll,oldm
           real(xdp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('XBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_xlasq1(n,d,e,work,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_xlamch('EPSILON')
           unfl = la_xlamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_xlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 work(i) = cs
                 work(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_xlasr('R','V','F',nru,n,work(1),work(n),u,ldu)

              if (ncc > 0) call la_xlasr('L','V','F',n,ncc,work(1),work(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=xdp))
              thresh = max(tol*sminoa,maxitr*(n*(n*unfl)))
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*(n*(n*unfl)))
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxitdivn = maxitr*n
           iterdivn = 0
           iter = -1
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter >= n) then
              iter = iter - n
              iterdivn = iterdivn + 1
              if (iterdivn >= maxitdivn) go to 200
           end if
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_xlasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_xrot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr,sinr &
                        )
              if (nru > 0) call la_xrot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)
              if (ncc > 0) call la_xrot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_xlas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_xlas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_xlartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_xlartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll + 1) = cs
                    work(i - ll + 1 + nm1) = sn
                    work(i - ll + 1 + nm12) = oldcs
                    work(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_xlasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_xlasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_xlasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_xlartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_xlartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll) = cs
                    work(i - ll + nm1) = -sn
                    work(i - ll + nm12) = oldcs
                    work(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_xlasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_xlasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_xlasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_xlartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_xlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    work(i - ll + 1) = cosr
                    work(i - ll + 1 + nm1) = sinr
                    work(i - ll + 1 + nm12) = cosl
                    work(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_xlasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_xlasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_xlasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_xlartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_xlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    work(i - ll) = cosr
                    work(i - ll + nm1) = -sinr
                    work(i - ll + nm12) = cosl
                    work(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_xlasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_xlasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_xlasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_xscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_xswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_xswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_xswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_xbdsqr
#endif
#ifdef LA_WITH_QP
     !> QBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**T
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**T*VT instead of
     !> P**T, for given real input matrices U and VT.  When U and VT are the
     !> orthogonal matrices that reduce a general matrix A to bidiagonal
     !> form:  A = U*B*VT, as computed by QGEBRD, then
     !> A = (U*Q) * S * (P**T*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
     !> for a given real input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_qbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,work, &
               info)
        use la_constants_qp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(qp),intent(inout) :: c(ldc,*),d(*),e(*),u(ldu,*),vt(ldvt,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: hndrth = 0.01_qp
           real(qp),parameter :: hndrd = 100.0_qp
           real(qp),parameter :: meigth = -0.125_qp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,iterdivn,j,ll,lll,m,maxitdivn,nm1,nm12, &
                     nm13,oldll,oldm
           real(qp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('QBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_qlasq1(n,d,e,work,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_qlamch('EPSILON')
           unfl = la_qlamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_qlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 work(i) = cs
                 work(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_qlasr('R','V','F',nru,n,work(1),work(n),u,ldu)

              if (ncc > 0) call la_qlasr('L','V','F',n,ncc,work(1),work(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=qp))
              thresh = max(tol*sminoa,maxitr*(n*(n*unfl)))
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*(n*(n*unfl)))
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxitdivn = maxitr*n
           iterdivn = 0
           iter = -1
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter >= n) then
              iter = iter - n
              iterdivn = iterdivn + 1
              if (iterdivn >= maxitdivn) go to 200
           end if
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_qlasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_qrot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr,sinr &
                        )
              if (nru > 0) call la_qrot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)
              if (ncc > 0) call la_qrot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_qlas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_qlas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_qlartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_qlartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll + 1) = cs
                    work(i - ll + 1 + nm1) = sn
                    work(i - ll + 1 + nm12) = oldcs
                    work(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_qlasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_qlasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_qlasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_qlartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_qlartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    work(i - ll) = cs
                    work(i - ll + nm1) = -sn
                    work(i - ll + nm12) = oldcs
                    work(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_qlasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_qlasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_qlasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_qlartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_qlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    work(i - ll + 1) = cosr
                    work(i - ll + 1 + nm1) = sinr
                    work(i - ll + 1 + nm12) = cosl
                    work(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_qlasr('L','V','F',m - ll + 1,ncvt,work(1),work(n), &
                           vt(ll,1),ldvt)
                 if (nru > 0) call la_qlasr('R','V','F',nru,m - ll + 1,work(nm12 + 1),work( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_qlasr('L','V','F',m - ll + 1,ncc,work(nm12 + 1),work( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_qlartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_qlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    work(i - ll) = cosr
                    work(i - ll + nm1) = -sinr
                    work(i - ll + nm12) = cosl
                    work(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_qlasr('L','V','B',m - ll + 1,ncvt,work(nm12 + 1),work( &
                           nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_qlasr('R','V','B',nru,m - ll + 1,work(1),work(n),u( &
                            1,ll),ldu)
                 if (ncc > 0) call la_qlasr('L','V','B',m - ll + 1,ncc,work(1),work(n),c( &
                            ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_qscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_qswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_qswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_qswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_qbdsqr
#endif

     !> SLASQ1: computes the singular values of a real N-by-N bidiagonal
     !> matrix with diagonal D and off-diagonal E. The singular values
     !> are computed to high relative accuracy, in the absence of
     !> denormalization, underflow and overflow. The algorithm was first
     !> presented in
     !> "Accurate singular values and differential qd algorithms" by K. V.
     !> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
     !> 1994,
     !> and the present implementation is described in "An implementation of
     !> the dqds Algorithm (Positive Case)", LAPACK Working Note.

     pure subroutine la_slasq1(n,d,e,work,info)
        use la_constants_sp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iinfo
           real(sp) :: eps,scale,safmin,sigmn,sigmx
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
              call la_xerbla('SLASQ1',-info)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              d(1) = abs(d(1))
              return
           else if (n == 2) then
              call la_slas2(d(1),e(1),d(2),sigmn,sigmx)
              d(1) = sigmx
              d(2) = sigmn
              return
           end if
           ! estimate the largest singular value.
           sigmx = zero
           do i = 1,n - 1
              d(i) = abs(d(i))
              sigmx = max(sigmx,abs(e(i)))
           end do
           d(n) = abs(d(n))
           ! early return if sigmx is zero (matrix is already diagonal).
           if (sigmx == zero) then
              call la_slasrt('D',n,d,iinfo)
              return
           end if
           do i = 1,n
              sigmx = max(sigmx,d(i))
           end do
           ! copy d and e into work (in the z format) and scale (squaring the
           ! input data makes scaling by a power of the radix pointless).
           eps = la_slamch('PRECISION')
           safmin = la_slamch('SAFE MINIMUM')
           scale = sqrt(eps/safmin)
           call la_scopy(n,d,1,work(1),2)
           call la_scopy(n - 1,e,1,work(2),2)
           call la_slascl('G',0,0,sigmx,scale,2*n - 1,1,work,2*n - 1,iinfo)
           ! compute the q's and e's.
           do i = 1,2*n - 1
              work(i) = work(i)**2
           end do
           work(2*n) = zero
           call la_slasq2(n,work,info)
           if (info == 0) then
              do i = 1,n
                 d(i) = sqrt(work(i))
              end do
              call la_slascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
           else if (info == 2) then
           ! maximum number of iterations exceeded.  move data from work
           ! into d and e so the calling subroutine can try to finish
              do i = 1,n
                 d(i) = sqrt(work(2*i - 1))
                 e(i) = sqrt(work(2*i))
              end do
              call la_slascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
              call la_slascl('G',0,0,scale,sigmx,n,1,e,n,iinfo)
           end if
           return
     end subroutine la_slasq1
     !> DLASQ1: computes the singular values of a real N-by-N bidiagonal
     !> matrix with diagonal D and off-diagonal E. The singular values
     !> are computed to high relative accuracy, in the absence of
     !> denormalization, underflow and overflow. The algorithm was first
     !> presented in
     !> "Accurate singular values and differential qd algorithms" by K. V.
     !> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
     !> 1994,
     !> and the present implementation is described in "An implementation of
     !> the dqds Algorithm (Positive Case)", LAPACK Working Note.

     pure subroutine la_dlasq1(n,d,e,work,info)
        use la_constants_dp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iinfo
           real(dp) :: eps,scale,safmin,sigmn,sigmx
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
              call la_xerbla('DLASQ1',-info)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              d(1) = abs(d(1))
              return
           else if (n == 2) then
              call la_dlas2(d(1),e(1),d(2),sigmn,sigmx)
              d(1) = sigmx
              d(2) = sigmn
              return
           end if
           ! estimate the largest singular value.
           sigmx = zero
           do i = 1,n - 1
              d(i) = abs(d(i))
              sigmx = max(sigmx,abs(e(i)))
           end do
           d(n) = abs(d(n))
           ! early return if sigmx is zero (matrix is already diagonal).
           if (sigmx == zero) then
              call la_dlasrt('D',n,d,iinfo)
              return
           end if
           do i = 1,n
              sigmx = max(sigmx,d(i))
           end do
           ! copy d and e into work (in the z format) and scale (squaring the
           ! input data makes scaling by a power of the radix pointless).
           eps = la_dlamch('PRECISION')
           safmin = la_dlamch('SAFE MINIMUM')
           scale = sqrt(eps/safmin)
           call la_dcopy(n,d,1,work(1),2)
           call la_dcopy(n - 1,e,1,work(2),2)
           call la_dlascl('G',0,0,sigmx,scale,2*n - 1,1,work,2*n - 1,iinfo)
           ! compute the q's and e's.
           do i = 1,2*n - 1
              work(i) = work(i)**2
           end do
           work(2*n) = zero
           call la_dlasq2(n,work,info)
           if (info == 0) then
              do i = 1,n
                 d(i) = sqrt(work(i))
              end do
              call la_dlascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
           else if (info == 2) then
           ! maximum number of iterations exceeded.  move data from work
           ! into d and e so the calling subroutine can try to finish
              do i = 1,n
                 d(i) = sqrt(work(2*i - 1))
                 e(i) = sqrt(work(2*i))
              end do
              call la_dlascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
              call la_dlascl('G',0,0,scale,sigmx,n,1,e,n,iinfo)
           end if
           return
     end subroutine la_dlasq1
#ifdef LA_WITH_XDP
     !> XLASQ1: computes the singular values of a real N-by-N bidiagonal
     !> matrix with diagonal D and off-diagonal E. The singular values
     !> are computed to high relative accuracy, in the absence of
     !> denormalization, underflow and overflow. The algorithm was first
     !> presented in
     !> "Accurate singular values and differential qd algorithms" by K. V.
     !> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
     !> 1994,
     !> and the present implementation is described in "An implementation of
     !> the dqds Algorithm (Positive Case)", LAPACK Working Note.

     pure subroutine la_xlasq1(n,d,e,work,info)
        use la_constants_xdp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(inout) :: d(*),e(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iinfo
           real(xdp) :: eps,scale,safmin,sigmn,sigmx
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
              call la_xerbla('XLASQ1',-info)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              d(1) = abs(d(1))
              return
           else if (n == 2) then
              call la_xlas2(d(1),e(1),d(2),sigmn,sigmx)
              d(1) = sigmx
              d(2) = sigmn
              return
           end if
           ! estimate the largest singular value.
           sigmx = zero
           do i = 1,n - 1
              d(i) = abs(d(i))
              sigmx = max(sigmx,abs(e(i)))
           end do
           d(n) = abs(d(n))
           ! early return if sigmx is zero (matrix is already diagonal).
           if (sigmx == zero) then
              call la_xlasrt('D',n,d,iinfo)
              return
           end if
           do i = 1,n
              sigmx = max(sigmx,d(i))
           end do
           ! copy d and e into work (in the z format) and scale (squaring the
           ! input data makes scaling by a power of the radix pointless).
           eps = la_xlamch('PRECISION')
           safmin = la_xlamch('SAFE MINIMUM')
           scale = sqrt(eps/safmin)
           call la_xcopy(n,d,1,work(1),2)
           call la_xcopy(n - 1,e,1,work(2),2)
           call la_xlascl('G',0,0,sigmx,scale,2*n - 1,1,work,2*n - 1,iinfo)
           ! compute the q's and e's.
           do i = 1,2*n - 1
              work(i) = work(i)**2
           end do
           work(2*n) = zero
           call la_xlasq2(n,work,info)
           if (info == 0) then
              do i = 1,n
                 d(i) = sqrt(work(i))
              end do
              call la_xlascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
           else if (info == 2) then
           ! maximum number of iterations exceeded.  move data from work
           ! into d and e so the calling subroutine can try to finish
              do i = 1,n
                 d(i) = sqrt(work(2*i - 1))
                 e(i) = sqrt(work(2*i))
              end do
              call la_xlascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
              call la_xlascl('G',0,0,scale,sigmx,n,1,e,n,iinfo)
           end if
           return
     end subroutine la_xlasq1
#endif
#ifdef LA_WITH_QP
     !> QLASQ1: computes the singular values of a real N-by-N bidiagonal
     !> matrix with diagonal D and off-diagonal E. The singular values
     !> are computed to high relative accuracy, in the absence of
     !> denormalization, underflow and overflow. The algorithm was first
     !> presented in
     !> "Accurate singular values and differential qd algorithms" by K. V.
     !> Fernando and B. N. Parlett, Numer. Math., Vol-67, No. 2, pp. 191-230,
     !> 1994,
     !> and the present implementation is described in "An implementation of
     !> the dqds Algorithm (Positive Case)", LAPACK Working Note.

     pure subroutine la_qlasq1(n,d,e,work,info)
        use la_constants_qp,only:zero
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i,iinfo
           real(qp) :: eps,scale,safmin,sigmn,sigmx
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           info = 0
           if (n < 0) then
              info = -1
              call la_xerbla('QLASQ1',-info)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              d(1) = abs(d(1))
              return
           else if (n == 2) then
              call la_qlas2(d(1),e(1),d(2),sigmn,sigmx)
              d(1) = sigmx
              d(2) = sigmn
              return
           end if
           ! estimate the largest singular value.
           sigmx = zero
           do i = 1,n - 1
              d(i) = abs(d(i))
              sigmx = max(sigmx,abs(e(i)))
           end do
           d(n) = abs(d(n))
           ! early return if sigmx is zero (matrix is already diagonal).
           if (sigmx == zero) then
              call la_qlasrt('D',n,d,iinfo)
              return
           end if
           do i = 1,n
              sigmx = max(sigmx,d(i))
           end do
           ! copy d and e into work (in the z format) and scale (squaring the
           ! input data makes scaling by a power of the radix pointless).
           eps = la_qlamch('PRECISION')
           safmin = la_qlamch('SAFE MINIMUM')
           scale = sqrt(eps/safmin)
           call la_qcopy(n,d,1,work(1),2)
           call la_qcopy(n - 1,e,1,work(2),2)
           call la_qlascl('G',0,0,sigmx,scale,2*n - 1,1,work,2*n - 1,iinfo)
           ! compute the q's and e's.
           do i = 1,2*n - 1
              work(i) = work(i)**2
           end do
           work(2*n) = zero
           call la_qlasq2(n,work,info)
           if (info == 0) then
              do i = 1,n
                 d(i) = sqrt(work(i))
              end do
              call la_qlascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
           else if (info == 2) then
           ! maximum number of iterations exceeded.  move data from work
           ! into d and e so the calling subroutine can try to finish
              do i = 1,n
                 d(i) = sqrt(work(2*i - 1))
                 e(i) = sqrt(work(2*i))
              end do
              call la_qlascl('G',0,0,scale,sigmx,n,1,d,n,iinfo)
              call la_qlascl('G',0,0,scale,sigmx,n,1,e,n,iinfo)
           end if
           return
     end subroutine la_qlasq1
#endif

     !> SLASQ2: computes all the eigenvalues of the symmetric positive
     !> definite tridiagonal matrix associated with the qd array Z to high
     !> relative accuracy are computed to high relative accuracy, in the
     !> absence of denormalization, underflow and overflow.
     !> To see the relation of Z to the tridiagonal matrix, let L be a
     !> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
     !> let U be an upper bidiagonal matrix with 1's above and diagonal
     !> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
     !> symmetric tridiagonal to which it is similar.
     !> Note : SLASQ2 defines a logical variable, IEEE, which is true
     !> on machines which follow ieee-754 floating-point standard in their
     !> handling of infinities and NaNs, and false otherwise. This variable
     !> is passed to SLASQ3.

     pure subroutine la_slasq2(n,z,info)
        use la_constants_sp,only:zero,half,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: cbias = 1.50_sp
           real(sp),parameter :: hundrd = 100.0_sp

           ! Local Scalars
           logical(lk) :: ieee
           integer(ilp) :: i0,i4,iinfo,ipn4,iter,iwhila,iwhilb,k,kmin,n0,nbig,ndiv, &
                     nfail,pp,splt,ttype,i1,n1
           real(sp) :: d,dee,deemin,desig,dmin,dmin1,dmin2,dn,dn1,dn2,e,emax,emin, &
           eps,g,oldemn,qmax,qmin,s,safmin,sigma,t,tau,temp,tol,tol2,trace,zmax, &
                     tempe,tempq
           ! Intrinsic Functions
           intrinsic :: abs,max,min,real,sqrt
           ! Executable Statements
           ! test the input arguments.
           ! (in case la_slasq2 is not called by la_slasq1)
           info = 0
           eps = la_slamch('PRECISION')
           safmin = la_slamch('SAFE MINIMUM')
           tol = eps*hundrd
           tol2 = tol**2
           if (n < 0) then
              info = -1
              call la_xerbla('SLASQ2',1)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              ! 1-by-1 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('SLASQ2',2)
              end if
              return
           else if (n == 2) then
              ! 2-by-2 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('SLASQ2',2)
                 return
              else if (z(2) < zero) then
                 info = -202
                 call la_xerbla('SLASQ2',2)
                 return
              else if (z(3) < zero) then
                info = -203
                call la_xerbla('SLASQ2',2)
                return
              else if (z(3) > z(1)) then
                 d = z(3)
                 z(3) = z(1)
                 z(1) = d
              end if
              z(5) = z(1) + z(2) + z(3)
              if (z(2) > z(3)*tol2) then
                 t = half*((z(1) - z(3)) + z(2))
                 s = z(3)*(z(2)/t)
                 if (s <= t) then
                    s = z(3)*(z(2)/(t*(one + sqrt(one + s/t))))
                 else
                    s = z(3)*(z(2)/(t + sqrt(t)*sqrt(t + s)))
                 end if
                 t = z(1) + (s + z(2))
                 z(3) = z(3)*(z(1)/t)
                 z(1) = t
              end if
              z(2) = z(3)
              z(6) = z(2) + z(1)
              return
           end if
           ! check for negative data and compute sums of q's and e's.
           z(2*n) = zero
           emin = z(2)
           qmax = zero
           zmax = zero
           d = zero
           e = zero
           do k = 1,2*(n - 1),2
              if (z(k) < zero) then
                 info = -(200 + k)
                 call la_xerbla('SLASQ2',2)
                 return
              else if (z(k + 1) < zero) then
                 info = -(200 + k + 1)
                 call la_xerbla('SLASQ2',2)
                 return
              end if
              d = d + z(k)
              e = e + z(k + 1)
              qmax = max(qmax,z(k))
              emin = min(emin,z(k + 1))
              zmax = max(qmax,zmax,z(k + 1))
           end do
           if (z(2*n - 1) < zero) then
              info = -(200 + 2*n - 1)
              call la_xerbla('SLASQ2',2)
              return
           end if
           d = d + z(2*n - 1)
           qmax = max(qmax,z(2*n - 1))
           zmax = max(qmax,zmax)
           ! check for diagonality.
           if (e == zero) then
              do k = 2,n
                 z(k) = z(2*k - 1)
              end do
              call la_slasrt('D',n,z,iinfo)
              z(2*n - 1) = d
              return
           end if
           trace = d + e
           ! check for zero data.
           if (trace == zero) then
              z(2*n - 1) = zero
              return
           end if
           ! check whether the machine is ieee conformable.
           ! ieee = ( la_ilaenv( 10, 'slasq2', 'n', 1, 2, 3, 4 )==1 )
           ! [11/15/2008] the case ieee=.true. has a problem in single precision with
           ! some the test matrices of type 16. the double precision code is fine.
           ieee = .false.
           ! rearrange data for locality: z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
           do k = 2*n,2,-2
              z(2*k) = zero
              z(2*k - 1) = z(k)
              z(2*k - 2) = zero
              z(2*k - 3) = z(k - 1)
           end do
           i0 = 1
           n0 = n
           ! reverse the qd-array, if warranted.
           if (cbias*z(4*i0 - 3) < z(4*n0 - 3)) then
              ipn4 = 4*(i0 + n0)
              do i4 = 4*i0,2*(i0 + n0 - 1),4
                 temp = z(i4 - 3)
                 z(i4 - 3) = z(ipn4 - i4 - 3)
                 z(ipn4 - i4 - 3) = temp
                 temp = z(i4 - 1)
                 z(i4 - 1) = z(ipn4 - i4 - 5)
                 z(ipn4 - i4 - 5) = temp
              end do
           end if
           ! initial split checking via dqd and li's test.
           pp = 0
           loop_80: do k = 1,2
              d = z(4*n0 + pp - 3)
              do i4 = 4*(n0 - 1) + pp,4*i0 + pp,-4
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    d = z(i4 - 3)
                 else
                    d = z(i4 - 3)*(d/(d + z(i4 - 1)))
                 end if
              end do
              ! dqd maps z to zz plus li's test.
              emin = z(4*i0 + pp + 1)
              d = z(4*i0 + pp - 3)
              do i4 = 4*i0 + pp,4*(n0 - 1) + pp,4
                 z(i4 - 2*pp - 2) = d + z(i4 - 1)
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    z(i4 - 2*pp - 2) = d
                    z(i4 - 2*pp) = zero
                    d = z(i4 + 1)
                 else if (safmin*z(i4 + 1) < z(i4 - 2*pp - 2) .and. safmin*z(i4 - 2*pp - 2) < z(i4 + 1)) &
                           then
                    temp = z(i4 + 1)/z(i4 - 2*pp - 2)
                    z(i4 - 2*pp) = z(i4 - 1)*temp
                    d = d*temp
                 else
                    z(i4 - 2*pp) = z(i4 + 1)*(z(i4 - 1)/z(i4 - 2*pp - 2))
                    d = z(i4 + 1)*(d/z(i4 - 2*pp - 2))
                 end if
                 emin = min(emin,z(i4 - 2*pp))
              end do
              z(4*n0 - pp - 2) = d
              ! now find qmax.
              qmax = z(4*i0 - pp - 2)
              do i4 = 4*i0 - pp + 2,4*n0 - pp - 2,4
                 qmax = max(qmax,z(i4))
              end do
              ! prepare for the next iteration on k.
              pp = 1 - pp
           end do loop_80
           ! initialise variables to pass to la_slasq3.
           ttype = 0
           dmin1 = zero
           dmin2 = zero
           dn = zero
           dn1 = zero
           dn2 = zero
           g = zero
           tau = zero
           iter = 2
           nfail = 0
           ndiv = 2*(n0 - i0)
           loop_160: do iwhila = 1,n + 1
              if (n0 < 1) go to 170
              ! while array unfinished do
              ! e(n0) holds the value of sigma when submatrix in i0:n0
              ! splits from the rest of the array, but is negated.
              desig = zero
              if (n0 == n) then
                 sigma = zero
              else
                 sigma = -z(4*n0 - 1)
              end if
              if (sigma < zero) then
                 info = 1
                 return
              end if
              ! find last unreduced submatrix's top index i0, find qmax and
              ! emin. find gershgorin-type bound if q's much greater than e's.
              emax = zero
              if (n0 > i0) then
                 emin = abs(z(4*n0 - 5))
              else
                 emin = zero
              end if
              qmin = z(4*n0 - 3)
              qmax = qmin
              do i4 = 4*n0,8,-4
                 if (z(i4 - 5) <= zero) go to 100
                 if (qmin >= four*emax) then
                    qmin = min(qmin,z(i4 - 3))
                    emax = max(emax,z(i4 - 5))
                 end if
                 qmax = max(qmax,z(i4 - 7) + z(i4 - 5))
                 emin = min(emin,z(i4 - 5))
              end do
              i4 = 4
              100 continue
              i0 = i4/4
              pp = 0
              if (n0 - i0 > 1) then
                 dee = z(4*i0 - 3)
                 deemin = dee
                 kmin = i0
                 do i4 = 4*i0 + 1,4*n0 - 3,4
                    dee = z(i4)*(dee/(dee + z(i4 - 2)))
                    if (dee <= deemin) then
                       deemin = dee
                       kmin = (i4 + 3)/4
                    end if
                 end do
                 if ((kmin - i0)*2 < n0 - kmin .and. deemin <= half*z(4*n0 - 3)) then
                    ipn4 = 4*(i0 + n0)
                    pp = 2
                    do i4 = 4*i0,2*(i0 + n0 - 1),4
                       temp = z(i4 - 3)
                       z(i4 - 3) = z(ipn4 - i4 - 3)
                       z(ipn4 - i4 - 3) = temp
                       temp = z(i4 - 2)
                       z(i4 - 2) = z(ipn4 - i4 - 2)
                       z(ipn4 - i4 - 2) = temp
                       temp = z(i4 - 1)
                       z(i4 - 1) = z(ipn4 - i4 - 5)
                       z(ipn4 - i4 - 5) = temp
                       temp = z(i4)
                       z(i4) = z(ipn4 - i4 - 4)
                       z(ipn4 - i4 - 4) = temp
                    end do
                 end if
              end if
              ! put -(initial shift) into dmin.
              dmin = -max(zero,qmin - two*sqrt(qmin)*sqrt(emax))
              ! now i0:n0 is unreduced.
              ! pp = 0 for ping, pp = 1 for pong.
              ! pp = 2 indicates that flipping was applied to the z array and
                     ! and that the tests for deflation upon entry in la_slasq3
                     ! should not be performed.
              nbig = 100*(n0 - i0 + 1)
              loop_140: do iwhilb = 1,nbig
                 if (i0 > n0) go to 150
                 ! while submatrix unfinished take a good dqds step.
                 call la_slasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
                           ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
                 pp = 1 - pp
                 ! when emin is very small check for splits.
                 if (pp == 0 .and. n0 - i0 >= 3) then
                    if (z(4*n0) <= tol2*qmax .or. z(4*n0 - 1) <= tol2*sigma) then
                       splt = i0 - 1
                       qmax = z(4*i0 - 3)
                       emin = z(4*i0 - 1)
                       oldemn = z(4*i0)
                       do i4 = 4*i0,4*(n0 - 3),4
                          if (z(i4) <= tol2*z(i4 - 3) .or. z(i4 - 1) <= tol2*sigma) then
                             z(i4 - 1) = -sigma
                             splt = i4/4
                             qmax = zero
                             emin = z(i4 + 3)
                             oldemn = z(i4 + 4)
                          else
                             qmax = max(qmax,z(i4 + 1))
                             emin = min(emin,z(i4 - 1))
                             oldemn = min(oldemn,z(i4))
                          end if
                       end do
                       z(4*n0 - 1) = emin
                       z(4*n0) = oldemn
                       i0 = splt + 1
                    end if
                 end if
              end do loop_140
              info = 2
              ! maximum number of iterations exceeded, restore the shift
              ! sigma and place the new d's and e's in a qd array.
              ! this might need to be done for several blocks
              i1 = i0
              n1 = n0
              145 continue
              tempq = z(4*i0 - 3)
              z(4*i0 - 3) = z(4*i0 - 3) + sigma
              do k = i0 + 1,n0
                 tempe = z(4*k - 5)
                 z(4*k - 5) = z(4*k - 5)*(tempq/z(4*k - 7))
                 tempq = z(4*k - 3)
                 z(4*k - 3) = z(4*k - 3) + sigma + tempe - z(4*k - 5)
              end do
              ! prepare to do this on the previous block if there is one
              if (i1 > 1) then
                 n1 = i1 - 1
                 do while ((i1 >= 2) .and. (z(4*i1 - 5) >= zero))
                    i1 = i1 - 1
                 end do
                 if (i1 >= 1) then
                 sigma = -z(4*n1 - 1)
                 go to 145
                 end if
              end if
              do k = 1,n
                 z(2*k - 1) = z(4*k - 3)
              ! only the block 1..n0 is unfinished.  the rest of the e's
              ! must be essentially zero, although sometimes other data
              ! has been stored in them.
                 if (k < n0) then
                    z(2*k) = z(4*k - 1)
                 else
                    z(2*k) = 0
                 end if
              end do
              return
              ! end iwhilb
              150 continue
           end do loop_160
           info = 3
           return
           ! end iwhila
           170 continue
           ! move q's to the front.
           do k = 2,n
              z(k) = z(4*k - 3)
           end do
           ! sort and compute sum of eigenvalues.
           call la_slasrt('D',n,z,iinfo)
           e = zero
           do k = n,1,-1
              e = e + z(k)
           end do
           ! store trace, sum(eigenvalues) and information on performance.
           z(2*n + 1) = trace
           z(2*n + 2) = e
           z(2*n + 3) = real(iter,KIND=sp)
           z(2*n + 4) = real(ndiv,KIND=sp)/real(n**2,KIND=sp)
           z(2*n + 5) = hundrd*nfail/real(iter,KIND=sp)
           return
     end subroutine la_slasq2
     !> DLASQ2: computes all the eigenvalues of the symmetric positive
     !> definite tridiagonal matrix associated with the qd array Z to high
     !> relative accuracy are computed to high relative accuracy, in the
     !> absence of denormalization, underflow and overflow.
     !> To see the relation of Z to the tridiagonal matrix, let L be a
     !> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
     !> let U be an upper bidiagonal matrix with 1's above and diagonal
     !> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
     !> symmetric tridiagonal to which it is similar.
     !> Note : DLASQ2 defines a logical variable, IEEE, which is true
     !> on machines which follow ieee-754 floating-point standard in their
     !> handling of infinities and NaNs, and false otherwise. This variable
     !> is passed to DLASQ3.

     pure subroutine la_dlasq2(n,z,info)
        use la_constants_dp,only:zero,half,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: cbias = 1.50_dp
           real(dp),parameter :: hundrd = 100.0_dp

           ! Local Scalars
           logical(lk) :: ieee
           integer(ilp) :: i0,i1,i4,iinfo,ipn4,iter,iwhila,iwhilb,k,kmin,n0,n1,nbig, &
                     ndiv,nfail,pp,splt,ttype
           real(dp) :: d,dee,deemin,desig,dmin,dmin1,dmin2,dn,dn1,dn2,e,emax,emin, &
           eps,g,oldemn,qmax,qmin,s,safmin,sigma,t,tau,temp,tol,tol2,trace,zmax, &
                     tempe,tempq
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sqrt
           ! Executable Statements
           ! test the input arguments.
           ! (in case la_dlasq2 is not called by la_dlasq1)
           info = 0
           eps = la_dlamch('PRECISION')
           safmin = la_dlamch('SAFE MINIMUM')
           tol = eps*hundrd
           tol2 = tol**2
           if (n < 0) then
              info = -1
              call la_xerbla('DLASQ2',1)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              ! 1-by-1 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('DLASQ2',2)
              end if
              return
           else if (n == 2) then
              ! 2-by-2 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('DLASQ2',2)
                 return
              else if (z(2) < zero) then
                 info = -202
                 call la_xerbla('DLASQ2',2)
                 return
              else if (z(3) < zero) then
                info = -203
                call la_xerbla('DLASQ2',2)
                return
              else if (z(3) > z(1)) then
                 d = z(3)
                 z(3) = z(1)
                 z(1) = d
              end if
              z(5) = z(1) + z(2) + z(3)
              if (z(2) > z(3)*tol2) then
                 t = half*((z(1) - z(3)) + z(2))
                 s = z(3)*(z(2)/t)
                 if (s <= t) then
                    s = z(3)*(z(2)/(t*(one + sqrt(one + s/t))))
                 else
                    s = z(3)*(z(2)/(t + sqrt(t)*sqrt(t + s)))
                 end if
                 t = z(1) + (s + z(2))
                 z(3) = z(3)*(z(1)/t)
                 z(1) = t
              end if
              z(2) = z(3)
              z(6) = z(2) + z(1)
              return
           end if
           ! check for negative data and compute sums of q's and e's.
           z(2*n) = zero
           emin = z(2)
           qmax = zero
           zmax = zero
           d = zero
           e = zero
           do k = 1,2*(n - 1),2
              if (z(k) < zero) then
                 info = -(200 + k)
                 call la_xerbla('DLASQ2',2)
                 return
              else if (z(k + 1) < zero) then
                 info = -(200 + k + 1)
                 call la_xerbla('DLASQ2',2)
                 return
              end if
              d = d + z(k)
              e = e + z(k + 1)
              qmax = max(qmax,z(k))
              emin = min(emin,z(k + 1))
              zmax = max(qmax,zmax,z(k + 1))
           end do
           if (z(2*n - 1) < zero) then
              info = -(200 + 2*n - 1)
              call la_xerbla('DLASQ2',2)
              return
           end if
           d = d + z(2*n - 1)
           qmax = max(qmax,z(2*n - 1))
           zmax = max(qmax,zmax)
           ! check for diagonality.
           if (e == zero) then
              do k = 2,n
                 z(k) = z(2*k - 1)
              end do
              call la_dlasrt('D',n,z,iinfo)
              z(2*n - 1) = d
              return
           end if
           trace = d + e
           ! check for zero data.
           if (trace == zero) then
              z(2*n - 1) = zero
              return
           end if
           ! check whether the machine is ieee conformable.
           ieee = (la_ilaenv(10,'DLASQ2','N',1,2,3,4) == 1)
           ! rearrange data for locality: z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
           do k = 2*n,2,-2
              z(2*k) = zero
              z(2*k - 1) = z(k)
              z(2*k - 2) = zero
              z(2*k - 3) = z(k - 1)
           end do
           i0 = 1
           n0 = n
           ! reverse the qd-array, if warranted.
           if (cbias*z(4*i0 - 3) < z(4*n0 - 3)) then
              ipn4 = 4*(i0 + n0)
              do i4 = 4*i0,2*(i0 + n0 - 1),4
                 temp = z(i4 - 3)
                 z(i4 - 3) = z(ipn4 - i4 - 3)
                 z(ipn4 - i4 - 3) = temp
                 temp = z(i4 - 1)
                 z(i4 - 1) = z(ipn4 - i4 - 5)
                 z(ipn4 - i4 - 5) = temp
              end do
           end if
           ! initial split checking via dqd and li's test.
           pp = 0
           loop_80: do k = 1,2
              d = z(4*n0 + pp - 3)
              do i4 = 4*(n0 - 1) + pp,4*i0 + pp,-4
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    d = z(i4 - 3)
                 else
                    d = z(i4 - 3)*(d/(d + z(i4 - 1)))
                 end if
              end do
              ! dqd maps z to zz plus li's test.
              emin = z(4*i0 + pp + 1)
              d = z(4*i0 + pp - 3)
              do i4 = 4*i0 + pp,4*(n0 - 1) + pp,4
                 z(i4 - 2*pp - 2) = d + z(i4 - 1)
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    z(i4 - 2*pp - 2) = d
                    z(i4 - 2*pp) = zero
                    d = z(i4 + 1)
                 else if (safmin*z(i4 + 1) < z(i4 - 2*pp - 2) .and. safmin*z(i4 - 2*pp - 2) < z(i4 + 1)) &
                           then
                    temp = z(i4 + 1)/z(i4 - 2*pp - 2)
                    z(i4 - 2*pp) = z(i4 - 1)*temp
                    d = d*temp
                 else
                    z(i4 - 2*pp) = z(i4 + 1)*(z(i4 - 1)/z(i4 - 2*pp - 2))
                    d = z(i4 + 1)*(d/z(i4 - 2*pp - 2))
                 end if
                 emin = min(emin,z(i4 - 2*pp))
              end do
              z(4*n0 - pp - 2) = d
              ! now find qmax.
              qmax = z(4*i0 - pp - 2)
              do i4 = 4*i0 - pp + 2,4*n0 - pp - 2,4
                 qmax = max(qmax,z(i4))
              end do
              ! prepare for the next iteration on k.
              pp = 1 - pp
           end do loop_80
           ! initialise variables to pass to la_dlasq3.
           ttype = 0
           dmin1 = zero
           dmin2 = zero
           dn = zero
           dn1 = zero
           dn2 = zero
           g = zero
           tau = zero
           iter = 2
           nfail = 0
           ndiv = 2*(n0 - i0)
           loop_160: do iwhila = 1,n + 1
              if (n0 < 1) go to 170
              ! while array unfinished do
              ! e(n0) holds the value of sigma when submatrix in i0:n0
              ! splits from the rest of the array, but is negated.
              desig = zero
              if (n0 == n) then
                 sigma = zero
              else
                 sigma = -z(4*n0 - 1)
              end if
              if (sigma < zero) then
                 info = 1
                 return
              end if
              ! find last unreduced submatrix's top index i0, find qmax and
              ! emin. find gershgorin-type bound if q's much greater than e's.
              emax = zero
              if (n0 > i0) then
                 emin = abs(z(4*n0 - 5))
              else
                 emin = zero
              end if
              qmin = z(4*n0 - 3)
              qmax = qmin
              do i4 = 4*n0,8,-4
                 if (z(i4 - 5) <= zero) go to 100
                 if (qmin >= four*emax) then
                    qmin = min(qmin,z(i4 - 3))
                    emax = max(emax,z(i4 - 5))
                 end if
                 qmax = max(qmax,z(i4 - 7) + z(i4 - 5))
                 emin = min(emin,z(i4 - 5))
              end do
              i4 = 4
              100 continue
              i0 = i4/4
              pp = 0
              if (n0 - i0 > 1) then
                 dee = z(4*i0 - 3)
                 deemin = dee
                 kmin = i0
                 do i4 = 4*i0 + 1,4*n0 - 3,4
                    dee = z(i4)*(dee/(dee + z(i4 - 2)))
                    if (dee <= deemin) then
                       deemin = dee
                       kmin = (i4 + 3)/4
                    end if
                 end do
                 if ((kmin - i0)*2 < n0 - kmin .and. deemin <= half*z(4*n0 - 3)) then
                    ipn4 = 4*(i0 + n0)
                    pp = 2
                    do i4 = 4*i0,2*(i0 + n0 - 1),4
                       temp = z(i4 - 3)
                       z(i4 - 3) = z(ipn4 - i4 - 3)
                       z(ipn4 - i4 - 3) = temp
                       temp = z(i4 - 2)
                       z(i4 - 2) = z(ipn4 - i4 - 2)
                       z(ipn4 - i4 - 2) = temp
                       temp = z(i4 - 1)
                       z(i4 - 1) = z(ipn4 - i4 - 5)
                       z(ipn4 - i4 - 5) = temp
                       temp = z(i4)
                       z(i4) = z(ipn4 - i4 - 4)
                       z(ipn4 - i4 - 4) = temp
                    end do
                 end if
              end if
              ! put -(initial shift) into dmin.
              dmin = -max(zero,qmin - two*sqrt(qmin)*sqrt(emax))
              ! now i0:n0 is unreduced.
              ! pp = 0 for ping, pp = 1 for pong.
              ! pp = 2 indicates that flipping was applied to the z array and
                     ! and that the tests for deflation upon entry in la_dlasq3
                     ! should not be performed.
              nbig = 100*(n0 - i0 + 1)
              loop_140: do iwhilb = 1,nbig
                 if (i0 > n0) go to 150
                 ! while submatrix unfinished take a good dqds step.
                 call la_dlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
                           ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
                 pp = 1 - pp
                 ! when emin is very small check for splits.
                 if (pp == 0 .and. n0 - i0 >= 3) then
                    if (z(4*n0) <= tol2*qmax .or. z(4*n0 - 1) <= tol2*sigma) then
                       splt = i0 - 1
                       qmax = z(4*i0 - 3)
                       emin = z(4*i0 - 1)
                       oldemn = z(4*i0)
                       do i4 = 4*i0,4*(n0 - 3),4
                          if (z(i4) <= tol2*z(i4 - 3) .or. z(i4 - 1) <= tol2*sigma) then
                             z(i4 - 1) = -sigma
                             splt = i4/4
                             qmax = zero
                             emin = z(i4 + 3)
                             oldemn = z(i4 + 4)
                          else
                             qmax = max(qmax,z(i4 + 1))
                             emin = min(emin,z(i4 - 1))
                             oldemn = min(oldemn,z(i4))
                          end if
                       end do
                       z(4*n0 - 1) = emin
                       z(4*n0) = oldemn
                       i0 = splt + 1
                    end if
                 end if
              end do loop_140
              info = 2
              ! maximum number of iterations exceeded, restore the shift
              ! sigma and place the new d's and e's in a qd array.
              ! this might need to be done for several blocks
              i1 = i0
              n1 = n0
              145 continue
              tempq = z(4*i0 - 3)
              z(4*i0 - 3) = z(4*i0 - 3) + sigma
              do k = i0 + 1,n0
                 tempe = z(4*k - 5)
                 z(4*k - 5) = z(4*k - 5)*(tempq/z(4*k - 7))
                 tempq = z(4*k - 3)
                 z(4*k - 3) = z(4*k - 3) + sigma + tempe - z(4*k - 5)
              end do
              ! prepare to do this on the previous block if there is one
              if (i1 > 1) then
                 n1 = i1 - 1
                 do while ((i1 >= 2) .and. (z(4*i1 - 5) >= zero))
                    i1 = i1 - 1
                 end do
                 sigma = -z(4*n1 - 1)
                 go to 145
              end if
              do k = 1,n
                 z(2*k - 1) = z(4*k - 3)
              ! only the block 1..n0 is unfinished.  the rest of the e's
              ! must be essentially zero, although sometimes other data
              ! has been stored in them.
                 if (k < n0) then
                    z(2*k) = z(4*k - 1)
                 else
                    z(2*k) = 0
                 end if
              end do
              return
              ! end iwhilb
              150 continue
           end do loop_160
           info = 3
           return
           ! end iwhila
           170 continue
           ! move q's to the front.
           do k = 2,n
              z(k) = z(4*k - 3)
           end do
           ! sort and compute sum of eigenvalues.
           call la_dlasrt('D',n,z,iinfo)
           e = zero
           do k = n,1,-1
              e = e + z(k)
           end do
           ! store trace, sum(eigenvalues) and information on performance.
           z(2*n + 1) = trace
           z(2*n + 2) = e
           z(2*n + 3) = real(iter,KIND=dp)
           z(2*n + 4) = real(ndiv,KIND=dp)/real(n**2,KIND=dp)
           z(2*n + 5) = hundrd*nfail/real(iter,KIND=dp)
           return
     end subroutine la_dlasq2
#ifdef LA_WITH_XDP
     !> XLASQ2: computes all the eigenvalues of the symmetric positive
     !> definite tridiagonal matrix associated with the qd array Z to high
     !> relative accuracy are computed to high relative accuracy, in the
     !> absence of denormalization, underflow and overflow.
     !> To see the relation of Z to the tridiagonal matrix, let L be a
     !> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
     !> let U be an upper bidiagonal matrix with 1's above and diagonal
     !> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
     !> symmetric tridiagonal to which it is similar.
     !> Note : XLASQ2 defines a logical variable, IEEE, which is true
     !> on machines which follow ieee-754 floating-point standard in their
     !> handling of infinities and NaNs, and false otherwise. This variable
     !> is passed to XLASQ3.

     pure subroutine la_xlasq2(n,z,info)
        use la_constants_xdp,only:zero,half,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: cbias = 1.50_xdp
           real(xdp),parameter :: hundrd = 100.0_xdp

           ! Local Scalars
           logical(lk) :: ieee
           integer(ilp) :: i0,i1,i4,iinfo,ipn4,iter,iwhila,iwhilb,k,kmin,n0,n1,nbig, &
                     ndiv,nfail,pp,splt,ttype
           real(xdp) :: d,dee,deemin,desig,dmin,dmin1,dmin2,dn,dn1,dn2,e,emax,emin, &
           eps,g,oldemn,qmax,qmin,s,safmin,sigma,t,tau,temp,tol,tol2,trace,zmax, &
                     tempe,tempq
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sqrt
           ! Executable Statements
           ! test the input arguments.
           ! (in case la_xlasq2 is not called by la_xlasq1)
           info = 0
           eps = la_xlamch('PRECISION')
           safmin = la_xlamch('SAFE MINIMUM')
           tol = eps*hundrd
           tol2 = tol**2
           if (n < 0) then
              info = -1
              call la_xerbla('XLASQ2',1)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              ! 1-by-1 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('XLASQ2',2)
              end if
              return
           else if (n == 2) then
              ! 2-by-2 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('XLASQ2',2)
                 return
              else if (z(2) < zero) then
                 info = -202
                 call la_xerbla('XLASQ2',2)
                 return
              else if (z(3) < zero) then
                info = -203
                call la_xerbla('XLASQ2',2)
                return
              else if (z(3) > z(1)) then
                 d = z(3)
                 z(3) = z(1)
                 z(1) = d
              end if
              z(5) = z(1) + z(2) + z(3)
              if (z(2) > z(3)*tol2) then
                 t = half*((z(1) - z(3)) + z(2))
                 s = z(3)*(z(2)/t)
                 if (s <= t) then
                    s = z(3)*(z(2)/(t*(one + sqrt(one + s/t))))
                 else
                    s = z(3)*(z(2)/(t + sqrt(t)*sqrt(t + s)))
                 end if
                 t = z(1) + (s + z(2))
                 z(3) = z(3)*(z(1)/t)
                 z(1) = t
              end if
              z(2) = z(3)
              z(6) = z(2) + z(1)
              return
           end if
           ! check for negative data and compute sums of q's and e's.
           z(2*n) = zero
           emin = z(2)
           qmax = zero
           zmax = zero
           d = zero
           e = zero
           do k = 1,2*(n - 1),2
              if (z(k) < zero) then
                 info = -(200 + k)
                 call la_xerbla('XLASQ2',2)
                 return
              else if (z(k + 1) < zero) then
                 info = -(200 + k + 1)
                 call la_xerbla('XLASQ2',2)
                 return
              end if
              d = d + z(k)
              e = e + z(k + 1)
              qmax = max(qmax,z(k))
              emin = min(emin,z(k + 1))
              zmax = max(qmax,zmax,z(k + 1))
           end do
           if (z(2*n - 1) < zero) then
              info = -(200 + 2*n - 1)
              call la_xerbla('XLASQ2',2)
              return
           end if
           d = d + z(2*n - 1)
           qmax = max(qmax,z(2*n - 1))
           zmax = max(qmax,zmax)
           ! check for diagonality.
           if (e == zero) then
              do k = 2,n
                 z(k) = z(2*k - 1)
              end do
              call la_xlasrt('D',n,z,iinfo)
              z(2*n - 1) = d
              return
           end if
           trace = d + e
           ! check for zero data.
           if (trace == zero) then
              z(2*n - 1) = zero
              return
           end if
           ! check whether the machine is ieee conformable.
           ieee = (la_ilaenv(10,'XLASQ2','N',1,2,3,4) == 1)
           ! rearrange data for locality: z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
           do k = 2*n,2,-2
              z(2*k) = zero
              z(2*k - 1) = z(k)
              z(2*k - 2) = zero
              z(2*k - 3) = z(k - 1)
           end do
           i0 = 1
           n0 = n
           ! reverse the qd-array, if warranted.
           if (cbias*z(4*i0 - 3) < z(4*n0 - 3)) then
              ipn4 = 4*(i0 + n0)
              do i4 = 4*i0,2*(i0 + n0 - 1),4
                 temp = z(i4 - 3)
                 z(i4 - 3) = z(ipn4 - i4 - 3)
                 z(ipn4 - i4 - 3) = temp
                 temp = z(i4 - 1)
                 z(i4 - 1) = z(ipn4 - i4 - 5)
                 z(ipn4 - i4 - 5) = temp
              end do
           end if
           ! initial split checking via dqd and li's test.
           pp = 0
           loop_80: do k = 1,2
              d = z(4*n0 + pp - 3)
              do i4 = 4*(n0 - 1) + pp,4*i0 + pp,-4
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    d = z(i4 - 3)
                 else
                    d = z(i4 - 3)*(d/(d + z(i4 - 1)))
                 end if
              end do
              ! dqd maps z to zz plus li's test.
              emin = z(4*i0 + pp + 1)
              d = z(4*i0 + pp - 3)
              do i4 = 4*i0 + pp,4*(n0 - 1) + pp,4
                 z(i4 - 2*pp - 2) = d + z(i4 - 1)
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    z(i4 - 2*pp - 2) = d
                    z(i4 - 2*pp) = zero
                    d = z(i4 + 1)
                 else if (safmin*z(i4 + 1) < z(i4 - 2*pp - 2) .and. safmin*z(i4 - 2*pp - 2) < z(i4 + 1)) &
                           then
                    temp = z(i4 + 1)/z(i4 - 2*pp - 2)
                    z(i4 - 2*pp) = z(i4 - 1)*temp
                    d = d*temp
                 else
                    z(i4 - 2*pp) = z(i4 + 1)*(z(i4 - 1)/z(i4 - 2*pp - 2))
                    d = z(i4 + 1)*(d/z(i4 - 2*pp - 2))
                 end if
                 emin = min(emin,z(i4 - 2*pp))
              end do
              z(4*n0 - pp - 2) = d
              ! now find qmax.
              qmax = z(4*i0 - pp - 2)
              do i4 = 4*i0 - pp + 2,4*n0 - pp - 2,4
                 qmax = max(qmax,z(i4))
              end do
              ! prepare for the next iteration on k.
              pp = 1 - pp
           end do loop_80
           ! initialise variables to pass to la_xlasq3.
           ttype = 0
           dmin1 = zero
           dmin2 = zero
           dn = zero
           dn1 = zero
           dn2 = zero
           g = zero
           tau = zero
           iter = 2
           nfail = 0
           ndiv = 2*(n0 - i0)
           loop_160: do iwhila = 1,n + 1
              if (n0 < 1) go to 170
              ! while array unfinished do
              ! e(n0) holds the value of sigma when submatrix in i0:n0
              ! splits from the rest of the array, but is negated.
              desig = zero
              if (n0 == n) then
                 sigma = zero
              else
                 sigma = -z(4*n0 - 1)
              end if
              if (sigma < zero) then
                 info = 1
                 return
              end if
              ! find last unreduced submatrix's top index i0, find qmax and
              ! emin. find gershgorin-type bound if q's much greater than e's.
              emax = zero
              if (n0 > i0) then
                 emin = abs(z(4*n0 - 5))
              else
                 emin = zero
              end if
              qmin = z(4*n0 - 3)
              qmax = qmin
              do i4 = 4*n0,8,-4
                 if (z(i4 - 5) <= zero) go to 100
                 if (qmin >= four*emax) then
                    qmin = min(qmin,z(i4 - 3))
                    emax = max(emax,z(i4 - 5))
                 end if
                 qmax = max(qmax,z(i4 - 7) + z(i4 - 5))
                 emin = min(emin,z(i4 - 5))
              end do
              i4 = 4
              100 continue
              i0 = i4/4
              pp = 0
              if (n0 - i0 > 1) then
                 dee = z(4*i0 - 3)
                 deemin = dee
                 kmin = i0
                 do i4 = 4*i0 + 1,4*n0 - 3,4
                    dee = z(i4)*(dee/(dee + z(i4 - 2)))
                    if (dee <= deemin) then
                       deemin = dee
                       kmin = (i4 + 3)/4
                    end if
                 end do
                 if ((kmin - i0)*2 < n0 - kmin .and. deemin <= half*z(4*n0 - 3)) then
                    ipn4 = 4*(i0 + n0)
                    pp = 2
                    do i4 = 4*i0,2*(i0 + n0 - 1),4
                       temp = z(i4 - 3)
                       z(i4 - 3) = z(ipn4 - i4 - 3)
                       z(ipn4 - i4 - 3) = temp
                       temp = z(i4 - 2)
                       z(i4 - 2) = z(ipn4 - i4 - 2)
                       z(ipn4 - i4 - 2) = temp
                       temp = z(i4 - 1)
                       z(i4 - 1) = z(ipn4 - i4 - 5)
                       z(ipn4 - i4 - 5) = temp
                       temp = z(i4)
                       z(i4) = z(ipn4 - i4 - 4)
                       z(ipn4 - i4 - 4) = temp
                    end do
                 end if
              end if
              ! put -(initial shift) into dmin.
              dmin = -max(zero,qmin - two*sqrt(qmin)*sqrt(emax))
              ! now i0:n0 is unreduced.
              ! pp = 0 for ping, pp = 1 for pong.
              ! pp = 2 indicates that flipping was applied to the z array and
                     ! and that the tests for deflation upon entry in la_xlasq3
                     ! should not be performed.
              nbig = 100*(n0 - i0 + 1)
              loop_140: do iwhilb = 1,nbig
                 if (i0 > n0) go to 150
                 ! while submatrix unfinished take a good dqds step.
                 call la_xlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
                           ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
                 pp = 1 - pp
                 ! when emin is very small check for splits.
                 if (pp == 0 .and. n0 - i0 >= 3) then
                    if (z(4*n0) <= tol2*qmax .or. z(4*n0 - 1) <= tol2*sigma) then
                       splt = i0 - 1
                       qmax = z(4*i0 - 3)
                       emin = z(4*i0 - 1)
                       oldemn = z(4*i0)
                       do i4 = 4*i0,4*(n0 - 3),4
                          if (z(i4) <= tol2*z(i4 - 3) .or. z(i4 - 1) <= tol2*sigma) then
                             z(i4 - 1) = -sigma
                             splt = i4/4
                             qmax = zero
                             emin = z(i4 + 3)
                             oldemn = z(i4 + 4)
                          else
                             qmax = max(qmax,z(i4 + 1))
                             emin = min(emin,z(i4 - 1))
                             oldemn = min(oldemn,z(i4))
                          end if
                       end do
                       z(4*n0 - 1) = emin
                       z(4*n0) = oldemn
                       i0 = splt + 1
                    end if
                 end if
              end do loop_140
              info = 2
              ! maximum number of iterations exceeded, restore the shift
              ! sigma and place the new d's and e's in a qd array.
              ! this might need to be done for several blocks
              i1 = i0
              n1 = n0
              145 continue
              tempq = z(4*i0 - 3)
              z(4*i0 - 3) = z(4*i0 - 3) + sigma
              do k = i0 + 1,n0
                 tempe = z(4*k - 5)
                 z(4*k - 5) = z(4*k - 5)*(tempq/z(4*k - 7))
                 tempq = z(4*k - 3)
                 z(4*k - 3) = z(4*k - 3) + sigma + tempe - z(4*k - 5)
              end do
              ! prepare to do this on the previous block if there is one
              if (i1 > 1) then
                 n1 = i1 - 1
                 do while ((i1 >= 2) .and. (z(4*i1 - 5) >= zero))
                    i1 = i1 - 1
                 end do
                 sigma = -z(4*n1 - 1)
                 go to 145
              end if
              do k = 1,n
                 z(2*k - 1) = z(4*k - 3)
              ! only the block 1..n0 is unfinished.  the rest of the e's
              ! must be essentially zero, although sometimes other data
              ! has been stored in them.
                 if (k < n0) then
                    z(2*k) = z(4*k - 1)
                 else
                    z(2*k) = 0
                 end if
              end do
              return
              ! end iwhilb
              150 continue
           end do loop_160
           info = 3
           return
           ! end iwhila
           170 continue
           ! move q's to the front.
           do k = 2,n
              z(k) = z(4*k - 3)
           end do
           ! sort and compute sum of eigenvalues.
           call la_xlasrt('D',n,z,iinfo)
           e = zero
           do k = n,1,-1
              e = e + z(k)
           end do
           ! store trace, sum(eigenvalues) and information on performance.
           z(2*n + 1) = trace
           z(2*n + 2) = e
           z(2*n + 3) = real(iter,KIND=xdp)
           z(2*n + 4) = real(ndiv,KIND=xdp)/real(n**2,KIND=xdp)
           z(2*n + 5) = hundrd*nfail/real(iter,KIND=xdp)
           return
     end subroutine la_xlasq2
#endif
#ifdef LA_WITH_QP
     !> QLASQ2: computes all the eigenvalues of the symmetric positive
     !> definite tridiagonal matrix associated with the qd array Z to high
     !> relative accuracy are computed to high relative accuracy, in the
     !> absence of denormalization, underflow and overflow.
     !> To see the relation of Z to the tridiagonal matrix, let L be a
     !> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
     !> let U be an upper bidiagonal matrix with 1's above and diagonal
     !> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
     !> symmetric tridiagonal to which it is similar.
     !> Note : QLASQ2 defines a logical variable, IEEE, which is true
     !> on machines which follow ieee-754 floating-point standard in their
     !> handling of infinities and NaNs, and false otherwise. This variable
     !> is passed to QLASQ3.

     pure subroutine la_qlasq2(n,z,info)
        use la_constants_qp,only:zero,half,one,two,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: z(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: cbias = 1.50_qp
           real(qp),parameter :: hundrd = 100.0_qp

           ! Local Scalars
           logical(lk) :: ieee
           integer(ilp) :: i0,i1,i4,iinfo,ipn4,iter,iwhila,iwhilb,k,kmin,n0,n1,nbig, &
                     ndiv,nfail,pp,splt,ttype
           real(qp) :: d,dee,deemin,desig,dmin,dmin1,dmin2,dn,dn1,dn2,e,emax,emin, &
           eps,g,oldemn,qmax,qmin,s,safmin,sigma,t,tau,temp,tol,tol2,trace,zmax, &
                     tempe,tempq
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sqrt
           ! Executable Statements
           ! test the input arguments.
           ! (in case la_qlasq2 is not called by la_qlasq1)
           info = 0
           eps = la_qlamch('PRECISION')
           safmin = la_qlamch('SAFE MINIMUM')
           tol = eps*hundrd
           tol2 = tol**2
           if (n < 0) then
              info = -1
              call la_xerbla('QLASQ2',1)
              return
           else if (n == 0) then
              return
           else if (n == 1) then
              ! 1-by-1 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('QLASQ2',2)
              end if
              return
           else if (n == 2) then
              ! 2-by-2 case.
              if (z(1) < zero) then
                 info = -201
                 call la_xerbla('QLASQ2',2)
                 return
              else if (z(2) < zero) then
                 info = -202
                 call la_xerbla('QLASQ2',2)
                 return
              else if (z(3) < zero) then
                info = -203
                call la_xerbla('QLASQ2',2)
                return
              else if (z(3) > z(1)) then
                 d = z(3)
                 z(3) = z(1)
                 z(1) = d
              end if
              z(5) = z(1) + z(2) + z(3)
              if (z(2) > z(3)*tol2) then
                 t = half*((z(1) - z(3)) + z(2))
                 s = z(3)*(z(2)/t)
                 if (s <= t) then
                    s = z(3)*(z(2)/(t*(one + sqrt(one + s/t))))
                 else
                    s = z(3)*(z(2)/(t + sqrt(t)*sqrt(t + s)))
                 end if
                 t = z(1) + (s + z(2))
                 z(3) = z(3)*(z(1)/t)
                 z(1) = t
              end if
              z(2) = z(3)
              z(6) = z(2) + z(1)
              return
           end if
           ! check for negative data and compute sums of q's and e's.
           z(2*n) = zero
           emin = z(2)
           qmax = zero
           zmax = zero
           d = zero
           e = zero
           do k = 1,2*(n - 1),2
              if (z(k) < zero) then
                 info = -(200 + k)
                 call la_xerbla('QLASQ2',2)
                 return
              else if (z(k + 1) < zero) then
                 info = -(200 + k + 1)
                 call la_xerbla('QLASQ2',2)
                 return
              end if
              d = d + z(k)
              e = e + z(k + 1)
              qmax = max(qmax,z(k))
              emin = min(emin,z(k + 1))
              zmax = max(qmax,zmax,z(k + 1))
           end do
           if (z(2*n - 1) < zero) then
              info = -(200 + 2*n - 1)
              call la_xerbla('QLASQ2',2)
              return
           end if
           d = d + z(2*n - 1)
           qmax = max(qmax,z(2*n - 1))
           zmax = max(qmax,zmax)
           ! check for diagonality.
           if (e == zero) then
              do k = 2,n
                 z(k) = z(2*k - 1)
              end do
              call la_qlasrt('D',n,z,iinfo)
              z(2*n - 1) = d
              return
           end if
           trace = d + e
           ! check for zero data.
           if (trace == zero) then
              z(2*n - 1) = zero
              return
           end if
           ! check whether the machine is ieee conformable.
           ieee = (la_ilaenv(10,'QLASQ2','N',1,2,3,4) == 1)
           ! rearrange data for locality: z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
           do k = 2*n,2,-2
              z(2*k) = zero
              z(2*k - 1) = z(k)
              z(2*k - 2) = zero
              z(2*k - 3) = z(k - 1)
           end do
           i0 = 1
           n0 = n
           ! reverse the qd-array, if warranted.
           if (cbias*z(4*i0 - 3) < z(4*n0 - 3)) then
              ipn4 = 4*(i0 + n0)
              do i4 = 4*i0,2*(i0 + n0 - 1),4
                 temp = z(i4 - 3)
                 z(i4 - 3) = z(ipn4 - i4 - 3)
                 z(ipn4 - i4 - 3) = temp
                 temp = z(i4 - 1)
                 z(i4 - 1) = z(ipn4 - i4 - 5)
                 z(ipn4 - i4 - 5) = temp
              end do
           end if
           ! initial split checking via dqd and li's test.
           pp = 0
           loop_80: do k = 1,2
              d = z(4*n0 + pp - 3)
              do i4 = 4*(n0 - 1) + pp,4*i0 + pp,-4
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    d = z(i4 - 3)
                 else
                    d = z(i4 - 3)*(d/(d + z(i4 - 1)))
                 end if
              end do
              ! dqd maps z to zz plus li's test.
              emin = z(4*i0 + pp + 1)
              d = z(4*i0 + pp - 3)
              do i4 = 4*i0 + pp,4*(n0 - 1) + pp,4
                 z(i4 - 2*pp - 2) = d + z(i4 - 1)
                 if (z(i4 - 1) <= tol2*d) then
                    z(i4 - 1) = -zero
                    z(i4 - 2*pp - 2) = d
                    z(i4 - 2*pp) = zero
                    d = z(i4 + 1)
                 else if (safmin*z(i4 + 1) < z(i4 - 2*pp - 2) .and. safmin*z(i4 - 2*pp - 2) < z(i4 + 1)) &
                           then
                    temp = z(i4 + 1)/z(i4 - 2*pp - 2)
                    z(i4 - 2*pp) = z(i4 - 1)*temp
                    d = d*temp
                 else
                    z(i4 - 2*pp) = z(i4 + 1)*(z(i4 - 1)/z(i4 - 2*pp - 2))
                    d = z(i4 + 1)*(d/z(i4 - 2*pp - 2))
                 end if
                 emin = min(emin,z(i4 - 2*pp))
              end do
              z(4*n0 - pp - 2) = d
              ! now find qmax.
              qmax = z(4*i0 - pp - 2)
              do i4 = 4*i0 - pp + 2,4*n0 - pp - 2,4
                 qmax = max(qmax,z(i4))
              end do
              ! prepare for the next iteration on k.
              pp = 1 - pp
           end do loop_80
           ! initialise variables to pass to la_qlasq3.
           ttype = 0
           dmin1 = zero
           dmin2 = zero
           dn = zero
           dn1 = zero
           dn2 = zero
           g = zero
           tau = zero
           iter = 2
           nfail = 0
           ndiv = 2*(n0 - i0)
           loop_160: do iwhila = 1,n + 1
              if (n0 < 1) go to 170
              ! while array unfinished do
              ! e(n0) holds the value of sigma when submatrix in i0:n0
              ! splits from the rest of the array, but is negated.
              desig = zero
              if (n0 == n) then
                 sigma = zero
              else
                 sigma = -z(4*n0 - 1)
              end if
              if (sigma < zero) then
                 info = 1
                 return
              end if
              ! find last unreduced submatrix's top index i0, find qmax and
              ! emin. find gershgorin-type bound if q's much greater than e's.
              emax = zero
              if (n0 > i0) then
                 emin = abs(z(4*n0 - 5))
              else
                 emin = zero
              end if
              qmin = z(4*n0 - 3)
              qmax = qmin
              do i4 = 4*n0,8,-4
                 if (z(i4 - 5) <= zero) go to 100
                 if (qmin >= four*emax) then
                    qmin = min(qmin,z(i4 - 3))
                    emax = max(emax,z(i4 - 5))
                 end if
                 qmax = max(qmax,z(i4 - 7) + z(i4 - 5))
                 emin = min(emin,z(i4 - 5))
              end do
              i4 = 4
              100 continue
              i0 = i4/4
              pp = 0
              if (n0 - i0 > 1) then
                 dee = z(4*i0 - 3)
                 deemin = dee
                 kmin = i0
                 do i4 = 4*i0 + 1,4*n0 - 3,4
                    dee = z(i4)*(dee/(dee + z(i4 - 2)))
                    if (dee <= deemin) then
                       deemin = dee
                       kmin = (i4 + 3)/4
                    end if
                 end do
                 if ((kmin - i0)*2 < n0 - kmin .and. deemin <= half*z(4*n0 - 3)) then
                    ipn4 = 4*(i0 + n0)
                    pp = 2
                    do i4 = 4*i0,2*(i0 + n0 - 1),4
                       temp = z(i4 - 3)
                       z(i4 - 3) = z(ipn4 - i4 - 3)
                       z(ipn4 - i4 - 3) = temp
                       temp = z(i4 - 2)
                       z(i4 - 2) = z(ipn4 - i4 - 2)
                       z(ipn4 - i4 - 2) = temp
                       temp = z(i4 - 1)
                       z(i4 - 1) = z(ipn4 - i4 - 5)
                       z(ipn4 - i4 - 5) = temp
                       temp = z(i4)
                       z(i4) = z(ipn4 - i4 - 4)
                       z(ipn4 - i4 - 4) = temp
                    end do
                 end if
              end if
              ! put -(initial shift) into dmin.
              dmin = -max(zero,qmin - two*sqrt(qmin)*sqrt(emax))
              ! now i0:n0 is unreduced.
              ! pp = 0 for ping, pp = 1 for pong.
              ! pp = 2 indicates that flipping was applied to the z array and
                     ! and that the tests for deflation upon entry in la_qlasq3
                     ! should not be performed.
              nbig = 100*(n0 - i0 + 1)
              loop_140: do iwhilb = 1,nbig
                 if (i0 > n0) go to 150
                 ! while submatrix unfinished take a good dqds step.
                 call la_qlasq3(i0,n0,z,pp,dmin,sigma,desig,qmax,nfail,iter,ndiv, &
                           ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
                 pp = 1 - pp
                 ! when emin is very small check for splits.
                 if (pp == 0 .and. n0 - i0 >= 3) then
                    if (z(4*n0) <= tol2*qmax .or. z(4*n0 - 1) <= tol2*sigma) then
                       splt = i0 - 1
                       qmax = z(4*i0 - 3)
                       emin = z(4*i0 - 1)
                       oldemn = z(4*i0)
                       do i4 = 4*i0,4*(n0 - 3),4
                          if (z(i4) <= tol2*z(i4 - 3) .or. z(i4 - 1) <= tol2*sigma) then
                             z(i4 - 1) = -sigma
                             splt = i4/4
                             qmax = zero
                             emin = z(i4 + 3)
                             oldemn = z(i4 + 4)
                          else
                             qmax = max(qmax,z(i4 + 1))
                             emin = min(emin,z(i4 - 1))
                             oldemn = min(oldemn,z(i4))
                          end if
                       end do
                       z(4*n0 - 1) = emin
                       z(4*n0) = oldemn
                       i0 = splt + 1
                    end if
                 end if
              end do loop_140
              info = 2
              ! maximum number of iterations exceeded, restore the shift
              ! sigma and place the new d's and e's in a qd array.
              ! this might need to be done for several blocks
              i1 = i0
              n1 = n0
              145 continue
              tempq = z(4*i0 - 3)
              z(4*i0 - 3) = z(4*i0 - 3) + sigma
              do k = i0 + 1,n0
                 tempe = z(4*k - 5)
                 z(4*k - 5) = z(4*k - 5)*(tempq/z(4*k - 7))
                 tempq = z(4*k - 3)
                 z(4*k - 3) = z(4*k - 3) + sigma + tempe - z(4*k - 5)
              end do
              ! prepare to do this on the previous block if there is one
              if (i1 > 1) then
                 n1 = i1 - 1
                 do while ((i1 >= 2) .and. (z(4*i1 - 5) >= zero))
                    i1 = i1 - 1
                 end do
                 sigma = -z(4*n1 - 1)
                 go to 145
              end if
              do k = 1,n
                 z(2*k - 1) = z(4*k - 3)
              ! only the block 1..n0 is unfinished.  the rest of the e's
              ! must be essentially zero, although sometimes other data
              ! has been stored in them.
                 if (k < n0) then
                    z(2*k) = z(4*k - 1)
                 else
                    z(2*k) = 0
                 end if
              end do
              return
              ! end iwhilb
              150 continue
           end do loop_160
           info = 3
           return
           ! end iwhila
           170 continue
           ! move q's to the front.
           do k = 2,n
              z(k) = z(4*k - 3)
           end do
           ! sort and compute sum of eigenvalues.
           call la_qlasrt('D',n,z,iinfo)
           e = zero
           do k = n,1,-1
              e = e + z(k)
           end do
           ! store trace, sum(eigenvalues) and information on performance.
           z(2*n + 1) = trace
           z(2*n + 2) = e
           z(2*n + 3) = real(iter,KIND=qp)
           z(2*n + 4) = real(ndiv,KIND=qp)/real(n**2,KIND=qp)
           z(2*n + 5) = hundrd*nfail/real(iter,KIND=qp)
           return
     end subroutine la_qlasq2
#endif

     !> CBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**H
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**H*VT instead of
     !> P**H, for given complex input matrices U and VT.  When U and VT are
     !> the unitary matrices that reduce a general matrix A to bidiagonal
     !> form: A = U*B*VT, as computed by CGEBRD, then
     !> A = (U*Q) * S * (P**H*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**H*C
     !> for a given complex input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_cbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,rwork, &
                info)
        use la_constants_sp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(inout) :: c(ldc,*),u(ldu,*),vt(ldvt,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: hndrth = 0.01_sp
           real(sp),parameter :: hndrd = 100.0_sp
           real(sp),parameter :: meigth = -0.125_sp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,j,ll,lll,m,maxit,nm1,nm12,nm13,oldll, &
                     oldm
           real(sp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,max,min,real,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('CBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_slasq1(n,d,e,rwork,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_slamch('EPSILON')
           unfl = la_slamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_slartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 rwork(i) = cs
                 rwork(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_clasr('R','V','F',nru,n,rwork(1),rwork(n),u,ldu)

              if (ncc > 0) call la_clasr('L','V','F',n,ncc,rwork(1),rwork(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=sp))
              thresh = max(tol*sminoa,maxitr*n*n*unfl)
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*n*n*unfl)
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxit = maxitr*n*n
           iter = 0
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter > maxit) go to 200
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_slasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_csrot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr, &
                        sinr)
              if (nru > 0) call la_csrot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)

              if (ncc > 0) call la_csrot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_slas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_slas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_slartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_slartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll + 1) = cs
                    rwork(i - ll + 1 + nm1) = sn
                    rwork(i - ll + 1 + nm12) = oldcs
                    rwork(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_clasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_clasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_clasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_slartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_slartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll) = cs
                    rwork(i - ll + nm1) = -sn
                    rwork(i - ll + nm12) = oldcs
                    rwork(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_clasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_clasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_clasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_slartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_slartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    rwork(i - ll + 1) = cosr
                    rwork(i - ll + 1 + nm1) = sinr
                    rwork(i - ll + 1 + nm12) = cosl
                    rwork(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_clasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_clasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_clasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_slartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_slartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    rwork(i - ll) = cosr
                    rwork(i - ll + nm1) = -sinr
                    rwork(i - ll + nm12) = cosl
                    rwork(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_clasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_clasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_clasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_csscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_cswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_cswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_cswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_cbdsqr
     !> ZBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**H
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**H*VT instead of
     !> P**H, for given complex input matrices U and VT.  When U and VT are
     !> the unitary matrices that reduce a general matrix A to bidiagonal
     !> form: A = U*B*VT, as computed by ZGEBRD, then
     !> A = (U*Q) * S * (P**H*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**H*C
     !> for a given complex input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_zbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,rwork, &
                info)
        use la_constants_dp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(inout) :: c(ldc,*),u(ldu,*),vt(ldvt,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: hndrth = 0.01_dp
           real(dp),parameter :: hndrd = 100.0_dp
           real(dp),parameter :: meigth = -0.125_dp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,j,ll,lll,m,maxit,nm1,nm12,nm13,oldll, &
                     oldm
           real(dp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('ZBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_dlasq1(n,d,e,rwork,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_dlamch('EPSILON')
           unfl = la_dlamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_dlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 rwork(i) = cs
                 rwork(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_zlasr('R','V','F',nru,n,rwork(1),rwork(n),u,ldu)

              if (ncc > 0) call la_zlasr('L','V','F',n,ncc,rwork(1),rwork(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=dp))
              thresh = max(tol*sminoa,maxitr*n*n*unfl)
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*n*n*unfl)
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxit = maxitr*n*n
           iter = 0
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter > maxit) go to 200
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_dlasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_zdrot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr, &
                        sinr)
              if (nru > 0) call la_zdrot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)

              if (ncc > 0) call la_zdrot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_dlas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_dlas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_dlartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_dlartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll + 1) = cs
                    rwork(i - ll + 1 + nm1) = sn
                    rwork(i - ll + 1 + nm12) = oldcs
                    rwork(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_zlasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_zlasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_zlasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_dlartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_dlartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll) = cs
                    rwork(i - ll + nm1) = -sn
                    rwork(i - ll + nm12) = oldcs
                    rwork(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_zlasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_zlasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_zlasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_dlartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_dlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    rwork(i - ll + 1) = cosr
                    rwork(i - ll + 1 + nm1) = sinr
                    rwork(i - ll + 1 + nm12) = cosl
                    rwork(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_zlasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_zlasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_zlasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_dlartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_dlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    rwork(i - ll) = cosr
                    rwork(i - ll + nm1) = -sinr
                    rwork(i - ll + nm12) = cosl
                    rwork(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_zlasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_zlasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_zlasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_zdscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_zswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_zswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_zswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_zbdsqr
#ifdef LA_WITH_XDP
     !> YBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**H
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**H*VT instead of
     !> P**H, for given complex input matrices U and VT.  When U and VT are
     !> the unitary matrices that reduce a general matrix A to bidiagonal
     !> form: A = U*B*VT, as computed by YGEBRD, then
     !> A = (U*Q) * S * (P**H*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**H*C
     !> for a given complex input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_ybdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,rwork, &
                info)
        use la_constants_xdp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(xdp),intent(inout) :: d(*),e(*)
           real(xdp),intent(out) :: rwork(*)
           complex(xdp),intent(inout) :: c(ldc,*),u(ldu,*),vt(ldvt,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: hndrth = 0.01_xdp
           real(xdp),parameter :: hndrd = 100.0_xdp
           real(xdp),parameter :: meigth = -0.125_xdp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,j,ll,lll,m,maxit,nm1,nm12,nm13,oldll, &
                     oldm
           real(xdp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('YBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_xlasq1(n,d,e,rwork,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_xlamch('EPSILON')
           unfl = la_xlamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_xlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 rwork(i) = cs
                 rwork(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_ylasr('R','V','F',nru,n,rwork(1),rwork(n),u,ldu)

              if (ncc > 0) call la_ylasr('L','V','F',n,ncc,rwork(1),rwork(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=xdp))
              thresh = max(tol*sminoa,maxitr*n*n*unfl)
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*n*n*unfl)
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxit = maxitr*n*n
           iter = 0
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter > maxit) go to 200
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_xlasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_yxrot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr, &
                        sinr)
              if (nru > 0) call la_yxrot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)

              if (ncc > 0) call la_yxrot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_xlas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_xlas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_xlartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_xlartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll + 1) = cs
                    rwork(i - ll + 1 + nm1) = sn
                    rwork(i - ll + 1 + nm12) = oldcs
                    rwork(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_ylasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_ylasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_ylasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_xlartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_xlartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll) = cs
                    rwork(i - ll + nm1) = -sn
                    rwork(i - ll + nm12) = oldcs
                    rwork(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_ylasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_ylasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_ylasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_xlartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_xlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    rwork(i - ll + 1) = cosr
                    rwork(i - ll + 1 + nm1) = sinr
                    rwork(i - ll + 1 + nm12) = cosl
                    rwork(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_ylasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_ylasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_ylasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_xlartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_xlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    rwork(i - ll) = cosr
                    rwork(i - ll + nm1) = -sinr
                    rwork(i - ll + nm12) = cosl
                    rwork(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_ylasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_ylasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_ylasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_yxscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_yswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_yswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_yswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_ybdsqr
#endif
#ifdef LA_WITH_QP
     !> WBDSQR: computes the singular values and, optionally, the right and/or
     !> left singular vectors from the singular value decomposition (SVD) of
     !> a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
     !> zero-shift QR algorithm.  The SVD of B has the form
     !> B = Q * S * P**H
     !> where S is the diagonal matrix of singular values, Q is an orthogonal
     !> matrix of left singular vectors, and P is an orthogonal matrix of
     !> right singular vectors.  If left singular vectors are requested, this
     !> subroutine actually returns U*Q instead of Q, and, if right singular
     !> vectors are requested, this subroutine returns P**H*VT instead of
     !> P**H, for given complex input matrices U and VT.  When U and VT are
     !> the unitary matrices that reduce a general matrix A to bidiagonal
     !> form: A = U*B*VT, as computed by WGEBRD, then
     !> A = (U*Q) * S * (P**H*VT)
     !> is the SVD of A.  Optionally, the subroutine may also compute Q**H*C
     !> for a given complex input matrix C.
     !> See "Computing  Small Singular Values of Bidiagonal Matrices With
     !> Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
     !> LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
     !> no. 5, pp. 873-912, Sept 1990) and
     !> "Accurate singular values and differential qd algorithms," by
     !> B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
     !> Department, University of California at Berkeley, July 1992
     !> for a detailed description of the algorithm.

     pure subroutine la_wbdsqr(uplo,n,ncvt,nru,ncc,d,e,vt,ldvt,u,ldu,c,ldc,rwork, &
                info)
        use la_constants_qp,only:negone,zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldc,ldu,ldvt,n,ncc,ncvt,nru
           ! Array Arguments
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(inout) :: c(ldc,*),u(ldu,*),vt(ldvt,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: hndrth = 0.01_qp
           real(qp),parameter :: hndrd = 100.0_qp
           real(qp),parameter :: meigth = -0.125_qp
           integer(ilp),parameter :: maxitr = 6

           ! Local Scalars
           logical(lk) :: lower,rotate
           integer(ilp) :: i,idir,isub,iter,j,ll,lll,m,maxit,nm1,nm12,nm13,oldll, &
                     oldm
           real(qp) :: abse,abss,cosl,cosr,cs,eps,f,g,h,mu,oldcs,oldsn,r,shift, &
           sigmn,sigmx,sinl,sinr,sll,smax,smin,sminl,sminoa,sn,thresh,tol,tolmul, &
                     unfl
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lower = la_lsame(uplo,'L')
           if (.not. la_lsame(uplo,'U') .and. .not. lower) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ncvt < 0) then
              info = -3
           else if (nru < 0) then
              info = -4
           else if (ncc < 0) then
              info = -5
           else if ((ncvt == 0 .and. ldvt < 1) .or. (ncvt > 0 .and. ldvt < max(1,n))) then
              info = -9
           else if (ldu < max(1,nru)) then
              info = -11
           else if ((ncc == 0 .and. ldc < 1) .or. (ncc > 0 .and. ldc < max(1,n))) then
              info = -13
           end if
           if (info /= 0) then
              call la_xerbla('WBDSQR',-info)
              return
           end if
           if (n == 0) return
           if (n == 1) go to 160
           ! rotate is true if any singular vectors desired, false otherwise
           rotate = (ncvt > 0) .or. (nru > 0) .or. (ncc > 0)
           ! if no singular vectors desired, use qd algorithm
           if (.not. rotate) then
              call la_qlasq1(n,d,e,rwork,info)
           ! if info equals 2, dqds didn't finish, try to finish
              if (info /= 2) return
              info = 0
           end if
           nm1 = n - 1
           nm12 = nm1 + nm1
           nm13 = nm12 + nm1
           idir = 0
           ! get machine constants
           eps = la_qlamch('EPSILON')
           unfl = la_qlamch('SAFE MINIMUM')
           ! if matrix lower bidiagonal, rotate to be upper bidiagonal
           ! by applying givens rotations on the left
           if (lower) then
              do i = 1,n - 1
                 call la_qlartg(d(i),e(i),cs,sn,r)
                 d(i) = r
                 e(i) = sn*d(i + 1)
                 d(i + 1) = cs*d(i + 1)
                 rwork(i) = cs
                 rwork(nm1 + i) = sn
              end do
              ! update singular vectors if desired
              if (nru > 0) call la_wlasr('R','V','F',nru,n,rwork(1),rwork(n),u,ldu)

              if (ncc > 0) call la_wlasr('L','V','F',n,ncc,rwork(1),rwork(n),c,ldc)

           end if
           ! compute singular values to relative accuracy tol
           ! (by setting tol to be negative, algorithm will compute
           ! singular values to absolute accuracy abs(tol)*norm(input matrix))
           tolmul = max(ten,min(hndrd,eps**meigth))
           tol = tolmul*eps
           ! compute approximate maximum, minimum singular values
           smax = zero
           do i = 1,n
              smax = max(smax,abs(d(i)))
           end do
           do i = 1,n - 1
              smax = max(smax,abs(e(i)))
           end do
           sminl = zero
           if (tol >= zero) then
              ! relative accuracy desired
              sminoa = abs(d(1))
              if (sminoa == zero) go to 50
              mu = sminoa
              do i = 2,n
                 mu = abs(d(i))*(mu/(mu + abs(e(i - 1))))
                 sminoa = min(sminoa,mu)
                 if (sminoa == zero) go to 50
              end do
              50 continue
              sminoa = sminoa/sqrt(real(n,KIND=qp))
              thresh = max(tol*sminoa,maxitr*n*n*unfl)
           else
              ! absolute accuracy desired
              thresh = max(abs(tol)*smax,maxitr*n*n*unfl)
           end if
           ! prepare for main iteration loop for the singular values
           ! (maxit is the maximum number of passes through the inner
           ! loop permitted before nonconvergence signalled.)
           maxit = maxitr*n*n
           iter = 0
           oldll = -1
           oldm = -1
           ! m points to last element of unconverged part of matrix
           m = n
           ! begin main iteration loop
           60 continue
           ! check for convergence or exceeding iteration count
           if (m <= 1) go to 160
           if (iter > maxit) go to 200
           ! find diagonal block of matrix to work on
           if (tol < zero .and. abs(d(m)) <= thresh) d(m) = zero
           smax = abs(d(m))
           smin = smax
           do lll = 1,m - 1
              ll = m - lll
              abss = abs(d(ll))
              abse = abs(e(ll))
              if (tol < zero .and. abss <= thresh) d(ll) = zero
              if (abse <= thresh) go to 80
              smin = min(smin,abss)
              smax = max(smax,abss,abse)
           end do
           ll = 0
           go to 90
           80 continue
           e(ll) = zero
           ! matrix splits since e(ll) = 0
           if (ll == m - 1) then
              ! convergence of bottom singular value, return to top of loop
              m = m - 1
              go to 60
           end if
           90 continue
           ll = ll + 1
           ! e(ll) through e(m-1) are nonzero, e(ll-1) is zero
           if (ll == m - 1) then
              ! 2 by 2 block, handle separately
              call la_qlasv2(d(m - 1),e(m - 1),d(m),sigmn,sigmx,sinr,cosr,sinl,cosl &
                        )
              d(m - 1) = sigmx
              e(m - 1) = zero
              d(m) = sigmn
              ! compute singular vectors, if desired
              if (ncvt > 0) call la_wqrot(ncvt,vt(m - 1,1),ldvt,vt(m,1),ldvt,cosr, &
                        sinr)
              if (nru > 0) call la_wqrot(nru,u(1,m - 1),1,u(1,m),1,cosl,sinl)

              if (ncc > 0) call la_wqrot(ncc,c(m - 1,1),ldc,c(m,1),ldc,cosl,sinl)

              m = m - 2
              go to 60
           end if
           ! if working on new submatrix, choose shift direction
           ! (from larger end diagonal element towards smaller)
           if (ll > oldm .or. m < oldll) then
              if (abs(d(ll)) >= abs(d(m))) then
                 ! chase bulge from top (big end) to bottom (small end)
                 idir = 1
              else
                 ! chase bulge from bottom (big end) to top (small end)
                 idir = 2
              end if
           end if
           ! apply convergence tests
           if (idir == 1) then
              ! run convergence test in forward direction
              ! first apply standard test to bottom of matrix
              if (abs(e(m - 1)) <= abs(tol)*abs(d(m)) .or. (tol < zero .and. abs(e(m - 1)) &
                        <= thresh)) then
                 e(m - 1) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion forward
                 mu = abs(d(ll))
                 sminl = mu
                 do lll = ll,m - 1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll + 1))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           else
              ! run convergence test in backward direction
              ! first apply standard test to top of matrix
              if (abs(e(ll)) <= abs(tol)*abs(d(ll)) .or. (tol < zero .and. abs(e(ll)) &
                        <= thresh)) then
                 e(ll) = zero
                 go to 60
              end if
              if (tol >= zero) then
                 ! if relative accuracy desired,
                 ! apply convergence criterion backward
                 mu = abs(d(m))
                 sminl = mu
                 do lll = m - 1,ll,-1
                    if (abs(e(lll)) <= tol*mu) then
                       e(lll) = zero
                       go to 60
                    end if
                    mu = abs(d(lll))*(mu/(mu + abs(e(lll))))
                    sminl = min(sminl,mu)
                 end do
              end if
           end if
           oldll = ll
           oldm = m
           ! compute shift.  first, test if shifting would ruin relative
           ! accuracy, and if so set the shift to zero.
           if (tol >= zero .and. n*tol*(sminl/smax) <= max(eps,hndrth*tol)) then
              ! use a zero shift to avoid loss of relative accuracy
              shift = zero
           else
              ! compute the shift from 2-by-2 block at end of matrix
              if (idir == 1) then
                 sll = abs(d(ll))
                 call la_qlas2(d(m - 1),e(m - 1),d(m),shift,r)
              else
                 sll = abs(d(m))
                 call la_qlas2(d(ll),e(ll),d(ll + 1),shift,r)
              end if
              ! test if shift negligible, and if so set to zero
              if (sll > zero) then
                 if ((shift/sll)**2 < eps) shift = zero
              end if
           end if
           ! increment iteration count
           iter = iter + m - ll
           ! if shift = 0, do simplified qr iteration
           if (shift == zero) then
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = ll,m - 1
                    call la_qlartg(d(i)*cs,e(i),cs,sn,r)
                    if (i > ll) e(i - 1) = oldsn*r
                    call la_qlartg(oldcs*r,d(i + 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll + 1) = cs
                    rwork(i - ll + 1 + nm1) = sn
                    rwork(i - ll + 1 + nm12) = oldcs
                    rwork(i - ll + 1 + nm13) = oldsn
                 end do
                 h = d(m)*cs
                 d(m) = h*oldcs
                 e(m - 1) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_wlasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_wlasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_wlasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 cs = one
                 oldcs = one
                 do i = m,ll + 1,-1
                    call la_qlartg(d(i)*cs,e(i - 1),cs,sn,r)
                    if (i < m) e(i) = oldsn*r
                    call la_qlartg(oldcs*r,d(i - 1)*sn,oldcs,oldsn,d(i))
                    rwork(i - ll) = cs
                    rwork(i - ll + nm1) = -sn
                    rwork(i - ll + nm12) = oldcs
                    rwork(i - ll + nm13) = -oldsn
                 end do
                 h = d(ll)*cs
                 d(ll) = h*oldcs
                 e(ll) = h*oldsn
                 ! update singular vectors
                 if (ncvt > 0) call la_wlasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_wlasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_wlasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
              end if
           else
              ! use nonzero shift
              if (idir == 1) then
                 ! chase bulge from top to bottom
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(ll)) - shift)*(sign(one,d(ll)) + shift/d(ll))
                 g = e(ll)
                 do i = ll,m - 1
                    call la_qlartg(f,g,cosr,sinr,r)
                    if (i > ll) e(i - 1) = r
                    f = cosr*d(i) + sinr*e(i)
                    e(i) = cosr*e(i) - sinr*d(i)
                    g = sinr*d(i + 1)
                    d(i + 1) = cosr*d(i + 1)
                    call la_qlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i) + sinl*d(i + 1)
                    d(i + 1) = cosl*d(i + 1) - sinl*e(i)
                    if (i < m - 1) then
                       g = sinl*e(i + 1)
                       e(i + 1) = cosl*e(i + 1)
                    end if
                    rwork(i - ll + 1) = cosr
                    rwork(i - ll + 1 + nm1) = sinr
                    rwork(i - ll + 1 + nm12) = cosl
                    rwork(i - ll + 1 + nm13) = sinl
                 end do
                 e(m - 1) = f
                 ! update singular vectors
                 if (ncvt > 0) call la_wlasr('L','V','F',m - ll + 1,ncvt,rwork(1),rwork(n) &
                           ,vt(ll,1),ldvt)
                 if (nru > 0) call la_wlasr('R','V','F',nru,m - ll + 1,rwork(nm12 + 1),rwork( &
                           nm13 + 1),u(1,ll),ldu)
                 if (ncc > 0) call la_wlasr('L','V','F',m - ll + 1,ncc,rwork(nm12 + 1),rwork( &
                           nm13 + 1),c(ll,1),ldc)
                 ! test convergence
                 if (abs(e(m - 1)) <= thresh) e(m - 1) = zero
              else
                 ! chase bulge from bottom to top
                 ! save cosines and sines for later singular vector updates
                 f = (abs(d(m)) - shift)*(sign(one,d(m)) + shift/d(m))
                 g = e(m - 1)
                 do i = m,ll + 1,-1
                    call la_qlartg(f,g,cosr,sinr,r)
                    if (i < m) e(i) = r
                    f = cosr*d(i) + sinr*e(i - 1)
                    e(i - 1) = cosr*e(i - 1) - sinr*d(i)
                    g = sinr*d(i - 1)
                    d(i - 1) = cosr*d(i - 1)
                    call la_qlartg(f,g,cosl,sinl,r)
                    d(i) = r
                    f = cosl*e(i - 1) + sinl*d(i - 1)
                    d(i - 1) = cosl*d(i - 1) - sinl*e(i - 1)
                    if (i > ll + 1) then
                       g = sinl*e(i - 2)
                       e(i - 2) = cosl*e(i - 2)
                    end if
                    rwork(i - ll) = cosr
                    rwork(i - ll + nm1) = -sinr
                    rwork(i - ll + nm12) = cosl
                    rwork(i - ll + nm13) = -sinl
                 end do
                 e(ll) = f
                 ! test convergence
                 if (abs(e(ll)) <= thresh) e(ll) = zero
                 ! update singular vectors if desired
                 if (ncvt > 0) call la_wlasr('L','V','B',m - ll + 1,ncvt,rwork(nm12 + 1), &
                           rwork(nm13 + 1),vt(ll,1),ldvt)
                 if (nru > 0) call la_wlasr('R','V','B',nru,m - ll + 1,rwork(1),rwork(n), &
                           u(1,ll),ldu)
                 if (ncc > 0) call la_wlasr('L','V','B',m - ll + 1,ncc,rwork(1),rwork(n), &
                           c(ll,1),ldc)
              end if
           end if
           ! qr iteration finished, go back and check convergence
           go to 60
           ! all singular values converged, so make them positive
           160 continue
           do i = 1,n
              if (d(i) < zero) then
                 d(i) = -d(i)
                 ! change sign of singular vectors, if desired
                 if (ncvt > 0) call la_wqscal(ncvt,negone,vt(i,1),ldvt)
              end if
           end do
           ! sort the singular values into decreasing order (insertion sort on
           ! singular values, but only one transposition per singular vector)
           do i = 1,n - 1
              ! scan for smallest d(i)
              isub = 1
              smin = d(1)
              do j = 2,n + 1 - i
                 if (d(j) <= smin) then
                    isub = j
                    smin = d(j)
                 end if
              end do
              if (isub /= n + 1 - i) then
                 ! swap singular values and vectors
                 d(isub) = d(n + 1 - i)
                 d(n + 1 - i) = smin
                 if (ncvt > 0) call la_wswap(ncvt,vt(isub,1),ldvt,vt(n + 1 - i,1),ldvt)

                 if (nru > 0) call la_wswap(nru,u(1,isub),1,u(1,n + 1 - i),1)
                 if (ncc > 0) call la_wswap(ncc,c(isub,1),ldc,c(n + 1 - i,1),ldc)

              end if
           end do
           go to 220
           ! maximum number of iterations exceeded, failure to converge
           200 continue
           info = 0
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           220 continue
           return
     end subroutine la_wbdsqr
#endif

end module la_lapack_svd_bidiag_qr
