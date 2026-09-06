!> Triangular systems: solve, inverse, condition estimation, refinement
module la_lapack_solve_tri_comp
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_tri
     use la_blas_level3_gen
     use la_blas_level3_sym
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_mnorm
     use la_lapack_blas_like_scalar
     use la_lapack_solve_aux
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slatbs
     public :: la_slatps
     public :: la_slatrs
     public :: la_slauu2
     public :: la_slauum
     public :: la_stbrfs
     public :: la_stbtrs
     public :: la_stprfs
     public :: la_stptri
     public :: la_stptrs
     public :: la_strrfs
     public :: la_strti2
     public :: la_strtri
     public :: la_strtrs
     public :: la_stbcon
     public :: la_stftri
     public :: la_stpcon
     public :: la_strcon
     public :: la_dlatbs
     public :: la_dlatps
     public :: la_dlatrs
     public :: la_dlauu2
     public :: la_dlauum
     public :: la_dtbrfs
     public :: la_dtbtrs
     public :: la_dtprfs
     public :: la_dtptri
     public :: la_dtptrs
     public :: la_dtrrfs
     public :: la_dtrti2
     public :: la_dtrtri
     public :: la_dtrtrs
     public :: la_dtbcon
     public :: la_dtftri
     public :: la_dtpcon
     public :: la_dtrcon
     public :: la_qlatbs
     public :: la_qlatps
     public :: la_qlatrs
     public :: la_qlauu2
     public :: la_qlauum
     public :: la_qtbrfs
     public :: la_qtbtrs
     public :: la_qtprfs
     public :: la_qtptri
     public :: la_qtptrs
     public :: la_qtrrfs
     public :: la_qtrti2
     public :: la_qtrtri
     public :: la_qtrtrs
     public :: la_qtbcon
     public :: la_qtftri
     public :: la_qtpcon
     public :: la_qtrcon
     public :: la_clatbs
     public :: la_clatps
     public :: la_clatrs
     public :: la_clauu2
     public :: la_clauum
     public :: la_ctbrfs
     public :: la_ctbtrs
     public :: la_ctprfs
     public :: la_ctptri
     public :: la_ctptrs
     public :: la_ctrrfs
     public :: la_ctrti2
     public :: la_ctrtri
     public :: la_ctrtrs
     public :: la_ctbcon
     public :: la_ctftri
     public :: la_ctpcon
     public :: la_ctrcon
     public :: la_zlatbs
     public :: la_zlatps
     public :: la_zlatrs
     public :: la_zlauu2
     public :: la_zlauum
     public :: la_ztbrfs
     public :: la_ztbtrs
     public :: la_ztprfs
     public :: la_ztptri
     public :: la_ztptrs
     public :: la_ztrrfs
     public :: la_ztrti2
     public :: la_ztrtri
     public :: la_ztrtrs
     public :: la_ztbcon
     public :: la_ztftri
     public :: la_ztpcon
     public :: la_ztrcon
     public :: la_wlatbs
     public :: la_wlatps
     public :: la_wlatrs
     public :: la_wlauu2
     public :: la_wlauum
     public :: la_wtbrfs
     public :: la_wtbtrs
     public :: la_wtprfs
     public :: la_wtptri
     public :: la_wtptrs
     public :: la_wtrrfs
     public :: la_wtrti2
     public :: la_wtrtri
     public :: la_wtrtrs
     public :: la_wtbcon
     public :: la_wtftri
     public :: la_wtpcon
     public :: la_wtrcon

     contains

     !> SLATBS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular band matrix.  Here A**T denotes the transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine STBSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_slatbs(uplo,trans,diag,normin,n,kd,ab,ldab,x,scale,cnorm, &
               info)
        use la_constants_sp,only:zero,half,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast,jlen,maind
           real(sp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (kd < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('SLATBS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    jlen = min(kd,j - 1)
                    cnorm(j) = la_sasum(jlen,ab(kd + 1 - jlen,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n
                    jlen = min(kd,n - j)
                    if (jlen > 0) then
                       cnorm(j) = la_sasum(jlen,ab(2,j),1)
                    else
                       cnorm(j) = zero
                    end if
                 end do
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_isamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_sscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_stbsv can be used.
           j = la_isamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = kd + 1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(ab(maind,j))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = kd + 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(ab(maind,j))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_stbsv(uplo,trans,diag,n,kd,ab,ldab,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_sscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_100: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = ab(maind,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 95
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_sscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    95 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_sscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                                   ! x(j)* a(max(1,j-kd):j-1,j)
                          jlen = min(kd,j - 1)
                          call la_saxpy(jlen,-x(j)*tscal,ab(kd + 1 - jlen,j),1,x(j - jlen &
                                    ),1)
                          i = la_isamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                    else if (j < n) then
                       ! compute the update
                          ! x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                                ! x(j) * a(j+1:min(j+kd,n),j)
                       jlen = min(kd,n - j)
                       if (jlen > 0) call la_saxpy(jlen,-x(j)*tscal,ab(2,j),1,x(j + 1), &
                                  1)
                       i = j + la_isamax(n - j,x(j + 1),1)
                       xmax = abs(x(i))
                    end if
                 end do loop_100
              else
                 ! solve a**t * x = b
                 loop_140: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_sdot to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          sumj = la_sdot(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)
                       else
                          jlen = min(kd,n - j)
                          if (jlen > 0) sumj = la_sdot(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             sumj = sumj + (ab(kd + i - jlen,j)*uscal)*x(j - jlen - 1 + i)
                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             sumj = sumj + (ab(i + 1,j)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 135
                       end if
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_sscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_sscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       135 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                 end do loop_140
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_sscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_slatbs
     !> DLATBS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular band matrix.  Here A**T denotes the transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine DTBSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_dlatbs(uplo,trans,diag,normin,n,kd,ab,ldab,x,scale,cnorm, &
               info)
        use la_constants_dp,only:zero,half,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast,jlen,maind
           real(dp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (kd < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DLATBS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    jlen = min(kd,j - 1)
                    cnorm(j) = la_dasum(jlen,ab(kd + 1 - jlen,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n
                    jlen = min(kd,n - j)
                    if (jlen > 0) then
                       cnorm(j) = la_dasum(jlen,ab(2,j),1)
                    else
                       cnorm(j) = zero
                    end if
                 end do
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_idamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_dscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_dtbsv can be used.
           j = la_idamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = kd + 1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(ab(maind,j))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = kd + 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(ab(maind,j))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_dtbsv(uplo,trans,diag,n,kd,ab,ldab,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_dscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = ab(maind,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 100
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_dscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    100 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_dscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                                   ! x(j)* a(max(1,j-kd):j-1,j)
                          jlen = min(kd,j - 1)
                          call la_daxpy(jlen,-x(j)*tscal,ab(kd + 1 - jlen,j),1,x(j - jlen &
                                    ),1)
                          i = la_idamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                    else if (j < n) then
                       ! compute the update
                          ! x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                                ! x(j) * a(j+1:min(j+kd,n),j)
                       jlen = min(kd,n - j)
                       if (jlen > 0) call la_daxpy(jlen,-x(j)*tscal,ab(2,j),1,x(j + 1), &
                                  1)
                       i = j + la_idamax(n - j,x(j + 1),1)
                       xmax = abs(x(i))
                    end if
                 end do loop_110
              else
                 ! solve a**t * x = b
                 loop_160: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_ddot to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          sumj = la_ddot(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)
                       else
                          jlen = min(kd,n - j)
                          if (jlen > 0) sumj = la_ddot(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             sumj = sumj + (ab(kd + i - jlen,j)*uscal)*x(j - jlen - 1 + i)
                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             sumj = sumj + (ab(i + 1,j)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 150
                       end if
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_dscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_dscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       150 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                 end do loop_160
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_dscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_dlatbs
     !> QLATBS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular band matrix.  Here A**T denotes the transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine QTBSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_qlatbs(uplo,trans,diag,normin,n,kd,ab,ldab,x,scale,cnorm, &
               info)
        use la_constants_qp,only:zero,half,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast,jlen,maind
           real(qp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (kd < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QLATBS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    jlen = min(kd,j - 1)
                    cnorm(j) = la_qasum(jlen,ab(kd + 1 - jlen,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n
                    jlen = min(kd,n - j)
                    if (jlen > 0) then
                       cnorm(j) = la_qasum(jlen,ab(2,j),1)
                    else
                       cnorm(j) = zero
                    end if
                 end do
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_iqamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_qscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_qtbsv can be used.
           j = la_iqamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = kd + 1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(ab(maind,j))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = kd + 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(ab(maind,j))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_qtbsv(uplo,trans,diag,n,kd,ab,ldab,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_qscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = ab(maind,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 100
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_qscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    100 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_qscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                                   ! x(j)* a(max(1,j-kd):j-1,j)
                          jlen = min(kd,j - 1)
                          call la_qaxpy(jlen,-x(j)*tscal,ab(kd + 1 - jlen,j),1,x(j - jlen &
                                    ),1)
                          i = la_iqamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                    else if (j < n) then
                       ! compute the update
                          ! x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                                ! x(j) * a(j+1:min(j+kd,n),j)
                       jlen = min(kd,n - j)
                       if (jlen > 0) call la_qaxpy(jlen,-x(j)*tscal,ab(2,j),1,x(j + 1), &
                                  1)
                       i = j + la_iqamax(n - j,x(j + 1),1)
                       xmax = abs(x(i))
                    end if
                 end do loop_110
              else
                 ! solve a**t * x = b
                 loop_160: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_qdot to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          sumj = la_qdot(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)
                       else
                          jlen = min(kd,n - j)
                          if (jlen > 0) sumj = la_qdot(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             sumj = sumj + (ab(kd + i - jlen,j)*uscal)*x(j - jlen - 1 + i)
                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             sumj = sumj + (ab(i + 1,j)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 150
                       end if
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_qscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_qscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       150 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                 end do loop_160
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_qscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_qlatbs

     !> SLATPS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular matrix stored in packed form.  Here A**T denotes the
     !> transpose of A, x and b are n-element vectors, and s is a scaling
     !> factor, usually less than or equal to 1, chosen so that the
     !> components of x will be less than the overflow threshold.  If the
     !> unscaled problem will not cause overflow, the Level 2 BLAS routine
     !> STPSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
     !> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.

     pure subroutine la_slatps(uplo,trans,diag,normin,n,ap,x,scale,cnorm,info)
        use la_constants_sp,only:zero,half,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,ip,j,jfirst,jinc,jlast,jlen
           real(sp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('SLATPS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 ip = 1
                 do j = 1,n
                    cnorm(j) = la_sasum(j - 1,ap(ip),1)
                    ip = ip + j
                 end do
              else
                 ! a is lower triangular.
                 ip = 1
                 do j = 1,n - 1
                    cnorm(j) = la_sasum(n - j,ap(ip + 1),1)
                    ip = ip + n - j + 1
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_isamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_sscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_stpsv can be used.
           j = la_isamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = n
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(ap(ip))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                    ip = ip + jinc*jlen
                    jlen = jlen - 1
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(ap(ip))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_stpsv(uplo,trans,diag,n,ap,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_sscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 ip = jfirst*(jfirst + 1)/2
                 loop_100: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = ap(ip)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 95
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_sscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    95 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_sscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_saxpy(j - 1,-x(j)*tscal,ap(ip - j + 1),1,x,1)
                          i = la_isamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                       ip = ip - j
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_saxpy(n - j,-x(j)*tscal,ap(ip + 1),1,x(j + 1),1)

                          i = j + la_isamax(n - j,x(j + 1),1)
                          xmax = abs(x(i))
                       end if
                       ip = ip + n - j + 1
                    end if
                 end do loop_100
              else
                 ! solve a**t * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_140: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_sdot to perform the dot product.
                       if (upper) then
                          sumj = la_sdot(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          sumj = la_sdot(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             sumj = sumj + (ap(ip - j + i)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             sumj = sumj + (ap(ip + i)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 135
                       end if
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_sscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_sscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       135 continue
                    else
                       ! compute x(j) := x(j) / a(j,j)  - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_140
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_sscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_slatps
     !> DLATPS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular matrix stored in packed form.  Here A**T denotes the
     !> transpose of A, x and b are n-element vectors, and s is a scaling
     !> factor, usually less than or equal to 1, chosen so that the
     !> components of x will be less than the overflow threshold.  If the
     !> unscaled problem will not cause overflow, the Level 2 BLAS routine
     !> DTPSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
     !> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.

     pure subroutine la_dlatps(uplo,trans,diag,normin,n,ap,x,scale,cnorm,info)
        use la_constants_dp,only:zero,half,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,ip,j,jfirst,jinc,jlast,jlen
           real(dp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DLATPS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 ip = 1
                 do j = 1,n
                    cnorm(j) = la_dasum(j - 1,ap(ip),1)
                    ip = ip + j
                 end do
              else
                 ! a is lower triangular.
                 ip = 1
                 do j = 1,n - 1
                    cnorm(j) = la_dasum(n - j,ap(ip + 1),1)
                    ip = ip + n - j + 1
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_idamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_dscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_dtpsv can be used.
           j = la_idamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = n
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(ap(ip))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                    ip = ip + jinc*jlen
                    jlen = jlen - 1
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(ap(ip))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_dtpsv(uplo,trans,diag,n,ap,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_dscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 ip = jfirst*(jfirst + 1)/2
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = ap(ip)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 100
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_dscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    100 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_dscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_daxpy(j - 1,-x(j)*tscal,ap(ip - j + 1),1,x,1)
                          i = la_idamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                       ip = ip - j
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_daxpy(n - j,-x(j)*tscal,ap(ip + 1),1,x(j + 1),1)

                          i = j + la_idamax(n - j,x(j + 1),1)
                          xmax = abs(x(i))
                       end if
                       ip = ip + n - j + 1
                    end if
                 end do loop_110
              else
                 ! solve a**t * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_160: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_ddot to perform the dot product.
                       if (upper) then
                          sumj = la_ddot(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          sumj = la_ddot(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             sumj = sumj + (ap(ip - j + i)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             sumj = sumj + (ap(ip + i)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 150
                       end if
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_dscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_dscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       150 continue
                    else
                       ! compute x(j) := x(j) / a(j,j)  - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_160
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_dscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_dlatps
     !> QLATPS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular matrix stored in packed form.  Here A**T denotes the
     !> transpose of A, x and b are n-element vectors, and s is a scaling
     !> factor, usually less than or equal to 1, chosen so that the
     !> components of x will be less than the overflow threshold.  If the
     !> unscaled problem will not cause overflow, the Level 2 BLAS routine
     !> QTPSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
     !> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.

     pure subroutine la_qlatps(uplo,trans,diag,normin,n,ap,x,scale,cnorm,info)
        use la_constants_qp,only:zero,half,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,ip,j,jfirst,jinc,jlast,jlen
           real(qp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QLATPS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 ip = 1
                 do j = 1,n
                    cnorm(j) = la_qasum(j - 1,ap(ip),1)
                    ip = ip + j
                 end do
              else
                 ! a is lower triangular.
                 ip = 1
                 do j = 1,n - 1
                    cnorm(j) = la_qasum(n - j,ap(ip + 1),1)
                    ip = ip + n - j + 1
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_iqamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_qscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_qtpsv can be used.
           j = la_iqamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = n
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(ap(ip))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                    ip = ip + jinc*jlen
                    jlen = jlen - 1
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(ap(ip))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_qtpsv(uplo,trans,diag,n,ap,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_qscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 ip = jfirst*(jfirst + 1)/2
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = ap(ip)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 100
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_qscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    100 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_qscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_qaxpy(j - 1,-x(j)*tscal,ap(ip - j + 1),1,x,1)
                          i = la_iqamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                       ip = ip - j
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_qaxpy(n - j,-x(j)*tscal,ap(ip + 1),1,x(j + 1),1)

                          i = j + la_iqamax(n - j,x(j + 1),1)
                          xmax = abs(x(i))
                       end if
                       ip = ip + n - j + 1
                    end if
                 end do loop_110
              else
                 ! solve a**t * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_160: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_qdot to perform the dot product.
                       if (upper) then
                          sumj = la_qdot(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          sumj = la_qdot(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             sumj = sumj + (ap(ip - j + i)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             sumj = sumj + (ap(ip + i)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 150
                       end if
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_qscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_qscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       150 continue
                    else
                       ! compute x(j) := x(j) / a(j,j)  - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_160
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_qscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_qlatps

     !> SLATRS: solves one of the triangular systems
     !> A *x = s*b  or  A**T*x = s*b
     !> with scaling to prevent overflow.  Here A is an upper or lower
     !> triangular matrix, A**T denotes the transpose of A, x and b are
     !> n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine STRSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_slatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
        use la_constants_sp,only:zero,half,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast
           real(sp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SLATRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    cnorm(j) = la_sasum(j - 1,a(1,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n - 1
                    cnorm(j) = la_sasum(n - j,a(j + 1,j),1)
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_isamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_sscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_strsv can be used.
           j = la_isamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(a(j,j))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(a(j,j))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_strsv(uplo,trans,diag,n,a,lda,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_sscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_100: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = a(j,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 95
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_sscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    95 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_sscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_saxpy(j - 1,-x(j)*tscal,a(1,j),1,x,1)
                          i = la_isamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_saxpy(n - j,-x(j)*tscal,a(j + 1,j),1,x(j + 1),1)

                          i = j + la_isamax(n - j,x(j + 1),1)
                          xmax = abs(x(i))
                       end if
                    end if
                 end do loop_100
              else
                 ! solve a**t * x = b
                 loop_140: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_sscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_sdot to perform the dot product.
                       if (upper) then
                          sumj = la_sdot(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          sumj = la_sdot(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             sumj = sumj + (a(i,j)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             sumj = sumj + (a(i,j)*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 135
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_sscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_sscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       135 continue
                    else
                       ! compute x(j) := x(j) / a(j,j)  - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                 end do loop_140
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_sscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_slatrs
     !> DLATRS: solves one of the triangular systems
     !> A *x = s*b  or  A**T *x = s*b
     !> with scaling to prevent overflow.  Here A is an upper or lower
     !> triangular matrix, A**T denotes the transpose of A, x and b are
     !> n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_dlatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
        use la_constants_dp,only:zero,half,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast
           real(dp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DLATRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    cnorm(j) = la_dasum(j - 1,a(1,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n - 1
                    cnorm(j) = la_dasum(n - j,a(j + 1,j),1)
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_idamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_dscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_dtrsv can be used.
           j = la_idamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(a(j,j))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(a(j,j))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_dtrsv(uplo,trans,diag,n,a,lda,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_dscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = a(j,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 100
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_dscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    100 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_dscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_daxpy(j - 1,-x(j)*tscal,a(1,j),1,x,1)
                          i = la_idamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_daxpy(n - j,-x(j)*tscal,a(j + 1,j),1,x(j + 1),1)

                          i = j + la_idamax(n - j,x(j + 1),1)
                          xmax = abs(x(i))
                       end if
                    end if
                 end do loop_110
              else
                 ! solve a**t * x = b
                 loop_160: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_dscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_ddot to perform the dot product.
                       if (upper) then
                          sumj = la_ddot(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          sumj = la_ddot(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             sumj = sumj + (a(i,j)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             sumj = sumj + (a(i,j)*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 150
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_dscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_dscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       150 continue
                    else
                       ! compute x(j) := x(j) / a(j,j)  - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                 end do loop_160
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_dscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_dlatrs
     !> QLATRS: solves one of the triangular systems
     !> A *x = s*b  or  A**T *x = s*b
     !> with scaling to prevent overflow.  Here A is an upper or lower
     !> triangular matrix, A**T denotes the transpose of A, x and b are
     !> n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine QTRSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_qlatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
        use la_constants_qp,only:zero,half,one

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: cnorm(*),x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast
           real(qp) :: bignum,grow,rec,smlnum,sumj,tjj,tjjs,tmax,tscal,uscal,xbnd,xj, &
                     xmax
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QLATRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    cnorm(j) = la_qasum(j - 1,a(1,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n - 1
                    cnorm(j) = la_qasum(n - j,a(j + 1,j),1)
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum.
           imax = la_iqamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum) then
              tscal = one
           else
              tscal = one/(smlnum*tmax)
              call la_qscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_qtrsv can be used.
           j = la_iqamax(n,x,1)
           xmax = abs(x(j))
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 50
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! m(j) = g(j-1) / abs(a(j,j))
                    tjj = abs(a(j,j))
                    xbnd = min(xbnd,min(one,tjj)*grow)
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 50
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              50 continue
           else
              ! compute the growth in a**t * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 80
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = one/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                    tjj = abs(a(j,j))
                    if (xj > tjj) xbnd = xbnd*(tjj/xj)
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,one/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 80
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              80 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_qtrsv(uplo,trans,diag,n,a,lda,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = bignum/xmax
                 call la_qscal(n,scale,x,1)
                 xmax = bignum
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = abs(x(j))
                    if (nounit) then
                       tjjs = a(j,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 100
                    end if
                    tjj = abs(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_qscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = x(j)/tjjs
                       xj = abs(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    100 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_qscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_qaxpy(j - 1,-x(j)*tscal,a(1,j),1,x,1)
                          i = la_iqamax(j - 1,x,1)
                          xmax = abs(x(i))
                       end if
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_qaxpy(n - j,-x(j)*tscal,a(j + 1,j),1,x(j + 1),1)

                          i = j + la_iqamax(n - j,x(j + 1),1)
                          xmax = abs(x(i))
                       end if
                    end if
                 end do loop_110
              else
                 ! solve a**t * x = b
                 loop_160: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = abs(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = abs(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = uscal/tjjs
                       end if
                       if (rec < one) then
                          call la_qscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    sumj = zero
                    if (uscal == one) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_qdot to perform the dot product.
                       if (upper) then
                          sumj = la_qdot(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          sumj = la_qdot(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             sumj = sumj + (a(i,j)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             sumj = sumj + (a(i,j)*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == tscal) then
                       ! compute x(j) := ( x(j) - sumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - sumj
                       xj = abs(x(j))
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 150
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = abs(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_qscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = x(j)/tjjs
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_qscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = x(j)/tjjs
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0, and compute a solution to a**t*x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       150 continue
                    else
                       ! compute x(j) := x(j) / a(j,j)  - sumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = x(j)/tjjs - sumj
                    end if
                    xmax = max(xmax,abs(x(j)))
                 end do loop_160
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_qscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_qlatrs

     !> SLAUU2: computes the product U * U**T or L**T * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the unblocked form of the algorithm, calling Level 2 BLAS.

     pure subroutine la_slauu2(uplo,n,a,lda,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SLAUU2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! compute the product u * u**t.
              do i = 1,n
                 aii = a(i,i)
                 if (i < n) then
                    a(i,i) = la_sdot(n - i + 1,a(i,i),lda,a(i,i),lda)
                    call la_sgemv('NO TRANSPOSE',i - 1,n - i,one,a(1,i + 1),lda,a(i,i + 1) &
                              ,lda,aii,a(1,i),1)
                 else
                    call la_sscal(i,aii,a(1,i),1)
                 end if
              end do
           else
              ! compute the product l**t * l.
              do i = 1,n
                 aii = a(i,i)
                 if (i < n) then
                    a(i,i) = la_sdot(n - i + 1,a(i,i),1,a(i,i),1)
                    call la_sgemv('TRANSPOSE',n - i,i - 1,one,a(i + 1,1),lda,a(i + 1,i), &
                              1,aii,a(i,1),lda)
                 else
                    call la_sscal(i,aii,a(i,1),lda)
                 end if
              end do
           end if
           return
     end subroutine la_slauu2
     !> DLAUU2: computes the product U * U**T or L**T * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the unblocked form of the algorithm, calling Level 2 BLAS.

     pure subroutine la_dlauu2(uplo,n,a,lda,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DLAUU2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! compute the product u * u**t.
              do i = 1,n
                 aii = a(i,i)
                 if (i < n) then
                    a(i,i) = la_ddot(n - i + 1,a(i,i),lda,a(i,i),lda)
                    call la_dgemv('NO TRANSPOSE',i - 1,n - i,one,a(1,i + 1),lda,a(i,i + 1) &
                              ,lda,aii,a(1,i),1)
                 else
                    call la_dscal(i,aii,a(1,i),1)
                 end if
              end do
           else
              ! compute the product l**t * l.
              do i = 1,n
                 aii = a(i,i)
                 if (i < n) then
                    a(i,i) = la_ddot(n - i + 1,a(i,i),1,a(i,i),1)
                    call la_dgemv('TRANSPOSE',n - i,i - 1,one,a(i + 1,1),lda,a(i + 1,i), &
                              1,aii,a(i,1),lda)
                 else
                    call la_dscal(i,aii,a(i,1),lda)
                 end if
              end do
           end if
           return
     end subroutine la_dlauu2
     !> QLAUU2: computes the product U * U**T or L**T * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the unblocked form of the algorithm, calling Level 2 BLAS.

     pure subroutine la_qlauu2(uplo,n,a,lda,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QLAUU2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! compute the product u * u**t.
              do i = 1,n
                 aii = a(i,i)
                 if (i < n) then
                    a(i,i) = la_qdot(n - i + 1,a(i,i),lda,a(i,i),lda)
                    call la_qgemv('NO TRANSPOSE',i - 1,n - i,one,a(1,i + 1),lda,a(i,i + 1) &
                              ,lda,aii,a(1,i),1)
                 else
                    call la_qscal(i,aii,a(1,i),1)
                 end if
              end do
           else
              ! compute the product l**t * l.
              do i = 1,n
                 aii = a(i,i)
                 if (i < n) then
                    a(i,i) = la_qdot(n - i + 1,a(i,i),1,a(i,i),1)
                    call la_qgemv('TRANSPOSE',n - i,i - 1,one,a(i + 1,1),lda,a(i + 1,i), &
                              1,aii,a(i,1),lda)
                 else
                    call la_qscal(i,aii,a(i,1),lda)
                 end if
              end do
           end if
           return
     end subroutine la_qlauu2

     !> SLAUUM: computes the product U * U**T or L**T * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the blocked form of the algorithm, calling Level 3 BLAS.

     pure subroutine la_slauum(uplo,n,a,lda,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ib,nb
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SLAUUM',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'SLAUUM',uplo,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_slauu2(uplo,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute the product u * u**t.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_strmm('RIGHT','UPPER','TRANSPOSE','NON-UNIT',i - 1,ib,one,a( &
                              i,i),lda,a(1,i),lda)
                    call la_slauu2('UPPER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_sgemm('NO TRANSPOSE','TRANSPOSE',i - 1,ib,n - i - ib + 1,one,a( &
                                 1,i + ib),lda,a(i,i + ib),lda,one,a(1,i),lda)
                       call la_ssyrk('UPPER','NO TRANSPOSE',ib,n - i - ib + 1,one,a(i,i + ib), &
                                  lda,one,a(i,i),lda)
                    end if
                 end do
              else
                 ! compute the product l**t * l.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_strmm('LEFT','LOWER','TRANSPOSE','NON-UNIT',ib,i - 1,one,a( &
                              i,i),lda,a(i,1),lda)
                    call la_slauu2('LOWER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_sgemm('TRANSPOSE','NO TRANSPOSE',ib,i - 1,n - i - ib + 1,one,a( &
                                 i + ib,i),lda,a(i + ib,1),lda,one,a(i,1),lda)
                       call la_ssyrk('LOWER','TRANSPOSE',ib,n - i - ib + 1,one,a(i + ib,i), &
                                 lda,one,a(i,i),lda)
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_slauum
     !> DLAUUM: computes the product U * U**T or L**T * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the blocked form of the algorithm, calling Level 3 BLAS.

     pure subroutine la_dlauum(uplo,n,a,lda,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ib,nb
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DLAUUM',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'DLAUUM',uplo,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_dlauu2(uplo,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute the product u * u**t.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_dtrmm('RIGHT','UPPER','TRANSPOSE','NON-UNIT',i - 1,ib,one,a( &
                              i,i),lda,a(1,i),lda)
                    call la_dlauu2('UPPER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_dgemm('NO TRANSPOSE','TRANSPOSE',i - 1,ib,n - i - ib + 1,one,a( &
                                 1,i + ib),lda,a(i,i + ib),lda,one,a(1,i),lda)
                       call la_dsyrk('UPPER','NO TRANSPOSE',ib,n - i - ib + 1,one,a(i,i + ib), &
                                  lda,one,a(i,i),lda)
                    end if
                 end do
              else
                 ! compute the product l**t * l.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_dtrmm('LEFT','LOWER','TRANSPOSE','NON-UNIT',ib,i - 1,one,a( &
                              i,i),lda,a(i,1),lda)
                    call la_dlauu2('LOWER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_dgemm('TRANSPOSE','NO TRANSPOSE',ib,i - 1,n - i - ib + 1,one,a( &
                                 i + ib,i),lda,a(i + ib,1),lda,one,a(i,1),lda)
                       call la_dsyrk('LOWER','TRANSPOSE',ib,n - i - ib + 1,one,a(i + ib,i), &
                                 lda,one,a(i,i),lda)
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_dlauum
     !> QLAUUM: computes the product U * U**T or L**T * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the blocked form of the algorithm, calling Level 3 BLAS.

     pure subroutine la_qlauum(uplo,n,a,lda,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ib,nb
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QLAUUM',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'QLAUUM',uplo,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_qlauu2(uplo,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute the product u * u**t.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_qtrmm('RIGHT','UPPER','TRANSPOSE','NON-UNIT',i - 1,ib,one,a( &
                              i,i),lda,a(1,i),lda)
                    call la_qlauu2('UPPER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_qgemm('NO TRANSPOSE','TRANSPOSE',i - 1,ib,n - i - ib + 1,one,a( &
                                 1,i + ib),lda,a(i,i + ib),lda,one,a(1,i),lda)
                       call la_qsyrk('UPPER','NO TRANSPOSE',ib,n - i - ib + 1,one,a(i,i + ib), &
                                  lda,one,a(i,i),lda)
                    end if
                 end do
              else
                 ! compute the product l**t * l.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_qtrmm('LEFT','LOWER','TRANSPOSE','NON-UNIT',ib,i - 1,one,a( &
                              i,i),lda,a(i,1),lda)
                    call la_qlauu2('LOWER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_qgemm('TRANSPOSE','NO TRANSPOSE',ib,i - 1,n - i - ib + 1,one,a( &
                                 i + ib,i),lda,a(i + ib,1),lda,one,a(i,1),lda)
                       call la_qsyrk('LOWER','TRANSPOSE',ib,n - i - ib + 1,one,a(i + ib,i), &
                                 lda,one,a(i,i),lda)
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_qlauum

     !> STBRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular band
     !> coefficient matrix.
     !> The solution matrix X must be computed by STBTRS or some other
     !> means before entering this routine.  STBRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_stbrfs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,x,ldx,ferr, &
                berr,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: ab(ldab,*),b(ldb,*),x(ldx,*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('STBRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = kd + 2
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_scopy(n,x(1,j),1,work(n + 1),1)
              call la_stbmv(uplo,trans,diag,n,kd,ab,ldab,work(n + 1),1)
              call la_saxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = max(1,k - kd),k
                             work(i) = work(i) + abs(ab(kd + 1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = max(1,k - kd),k - 1
                             work(i) = work(i) + abs(ab(kd + 1 + i - k,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,min(n,k + kd)
                             work(i) = work(i) + abs(ab(1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             work(i) = work(i) + abs(ab(1 + i - k,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = max(1,k - kd),k
                             s = s + abs(ab(kd + 1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = max(1,k - kd),k - 1
                             s = s + abs(ab(kd + 1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,min(n,k + kd)
                             s = s + abs(ab(1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             s = s + abs(ab(1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_slacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_slacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_stbsv(uplo,transt,diag,n,kd,ab,ldab,work(n + 1),1)

                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_stbsv(uplo,trans,diag,n,kd,ab,ldab,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_stbrfs
     !> DTBRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular band
     !> coefficient matrix.
     !> The solution matrix X must be computed by DTBTRS or some other
     !> means before entering this routine.  DTBRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_dtbrfs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,x,ldx,ferr, &
                berr,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: ab(ldab,*),b(ldb,*),x(ldx,*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('DTBRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = kd + 2
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_dcopy(n,x(1,j),1,work(n + 1),1)
              call la_dtbmv(uplo,trans,diag,n,kd,ab,ldab,work(n + 1),1)
              call la_daxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = max(1,k - kd),k
                             work(i) = work(i) + abs(ab(kd + 1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = max(1,k - kd),k - 1
                             work(i) = work(i) + abs(ab(kd + 1 + i - k,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,min(n,k + kd)
                             work(i) = work(i) + abs(ab(1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             work(i) = work(i) + abs(ab(1 + i - k,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = max(1,k - kd),k
                             s = s + abs(ab(kd + 1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = max(1,k - kd),k - 1
                             s = s + abs(ab(kd + 1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,min(n,k + kd)
                             s = s + abs(ab(1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             s = s + abs(ab(1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_dlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_dlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_dtbsv(uplo,transt,diag,n,kd,ab,ldab,work(n + 1),1)

                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_dtbsv(uplo,trans,diag,n,kd,ab,ldab,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_dtbrfs
     !> QTBRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular band
     !> coefficient matrix.
     !> The solution matrix X must be computed by QTBTRS or some other
     !> means before entering this routine.  QTBRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_qtbrfs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,x,ldx,ferr, &
                berr,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: ab(ldab,*),b(ldb,*),x(ldx,*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('QTBRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = kd + 2
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_qcopy(n,x(1,j),1,work(n + 1),1)
              call la_qtbmv(uplo,trans,diag,n,kd,ab,ldab,work(n + 1),1)
              call la_qaxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = max(1,k - kd),k
                             work(i) = work(i) + abs(ab(kd + 1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = max(1,k - kd),k - 1
                             work(i) = work(i) + abs(ab(kd + 1 + i - k,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,min(n,k + kd)
                             work(i) = work(i) + abs(ab(1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             work(i) = work(i) + abs(ab(1 + i - k,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = max(1,k - kd),k
                             s = s + abs(ab(kd + 1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = max(1,k - kd),k - 1
                             s = s + abs(ab(kd + 1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,min(n,k + kd)
                             s = s + abs(ab(1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             s = s + abs(ab(1 + i - k,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_qlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_qlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_qtbsv(uplo,transt,diag,n,kd,ab,ldab,work(n + 1),1)

                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_qtbsv(uplo,trans,diag,n,kd,ab,ldab,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_qtbrfs

     !> STBTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular band matrix of order N, and B is an
     !> N-by NRHS matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_stbtrs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('STBTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 do info = 1,n
                    if (ab(kd + 1,info) == zero) return
                 end do
              else
                 do info = 1,n
                    if (ab(1,info) == zero) return
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           do j = 1,nrhs
              call la_stbsv(uplo,trans,diag,n,kd,ab,ldab,b(1,j),1)
           end do
           return
     end subroutine la_stbtrs
     !> DTBTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular band matrix of order N, and B is an
     !> N-by NRHS matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_dtbtrs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('DTBTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 do info = 1,n
                    if (ab(kd + 1,info) == zero) return
                 end do
              else
                 do info = 1,n
                    if (ab(1,info) == zero) return
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           do j = 1,nrhs
              call la_dtbsv(uplo,trans,diag,n,kd,ab,ldab,b(1,j),1)
           end do
           return
     end subroutine la_dtbtrs
     !> QTBTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular band matrix of order N, and B is an
     !> N-by NRHS matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_qtbtrs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('QTBTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 do info = 1,n
                    if (ab(kd + 1,info) == zero) return
                 end do
              else
                 do info = 1,n
                    if (ab(1,info) == zero) return
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           do j = 1,nrhs
              call la_qtbsv(uplo,trans,diag,n,kd,ab,ldab,b(1,j),1)
           end do
           return
     end subroutine la_qtbtrs

     !> STPRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular packed
     !> coefficient matrix.
     !> The solution matrix X must be computed by STPTRS or some other
     !> means before entering this routine.  STPRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_stprfs(uplo,trans,diag,n,nrhs,ap,b,ldb,x,ldx,ferr,berr, &
               work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: ap(*),b(ldb,*),x(ldx,*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,kc,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('STPRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_scopy(n,x(1,j),1,work(n + 1),1)
              call la_stpmv(uplo,trans,diag,n,ap,work(n + 1),1)
              call la_saxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k
                             work(i) = work(i) + abs(ap(kc + i - 1))*xk
                          end do
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k - 1
                             work(i) = work(i) + abs(ap(kc + i - 1))*xk
                          end do
                          work(k) = work(k) + xk
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,n
                             work(i) = work(i) + abs(ap(kc + i - k))*xk
                          end do
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,n
                             work(i) = work(i) + abs(ap(kc + i - k))*xk
                          end do
                          work(k) = work(k) + xk
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + abs(ap(kc + i - 1))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = 1,k - 1
                             s = s + abs(ap(kc + i - 1))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + abs(ap(kc + i - k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,n
                             s = s + abs(ap(kc + i - k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_slacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_slacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_stpsv(uplo,transt,diag,n,ap,work(n + 1),1)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_stpsv(uplo,trans,diag,n,ap,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_stprfs
     !> DTPRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular packed
     !> coefficient matrix.
     !> The solution matrix X must be computed by DTPTRS or some other
     !> means before entering this routine.  DTPRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_dtprfs(uplo,trans,diag,n,nrhs,ap,b,ldb,x,ldx,ferr,berr, &
               work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: ap(*),b(ldb,*),x(ldx,*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,kc,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('DTPRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_dcopy(n,x(1,j),1,work(n + 1),1)
              call la_dtpmv(uplo,trans,diag,n,ap,work(n + 1),1)
              call la_daxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k
                             work(i) = work(i) + abs(ap(kc + i - 1))*xk
                          end do
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k - 1
                             work(i) = work(i) + abs(ap(kc + i - 1))*xk
                          end do
                          work(k) = work(k) + xk
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,n
                             work(i) = work(i) + abs(ap(kc + i - k))*xk
                          end do
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,n
                             work(i) = work(i) + abs(ap(kc + i - k))*xk
                          end do
                          work(k) = work(k) + xk
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + abs(ap(kc + i - 1))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = 1,k - 1
                             s = s + abs(ap(kc + i - 1))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + abs(ap(kc + i - k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,n
                             s = s + abs(ap(kc + i - k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_dlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_dlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_dtpsv(uplo,transt,diag,n,ap,work(n + 1),1)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_dtpsv(uplo,trans,diag,n,ap,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_dtprfs
     !> QTPRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular packed
     !> coefficient matrix.
     !> The solution matrix X must be computed by QTPTRS or some other
     !> means before entering this routine.  QTPRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_qtprfs(uplo,trans,diag,n,nrhs,ap,b,ldb,x,ldx,ferr,berr, &
               work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: ap(*),b(ldb,*),x(ldx,*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,kc,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('QTPRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_qcopy(n,x(1,j),1,work(n + 1),1)
              call la_qtpmv(uplo,trans,diag,n,ap,work(n + 1),1)
              call la_qaxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k
                             work(i) = work(i) + abs(ap(kc + i - 1))*xk
                          end do
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k - 1
                             work(i) = work(i) + abs(ap(kc + i - 1))*xk
                          end do
                          work(k) = work(k) + xk
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,n
                             work(i) = work(i) + abs(ap(kc + i - k))*xk
                          end do
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,n
                             work(i) = work(i) + abs(ap(kc + i - k))*xk
                          end do
                          work(k) = work(k) + xk
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + abs(ap(kc + i - 1))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = 1,k - 1
                             s = s + abs(ap(kc + i - 1))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + abs(ap(kc + i - k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,n
                             s = s + abs(ap(kc + i - k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_qlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_qlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_qtpsv(uplo,transt,diag,n,ap,work(n + 1),1)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_qtpsv(uplo,trans,diag,n,ap,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_qtprfs

     !> STPTRI: computes the inverse of a real upper or lower triangular
     !> matrix A stored in packed format.

     pure subroutine la_stptri(uplo,diag,n,ap,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: ap(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc,jclast,jj
           real(sp) :: ajj
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('STPTRI',-info)
              return
           end if
           ! check for singularity if non-unit.
           if (nounit) then
              if (upper) then
                 jj = 0
                 do info = 1,n
                    jj = jj + info
                    if (ap(jj) == zero) return
                 end do
              else
                 jj = 1
                 do info = 1,n
                    if (ap(jj) == zero) return
                    jj = jj + n - info + 1
                 end do
              end if
              info = 0
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              jc = 1
              do j = 1,n
                 if (nounit) then
                    ap(jc + j - 1) = one/ap(jc + j - 1)
                    ajj = -ap(jc + j - 1)
                 else
                    ajj = -one
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_stpmv('UPPER','NO TRANSPOSE',diag,j - 1,ap,ap(jc),1)
                 call la_sscal(j - 1,ajj,ap(jc),1)
                 jc = jc + j
              end do
           else
              ! compute inverse of lower triangular matrix.
              jc = n*(n + 1)/2
              do j = n,1,-1
                 if (nounit) then
                    ap(jc) = one/ap(jc)
                    ajj = -ap(jc)
                 else
                    ajj = -one
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_stpmv('LOWER','NO TRANSPOSE',diag,n - j,ap(jclast),ap(jc + 1) &
                              ,1)
                    call la_sscal(n - j,ajj,ap(jc + 1),1)
                 end if
                 jclast = jc
                 jc = jc - n + j - 2
              end do
           end if
           return
     end subroutine la_stptri
     !> DTPTRI: computes the inverse of a real upper or lower triangular
     !> matrix A stored in packed format.

     pure subroutine la_dtptri(uplo,diag,n,ap,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: ap(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc,jclast,jj
           real(dp) :: ajj
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('DTPTRI',-info)
              return
           end if
           ! check for singularity if non-unit.
           if (nounit) then
              if (upper) then
                 jj = 0
                 do info = 1,n
                    jj = jj + info
                    if (ap(jj) == zero) return
                 end do
              else
                 jj = 1
                 do info = 1,n
                    if (ap(jj) == zero) return
                    jj = jj + n - info + 1
                 end do
              end if
              info = 0
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              jc = 1
              do j = 1,n
                 if (nounit) then
                    ap(jc + j - 1) = one/ap(jc + j - 1)
                    ajj = -ap(jc + j - 1)
                 else
                    ajj = -one
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_dtpmv('UPPER','NO TRANSPOSE',diag,j - 1,ap,ap(jc),1)
                 call la_dscal(j - 1,ajj,ap(jc),1)
                 jc = jc + j
              end do
           else
              ! compute inverse of lower triangular matrix.
              jc = n*(n + 1)/2
              do j = n,1,-1
                 if (nounit) then
                    ap(jc) = one/ap(jc)
                    ajj = -ap(jc)
                 else
                    ajj = -one
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_dtpmv('LOWER','NO TRANSPOSE',diag,n - j,ap(jclast),ap(jc + 1) &
                              ,1)
                    call la_dscal(n - j,ajj,ap(jc + 1),1)
                 end if
                 jclast = jc
                 jc = jc - n + j - 2
              end do
           end if
           return
     end subroutine la_dtptri
     !> QTPTRI: computes the inverse of a real upper or lower triangular
     !> matrix A stored in packed format.

     pure subroutine la_qtptri(uplo,diag,n,ap,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: ap(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc,jclast,jj
           real(qp) :: ajj
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('QTPTRI',-info)
              return
           end if
           ! check for singularity if non-unit.
           if (nounit) then
              if (upper) then
                 jj = 0
                 do info = 1,n
                    jj = jj + info
                    if (ap(jj) == zero) return
                 end do
              else
                 jj = 1
                 do info = 1,n
                    if (ap(jj) == zero) return
                    jj = jj + n - info + 1
                 end do
              end if
              info = 0
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              jc = 1
              do j = 1,n
                 if (nounit) then
                    ap(jc + j - 1) = one/ap(jc + j - 1)
                    ajj = -ap(jc + j - 1)
                 else
                    ajj = -one
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_qtpmv('UPPER','NO TRANSPOSE',diag,j - 1,ap,ap(jc),1)
                 call la_qscal(j - 1,ajj,ap(jc),1)
                 jc = jc + j
              end do
           else
              ! compute inverse of lower triangular matrix.
              jc = n*(n + 1)/2
              do j = n,1,-1
                 if (nounit) then
                    ap(jc) = one/ap(jc)
                    ajj = -ap(jc)
                 else
                    ajj = -one
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_qtpmv('LOWER','NO TRANSPOSE',diag,n - j,ap(jclast),ap(jc + 1) &
                              ,1)
                    call la_qscal(n - j,ajj,ap(jc + 1),1)
                 end if
                 jclast = jc
                 jc = jc - n + j - 2
              end do
           end if
           return
     end subroutine la_qtptri

     !> STPTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular matrix of order N stored in packed format,
     !> and B is an N-by-NRHS matrix.  A check is made to verify that A is
     !> nonsingular.

     pure subroutine la_stptrs(uplo,trans,diag,n,nrhs,ap,b,ldb,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(in) :: ap(*)
           real(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('STPTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 jc = 1
                 do info = 1,n
                    if (ap(jc + info - 1) == zero) return
                    jc = jc + info
                 end do
              else
                 jc = 1
                 do info = 1,n
                    if (ap(jc) == zero) return
                    jc = jc + n - info + 1
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           do j = 1,nrhs
              call la_stpsv(uplo,trans,diag,n,ap,b(1,j),1)
           end do
           return
     end subroutine la_stptrs
     !> DTPTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular matrix of order N stored in packed format,
     !> and B is an N-by-NRHS matrix.  A check is made to verify that A is
     !> nonsingular.

     pure subroutine la_dtptrs(uplo,trans,diag,n,nrhs,ap,b,ldb,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(in) :: ap(*)
           real(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('DTPTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 jc = 1
                 do info = 1,n
                    if (ap(jc + info - 1) == zero) return
                    jc = jc + info
                 end do
              else
                 jc = 1
                 do info = 1,n
                    if (ap(jc) == zero) return
                    jc = jc + n - info + 1
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           do j = 1,nrhs
              call la_dtpsv(uplo,trans,diag,n,ap,b(1,j),1)
           end do
           return
     end subroutine la_dtptrs
     !> QTPTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular matrix of order N stored in packed format,
     !> and B is an N-by-NRHS matrix.  A check is made to verify that A is
     !> nonsingular.

     pure subroutine la_qtptrs(uplo,trans,diag,n,nrhs,ap,b,ldb,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(in) :: ap(*)
           real(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('QTPTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 jc = 1
                 do info = 1,n
                    if (ap(jc + info - 1) == zero) return
                    jc = jc + info
                 end do
              else
                 jc = 1
                 do info = 1,n
                    if (ap(jc) == zero) return
                    jc = jc + n - info + 1
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           do j = 1,nrhs
              call la_qtpsv(uplo,trans,diag,n,ap,b(1,j),1)
           end do
           return
     end subroutine la_qtptrs

     !> STRRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular
     !> coefficient matrix.
     !> The solution matrix X must be computed by STRTRS or some other
     !> means before entering this routine.  STRRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_strrfs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,x,ldx,ferr,berr, &
                work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*),b(ldb,*),x(ldx,*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('STRRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_scopy(n,x(1,j),1,work(n + 1),1)
              call la_strmv(uplo,trans,diag,n,a,lda,work(n + 1),1)
              call la_saxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k - 1
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,n
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,n
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = 1,k - 1
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,n
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_slacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_slacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_strsv(uplo,transt,diag,n,a,lda,work(n + 1),1)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_strsv(uplo,trans,diag,n,a,lda,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_strrfs
     !> DTRRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular
     !> coefficient matrix.
     !> The solution matrix X must be computed by DTRTRS or some other
     !> means before entering this routine.  DTRRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_dtrrfs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,x,ldx,ferr,berr, &
                work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*),b(ldb,*),x(ldx,*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('DTRRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_dcopy(n,x(1,j),1,work(n + 1),1)
              call la_dtrmv(uplo,trans,diag,n,a,lda,work(n + 1),1)
              call la_daxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k - 1
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,n
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,n
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = 1,k - 1
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,n
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_dlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_dlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_dtrsv(uplo,transt,diag,n,a,lda,work(n + 1),1)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_dtrsv(uplo,trans,diag,n,a,lda,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_dtrrfs
     !> QTRRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular
     !> coefficient matrix.
     !> The solution matrix X must be computed by QTRTRS or some other
     !> means before entering this routine.  QTRRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_qtrrfs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,x,ldx,ferr,berr, &
                work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*),b(ldb,*),x(ldx,*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transt
           integer(ilp) :: i,j,k,kase,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('QTRRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transt = 'T'
           else
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a or a**t, depending on trans.
              call la_qcopy(n,x(1,j),1,work(n + 1),1)
              call la_qtrmv(uplo,trans,diag,n,a,lda,work(n + 1),1)
              call la_qaxpy(n,-one,b(1,j),1,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = 1,k - 1
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k,n
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = abs(x(k,j))
                          do i = k + 1,n
                             work(i) = work(i) + abs(a(i,k))*xk
                          end do
                          work(k) = work(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**t)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = 1,k - 1
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    else
                       do k = 1,n
                          s = abs(x(k,j))
                          do i = k + 1,n
                             s = s + abs(a(i,k))*abs(x(i,j))
                          end do
                          work(k) = work(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (work(i) > safe2) then
                    s = max(s,abs(work(n + i))/work(i))
                 else
                    s = max(s, (abs(work(n + i)) + safe1)/(work(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_qlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_qlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**t).
                    call la_qtrsv(uplo,transt,diag,n,a,lda,work(n + 1),1)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_qtrsv(uplo,trans,diag,n,a,lda,work(n + 1),1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_qtrrfs

     !> STRTI2: computes the inverse of a real upper or lower triangular
     !> matrix.
     !> This is the Level 2 BLAS version of the algorithm.

     pure subroutine la_strti2(uplo,diag,n,a,lda,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           real(sp) :: ajj
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('STRTI2',-info)
              return
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              do j = 1,n
                 if (nounit) then
                    a(j,j) = one/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -one
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_strmv('UPPER','NO TRANSPOSE',diag,j - 1,a,lda,a(1,j),1)

                 call la_sscal(j - 1,ajj,a(1,j),1)
              end do
           else
              ! compute inverse of lower triangular matrix.
              do j = n,1,-1
                 if (nounit) then
                    a(j,j) = one/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -one
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_strmv('LOWER','NO TRANSPOSE',diag,n - j,a(j + 1,j + 1),lda,a( &
                              j + 1,j),1)
                    call la_sscal(n - j,ajj,a(j + 1,j),1)
                 end if
              end do
           end if
           return
     end subroutine la_strti2
     !> DTRTI2: computes the inverse of a real upper or lower triangular
     !> matrix.
     !> This is the Level 2 BLAS version of the algorithm.

     pure subroutine la_dtrti2(uplo,diag,n,a,lda,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           real(dp) :: ajj
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DTRTI2',-info)
              return
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              do j = 1,n
                 if (nounit) then
                    a(j,j) = one/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -one
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_dtrmv('UPPER','NO TRANSPOSE',diag,j - 1,a,lda,a(1,j),1)

                 call la_dscal(j - 1,ajj,a(1,j),1)
              end do
           else
              ! compute inverse of lower triangular matrix.
              do j = n,1,-1
                 if (nounit) then
                    a(j,j) = one/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -one
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_dtrmv('LOWER','NO TRANSPOSE',diag,n - j,a(j + 1,j + 1),lda,a( &
                              j + 1,j),1)
                    call la_dscal(n - j,ajj,a(j + 1,j),1)
                 end if
              end do
           end if
           return
     end subroutine la_dtrti2
     !> QTRTI2: computes the inverse of a real upper or lower triangular
     !> matrix.
     !> This is the Level 2 BLAS version of the algorithm.

     pure subroutine la_qtrti2(uplo,diag,n,a,lda,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           real(qp) :: ajj
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QTRTI2',-info)
              return
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              do j = 1,n
                 if (nounit) then
                    a(j,j) = one/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -one
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_qtrmv('UPPER','NO TRANSPOSE',diag,j - 1,a,lda,a(1,j),1)

                 call la_qscal(j - 1,ajj,a(1,j),1)
              end do
           else
              ! compute inverse of lower triangular matrix.
              do j = n,1,-1
                 if (nounit) then
                    a(j,j) = one/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -one
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_qtrmv('LOWER','NO TRANSPOSE',diag,n - j,a(j + 1,j + 1),lda,a( &
                              j + 1,j),1)
                    call la_qscal(n - j,ajj,a(j + 1,j),1)
                 end if
              end do
           end if
           return
     end subroutine la_qtrti2

     !> STRTRI: computes the inverse of a real upper or lower triangular
     !> matrix A.
     !> This is the Level 3 BLAS version of the algorithm.

     pure subroutine la_strtri(uplo,diag,n,a,lda,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jb,nb,nn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('STRTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity if non-unit.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == zero) return
              end do
              info = 0
           end if
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'STRTRI',uplo//diag,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_strti2(uplo,diag,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute inverse of upper triangular matrix
                 do j = 1,n,nb
                    jb = min(nb,n - j + 1)
                    ! compute rows 1:j-1 of current block column
                    call la_strmm('LEFT','UPPER','NO TRANSPOSE',diag,j - 1,jb,one,a,lda, &
                               a(1,j),lda)
                    call la_strsm('RIGHT','UPPER','NO TRANSPOSE',diag,j - 1,jb,-one,a(j, &
                               j),lda,a(1,j),lda)
                    ! compute inverse of current diagonal block
                    call la_strti2('UPPER',diag,jb,a(j,j),lda,info)
                 end do
              else
                 ! compute inverse of lower triangular matrix
                 nn = ((n - 1)/nb)*nb + 1
                 do j = nn,1,-nb
                    jb = min(nb,n - j + 1)
                    if (j + jb <= n) then
                       ! compute rows j+jb:n of current block column
                       call la_strmm('LEFT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,one, &
                                  a(j + jb,j + jb),lda,a(j + jb,j),lda)
                       call la_strsm('RIGHT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,- &
                                 one,a(j,j),lda,a(j + jb,j),lda)
                    end if
                    ! compute inverse of current diagonal block
                    call la_strti2('LOWER',diag,jb,a(j,j),lda,info)
                 end do
              end if
           end if
           return
     end subroutine la_strtri
     !> DTRTRI: computes the inverse of a real upper or lower triangular
     !> matrix A.
     !> This is the Level 3 BLAS version of the algorithm.

     pure subroutine la_dtrtri(uplo,diag,n,a,lda,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jb,nb,nn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DTRTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity if non-unit.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == zero) return
              end do
              info = 0
           end if
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'DTRTRI',uplo//diag,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_dtrti2(uplo,diag,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute inverse of upper triangular matrix
                 do j = 1,n,nb
                    jb = min(nb,n - j + 1)
                    ! compute rows 1:j-1 of current block column
                    call la_dtrmm('LEFT','UPPER','NO TRANSPOSE',diag,j - 1,jb,one,a,lda, &
                               a(1,j),lda)
                    call la_dtrsm('RIGHT','UPPER','NO TRANSPOSE',diag,j - 1,jb,-one,a(j, &
                               j),lda,a(1,j),lda)
                    ! compute inverse of current diagonal block
                    call la_dtrti2('UPPER',diag,jb,a(j,j),lda,info)
                 end do
              else
                 ! compute inverse of lower triangular matrix
                 nn = ((n - 1)/nb)*nb + 1
                 do j = nn,1,-nb
                    jb = min(nb,n - j + 1)
                    if (j + jb <= n) then
                       ! compute rows j+jb:n of current block column
                       call la_dtrmm('LEFT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,one, &
                                  a(j + jb,j + jb),lda,a(j + jb,j),lda)
                       call la_dtrsm('RIGHT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,- &
                                 one,a(j,j),lda,a(j + jb,j),lda)
                    end if
                    ! compute inverse of current diagonal block
                    call la_dtrti2('LOWER',diag,jb,a(j,j),lda,info)
                 end do
              end if
           end if
           return
     end subroutine la_dtrtri
     !> QTRTRI: computes the inverse of a real upper or lower triangular
     !> matrix A.
     !> This is the Level 3 BLAS version of the algorithm.

     pure subroutine la_qtrtri(uplo,diag,n,a,lda,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jb,nb,nn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QTRTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity if non-unit.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == zero) return
              end do
              info = 0
           end if
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'QTRTRI',uplo//diag,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_qtrti2(uplo,diag,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute inverse of upper triangular matrix
                 do j = 1,n,nb
                    jb = min(nb,n - j + 1)
                    ! compute rows 1:j-1 of current block column
                    call la_qtrmm('LEFT','UPPER','NO TRANSPOSE',diag,j - 1,jb,one,a,lda, &
                               a(1,j),lda)
                    call la_qtrsm('RIGHT','UPPER','NO TRANSPOSE',diag,j - 1,jb,-one,a(j, &
                               j),lda,a(1,j),lda)
                    ! compute inverse of current diagonal block
                    call la_qtrti2('UPPER',diag,jb,a(j,j),lda,info)
                 end do
              else
                 ! compute inverse of lower triangular matrix
                 nn = ((n - 1)/nb)*nb + 1
                 do j = nn,1,-nb
                    jb = min(nb,n - j + 1)
                    if (j + jb <= n) then
                       ! compute rows j+jb:n of current block column
                       call la_qtrmm('LEFT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,one, &
                                  a(j + jb,j + jb),lda,a(j + jb,j),lda)
                       call la_qtrsm('RIGHT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,- &
                                 one,a(j,j),lda,a(j + jb,j),lda)
                    end if
                    ! compute inverse of current diagonal block
                    call la_qtrti2('LOWER',diag,jb,a(j,j),lda,info)
                 end do
              end if
           end if
           return
     end subroutine la_qtrtri

     !> STRTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular matrix of order N, and B is an N-by-NRHS
     !> matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_strtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('STRTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == zero) return
              end do
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           call la_strsm('LEFT',uplo,trans,diag,n,nrhs,one,a,lda,b,ldb)
           return
     end subroutine la_strtrs
     !> DTRTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular matrix of order N, and B is an N-by-NRHS
     !> matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_dtrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DTRTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == zero) return
              end do
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           call la_dtrsm('LEFT',uplo,trans,diag,n,nrhs,one,a,lda,b,ldb)
           return
     end subroutine la_dtrtrs
     !> QTRTRS: solves a triangular system of the form
     !> A * X = B  or  A**T * X = B,
     !> where A is a triangular matrix of order N, and B is an N-by-NRHS
     !> matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_qtrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QTRTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == zero) return
              end do
           end if
           info = 0
           ! solve a * x = b  or  a**t * x = b.
           call la_qtrsm('LEFT',uplo,trans,diag,n,nrhs,one,a,lda,b,ldb)
           return
     end subroutine la_qtrtrs

     !> STBCON: estimates the reciprocal of the condition number of a
     !> triangular band matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_stbcon(norm,uplo,diag,n,kd,ab,ldab,rcond,work,iwork,info)
        use la_constants_sp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: ab(ldab,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(sp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('STBCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_slamch('SAFE MINIMUM')*real(max(1,n),KIND=sp)
           ! compute the norm of the triangular matrix a.
           anorm = la_slantb(norm,uplo,diag,n,kd,ab,ldab,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_slacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_slatbs(uplo,'NO TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_slatbs(uplo,'TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,work(2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_isamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_srscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_stbcon
     !> DTBCON: estimates the reciprocal of the condition number of a
     !> triangular band matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_dtbcon(norm,uplo,diag,n,kd,ab,ldab,rcond,work,iwork,info)
        use la_constants_dp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: ab(ldab,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(dp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DTBCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_dlamch('SAFE MINIMUM')*real(max(1,n),KIND=dp)
           ! compute the norm of the triangular matrix a.
           anorm = la_dlantb(norm,uplo,diag,n,kd,ab,ldab,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_dlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_dlatbs(uplo,'NO TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_dlatbs(uplo,'TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,work(2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_idamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_drscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_dtbcon
     !> QTBCON: estimates the reciprocal of the condition number of a
     !> triangular band matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_qtbcon(norm,uplo,diag,n,kd,ab,ldab,rcond,work,iwork,info)
        use la_constants_qp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: ab(ldab,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(qp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QTBCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_qlamch('SAFE MINIMUM')*real(max(1,n),KIND=qp)
           ! compute the norm of the triangular matrix a.
           anorm = la_qlantb(norm,uplo,diag,n,kd,ab,ldab,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_qlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_qlatbs(uplo,'NO TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_qlatbs(uplo,'TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,work(2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_iqamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_qrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_qtbcon

     !> STFTRI: computes the inverse of a triangular matrix A stored in RFP
     !> format.
     !> This is a Level 3 BLAS version of the algorithm.

     pure subroutine la_stftri(transr,uplo,diag,n,a,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo,diag
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: a(0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('STFTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
           else
              nisodd = .true.
           end if
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1)
                    call la_strtri('L',diag,n1,a(0),n,info)
                    if (info > 0) return
                    call la_strmm('R','L','N',diag,n2,n1,-one,a(0),n,a(n1),n)

                    call la_strtri('U',diag,n2,a(n),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_strmm('L','U','T',diag,n2,n1,one,a(n),n,a(n1),n)

                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    call la_strtri('L',diag,n1,a(n2),n,info)
                    if (info > 0) return
                    call la_strmm('L','L','T',diag,n1,n2,-one,a(n2),n,a(0),n)

                    call la_strtri('U',diag,n2,a(n1),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_strmm('R','U','N',diag,n1,n2,one,a(n1),n,a(0),n)

                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0), t2 -> a(1), s -> a(0+n1*n1)
                    call la_strtri('U',diag,n1,a(0),n1,info)
                    if (info > 0) return
                    call la_strmm('L','U','N',diag,n1,n2,-one,a(0),n1,a(n1*n1), &
                              n1)
                    call la_strtri('L',diag,n2,a(1),n1,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_strmm('R','L','T',diag,n1,n2,one,a(1),n1,a(n1*n1), &
                              n1)
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0+n2*n2), t2 -> a(0+n1*n2), s -> a(0)
                    call la_strtri('U',diag,n1,a(n2*n2),n2,info)
                    if (info > 0) return
                    call la_strmm('R','U','T',diag,n2,n1,-one,a(n2*n2),n2,a(0), &
                              n2)
                    call la_strtri('L',diag,n2,a(n1*n2),n2,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_strmm('L','L','N',diag,n2,n1,one,a(n1*n2),n2,a(0), &
                              n2)
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    call la_strtri('L',diag,k,a(1),n + 1,info)
                    if (info > 0) return
                    call la_strmm('R','L','N',diag,k,k,-one,a(1),n + 1,a(k + 1),n + 1 &
                              )
                    call la_strtri('U',diag,k,a(0),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_strmm('L','U','T',diag,k,k,one,a(0),n + 1,a(k + 1),n + 1)

                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    call la_strtri('L',diag,k,a(k + 1),n + 1,info)
                    if (info > 0) return
                    call la_strmm('L','L','T',diag,k,k,-one,a(k + 1),n + 1,a(0),n + 1 &
                              )
                    call la_strtri('U',diag,k,a(k),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_strmm('R','U','N',diag,k,k,one,a(k),n + 1,a(0),n + 1)

                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    call la_strtri('U',diag,k,a(k),k,info)
                    if (info > 0) return
                    call la_strmm('L','U','N',diag,k,k,-one,a(k),k,a(k*(k + 1)), &
                              k)
                    call la_strtri('L',diag,k,a(0),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_strmm('R','L','T',diag,k,k,one,a(0),k,a(k*(k + 1)), &
                              k)
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    call la_strtri('U',diag,k,a(k*(k + 1)),k,info)
                    if (info > 0) return
                    call la_strmm('R','U','T',diag,k,k,-one,a(k*(k + 1)),k,a(0), &
                              k)
                    call la_strtri('L',diag,k,a(k*k),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_strmm('L','L','N',diag,k,k,one,a(k*k),k,a(0),k)

                 end if
              end if
           end if
           return
     end subroutine la_stftri
     !> DTFTRI: computes the inverse of a triangular matrix A stored in RFP
     !> format.
     !> This is a Level 3 BLAS version of the algorithm.

     pure subroutine la_dtftri(transr,uplo,diag,n,a,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo,diag
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: a(0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DTFTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
           else
              nisodd = .true.
           end if
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1)
                    call la_dtrtri('L',diag,n1,a(0),n,info)
                    if (info > 0) return
                    call la_dtrmm('R','L','N',diag,n2,n1,-one,a(0),n,a(n1),n)

                    call la_dtrtri('U',diag,n2,a(n),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_dtrmm('L','U','T',diag,n2,n1,one,a(n),n,a(n1),n)

                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    call la_dtrtri('L',diag,n1,a(n2),n,info)
                    if (info > 0) return
                    call la_dtrmm('L','L','T',diag,n1,n2,-one,a(n2),n,a(0),n)

                    call la_dtrtri('U',diag,n2,a(n1),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_dtrmm('R','U','N',diag,n1,n2,one,a(n1),n,a(0),n)

                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0), t2 -> a(1), s -> a(0+n1*n1)
                    call la_dtrtri('U',diag,n1,a(0),n1,info)
                    if (info > 0) return
                    call la_dtrmm('L','U','N',diag,n1,n2,-one,a(0),n1,a(n1*n1), &
                              n1)
                    call la_dtrtri('L',diag,n2,a(1),n1,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_dtrmm('R','L','T',diag,n1,n2,one,a(1),n1,a(n1*n1), &
                              n1)
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0+n2*n2), t2 -> a(0+n1*n2), s -> a(0)
                    call la_dtrtri('U',diag,n1,a(n2*n2),n2,info)
                    if (info > 0) return
                    call la_dtrmm('R','U','T',diag,n2,n1,-one,a(n2*n2),n2,a(0), &
                              n2)
                    call la_dtrtri('L',diag,n2,a(n1*n2),n2,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_dtrmm('L','L','N',diag,n2,n1,one,a(n1*n2),n2,a(0), &
                              n2)
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    call la_dtrtri('L',diag,k,a(1),n + 1,info)
                    if (info > 0) return
                    call la_dtrmm('R','L','N',diag,k,k,-one,a(1),n + 1,a(k + 1),n + 1 &
                              )
                    call la_dtrtri('U',diag,k,a(0),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_dtrmm('L','U','T',diag,k,k,one,a(0),n + 1,a(k + 1),n + 1)

                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    call la_dtrtri('L',diag,k,a(k + 1),n + 1,info)
                    if (info > 0) return
                    call la_dtrmm('L','L','T',diag,k,k,-one,a(k + 1),n + 1,a(0),n + 1 &
                              )
                    call la_dtrtri('U',diag,k,a(k),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_dtrmm('R','U','N',diag,k,k,one,a(k),n + 1,a(0),n + 1)

                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    call la_dtrtri('U',diag,k,a(k),k,info)
                    if (info > 0) return
                    call la_dtrmm('L','U','N',diag,k,k,-one,a(k),k,a(k*(k + 1)), &
                              k)
                    call la_dtrtri('L',diag,k,a(0),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_dtrmm('R','L','T',diag,k,k,one,a(0),k,a(k*(k + 1)), &
                              k)
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    call la_dtrtri('U',diag,k,a(k*(k + 1)),k,info)
                    if (info > 0) return
                    call la_dtrmm('R','U','T',diag,k,k,-one,a(k*(k + 1)),k,a(0), &
                              k)
                    call la_dtrtri('L',diag,k,a(k*k),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_dtrmm('L','L','N',diag,k,k,one,a(k*k),k,a(0),k)

                 end if
              end if
           end if
           return
     end subroutine la_dtftri
     !> QTFTRI: computes the inverse of a triangular matrix A stored in RFP
     !> format.
     !> This is a Level 3 BLAS version of the algorithm.

     pure subroutine la_qtftri(transr,uplo,diag,n,a,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo,diag
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: a(0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'T')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QTFTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
           else
              nisodd = .true.
           end if
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1)
                    call la_qtrtri('L',diag,n1,a(0),n,info)
                    if (info > 0) return
                    call la_qtrmm('R','L','N',diag,n2,n1,-one,a(0),n,a(n1),n)

                    call la_qtrtri('U',diag,n2,a(n),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_qtrmm('L','U','T',diag,n2,n1,one,a(n),n,a(n1),n)

                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    call la_qtrtri('L',diag,n1,a(n2),n,info)
                    if (info > 0) return
                    call la_qtrmm('L','L','T',diag,n1,n2,-one,a(n2),n,a(0),n)

                    call la_qtrtri('U',diag,n2,a(n1),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_qtrmm('R','U','N',diag,n1,n2,one,a(n1),n,a(0),n)

                 end if
              else
                 ! n is odd and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0), t2 -> a(1), s -> a(0+n1*n1)
                    call la_qtrtri('U',diag,n1,a(0),n1,info)
                    if (info > 0) return
                    call la_qtrmm('L','U','N',diag,n1,n2,-one,a(0),n1,a(n1*n1), &
                              n1)
                    call la_qtrtri('L',diag,n2,a(1),n1,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_qtrmm('R','L','T',diag,n1,n2,one,a(1),n1,a(n1*n1), &
                              n1)
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0+n2*n2), t2 -> a(0+n1*n2), s -> a(0)
                    call la_qtrtri('U',diag,n1,a(n2*n2),n2,info)
                    if (info > 0) return
                    call la_qtrmm('R','U','T',diag,n2,n1,-one,a(n2*n2),n2,a(0), &
                              n2)
                    call la_qtrtri('L',diag,n2,a(n1*n2),n2,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_qtrmm('L','L','N',diag,n2,n1,one,a(n1*n2),n2,a(0), &
                              n2)
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    call la_qtrtri('L',diag,k,a(1),n + 1,info)
                    if (info > 0) return
                    call la_qtrmm('R','L','N',diag,k,k,-one,a(1),n + 1,a(k + 1),n + 1 &
                              )
                    call la_qtrtri('U',diag,k,a(0),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_qtrmm('L','U','T',diag,k,k,one,a(0),n + 1,a(k + 1),n + 1)

                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    call la_qtrtri('L',diag,k,a(k + 1),n + 1,info)
                    if (info > 0) return
                    call la_qtrmm('L','L','T',diag,k,k,-one,a(k + 1),n + 1,a(0),n + 1 &
                              )
                    call la_qtrtri('U',diag,k,a(k),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_qtrmm('R','U','N',diag,k,k,one,a(k),n + 1,a(0),n + 1)

                 end if
              else
                 ! n is even and transr = 't'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    call la_qtrtri('U',diag,k,a(k),k,info)
                    if (info > 0) return
                    call la_qtrmm('L','U','N',diag,k,k,-one,a(k),k,a(k*(k + 1)), &
                              k)
                    call la_qtrtri('L',diag,k,a(0),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_qtrmm('R','L','T',diag,k,k,one,a(0),k,a(k*(k + 1)), &
                              k)
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    call la_qtrtri('U',diag,k,a(k*(k + 1)),k,info)
                    if (info > 0) return
                    call la_qtrmm('R','U','T',diag,k,k,-one,a(k*(k + 1)),k,a(0), &
                              k)
                    call la_qtrtri('L',diag,k,a(k*k),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_qtrmm('L','L','N',diag,k,k,one,a(k*k),k,a(0),k)

                 end if
              end if
           end if
           return
     end subroutine la_qtftri

     !> STPCON: estimates the reciprocal of the condition number of a packed
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_stpcon(norm,uplo,diag,n,ap,rcond,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: ap(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(sp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('STPCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_slamch('SAFE MINIMUM')*real(max(1,n),KIND=sp)
           ! compute the norm of the triangular matrix a.
           anorm = la_slantp(norm,uplo,diag,n,ap,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_slacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_slatps(uplo,'NO TRANSPOSE',diag,normin,n,ap,work,scale, &
                              work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_slatps(uplo,'TRANSPOSE',diag,normin,n,ap,work,scale,work( &
                              2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_isamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_srscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_stpcon
     !> DTPCON: estimates the reciprocal of the condition number of a packed
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_dtpcon(norm,uplo,diag,n,ap,rcond,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: ap(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(dp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DTPCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_dlamch('SAFE MINIMUM')*real(max(1,n),KIND=dp)
           ! compute the norm of the triangular matrix a.
           anorm = la_dlantp(norm,uplo,diag,n,ap,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_dlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_dlatps(uplo,'NO TRANSPOSE',diag,normin,n,ap,work,scale, &
                              work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_dlatps(uplo,'TRANSPOSE',diag,normin,n,ap,work,scale,work( &
                              2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_idamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_drscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_dtpcon
     !> QTPCON: estimates the reciprocal of the condition number of a packed
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_qtpcon(norm,uplo,diag,n,ap,rcond,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: ap(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(qp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QTPCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_qlamch('SAFE MINIMUM')*real(max(1,n),KIND=qp)
           ! compute the norm of the triangular matrix a.
           anorm = la_qlantp(norm,uplo,diag,n,ap,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_qlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_qlatps(uplo,'NO TRANSPOSE',diag,normin,n,ap,work,scale, &
                              work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_qlatps(uplo,'TRANSPOSE',diag,normin,n,ap,work,scale,work( &
                              2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_iqamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_qrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_qtpcon

     !> STRCON: estimates the reciprocal of the condition number of a
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_strcon(norm,uplo,diag,n,a,lda,rcond,work,iwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(sp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max,real
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('STRCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_slamch('SAFE MINIMUM')*real(max(1,n),KIND=sp)
           ! compute the norm of the triangular matrix a.
           anorm = la_slantr(norm,uplo,diag,n,n,a,lda,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_slacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_slatrs(uplo,'NO TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                               work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_slatrs(uplo,'TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                              work(2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_isamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_srscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_strcon
     !> DTRCON: estimates the reciprocal of the condition number of a
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_dtrcon(norm,uplo,diag,n,a,lda,rcond,work,iwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(dp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DTRCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_dlamch('SAFE MINIMUM')*real(max(1,n),KIND=dp)
           ! compute the norm of the triangular matrix a.
           anorm = la_dlantr(norm,uplo,diag,n,n,a,lda,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_dlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_dlatrs(uplo,'NO TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                               work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_dlatrs(uplo,'TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                              work(2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_idamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_drscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_dtrcon
     !> QTRCON: estimates the reciprocal of the condition number of a
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_qtrcon(norm,uplo,diag,n,a,lda,rcond,work,iwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(qp) :: ainvnm,anorm,scale,smlnum,xnorm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QTRCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_qlamch('SAFE MINIMUM')*real(max(1,n),KIND=qp)
           ! compute the norm of the triangular matrix a.
           anorm = la_qlantr(norm,uplo,diag,n,n,a,lda,work)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_qlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_qlatrs(uplo,'NO TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                               work(2*n + 1),info)
                 else
                    ! multiply by inv(a**t).
                    call la_qlatrs(uplo,'TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                              work(2*n + 1),info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_iqamax(n,work,1)
                    xnorm = abs(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_qrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_qtrcon

     !> CLATBS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular band matrix.  Here A**T denotes the transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine CTBSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_clatbs(uplo,trans,diag,normin,n,kd,ab,ldab,x,scale,cnorm, &
               info)
        use la_constants_sp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(inout) :: cnorm(*)
           complex(sp),intent(in) :: ab(ldab,*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast,jlen,maind
           real(sp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(sp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,min,real
           ! Statement Functions
           real(sp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=sp)/2.) + abs(aimag(zdum)/2.)
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (kd < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CLATBS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_slamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = smlnum/la_slamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    jlen = min(kd,j - 1)
                    cnorm(j) = la_scasum(jlen,ab(kd + 1 - jlen,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n
                    jlen = min(kd,n - j)
                    if (jlen > 0) then
                       cnorm(j) = la_scasum(jlen,ab(2,j),1)
                    else
                       cnorm(j) = zero
                    end if
                 end do
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_isamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_sscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_ctbsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = kd + 1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = ab(maind,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = kd + 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = ab(maind,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_ctbsv(uplo,trans,diag,n,kd,ab,ldab,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_csscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = ab(maind,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 105
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_cladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_cladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    105 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_csscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                                   ! x(j)* a(max(1,j-kd):j-1,j)
                          jlen = min(kd,j - 1)
                          call la_caxpy(jlen,-x(j)*tscal,ab(kd + 1 - jlen,j),1,x(j - jlen &
                                    ),1)
                          i = la_icamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                    else if (j < n) then
                       ! compute the update
                          ! x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                                ! x(j) * a(j+1:min(j+kd,n),j)
                       jlen = min(kd,n - j)
                       if (jlen > 0) call la_caxpy(jlen,-x(j)*tscal,ab(2,j),1,x(j + 1), &
                                  1)
                       i = j + la_icamax(n - j,x(j + 1),1)
                       xmax = cabs1(x(i))
                    end if
                 end do loop_110
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 loop_150: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_cladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=sp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_cdotu to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          csumj = la_cdotu(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)

                       else
                          jlen = min(kd,n - j)
                          if (jlen > 1) csumj = la_cdotu(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             csumj = csumj + (ab(kd + i - jlen,j)*uscal)*x(j - jlen - 1 + i)
                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             csumj = csumj + (ab(i + 1,j)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=sp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 145
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_csscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       145 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_cladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_150
              else
                 ! solve a**h * x = b
                 loop_190: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(ab(maind,j))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_cladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=sp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_cdotc to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          csumj = la_cdotc(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)

                       else
                          jlen = min(kd,n - j)
                          if (jlen > 1) csumj = la_cdotc(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             csumj = csumj + (conjg(ab(kd + i - jlen,j))*uscal)*x(j - jlen - 1 + i)

                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             csumj = csumj + (conjg(ab(i + 1,j))*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=sp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = conjg(ab(maind,j))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 185
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_csscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       185 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_cladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_190
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_sscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_clatbs
     !> ZLATBS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular band matrix.  Here A**T denotes the transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine ZTBSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_zlatbs(uplo,trans,diag,normin,n,kd,ab,ldab,x,scale,cnorm, &
               info)
        use la_constants_dp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(inout) :: cnorm(*)
           complex(dp),intent(in) :: ab(ldab,*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast,jlen,maind
           real(dp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(dp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(dp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=dp)/2._dp) + abs(aimag(zdum)/2._dp)

           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (kd < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZLATBS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_dlamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = smlnum/la_dlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    jlen = min(kd,j - 1)
                    cnorm(j) = la_dzasum(jlen,ab(kd + 1 - jlen,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n
                    jlen = min(kd,n - j)
                    if (jlen > 0) then
                       cnorm(j) = la_dzasum(jlen,ab(2,j),1)
                    else
                       cnorm(j) = zero
                    end if
                 end do
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_idamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_dscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_ztbsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = kd + 1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = ab(maind,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = kd + 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = ab(maind,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_ztbsv(uplo,trans,diag,n,kd,ab,ldab,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_zdscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_120: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = ab(maind,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 110
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_zladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_zladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    110 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_zdscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                                   ! x(j)* a(max(1,j-kd):j-1,j)
                          jlen = min(kd,j - 1)
                          call la_zaxpy(jlen,-x(j)*tscal,ab(kd + 1 - jlen,j),1,x(j - jlen &
                                    ),1)
                          i = la_izamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                    else if (j < n) then
                       ! compute the update
                          ! x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                                ! x(j) * a(j+1:min(j+kd,n),j)
                       jlen = min(kd,n - j)
                       if (jlen > 0) call la_zaxpy(jlen,-x(j)*tscal,ab(2,j),1,x(j + 1), &
                                  1)
                       i = j + la_izamax(n - j,x(j + 1),1)
                       xmax = cabs1(x(i))
                    end if
                 end do loop_120
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 loop_170: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_zladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=dp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_zdotu to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          csumj = la_zdotu(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)

                       else
                          jlen = min(kd,n - j)
                          if (jlen > 1) csumj = la_zdotu(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             csumj = csumj + (ab(kd + i - jlen,j)*uscal)*x(j - jlen - 1 + i)
                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             csumj = csumj + (ab(i + 1,j)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=dp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 160
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_zdscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       160 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_zladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_170
              else
                 ! solve a**h * x = b
                 loop_220: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(ab(maind,j))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_zladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=dp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_zdotc to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          csumj = la_zdotc(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)

                       else
                          jlen = min(kd,n - j)
                          if (jlen > 1) csumj = la_zdotc(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             csumj = csumj + (conjg(ab(kd + i - jlen,j))*uscal)*x(j - jlen - 1 + i)

                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             csumj = csumj + (conjg(ab(i + 1,j))*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=dp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = conjg(ab(maind,j))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 210
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_zdscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       210 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_zladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_220
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_dscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_zlatbs
     !> WLATBS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular band matrix.  Here A**T denotes the transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine WTBSV is called.  If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_wlatbs(uplo,trans,diag,normin,n,kd,ab,ldab,x,scale,cnorm, &
               info)
        use la_constants_qp,only:zero,half,one,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(inout) :: cnorm(*)
           complex(qp),intent(in) :: ab(ldab,*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast,jlen,maind
           real(qp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(qp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(qp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=qp)/2._qp) + abs(aimag(zdum)/2._qp)

           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (kd < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WLATBS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_qlamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = smlnum/la_qlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    jlen = min(kd,j - 1)
                    cnorm(j) = la_qwasum(jlen,ab(kd + 1 - jlen,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n
                    jlen = min(kd,n - j)
                    if (jlen > 0) then
                       cnorm(j) = la_qwasum(jlen,ab(2,j),1)
                    else
                       cnorm(j) = zero
                    end if
                 end do
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_iqamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_qscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_wtbsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = kd + 1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = ab(maind,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
                 maind = kd + 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
                 maind = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = ab(maind,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_wtbsv(uplo,trans,diag,n,kd,ab,ldab,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_wqscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_120: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = ab(maind,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 110
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_wladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_wladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    110 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_wqscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(max(1,j-kd):j-1) := x(max(1,j-kd):j-1) -
                                                   ! x(j)* a(max(1,j-kd):j-1,j)
                          jlen = min(kd,j - 1)
                          call la_waxpy(jlen,-x(j)*tscal,ab(kd + 1 - jlen,j),1,x(j - jlen &
                                    ),1)
                          i = la_iwamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                    else if (j < n) then
                       ! compute the update
                          ! x(j+1:min(j+kd,n)) := x(j+1:min(j+kd,n)) -
                                                ! x(j) * a(j+1:min(j+kd,n),j)
                       jlen = min(kd,n - j)
                       if (jlen > 0) call la_waxpy(jlen,-x(j)*tscal,ab(2,j),1,x(j + 1), &
                                  1)
                       i = j + la_iwamax(n - j,x(j + 1),1)
                       xmax = cabs1(x(i))
                    end if
                 end do loop_120
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 loop_170: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_wladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=qp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_wdotu to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          csumj = la_wdotu(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)

                       else
                          jlen = min(kd,n - j)
                          if (jlen > 1) csumj = la_wdotu(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             csumj = csumj + (ab(kd + i - jlen,j)*uscal)*x(j - jlen - 1 + i)
                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             csumj = csumj + (ab(i + 1,j)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=qp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ab(maind,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 160
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_wqscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       160 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_wladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_170
              else
                 ! solve a**h * x = b
                 loop_220: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(ab(maind,j))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_wladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=qp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_wdotc to perform the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          csumj = la_wdotc(jlen,ab(kd + 1 - jlen,j),1,x(j - jlen),1)

                       else
                          jlen = min(kd,n - j)
                          if (jlen > 1) csumj = la_wdotc(jlen,ab(2,j),1,x(j + 1),1)

                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          jlen = min(kd,j - 1)
                          do i = 1,jlen
                             csumj = csumj + (conjg(ab(kd + i - jlen,j))*uscal)*x(j - jlen - 1 + i)

                          end do
                       else
                          jlen = min(kd,n - j)
                          do i = 1,jlen
                             csumj = csumj + (conjg(ab(i + 1,j))*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=qp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = conjg(ab(maind,j))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 210
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_wqscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       210 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_wladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_220
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_qscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_wlatbs

     !> CLATPS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular matrix stored in packed form.  Here A**T denotes the
     !> transpose of A, A**H denotes the conjugate transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine CTPSV is called. If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_clatps(uplo,trans,diag,normin,n,ap,x,scale,cnorm,info)
        use la_constants_sp,only:zero,half,one,two

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(inout) :: cnorm(*)
           complex(sp),intent(in) :: ap(*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,ip,j,jfirst,jinc,jlast,jlen
           real(sp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(sp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,min,real
           ! Statement Functions
           real(sp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=sp)/2.) + abs(aimag(zdum)/2.)
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CLATPS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_slamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = smlnum/la_slamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 ip = 1
                 do j = 1,n
                    cnorm(j) = la_scasum(j - 1,ap(ip),1)
                    ip = ip + j
                 end do
              else
                 ! a is lower triangular.
                 ip = 1
                 do j = 1,n - 1
                    cnorm(j) = la_scasum(n - j,ap(ip + 1),1)
                    ip = ip + n - j + 1
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_isamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_sscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_ctpsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = n
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = ap(ip)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                    ip = ip + jinc*jlen
                    jlen = jlen - 1
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = ap(ip)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_ctpsv(uplo,trans,diag,n,ap,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_csscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 ip = jfirst*(jfirst + 1)/2
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = ap(ip)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 105
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_cladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_cladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    105 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_csscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_caxpy(j - 1,-x(j)*tscal,ap(ip - j + 1),1,x,1)
                          i = la_icamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                       ip = ip - j
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_caxpy(n - j,-x(j)*tscal,ap(ip + 1),1,x(j + 1),1)

                          i = j + la_icamax(n - j,x(j + 1),1)
                          xmax = cabs1(x(i))
                       end if
                       ip = ip + n - j + 1
                    end if
                 end do loop_110
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_150: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_cladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=sp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_cdotu to perform the dot product.
                       if (upper) then
                          csumj = la_cdotu(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          csumj = la_cdotu(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (ap(ip - j + i)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             csumj = csumj + (ap(ip + i)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=sp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 145
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_csscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       145 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_cladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_150
              else
                 ! solve a**h * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_190: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(ap(ip))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_cladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=sp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_cdotc to perform the dot product.
                       if (upper) then
                          csumj = la_cdotc(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          csumj = la_cdotc(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (conjg(ap(ip - j + i))*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             csumj = csumj + (conjg(ap(ip + i))*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=sp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = conjg(ap(ip))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 185
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_csscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       185 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_cladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_190
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_sscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_clatps
     !> ZLATPS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular matrix stored in packed form.  Here A**T denotes the
     !> transpose of A, A**H denotes the conjugate transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine ZTPSV is called. If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_zlatps(uplo,trans,diag,normin,n,ap,x,scale,cnorm,info)
        use la_constants_dp,only:zero,half,one,two

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(inout) :: cnorm(*)
           complex(dp),intent(in) :: ap(*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,ip,j,jfirst,jinc,jlast,jlen
           real(dp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(dp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(dp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=dp)/2._dp) + abs(aimag(zdum)/2._dp)

           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZLATPS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_dlamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = smlnum/la_dlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 ip = 1
                 do j = 1,n
                    cnorm(j) = la_dzasum(j - 1,ap(ip),1)
                    ip = ip + j
                 end do
              else
                 ! a is lower triangular.
                 ip = 1
                 do j = 1,n - 1
                    cnorm(j) = la_dzasum(n - j,ap(ip + 1),1)
                    ip = ip + n - j + 1
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_idamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_dscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_ztpsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = n
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = ap(ip)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                    ip = ip + jinc*jlen
                    jlen = jlen - 1
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = ap(ip)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_ztpsv(uplo,trans,diag,n,ap,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_zdscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 ip = jfirst*(jfirst + 1)/2
                 loop_120: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = ap(ip)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 110
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_zladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_zladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    110 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_zdscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_zaxpy(j - 1,-x(j)*tscal,ap(ip - j + 1),1,x,1)
                          i = la_izamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                       ip = ip - j
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_zaxpy(n - j,-x(j)*tscal,ap(ip + 1),1,x(j + 1),1)

                          i = j + la_izamax(n - j,x(j + 1),1)
                          xmax = cabs1(x(i))
                       end if
                       ip = ip + n - j + 1
                    end if
                 end do loop_120
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_170: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_zladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=dp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_zdotu to perform the dot product.
                       if (upper) then
                          csumj = la_zdotu(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          csumj = la_zdotu(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (ap(ip - j + i)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             csumj = csumj + (ap(ip + i)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=dp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 160
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_zdscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       160 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_zladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_170
              else
                 ! solve a**h * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_220: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(ap(ip))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_zladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=dp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_zdotc to perform the dot product.
                       if (upper) then
                          csumj = la_zdotc(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          csumj = la_zdotc(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (conjg(ap(ip - j + i))*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             csumj = csumj + (conjg(ap(ip + i))*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=dp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = conjg(ap(ip))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 210
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_zdscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       210 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_zladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_220
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_dscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_zlatps
     !> WLATPS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow, where A is an upper or lower
     !> triangular matrix stored in packed form.  Here A**T denotes the
     !> transpose of A, A**H denotes the conjugate transpose of A, x and b
     !> are n-element vectors, and s is a scaling factor, usually less than
     !> or equal to 1, chosen so that the components of x will be less than
     !> the overflow threshold.  If the unscaled problem will not cause
     !> overflow, the Level 2 BLAS routine WTPSV is called. If the matrix A
     !> is singular (A(j,j) = 0 for some j), then s is set to 0 and a
     !> non-trivial solution to A*x = 0 is returned.

     pure subroutine la_wlatps(uplo,trans,diag,normin,n,ap,x,scale,cnorm,info)
        use la_constants_qp,only:zero,half,one,two

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(inout) :: cnorm(*)
           complex(qp),intent(in) :: ap(*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,ip,j,jfirst,jinc,jlast,jlen
           real(qp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(qp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(qp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=qp)/2._qp) + abs(aimag(zdum)/2._qp)

           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WLATPS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_qlamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = smlnum/la_qlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 ip = 1
                 do j = 1,n
                    cnorm(j) = la_qwasum(j - 1,ap(ip),1)
                    ip = ip + j
                 end do
              else
                 ! a is lower triangular.
                 ip = 1
                 do j = 1,n - 1
                    cnorm(j) = la_qwasum(n - j,ap(ip + 1),1)
                    ip = ip + n - j + 1
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_iqamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_qscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_wtpsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = n
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = ap(ip)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                    ip = ip + jinc*jlen
                    jlen = jlen - 1
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = ap(ip)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_wtpsv(uplo,trans,diag,n,ap,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_wqscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 ip = jfirst*(jfirst + 1)/2
                 loop_120: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = ap(ip)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 110
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_wladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_wladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    110 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_wqscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_waxpy(j - 1,-x(j)*tscal,ap(ip - j + 1),1,x,1)
                          i = la_iwamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                       ip = ip - j
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_waxpy(n - j,-x(j)*tscal,ap(ip + 1),1,x(j + 1),1)

                          i = j + la_iwamax(n - j,x(j + 1),1)
                          xmax = cabs1(x(i))
                       end if
                       ip = ip + n - j + 1
                    end if
                 end do loop_120
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_170: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_wladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=qp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_wdotu to perform the dot product.
                       if (upper) then
                          csumj = la_wdotu(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          csumj = la_wdotu(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (ap(ip - j + i)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             csumj = csumj + (ap(ip + i)*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=qp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = ap(ip)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 160
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_wqscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       160 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_wladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_170
              else
                 ! solve a**h * x = b
                 ip = jfirst*(jfirst + 1)/2
                 jlen = 1
                 loop_220: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(ap(ip))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_wladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=qp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_wdotc to perform the dot product.
                       if (upper) then
                          csumj = la_wdotc(j - 1,ap(ip - j + 1),1,x,1)
                       else if (j < n) then
                          csumj = la_wdotc(n - j,ap(ip + 1),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (conjg(ap(ip - j + i))*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = 1,n - j
                             csumj = csumj + (conjg(ap(ip + i))*uscal)*x(j + i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=qp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                          tjjs = conjg(ap(ip))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 210
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_wqscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       210 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_wladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                    jlen = jlen + 1
                    ip = ip + jinc*jlen
                 end do loop_220
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_qscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_wlatps

     !> CLATRS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow.  Here A is an upper or lower
     !> triangular matrix, A**T denotes the transpose of A, A**H denotes the
     !> conjugate transpose of A, x and b are n-element vectors, and s is a
     !> scaling factor, usually less than or equal to 1, chosen so that the
     !> components of x will be less than the overflow threshold.  If the
     !> unscaled problem will not cause overflow, the Level 2 BLAS routine
     !> CTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
     !> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.

     pure subroutine la_clatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
        use la_constants_sp,only:zero,half,one,two

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(out) :: scale
           ! Array Arguments
           real(sp),intent(inout) :: cnorm(*)
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast
           real(sp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(sp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,min,real
           ! Statement Functions
           real(sp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=sp)/2.) + abs(aimag(zdum)/2.)
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CLATRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_slamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           smlnum = smlnum/la_slamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    cnorm(j) = la_scasum(j - 1,a(1,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n - 1
                    cnorm(j) = la_scasum(n - j,a(j + 1,j),1)
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_isamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_sscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_ctrsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = a(j,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = a(j,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_ctrsv(uplo,trans,diag,n,a,lda,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_csscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_110: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = a(j,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 105
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_cladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_cladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    105 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_csscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_caxpy(j - 1,-x(j)*tscal,a(1,j),1,x,1)
                          i = la_icamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_caxpy(n - j,-x(j)*tscal,a(j + 1,j),1,x(j + 1),1)

                          i = j + la_icamax(n - j,x(j + 1),1)
                          xmax = cabs1(x(i))
                       end if
                    end if
                 end do loop_110
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 loop_150: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_cladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=sp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_cdotu to perform the dot product.
                       if (upper) then
                          csumj = la_cdotu(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          csumj = la_cdotu(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (a(i,j)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             csumj = csumj + (a(i,j)*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=sp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 145
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_csscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       145 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_cladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_150
              else
                 ! solve a**h * x = b
                 loop_190: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(a(j,j))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_cladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_csscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=sp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_cdotc to perform the dot product.
                       if (upper) then
                          csumj = la_cdotc(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          csumj = la_cdotc(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (conjg(a(i,j))*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             csumj = csumj + (conjg(a(i,j))*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=sp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          tjjs = conjg(a(j,j))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 185
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_csscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_csscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_cladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       185 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_cladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_190
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_sscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_clatrs
     !> ZLATRS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow.  Here A is an upper or lower
     !> triangular matrix, A**T denotes the transpose of A, A**H denotes the
     !> conjugate transpose of A, x and b are n-element vectors, and s is a
     !> scaling factor, usually less than or equal to 1, chosen so that the
     !> components of x will be less than the overflow threshold.  If the
     !> unscaled problem will not cause overflow, the Level 2 BLAS routine
     !> ZTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
     !> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.

     pure subroutine la_zlatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
        use la_constants_dp,only:zero,half,one,two

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(out) :: scale
           ! Array Arguments
           real(dp),intent(inout) :: cnorm(*)
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast
           real(dp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(dp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(dp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=dp)/2._dp) + abs(aimag(zdum)/2._dp)

           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZLATRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_dlamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           smlnum = smlnum/la_dlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    cnorm(j) = la_dzasum(j - 1,a(1,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n - 1
                    cnorm(j) = la_dzasum(n - j,a(j + 1,j),1)
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_idamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_dscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_ztrsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = a(j,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = a(j,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_ztrsv(uplo,trans,diag,n,a,lda,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_zdscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_120: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = a(j,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 110
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_zladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_zladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    110 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_zdscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_zaxpy(j - 1,-x(j)*tscal,a(1,j),1,x,1)
                          i = la_izamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_zaxpy(n - j,-x(j)*tscal,a(j + 1,j),1,x(j + 1),1)

                          i = j + la_izamax(n - j,x(j + 1),1)
                          xmax = cabs1(x(i))
                       end if
                    end if
                 end do loop_120
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 loop_170: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_zladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=dp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_zdotu to perform the dot product.
                       if (upper) then
                          csumj = la_zdotu(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          csumj = la_zdotu(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (a(i,j)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             csumj = csumj + (a(i,j)*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=dp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 160
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_zdscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       160 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_zladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_170
              else
                 ! solve a**h * x = b
                 loop_220: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(a(j,j))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_zladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_zdscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=dp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_zdotc to perform the dot product.
                       if (upper) then
                          csumj = la_zdotc(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          csumj = la_zdotc(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (conjg(a(i,j))*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             csumj = csumj + (conjg(a(i,j))*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=dp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          tjjs = conjg(a(j,j))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 210
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_zdscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_zdscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_zladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       210 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_zladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_220
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_dscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_zlatrs
     !> WLATRS: solves one of the triangular systems
     !> A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
     !> with scaling to prevent overflow.  Here A is an upper or lower
     !> triangular matrix, A**T denotes the transpose of A, A**H denotes the
     !> conjugate transpose of A, x and b are n-element vectors, and s is a
     !> scaling factor, usually less than or equal to 1, chosen so that the
     !> components of x will be less than the overflow threshold.  If the
     !> unscaled problem will not cause overflow, the Level 2 BLAS routine
     !> WTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
     !> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.

     pure subroutine la_wlatrs(uplo,trans,diag,normin,n,a,lda,x,scale,cnorm,info)
        use la_constants_qp,only:zero,half,one,two

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,normin,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(out) :: scale
           ! Array Arguments
           real(qp),intent(inout) :: cnorm(*)
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: x(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           integer(ilp) :: i,imax,j,jfirst,jinc,jlast
           real(qp) :: bignum,grow,rec,smlnum,tjj,tmax,tscal,xbnd,xj,xmax
           complex(qp) :: csumj,tjjs,uscal,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,min
           ! Statement Functions
           real(qp) :: cabs1,cabs2
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           cabs2(zdum) = abs(real(zdum,KIND=qp)/2._qp) + abs(aimag(zdum)/2._qp)

           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           ! test the input parameters.
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (.not. la_lsame(normin,'Y') .and. .not. la_lsame(normin,'N')) &
                     then
              info = -4
           else if (n < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WLATRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine machine dependent parameters to control overflow.
           smlnum = la_qlamch('SAFE MINIMUM')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           smlnum = smlnum/la_qlamch('PRECISION')
           bignum = one/smlnum
           scale = one
           if (la_lsame(normin,'N')) then
              ! compute the 1-norm of each column, not including the diagonal.
              if (upper) then
                 ! a is upper triangular.
                 do j = 1,n
                    cnorm(j) = la_qwasum(j - 1,a(1,j),1)
                 end do
              else
                 ! a is lower triangular.
                 do j = 1,n - 1
                    cnorm(j) = la_qwasum(n - j,a(j + 1,j),1)
                 end do
                 cnorm(n) = zero
              end if
           end if
           ! scale the column norms by tscal if the maximum element in cnorm is
           ! greater than bignum/2.
           imax = la_iqamax(n,cnorm,1)
           tmax = cnorm(imax)
           if (tmax <= bignum*half) then
              tscal = one
           else
              tscal = half/(smlnum*tmax)
              call la_qscal(n,tscal,cnorm,1)
           end if
           ! compute a bound on the computed solution vector to see if the
           ! level 2 blas routine la_wtrsv can be used.
           xmax = zero
           do j = 1,n
              xmax = max(xmax,cabs2(x(j)))
           end do
           xbnd = xmax
           if (notran) then
              ! compute the growth in a * x = b.
              if (upper) then
                 jfirst = n
                 jlast = 1
                 jinc = -1
              else
                 jfirst = 1
                 jlast = n
                 jinc = 1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 60
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, g(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    tjjs = a(j,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = g(j-1) / abs(a(j,j))
                       xbnd = min(xbnd,min(one,tjj)*grow)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                    if (tjj + cnorm(j) >= smlnum) then
                       ! g(j) = g(j-1)*( 1 + cnorm(j) / abs(a(j,j)) )
                       grow = grow*(tjj/(tjj + cnorm(j)))
                    else
                       ! g(j) could overflow, set grow to 0.
                       grow = zero
                    end if
                 end do
                 grow = xbnd
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 60
                    ! g(j) = g(j-1)*( 1 + cnorm(j) )
                    grow = grow*(one/(one + cnorm(j)))
                 end do
              end if
              60 continue
           else
              ! compute the growth in a**t * x = b  or  a**h * x = b.
              if (upper) then
                 jfirst = 1
                 jlast = n
                 jinc = 1
              else
                 jfirst = n
                 jlast = 1
                 jinc = -1
              end if
              if (tscal /= one) then
                 grow = zero
                 go to 90
              end if
              if (nounit) then
                 ! a is non-unit triangular.
                 ! compute grow = 1/g(j) and xbnd = 1/m(j).
                 ! initially, m(0) = max{x(i), i=1,...,n}.
                 grow = half/max(xbnd,smlnum)
                 xbnd = grow
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = max( g(j-1), m(j-1)*( 1 + cnorm(j) ) )
                    xj = one + cnorm(j)
                    grow = min(grow,xbnd/xj)
                    tjjs = a(j,j)
                    tjj = cabs1(tjjs)
                    if (tjj >= smlnum) then
                       ! m(j) = m(j-1)*( 1 + cnorm(j) ) / abs(a(j,j))
                       if (xj > tjj) xbnd = xbnd*(tjj/xj)
                    else
                       ! m(j) could overflow, set xbnd to 0.
                       xbnd = zero
                    end if
                 end do
                 grow = min(grow,xbnd)
              else
                 ! a is unit triangular.
                 ! compute grow = 1/g(j), where g(0) = max{x(i), i=1,...,n}.
                 grow = min(one,half/max(xbnd,smlnum))
                 do j = jfirst,jlast,jinc
                    ! exit the loop if the growth factor is too small.
                    if (grow <= smlnum) go to 90
                    ! g(j) = ( 1 + cnorm(j) )*g(j-1)
                    xj = one + cnorm(j)
                    grow = grow/xj
                 end do
              end if
              90 continue
           end if
           if ((grow*tscal) > smlnum) then
              ! use the level 2 blas solve if the reciprocal of the bound on
              ! elements of x is not too small.
              call la_wtrsv(uplo,trans,diag,n,a,lda,x,1)
           else
              ! use a level 1 blas solve, scaling intermediate results.
              if (xmax > bignum*half) then
                 ! scale x so that its components are less than or equal to
                 ! bignum in absolute value.
                 scale = (bignum*half)/xmax
                 call la_wqscal(n,scale,x,1)
                 xmax = bignum
              else
                 xmax = xmax*two
              end if
              if (notran) then
                 ! solve a * x = b
                 loop_120: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) / a(j,j), scaling x if necessary.
                    xj = cabs1(x(j))
                    if (nounit) then
                       tjjs = a(j,j)*tscal
                    else
                       tjjs = tscal
                       if (tscal == one) go to 110
                    end if
                    tjj = cabs1(tjjs)
                    if (tjj > smlnum) then
                          ! abs(a(j,j)) > smlnum:
                       if (tjj < one) then
                          if (xj > tjj*bignum) then
                                ! scale x by 1/b(j).
                             rec = one/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                       end if
                       x(j) = la_wladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else if (tjj > zero) then
                          ! 0 < abs(a(j,j)) <= smlnum:
                       if (xj > tjj*bignum) then
                             ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum
                             ! to avoid overflow when dividing by a(j,j).
                          rec = (tjj*bignum)/xj
                          if (cnorm(j) > one) then
                                ! scale by 1/cnorm(j) to avoid overflow when
                                ! multiplying x(j) times column j.
                             rec = rec/cnorm(j)
                          end if
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                       x(j) = la_wladiv(x(j),tjjs)
                       xj = cabs1(x(j))
                    else
                          ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                          ! scale = 0, and compute a solution to a*x = 0.
                       do i = 1,n
                          x(i) = zero
                       end do
                       x(j) = one
                       xj = one
                       scale = zero
                       xmax = zero
                    end if
                    110 continue
                    ! scale x if necessary to avoid overflow when adding a
                    ! multiple of column j of a.
                    if (xj > one) then
                       rec = one/xj
                       if (cnorm(j) > (bignum - xmax)*rec) then
                          ! scale x by 1/(2*abs(x(j))).
                          rec = rec*half
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                       end if
                    else if (xj*cnorm(j) > (bignum - xmax)) then
                       ! scale x by 1/2.
                       call la_wqscal(n,half,x,1)
                       scale = scale*half
                    end if
                    if (upper) then
                       if (j > 1) then
                          ! compute the update
                             ! x(1:j-1) := x(1:j-1) - x(j) * a(1:j-1,j)
                          call la_waxpy(j - 1,-x(j)*tscal,a(1,j),1,x,1)
                          i = la_iwamax(j - 1,x,1)
                          xmax = cabs1(x(i))
                       end if
                    else
                       if (j < n) then
                          ! compute the update
                             ! x(j+1:n) := x(j+1:n) - x(j) * a(j+1:n,j)
                          call la_waxpy(n - j,-x(j)*tscal,a(j + 1,j),1,x(j + 1),1)

                          i = j + la_iwamax(n - j,x(j + 1),1)
                          xmax = cabs1(x(i))
                       end if
                    end if
                 end do loop_120
              else if (la_lsame(trans,'T')) then
                 ! solve a**t * x = b
                 loop_170: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_wladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=qp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_wdotu to perform the dot product.
                       if (upper) then
                          csumj = la_wdotu(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          csumj = la_wdotu(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (a(i,j)*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             csumj = csumj + (a(i,j)*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=qp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          tjjs = a(j,j)*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 160
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_wqscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**t *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       160 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_wladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_170
              else
                 ! solve a**h * x = b
                 loop_220: do j = jfirst,jlast,jinc
                    ! compute x(j) = b(j) - sum a(k,j)*x(k).
                                          ! k<>j
                    xj = cabs1(x(j))
                    uscal = tscal
                    rec = one/max(xmax,one)
                    if (cnorm(j) > (bignum - xj)*rec) then
                       ! if x(j) could overflow, scale x by 1/(2*xmax).
                       rec = rec*half
                       if (nounit) then
                          tjjs = conjg(a(j,j))*tscal
                       else
                          tjjs = tscal
                       end if
                       tjj = cabs1(tjjs)
                       if (tjj > one) then
                             ! divide by a(j,j) when scaling x if a(j,j) > 1.
                          rec = min(one,rec*tjj)
                          uscal = la_wladiv(uscal,tjjs)
                       end if
                       if (rec < one) then
                          call la_wqscal(n,rec,x,1)
                          scale = scale*rec
                          xmax = xmax*rec
                       end if
                    end if
                    csumj = zero
                    if (uscal == cmplx(one,KIND=qp)) then
                       ! if the scaling needed for a in the dot product is 1,
                       ! call la_wdotc to perform the dot product.
                       if (upper) then
                          csumj = la_wdotc(j - 1,a(1,j),1,x,1)
                       else if (j < n) then
                          csumj = la_wdotc(n - j,a(j + 1,j),1,x(j + 1),1)
                       end if
                    else
                       ! otherwise, use in-line code for the dot product.
                       if (upper) then
                          do i = 1,j - 1
                             csumj = csumj + (conjg(a(i,j))*uscal)*x(i)
                          end do
                       else if (j < n) then
                          do i = j + 1,n
                             csumj = csumj + (conjg(a(i,j))*uscal)*x(i)
                          end do
                       end if
                    end if
                    if (uscal == cmplx(tscal,KIND=qp)) then
                       ! compute x(j) := ( x(j) - csumj ) / a(j,j) if 1/a(j,j)
                       ! was not used to scale the dotproduct.
                       x(j) = x(j) - csumj
                       xj = cabs1(x(j))
                       if (nounit) then
                          tjjs = conjg(a(j,j))*tscal
                       else
                          tjjs = tscal
                          if (tscal == one) go to 210
                       end if
                          ! compute x(j) = x(j) / a(j,j), scaling if necessary.
                       tjj = cabs1(tjjs)
                       if (tjj > smlnum) then
                             ! abs(a(j,j)) > smlnum:
                          if (tjj < one) then
                             if (xj > tjj*bignum) then
                                   ! scale x by 1/abs(x(j)).
                                rec = one/xj
                                call la_wqscal(n,rec,x,1)
                                scale = scale*rec
                                xmax = xmax*rec
                             end if
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else if (tjj > zero) then
                             ! 0 < abs(a(j,j)) <= smlnum:
                          if (xj > tjj*bignum) then
                                ! scale x by (1/abs(x(j)))*abs(a(j,j))*bignum.
                             rec = (tjj*bignum)/xj
                             call la_wqscal(n,rec,x,1)
                             scale = scale*rec
                             xmax = xmax*rec
                          end if
                          x(j) = la_wladiv(x(j),tjjs)
                       else
                             ! a(j,j) = 0:  set x(1:n) = 0, x(j) = 1, and
                             ! scale = 0 and compute a solution to a**h *x = 0.
                          do i = 1,n
                             x(i) = zero
                          end do
                          x(j) = one
                          scale = zero
                          xmax = zero
                       end if
                       210 continue
                    else
                       ! compute x(j) := x(j) / a(j,j) - csumj if the dot
                       ! product has already been divided by 1/a(j,j).
                       x(j) = la_wladiv(x(j),tjjs) - csumj
                    end if
                    xmax = max(xmax,cabs1(x(j)))
                 end do loop_220
              end if
              scale = scale/tscal
           end if
           ! scale the column norms by 1/tscal for return.
           if (tscal /= one) then
              call la_qscal(n,one/tscal,cnorm,1)
           end if
           return
     end subroutine la_wlatrs

     !> CLAUU2: computes the product U * U**H or L**H * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the unblocked form of the algorithm, calling Level 2 BLAS.

     pure subroutine la_clauu2(uplo,n,a,lda,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(sp) :: aii
           ! Intrinsic Functions
           intrinsic :: cmplx,max,real
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CLAUU2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! compute the product u * u**h.
              do i = 1,n
                 aii = real(a(i,i),KIND=sp)
                 if (i < n) then
                    a(i,i) = aii*aii + real(la_cdotc(n - i,a(i,i + 1),lda,a(i,i + 1), &
                              lda),KIND=sp)
                    call la_clacgv(n - i,a(i,i + 1),lda)
                    call la_cgemv('NO TRANSPOSE',i - 1,n - i,cone,a(1,i + 1),lda,a(i,i + 1 &
                              ),lda,cmplx(aii,KIND=sp),a(1,i),1)
                    call la_clacgv(n - i,a(i,i + 1),lda)
                 else
                    call la_csscal(i,aii,a(1,i),1)
                 end if
              end do
           else
              ! compute the product l**h * l.
              do i = 1,n
                 aii = real(a(i,i),KIND=sp)
                 if (i < n) then
                    a(i,i) = aii*aii + real(la_cdotc(n - i,a(i + 1,i),1,a(i + 1,i),1) &
                              ,KIND=sp)
                    call la_clacgv(i - 1,a(i,1),lda)
                    call la_cgemv('CONJUGATE TRANSPOSE',n - i,i - 1,cone,a(i + 1,1),lda,a( &
                              i + 1,i),1,cmplx(aii,KIND=sp),a(i,1),lda)
                    call la_clacgv(i - 1,a(i,1),lda)
                 else
                    call la_csscal(i,aii,a(i,1),lda)
                 end if
              end do
           end if
           return
     end subroutine la_clauu2
     !> ZLAUU2: computes the product U * U**H or L**H * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the unblocked form of the algorithm, calling Level 2 BLAS.

     pure subroutine la_zlauu2(uplo,n,a,lda,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(dp) :: aii
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZLAUU2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! compute the product u * u**h.
              do i = 1,n
                 aii = real(a(i,i),KIND=dp)
                 if (i < n) then
                    a(i,i) = aii*aii + real(la_zdotc(n - i,a(i,i + 1),lda,a(i,i + 1), &
                              lda),KIND=dp)
                    call la_zlacgv(n - i,a(i,i + 1),lda)
                    call la_zgemv('NO TRANSPOSE',i - 1,n - i,cone,a(1,i + 1),lda,a(i,i + 1 &
                              ),lda,cmplx(aii,KIND=dp),a(1,i),1)
                    call la_zlacgv(n - i,a(i,i + 1),lda)
                 else
                    call la_zdscal(i,aii,a(1,i),1)
                 end if
              end do
           else
              ! compute the product l**h * l.
              do i = 1,n
                 aii = real(a(i,i),KIND=dp)
                 if (i < n) then
                    a(i,i) = aii*aii + real(la_zdotc(n - i,a(i + 1,i),1,a(i + 1,i),1) &
                              ,KIND=dp)
                    call la_zlacgv(i - 1,a(i,1),lda)
                    call la_zgemv('CONJUGATE TRANSPOSE',n - i,i - 1,cone,a(i + 1,1),lda,a( &
                              i + 1,i),1,cmplx(aii,KIND=dp),a(i,1),lda)
                    call la_zlacgv(i - 1,a(i,1),lda)
                 else
                    call la_zdscal(i,aii,a(i,1),lda)
                 end if
              end do
           end if
           return
     end subroutine la_zlauu2
     !> WLAUU2: computes the product U * U**H or L**H * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the unblocked form of the algorithm, calling Level 2 BLAS.

     pure subroutine la_wlauu2(uplo,n,a,lda,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(qp) :: aii
           ! Intrinsic Functions
           intrinsic :: real,cmplx,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WLAUU2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! compute the product u * u**h.
              do i = 1,n
                 aii = real(a(i,i),KIND=qp)
                 if (i < n) then
                    a(i,i) = aii*aii + real(la_wdotc(n - i,a(i,i + 1),lda,a(i,i + 1), &
                              lda),KIND=qp)
                    call la_wlacgv(n - i,a(i,i + 1),lda)
                    call la_wgemv('NO TRANSPOSE',i - 1,n - i,cone,a(1,i + 1),lda,a(i,i + 1 &
                              ),lda,cmplx(aii,KIND=qp),a(1,i),1)
                    call la_wlacgv(n - i,a(i,i + 1),lda)
                 else
                    call la_wqscal(i,aii,a(1,i),1)
                 end if
              end do
           else
              ! compute the product l**h * l.
              do i = 1,n
                 aii = real(a(i,i),KIND=qp)
                 if (i < n) then
                    a(i,i) = aii*aii + real(la_wdotc(n - i,a(i + 1,i),1,a(i + 1,i),1) &
                              ,KIND=qp)
                    call la_wlacgv(i - 1,a(i,1),lda)
                    call la_wgemv('CONJUGATE TRANSPOSE',n - i,i - 1,cone,a(i + 1,1),lda,a( &
                              i + 1,i),1,cmplx(aii,KIND=qp),a(i,1),lda)
                    call la_wlacgv(i - 1,a(i,1),lda)
                 else
                    call la_wqscal(i,aii,a(i,1),lda)
                 end if
              end do
           end if
           return
     end subroutine la_wlauu2

     !> CLAUUM: computes the product U * U**H or L**H * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the blocked form of the algorithm, calling Level 3 BLAS.

     pure subroutine la_clauum(uplo,n,a,lda,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ib,nb
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CLAUUM',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'CLAUUM',uplo,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_clauu2(uplo,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute the product u * u**h.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_ctrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1, &
                              ib,cone,a(i,i),lda,a(1,i),lda)
                    call la_clauu2('UPPER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',i - 1,ib,n - i - ib + 1, &
                                  cone,a(1,i + ib),lda,a(i,i + ib),lda,cone,a(1,i),lda)
                       call la_cherk('UPPER','NO TRANSPOSE',ib,n - i - ib + 1,one,a(i,i + ib), &
                                  lda,one,a(i,i),lda)
                    end if
                 end do
              else
                 ! compute the product l**h * l.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_ctrmm('LEFT','LOWER','CONJUGATE TRANSPOSE','NON-UNIT',ib,i - 1, &
                               cone,a(i,i),lda,a(i,1),lda)
                    call la_clauu2('LOWER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_cgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',ib,i - 1,n - i - ib + 1, &
                                  cone,a(i + ib,i),lda,a(i + ib,1),lda,cone,a(i,1),lda)
                       call la_cherk('LOWER','CONJUGATE TRANSPOSE',ib,n - i - ib + 1,one,a(i + &
                                 ib,i),lda,one,a(i,i),lda)
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_clauum
     !> ZLAUUM: computes the product U * U**H or L**H * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the blocked form of the algorithm, calling Level 3 BLAS.

     pure subroutine la_zlauum(uplo,n,a,lda,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ib,nb
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZLAUUM',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'ZLAUUM',uplo,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_zlauu2(uplo,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute the product u * u**h.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_ztrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1, &
                              ib,cone,a(i,i),lda,a(1,i),lda)
                    call la_zlauu2('UPPER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',i - 1,ib,n - i - ib + 1, &
                                  cone,a(1,i + ib),lda,a(i,i + ib),lda,cone,a(1,i),lda)
                       call la_zherk('UPPER','NO TRANSPOSE',ib,n - i - ib + 1,one,a(i,i + ib), &
                                  lda,one,a(i,i),lda)
                    end if
                 end do
              else
                 ! compute the product l**h * l.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_ztrmm('LEFT','LOWER','CONJUGATE TRANSPOSE','NON-UNIT',ib,i - 1, &
                               cone,a(i,i),lda,a(i,1),lda)
                    call la_zlauu2('LOWER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_zgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',ib,i - 1,n - i - ib + 1, &
                                  cone,a(i + ib,i),lda,a(i + ib,1),lda,cone,a(i,1),lda)
                       call la_zherk('LOWER','CONJUGATE TRANSPOSE',ib,n - i - ib + 1,one,a(i + &
                                 ib,i),lda,one,a(i,i),lda)
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_zlauum
     !> WLAUUM: computes the product U * U**H or L**H * L, where the triangular
     !> factor U or L is stored in the upper or lower triangular part of
     !> the array A.
     !> If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
     !> overwriting the factor U in A.
     !> If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
     !> overwriting the factor L in A.
     !> This is the blocked form of the algorithm, calling Level 3 BLAS.

     pure subroutine la_wlauum(uplo,n,a,lda,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ib,nb
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WLAUUM',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'WLAUUM',uplo,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_wlauu2(uplo,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute the product u * u**h.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_wtrmm('RIGHT','UPPER','CONJUGATE TRANSPOSE','NON-UNIT',i - 1, &
                              ib,cone,a(i,i),lda,a(1,i),lda)
                    call la_wlauu2('UPPER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',i - 1,ib,n - i - ib + 1, &
                                  cone,a(1,i + ib),lda,a(i,i + ib),lda,cone,a(1,i),lda)
                       call la_wherk('UPPER','NO TRANSPOSE',ib,n - i - ib + 1,one,a(i,i + ib), &
                                  lda,one,a(i,i),lda)
                    end if
                 end do
              else
                 ! compute the product l**h * l.
                 do i = 1,n,nb
                    ib = min(nb,n - i + 1)
                    call la_wtrmm('LEFT','LOWER','CONJUGATE TRANSPOSE','NON-UNIT',ib,i - 1, &
                               cone,a(i,i),lda,a(i,1),lda)
                    call la_wlauu2('LOWER',ib,a(i,i),lda,info)
                    if (i + ib <= n) then
                       call la_wgemm('CONJUGATE TRANSPOSE','NO TRANSPOSE',ib,i - 1,n - i - ib + 1, &
                                  cone,a(i + ib,i),lda,a(i + ib,1),lda,cone,a(i,1),lda)
                       call la_wherk('LOWER','CONJUGATE TRANSPOSE',ib,n - i - ib + 1,one,a(i + &
                                 ib,i),lda,one,a(i,i),lda)
                    end if
                 end do
              end if
           end if
           return
     end subroutine la_wlauum

     !> CTBRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular band
     !> coefficient matrix.
     !> The solution matrix X must be computed by CTBTRS or some other
     !> means before entering this routine.  CTBRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_ctbrfs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,x,ldx,ferr, &
                berr,work,rwork,info)
        use la_constants_sp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,ldx,n,nrhs
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: ab(ldab,*),b(ldb,*),x(ldx,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,min,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('CTBRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = kd + 2
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_ccopy(n,x(1,j),1,work,1)
              call la_ctbmv(uplo,trans,diag,n,kd,ab,ldab,work,1)
              call la_caxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = max(1,k - kd),k
                             rwork(i) = rwork(i) + cabs1(ab(kd + 1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = max(1,k - kd),k - 1
                             rwork(i) = rwork(i) + cabs1(ab(kd + 1 + i - k,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,min(n,k + kd)
                             rwork(i) = rwork(i) + cabs1(ab(1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             rwork(i) = rwork(i) + cabs1(ab(1 + i - k,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = max(1,k - kd),k
                             s = s + cabs1(ab(kd + 1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = max(1,k - kd),k - 1
                             s = s + cabs1(ab(kd + 1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,min(n,k + kd)
                             s = s + cabs1(ab(1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             s = s + cabs1(ab(1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_clacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_clacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_ctbsv(uplo,transt,diag,n,kd,ab,ldab,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ctbsv(uplo,transn,diag,n,kd,ab,ldab,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_ctbrfs
     !> ZTBRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular band
     !> coefficient matrix.
     !> The solution matrix X must be computed by ZTBTRS or some other
     !> means before entering this routine.  ZTBRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_ztbrfs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,x,ldx,ferr, &
                berr,work,rwork,info)
        use la_constants_dp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,ldx,n,nrhs
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: ab(ldab,*),b(ldb,*),x(ldx,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('ZTBRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = kd + 2
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_zcopy(n,x(1,j),1,work,1)
              call la_ztbmv(uplo,trans,diag,n,kd,ab,ldab,work,1)
              call la_zaxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = max(1,k - kd),k
                             rwork(i) = rwork(i) + cabs1(ab(kd + 1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = max(1,k - kd),k - 1
                             rwork(i) = rwork(i) + cabs1(ab(kd + 1 + i - k,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,min(n,k + kd)
                             rwork(i) = rwork(i) + cabs1(ab(1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             rwork(i) = rwork(i) + cabs1(ab(1 + i - k,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = max(1,k - kd),k
                             s = s + cabs1(ab(kd + 1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = max(1,k - kd),k - 1
                             s = s + cabs1(ab(kd + 1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,min(n,k + kd)
                             s = s + cabs1(ab(1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             s = s + cabs1(ab(1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_zlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_zlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_ztbsv(uplo,transt,diag,n,kd,ab,ldab,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ztbsv(uplo,transn,diag,n,kd,ab,ldab,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_ztbrfs
     !> WTBRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular band
     !> coefficient matrix.
     !> The solution matrix X must be computed by WTBTRS or some other
     !> means before entering this routine.  WTBRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_wtbrfs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,x,ldx,ferr, &
                berr,work,rwork,info)
        use la_constants_qp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,ldx,n,nrhs
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: ab(ldab,*),b(ldb,*),x(ldx,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('WTBRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = kd + 2
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_wcopy(n,x(1,j),1,work,1)
              call la_wtbmv(uplo,trans,diag,n,kd,ab,ldab,work,1)
              call la_waxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = max(1,k - kd),k
                             rwork(i) = rwork(i) + cabs1(ab(kd + 1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = max(1,k - kd),k - 1
                             rwork(i) = rwork(i) + cabs1(ab(kd + 1 + i - k,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,min(n,k + kd)
                             rwork(i) = rwork(i) + cabs1(ab(1 + i - k,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             rwork(i) = rwork(i) + cabs1(ab(1 + i - k,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = max(1,k - kd),k
                             s = s + cabs1(ab(kd + 1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = max(1,k - kd),k - 1
                             s = s + cabs1(ab(kd + 1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,min(n,k + kd)
                             s = s + cabs1(ab(1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,min(n,k + kd)
                             s = s + cabs1(ab(1 + i - k,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_wlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_wlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_wtbsv(uplo,transt,diag,n,kd,ab,ldab,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_wtbsv(uplo,transn,diag,n,kd,ab,ldab,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_wtbrfs

     !> CTBTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular band matrix of order N, and B is an
     !> N-by-NRHS matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_ctbtrs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(in) :: ab(ldab,*)
           complex(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CTBTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 do info = 1,n
                    if (ab(kd + 1,info) == czero) return
                 end do
              else
                 do info = 1,n
                    if (ab(1,info) == czero) return
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b,  a**t * x = b,  or  a**h * x = b.
           do j = 1,nrhs
              call la_ctbsv(uplo,trans,diag,n,kd,ab,ldab,b(1,j),1)
           end do
           return
     end subroutine la_ctbtrs
     !> ZTBTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular band matrix of order N, and B is an
     !> N-by-NRHS matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_ztbtrs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(in) :: ab(ldab,*)
           complex(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZTBTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 do info = 1,n
                    if (ab(kd + 1,info) == czero) return
                 end do
              else
                 do info = 1,n
                    if (ab(1,info) == czero) return
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b,  a**t * x = b,  or  a**h * x = b.
           do j = 1,nrhs
              call la_ztbsv(uplo,trans,diag,n,kd,ab,ldab,b(1,j),1)
           end do
           return
     end subroutine la_ztbtrs
     !> WTBTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular band matrix of order N, and B is an
     !> N-by-NRHS matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_wtbtrs(uplo,trans,diag,n,kd,nrhs,ab,ldab,b,ldb,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(in) :: ab(ldab,*)
           complex(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (nrhs < 0) then
              info = -6
           else if (ldab < kd + 1) then
              info = -8
           else if (ldb < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WTBTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 do info = 1,n
                    if (ab(kd + 1,info) == czero) return
                 end do
              else
                 do info = 1,n
                    if (ab(1,info) == czero) return
                 end do
              end if
           end if
           info = 0
           ! solve a * x = b,  a**t * x = b,  or  a**h * x = b.
           do j = 1,nrhs
              call la_wtbsv(uplo,trans,diag,n,kd,ab,ldab,b(1,j),1)
           end do
           return
     end subroutine la_wtbtrs

     !> CTPRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular packed
     !> coefficient matrix.
     !> The solution matrix X must be computed by CTPTRS or some other
     !> means before entering this routine.  CTPRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_ctprfs(uplo,trans,diag,n,nrhs,ap,b,ldb,x,ldx,ferr,berr, &
               work,rwork,info)
        use la_constants_sp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: ap(*),b(ldb,*),x(ldx,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,kc,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CTPRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_ccopy(n,x(1,j),1,work,1)
              call la_ctpmv(uplo,trans,diag,n,ap,work,1)
              call la_caxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - 1))*xk
                          end do
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k - 1
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - 1))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,n
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - k))*xk
                          end do
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,n
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + cabs1(ap(kc + i - 1))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = 1,k - 1
                             s = s + cabs1(ap(kc + i - 1))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + cabs1(ap(kc + i - k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,n
                             s = s + cabs1(ap(kc + i - k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_clacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_clacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_ctpsv(uplo,transt,diag,n,ap,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ctpsv(uplo,transn,diag,n,ap,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_ctprfs
     !> ZTPRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular packed
     !> coefficient matrix.
     !> The solution matrix X must be computed by ZTPTRS or some other
     !> means before entering this routine.  ZTPRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_ztprfs(uplo,trans,diag,n,nrhs,ap,b,ldb,x,ldx,ferr,berr, &
               work,rwork,info)
        use la_constants_dp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: ap(*),b(ldb,*),x(ldx,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,kc,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZTPRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_zcopy(n,x(1,j),1,work,1)
              call la_ztpmv(uplo,trans,diag,n,ap,work,1)
              call la_zaxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - 1))*xk
                          end do
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k - 1
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - 1))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,n
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - k))*xk
                          end do
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,n
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + cabs1(ap(kc + i - 1))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = 1,k - 1
                             s = s + cabs1(ap(kc + i - 1))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + cabs1(ap(kc + i - k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,n
                             s = s + cabs1(ap(kc + i - k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_zlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_zlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_ztpsv(uplo,transt,diag,n,ap,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ztpsv(uplo,transn,diag,n,ap,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_ztprfs
     !> WTPRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular packed
     !> coefficient matrix.
     !> The solution matrix X must be computed by WTPTRS or some other
     !> means before entering this routine.  WTPRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_wtprfs(uplo,trans,diag,n,nrhs,ap,b,ldb,x,ldx,ferr,berr, &
               work,rwork,info)
        use la_constants_qp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: ap(*),b(ldb,*),x(ldx,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,kc,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WTPRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_wcopy(n,x(1,j),1,work,1)
              call la_wtpmv(uplo,trans,diag,n,ap,work,1)
              call la_waxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - 1))*xk
                          end do
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k - 1
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - 1))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,n
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - k))*xk
                          end do
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,n
                             rwork(i) = rwork(i) + cabs1(ap(kc + i - k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + cabs1(ap(kc + i - 1))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + k
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = 1,k - 1
                             s = s + cabs1(ap(kc + i - 1))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + k
                       end do
                    end if
                 else
                    kc = 1
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + cabs1(ap(kc + i - k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + n - k + 1
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,n
                             s = s + cabs1(ap(kc + i - k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                          kc = kc + n - k + 1
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_wlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_wlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_wtpsv(uplo,transt,diag,n,ap,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_wtpsv(uplo,transn,diag,n,ap,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_wtprfs

     !> CTPTRI: computes the inverse of a complex upper or lower triangular
     !> matrix A stored in packed format.

     pure subroutine la_ctptri(uplo,diag,n,ap,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(sp),intent(inout) :: ap(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc,jclast,jj
           complex(sp) :: ajj
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('CTPTRI',-info)
              return
           end if
           ! check for singularity if non-unit.
           if (nounit) then
              if (upper) then
                 jj = 0
                 do info = 1,n
                    jj = jj + info
                    if (ap(jj) == czero) return
                 end do
              else
                 jj = 1
                 do info = 1,n
                    if (ap(jj) == czero) return
                    jj = jj + n - info + 1
                 end do
              end if
              info = 0
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              jc = 1
              do j = 1,n
                 if (nounit) then
                    ap(jc + j - 1) = cone/ap(jc + j - 1)
                    ajj = -ap(jc + j - 1)
                 else
                    ajj = -cone
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_ctpmv('UPPER','NO TRANSPOSE',diag,j - 1,ap,ap(jc),1)
                 call la_cscal(j - 1,ajj,ap(jc),1)
                 jc = jc + j
              end do
           else
              ! compute inverse of lower triangular matrix.
              jc = n*(n + 1)/2
              do j = n,1,-1
                 if (nounit) then
                    ap(jc) = cone/ap(jc)
                    ajj = -ap(jc)
                 else
                    ajj = -cone
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_ctpmv('LOWER','NO TRANSPOSE',diag,n - j,ap(jclast),ap(jc + 1) &
                              ,1)
                    call la_cscal(n - j,ajj,ap(jc + 1),1)
                 end if
                 jclast = jc
                 jc = jc - n + j - 2
              end do
           end if
           return
     end subroutine la_ctptri
     !> ZTPTRI: computes the inverse of a complex upper or lower triangular
     !> matrix A stored in packed format.

     pure subroutine la_ztptri(uplo,diag,n,ap,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(dp),intent(inout) :: ap(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc,jclast,jj
           complex(dp) :: ajj
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('ZTPTRI',-info)
              return
           end if
           ! check for singularity if non-unit.
           if (nounit) then
              if (upper) then
                 jj = 0
                 do info = 1,n
                    jj = jj + info
                    if (ap(jj) == czero) return
                 end do
              else
                 jj = 1
                 do info = 1,n
                    if (ap(jj) == czero) return
                    jj = jj + n - info + 1
                 end do
              end if
              info = 0
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              jc = 1
              do j = 1,n
                 if (nounit) then
                    ap(jc + j - 1) = cone/ap(jc + j - 1)
                    ajj = -ap(jc + j - 1)
                 else
                    ajj = -cone
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_ztpmv('UPPER','NO TRANSPOSE',diag,j - 1,ap,ap(jc),1)
                 call la_zscal(j - 1,ajj,ap(jc),1)
                 jc = jc + j
              end do
           else
              ! compute inverse of lower triangular matrix.
              jc = n*(n + 1)/2
              do j = n,1,-1
                 if (nounit) then
                    ap(jc) = cone/ap(jc)
                    ajj = -ap(jc)
                 else
                    ajj = -cone
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_ztpmv('LOWER','NO TRANSPOSE',diag,n - j,ap(jclast),ap(jc + 1) &
                              ,1)
                    call la_zscal(n - j,ajj,ap(jc + 1),1)
                 end if
                 jclast = jc
                 jc = jc - n + j - 2
              end do
           end if
           return
     end subroutine la_ztptri
     !> WTPTRI: computes the inverse of a complex upper or lower triangular
     !> matrix A stored in packed format.

     pure subroutine la_wtptri(uplo,diag,n,ap,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(qp),intent(inout) :: ap(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc,jclast,jj
           complex(qp) :: ajj
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           end if
           if (info /= 0) then
              call la_xerbla('WTPTRI',-info)
              return
           end if
           ! check for singularity if non-unit.
           if (nounit) then
              if (upper) then
                 jj = 0
                 do info = 1,n
                    jj = jj + info
                    if (ap(jj) == czero) return
                 end do
              else
                 jj = 1
                 do info = 1,n
                    if (ap(jj) == czero) return
                    jj = jj + n - info + 1
                 end do
              end if
              info = 0
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              jc = 1
              do j = 1,n
                 if (nounit) then
                    ap(jc + j - 1) = cone/ap(jc + j - 1)
                    ajj = -ap(jc + j - 1)
                 else
                    ajj = -cone
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_wtpmv('UPPER','NO TRANSPOSE',diag,j - 1,ap,ap(jc),1)
                 call la_wscal(j - 1,ajj,ap(jc),1)
                 jc = jc + j
              end do
           else
              ! compute inverse of lower triangular matrix.
              jc = n*(n + 1)/2
              do j = n,1,-1
                 if (nounit) then
                    ap(jc) = cone/ap(jc)
                    ajj = -ap(jc)
                 else
                    ajj = -cone
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_wtpmv('LOWER','NO TRANSPOSE',diag,n - j,ap(jclast),ap(jc + 1) &
                              ,1)
                    call la_wscal(n - j,ajj,ap(jc + 1),1)
                 end if
                 jclast = jc
                 jc = jc - n + j - 2
              end do
           end if
           return
     end subroutine la_wtptri

     !> CTPTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular matrix of order N stored in packed format,
     !> and B is an N-by-NRHS matrix.  A check is made to verify that A is
     !> nonsingular.

     pure subroutine la_ctptrs(uplo,trans,diag,n,nrhs,ap,b,ldb,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(in) :: ap(*)
           complex(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CTPTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 jc = 1
                 do info = 1,n
                    if (ap(jc + info - 1) == czero) return
                    jc = jc + info
                 end do
              else
                 jc = 1
                 do info = 1,n
                    if (ap(jc) == czero) return
                    jc = jc + n - info + 1
                 end do
              end if
           end if
           info = 0
           ! solve  a * x = b,  a**t * x = b,  or  a**h * x = b.
           do j = 1,nrhs
              call la_ctpsv(uplo,trans,diag,n,ap,b(1,j),1)
           end do
           return
     end subroutine la_ctptrs
     !> ZTPTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular matrix of order N stored in packed format,
     !> and B is an N-by-NRHS matrix.  A check is made to verify that A is
     !> nonsingular.

     pure subroutine la_ztptrs(uplo,trans,diag,n,nrhs,ap,b,ldb,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(in) :: ap(*)
           complex(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZTPTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 jc = 1
                 do info = 1,n
                    if (ap(jc + info - 1) == czero) return
                    jc = jc + info
                 end do
              else
                 jc = 1
                 do info = 1,n
                    if (ap(jc) == czero) return
                    jc = jc + n - info + 1
                 end do
              end if
           end if
           info = 0
           ! solve  a * x = b,  a**t * x = b,  or  a**h * x = b.
           do j = 1,nrhs
              call la_ztpsv(uplo,trans,diag,n,ap,b(1,j),1)
           end do
           return
     end subroutine la_ztptrs
     !> WTPTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular matrix of order N stored in packed format,
     !> and B is an N-by-NRHS matrix.  A check is made to verify that A is
     !> nonsingular.

     pure subroutine la_wtptrs(uplo,trans,diag,n,nrhs,ap,b,ldb,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(in) :: ap(*)
           complex(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jc
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WTPTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              if (upper) then
                 jc = 1
                 do info = 1,n
                    if (ap(jc + info - 1) == czero) return
                    jc = jc + info
                 end do
              else
                 jc = 1
                 do info = 1,n
                    if (ap(jc) == czero) return
                    jc = jc + n - info + 1
                 end do
              end if
           end if
           info = 0
           ! solve  a * x = b,  a**t * x = b,  or  a**h * x = b.
           do j = 1,nrhs
              call la_wtpsv(uplo,trans,diag,n,ap,b(1,j),1)
           end do
           return
     end subroutine la_wtptrs

     !> CTRRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular
     !> coefficient matrix.
     !> The solution matrix X must be computed by CTRTRS or some other
     !> means before entering this routine.  CTRRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_ctrrfs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,x,ldx,ferr,berr, &
                work,rwork,info)
        use la_constants_sp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: a(lda,*),b(ldb,*),x(ldx,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('CTRRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_ccopy(n,x(1,j),1,work,1)
              call la_ctrmv(uplo,trans,diag,n,a,lda,work,1)
              call la_caxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k - 1
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,n
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,n
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = 1,k - 1
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,n
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_clacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_clacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_ctrsv(uplo,transt,diag,n,a,lda,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ctrsv(uplo,transn,diag,n,a,lda,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_ctrrfs
     !> ZTRRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular
     !> coefficient matrix.
     !> The solution matrix X must be computed by ZTRTRS or some other
     !> means before entering this routine.  ZTRRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_ztrrfs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,x,ldx,ferr,berr, &
                work,rwork,info)
        use la_constants_dp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: a(lda,*),b(ldb,*),x(ldx,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('ZTRRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_zcopy(n,x(1,j),1,work,1)
              call la_ztrmv(uplo,trans,diag,n,a,lda,work,1)
              call la_zaxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k - 1
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,n
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,n
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = 1,k - 1
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,n
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_zlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_zlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_ztrsv(uplo,transt,diag,n,a,lda,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ztrsv(uplo,transn,diag,n,a,lda,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_ztrrfs
     !> WTRRFS: provides error bounds and backward error estimates for the
     !> solution to a system of linear equations with a triangular
     !> coefficient matrix.
     !> The solution matrix X must be computed by WTRTRS or some other
     !> means before entering this routine.  WTRRFS does not do iterative
     !> refinement because doing so cannot improve the backward error.

     pure subroutine la_wtrrfs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,x,ldx,ferr,berr, &
                work,rwork,info)
        use la_constants_qp,only:zero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,ldx,n,nrhs
           ! Array Arguments
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: a(lda,*),b(ldb,*),x(ldx,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: notran,nounit,upper
           character :: transn,transt
           integer(ilp) :: i,j,k,kase,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           notran = la_lsame(trans,'N')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. notran .and. .not. la_lsame(trans,'T') .and. .not. la_lsame( &
                     trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           else if (ldx < max(1,n)) then
              info = -11
           end if
           if (info /= 0) then
              call la_xerbla('WTRRFS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) then
              do j = 1,nrhs
                 ferr(j) = zero
                 berr(j) = zero
              end do
              return
           end if
           if (notran) then
              transn = 'N'
              transt = 'C'
           else
              transn = 'C'
              transt = 'N'
           end if
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_250: do j = 1,nrhs
              ! compute residual r = b - op(a) * x,
              ! where op(a) = a, a**t, or a**h, depending on trans.
              call la_wcopy(n,x(1,j),1,work,1)
              call la_wtrmv(uplo,trans,diag,n,a,lda,work,1)
              call la_waxpy(n,-cone,b(1,j),1,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(op(a))*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              if (notran) then
                 ! compute abs(a)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = 1,k - 1
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k,n
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                       end do
                    else
                       do k = 1,n
                          xk = cabs1(x(k,j))
                          do i = k + 1,n
                             rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                          end do
                          rwork(k) = rwork(k) + xk
                       end do
                    end if
                 end if
              else
                 ! compute abs(a**h)*abs(x) + abs(b).
                 if (upper) then
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = 1,k
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = 1,k - 1
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 else
                    if (nounit) then
                       do k = 1,n
                          s = zero
                          do i = k,n
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    else
                       do k = 1,n
                          s = cabs1(x(k,j))
                          do i = k + 1,n
                             s = s + cabs1(a(i,k))*cabs1(x(i,j))
                          end do
                          rwork(k) = rwork(k) + s
                       end do
                    end if
                 end if
              end if
              s = zero
              do i = 1,n
                 if (rwork(i) > safe2) then
                    s = max(s,cabs1(work(i))/rwork(i))
                 else
                    s = max(s, (cabs1(work(i)) + safe1)/(rwork(i) + safe1))
                 end if
              end do
              berr(j) = s
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(op(a)))*
                 ! ( abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(op(a)) is the inverse of op(a)
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(op(a))*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(op(a))*abs(x) + abs(b) is less than safe2.
              ! use la_wlacn2 to estimate the infinity-norm of the matrix
                 ! inv(op(a)) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(op(a))*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              210 continue
              call la_wlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(op(a)**h).
                    call la_wtrsv(uplo,transt,diag,n,a,lda,work,1)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else
                    ! multiply by inv(op(a))*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_wtrsv(uplo,transn,diag,n,a,lda,work,1)
                 end if
                 go to 210
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_250
           return
     end subroutine la_wtrrfs

     !> CTRTI2: computes the inverse of a complex upper or lower triangular
     !> matrix.
     !> This is the Level 2 BLAS version of the algorithm.

     pure subroutine la_ctrti2(uplo,diag,n,a,lda,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           complex(sp) :: ajj
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CTRTI2',-info)
              return
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              do j = 1,n
                 if (nounit) then
                    a(j,j) = cone/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -cone
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_ctrmv('UPPER','NO TRANSPOSE',diag,j - 1,a,lda,a(1,j),1)

                 call la_cscal(j - 1,ajj,a(1,j),1)
              end do
           else
              ! compute inverse of lower triangular matrix.
              do j = n,1,-1
                 if (nounit) then
                    a(j,j) = cone/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -cone
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_ctrmv('LOWER','NO TRANSPOSE',diag,n - j,a(j + 1,j + 1),lda,a( &
                              j + 1,j),1)
                    call la_cscal(n - j,ajj,a(j + 1,j),1)
                 end if
              end do
           end if
           return
     end subroutine la_ctrti2
     !> ZTRTI2: computes the inverse of a complex upper or lower triangular
     !> matrix.
     !> This is the Level 2 BLAS version of the algorithm.

     pure subroutine la_ztrti2(uplo,diag,n,a,lda,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           complex(dp) :: ajj
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZTRTI2',-info)
              return
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              do j = 1,n
                 if (nounit) then
                    a(j,j) = cone/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -cone
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_ztrmv('UPPER','NO TRANSPOSE',diag,j - 1,a,lda,a(1,j),1)

                 call la_zscal(j - 1,ajj,a(1,j),1)
              end do
           else
              ! compute inverse of lower triangular matrix.
              do j = n,1,-1
                 if (nounit) then
                    a(j,j) = cone/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -cone
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_ztrmv('LOWER','NO TRANSPOSE',diag,n - j,a(j + 1,j + 1),lda,a( &
                              j + 1,j),1)
                    call la_zscal(n - j,ajj,a(j + 1,j),1)
                 end if
              end do
           end if
           return
     end subroutine la_ztrti2
     !> WTRTI2: computes the inverse of a complex upper or lower triangular
     !> matrix.
     !> This is the Level 2 BLAS version of the algorithm.

     pure subroutine la_wtrti2(uplo,diag,n,a,lda,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j
           complex(qp) :: ajj
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WTRTI2',-info)
              return
           end if
           if (upper) then
              ! compute inverse of upper triangular matrix.
              do j = 1,n
                 if (nounit) then
                    a(j,j) = cone/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -cone
                 end if
                 ! compute elements 1:j-1 of j-th column.
                 call la_wtrmv('UPPER','NO TRANSPOSE',diag,j - 1,a,lda,a(1,j),1)

                 call la_wscal(j - 1,ajj,a(1,j),1)
              end do
           else
              ! compute inverse of lower triangular matrix.
              do j = n,1,-1
                 if (nounit) then
                    a(j,j) = cone/a(j,j)
                    ajj = -a(j,j)
                 else
                    ajj = -cone
                 end if
                 if (j < n) then
                    ! compute elements j+1:n of j-th column.
                    call la_wtrmv('LOWER','NO TRANSPOSE',diag,n - j,a(j + 1,j + 1),lda,a( &
                              j + 1,j),1)
                    call la_wscal(n - j,ajj,a(j + 1,j),1)
                 end if
              end do
           end if
           return
     end subroutine la_wtrti2

     !> CTRTRI: computes the inverse of a complex upper or lower triangular
     !> matrix A.
     !> This is the Level 3 BLAS version of the algorithm.

     pure subroutine la_ctrtri(uplo,diag,n,a,lda,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jb,nb,nn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CTRTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity if non-unit.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == czero) return
              end do
              info = 0
           end if
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'CTRTRI',uplo//diag,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_ctrti2(uplo,diag,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute inverse of upper triangular matrix
                 do j = 1,n,nb
                    jb = min(nb,n - j + 1)
                    ! compute rows 1:j-1 of current block column
                    call la_ctrmm('LEFT','UPPER','NO TRANSPOSE',diag,j - 1,jb,cone,a, &
                              lda,a(1,j),lda)
                    call la_ctrsm('RIGHT','UPPER','NO TRANSPOSE',diag,j - 1,jb,-cone,a( &
                              j,j),lda,a(1,j),lda)
                    ! compute inverse of current diagonal block
                    call la_ctrti2('UPPER',diag,jb,a(j,j),lda,info)
                 end do
              else
                 ! compute inverse of lower triangular matrix
                 nn = ((n - 1)/nb)*nb + 1
                 do j = nn,1,-nb
                    jb = min(nb,n - j + 1)
                    if (j + jb <= n) then
                       ! compute rows j+jb:n of current block column
                       call la_ctrmm('LEFT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb, &
                                 cone,a(j + jb,j + jb),lda,a(j + jb,j),lda)
                       call la_ctrsm('RIGHT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,- &
                                 cone,a(j,j),lda,a(j + jb,j),lda)
                    end if
                    ! compute inverse of current diagonal block
                    call la_ctrti2('LOWER',diag,jb,a(j,j),lda,info)
                 end do
              end if
           end if
           return
     end subroutine la_ctrtri
     !> ZTRTRI: computes the inverse of a complex upper or lower triangular
     !> matrix A.
     !> This is the Level 3 BLAS version of the algorithm.

     pure subroutine la_ztrtri(uplo,diag,n,a,lda,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jb,nb,nn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZTRTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity if non-unit.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == czero) return
              end do
              info = 0
           end if
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'ZTRTRI',uplo//diag,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_ztrti2(uplo,diag,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute inverse of upper triangular matrix
                 do j = 1,n,nb
                    jb = min(nb,n - j + 1)
                    ! compute rows 1:j-1 of current block column
                    call la_ztrmm('LEFT','UPPER','NO TRANSPOSE',diag,j - 1,jb,cone,a, &
                              lda,a(1,j),lda)
                    call la_ztrsm('RIGHT','UPPER','NO TRANSPOSE',diag,j - 1,jb,-cone,a( &
                              j,j),lda,a(1,j),lda)
                    ! compute inverse of current diagonal block
                    call la_ztrti2('UPPER',diag,jb,a(j,j),lda,info)
                 end do
              else
                 ! compute inverse of lower triangular matrix
                 nn = ((n - 1)/nb)*nb + 1
                 do j = nn,1,-nb
                    jb = min(nb,n - j + 1)
                    if (j + jb <= n) then
                       ! compute rows j+jb:n of current block column
                       call la_ztrmm('LEFT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb, &
                                 cone,a(j + jb,j + jb),lda,a(j + jb,j),lda)
                       call la_ztrsm('RIGHT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,- &
                                 cone,a(j,j),lda,a(j + jb,j),lda)
                    end if
                    ! compute inverse of current diagonal block
                    call la_ztrti2('LOWER',diag,jb,a(j,j),lda,info)
                 end do
              end if
           end if
           return
     end subroutine la_ztrtri
     !> WTRTRI: computes the inverse of a complex upper or lower triangular
     !> matrix A.
     !> This is the Level 3 BLAS version of the algorithm.

     pure subroutine la_wtrtri(uplo,diag,n,a,lda,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,upper
           integer(ilp) :: j,jb,nb,nn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           nounit = la_lsame(diag,'N')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WTRTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity if non-unit.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == czero) return
              end do
              info = 0
           end if
           ! determine the block size for this environment.
           nb = la_ilaenv(1,'WTRTRI',uplo//diag,n,-1,-1,-1)
           if (nb <= 1 .or. nb >= n) then
              ! use unblocked code
              call la_wtrti2(uplo,diag,n,a,lda,info)
           else
              ! use blocked code
              if (upper) then
                 ! compute inverse of upper triangular matrix
                 do j = 1,n,nb
                    jb = min(nb,n - j + 1)
                    ! compute rows 1:j-1 of current block column
                    call la_wtrmm('LEFT','UPPER','NO TRANSPOSE',diag,j - 1,jb,cone,a, &
                              lda,a(1,j),lda)
                    call la_wtrsm('RIGHT','UPPER','NO TRANSPOSE',diag,j - 1,jb,-cone,a( &
                              j,j),lda,a(1,j),lda)
                    ! compute inverse of current diagonal block
                    call la_wtrti2('UPPER',diag,jb,a(j,j),lda,info)
                 end do
              else
                 ! compute inverse of lower triangular matrix
                 nn = ((n - 1)/nb)*nb + 1
                 do j = nn,1,-nb
                    jb = min(nb,n - j + 1)
                    if (j + jb <= n) then
                       ! compute rows j+jb:n of current block column
                       call la_wtrmm('LEFT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb, &
                                 cone,a(j + jb,j + jb),lda,a(j + jb,j),lda)
                       call la_wtrsm('RIGHT','LOWER','NO TRANSPOSE',diag,n - j - jb + 1,jb,- &
                                 cone,a(j,j),lda,a(j + jb,j),lda)
                    end if
                    ! compute inverse of current diagonal block
                    call la_wtrti2('LOWER',diag,jb,a(j,j),lda,info)
                 end do
              end if
           end if
           return
     end subroutine la_wtrtri

     !> CTRTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular matrix of order N, and B is an N-by-NRHS
     !> matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_ctrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('CTRTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == czero) return
              end do
           end if
           info = 0
           ! solve a * x = b,  a**t * x = b,  or  a**h * x = b.
           call la_ctrsm('LEFT',uplo,trans,diag,n,nrhs,cone,a,lda,b,ldb)
           return
     end subroutine la_ctrtrs
     !> ZTRTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular matrix of order N, and B is an N-by-NRHS
     !> matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_ztrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZTRTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == czero) return
              end do
           end if
           info = 0
           ! solve a * x = b,  a**t * x = b,  or  a**h * x = b.
           call la_ztrsm('LEFT',uplo,trans,diag,n,nrhs,cone,a,lda,b,ldb)
           return
     end subroutine la_ztrtrs
     !> WTRTRS: solves a triangular system of the form
     !> A * X = B,  A**T * X = B,  or  A**H * X = B,
     !> where A is a triangular matrix of order N, and B is an N-by-NRHS
     !> matrix.  A check is made to verify that A is nonsingular.

     pure subroutine la_wtrtrs(uplo,trans,diag,n,nrhs,a,lda,b,ldb,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,trans,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           nounit = la_lsame(diag,'N')
           if (.not. la_lsame(uplo,'U') .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. la_lsame(trans,'N') .and. .not. la_lsame(trans,'T') .and. &
                     .not. la_lsame(trans,'C')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (nrhs < 0) then
              info = -5
           else if (lda < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WTRTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check for singularity.
           if (nounit) then
              do info = 1,n
                 if (a(info,info) == czero) return
              end do
           end if
           info = 0
           ! solve a * x = b,  a**t * x = b,  or  a**h * x = b.
           call la_wtrsm('LEFT',uplo,trans,diag,n,nrhs,cone,a,lda,b,ldb)
           return
     end subroutine la_wtrtrs

     !> CTBCON: estimates the reciprocal of the condition number of a
     !> triangular band matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_ctbcon(norm,uplo,diag,n,kd,ab,ldab,rcond,work,rwork,info)
        use la_constants_sp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: ab(ldab,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(sp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('CTBCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_slamch('SAFE MINIMUM')*real(max(n,1),KIND=sp)
           ! compute the 1-norm of the triangular matrix a or a**h.
           anorm = la_clantb(norm,uplo,diag,n,kd,ab,ldab,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the 1-norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_clatbs(uplo,'NO TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_clatbs(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,kd,ab,ldab, &
                               work,scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_icamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_csrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_ctbcon
     !> ZTBCON: estimates the reciprocal of the condition number of a
     !> triangular band matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_ztbcon(norm,uplo,diag,n,kd,ab,ldab,rcond,work,rwork,info)
        use la_constants_dp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: ab(ldab,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(dp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('ZTBCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_dlamch('SAFE MINIMUM')*real(max(n,1),KIND=dp)
           ! compute the 1-norm of the triangular matrix a or a**h.
           anorm = la_zlantb(norm,uplo,diag,n,kd,ab,ldab,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the 1-norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_zlatbs(uplo,'NO TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_zlatbs(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,kd,ab,ldab, &
                               work,scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_izamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_zdrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_ztbcon
     !> WTBCON: estimates the reciprocal of the condition number of a
     !> triangular band matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_wtbcon(norm,uplo,diag,n,kd,ab,ldab,rcond,work,rwork,info)
        use la_constants_qp,only:zero,one

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kd,ldab,n
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: ab(ldab,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(qp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (kd < 0) then
              info = -5
           else if (ldab < kd + 1) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('WTBCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_qlamch('SAFE MINIMUM')*real(max(n,1),KIND=qp)
           ! compute the 1-norm of the triangular matrix a or a**h.
           anorm = la_wlantb(norm,uplo,diag,n,kd,ab,ldab,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the 1-norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_wlatbs(uplo,'NO TRANSPOSE',diag,normin,n,kd,ab,ldab,work, &
                              scale,rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_wlatbs(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,kd,ab,ldab, &
                               work,scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_iwamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_wqrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_wtbcon

     !> CTFTRI: computes the inverse of a triangular matrix A stored in RFP
     !> format.
     !> This is a Level 3 BLAS version of the algorithm.

     pure subroutine la_ctftri(transr,uplo,diag,n,a,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo,diag
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(sp),intent(inout) :: a(0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CTFTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
           else
              nisodd = .true.
           end if
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1)
                    call la_ctrtri('L',diag,n1,a(0),n,info)
                    if (info > 0) return
                    call la_ctrmm('R','L','N',diag,n2,n1,-cone,a(0),n,a(n1),n)

                    call la_ctrtri('U',diag,n2,a(n),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ctrmm('L','U','C',diag,n2,n1,cone,a(n),n,a(n1),n)

                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    call la_ctrtri('L',diag,n1,a(n2),n,info)
                    if (info > 0) return
                    call la_ctrmm('L','L','C',diag,n1,n2,-cone,a(n2),n,a(0),n)

                    call la_ctrtri('U',diag,n2,a(n1),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ctrmm('R','U','N',diag,n1,n2,cone,a(n1),n,a(0),n)

                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0), t2 -> a(1), s -> a(0+n1*n1)
                    call la_ctrtri('U',diag,n1,a(0),n1,info)
                    if (info > 0) return
                    call la_ctrmm('L','U','N',diag,n1,n2,-cone,a(0),n1,a(n1*n1), &
                              n1)
                    call la_ctrtri('L',diag,n2,a(1),n1,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ctrmm('R','L','C',diag,n1,n2,cone,a(1),n1,a(n1*n1), &
                              n1)
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0+n2*n2), t2 -> a(0+n1*n2), s -> a(0)
                    call la_ctrtri('U',diag,n1,a(n2*n2),n2,info)
                    if (info > 0) return
                    call la_ctrmm('R','U','C',diag,n2,n1,-cone,a(n2*n2),n2,a(0), &
                              n2)
                    call la_ctrtri('L',diag,n2,a(n1*n2),n2,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ctrmm('L','L','N',diag,n2,n1,cone,a(n1*n2),n2,a(0), &
                              n2)
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    call la_ctrtri('L',diag,k,a(1),n + 1,info)
                    if (info > 0) return
                    call la_ctrmm('R','L','N',diag,k,k,-cone,a(1),n + 1,a(k + 1),n + &
                              1)
                    call la_ctrtri('U',diag,k,a(0),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ctrmm('L','U','C',diag,k,k,cone,a(0),n + 1,a(k + 1),n + 1 &
                              )
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    call la_ctrtri('L',diag,k,a(k + 1),n + 1,info)
                    if (info > 0) return
                    call la_ctrmm('L','L','C',diag,k,k,-cone,a(k + 1),n + 1,a(0),n + &
                              1)
                    call la_ctrtri('U',diag,k,a(k),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ctrmm('R','U','N',diag,k,k,cone,a(k),n + 1,a(0),n + 1)

                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    call la_ctrtri('U',diag,k,a(k),k,info)
                    if (info > 0) return
                    call la_ctrmm('L','U','N',diag,k,k,-cone,a(k),k,a(k*(k + 1)), &
                               k)
                    call la_ctrtri('L',diag,k,a(0),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ctrmm('R','L','C',diag,k,k,cone,a(0),k,a(k*(k + 1)), &
                              k)
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    call la_ctrtri('U',diag,k,a(k*(k + 1)),k,info)
                    if (info > 0) return
                    call la_ctrmm('R','U','C',diag,k,k,-cone,a(k*(k + 1)),k,a(0), &
                               k)
                    call la_ctrtri('L',diag,k,a(k*k),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ctrmm('L','L','N',diag,k,k,cone,a(k*k),k,a(0),k)

                 end if
              end if
           end if
           return
     end subroutine la_ctftri
     !> ZTFTRI: computes the inverse of a triangular matrix A stored in RFP
     !> format.
     !> This is a Level 3 BLAS version of the algorithm.

     pure subroutine la_ztftri(transr,uplo,diag,n,a,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo,diag
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(dp),intent(inout) :: a(0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZTFTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
           else
              nisodd = .true.
           end if
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1)
                    call la_ztrtri('L',diag,n1,a(0),n,info)
                    if (info > 0) return
                    call la_ztrmm('R','L','N',diag,n2,n1,-cone,a(0),n,a(n1),n)

                    call la_ztrtri('U',diag,n2,a(n),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ztrmm('L','U','C',diag,n2,n1,cone,a(n),n,a(n1),n)

                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    call la_ztrtri('L',diag,n1,a(n2),n,info)
                    if (info > 0) return
                    call la_ztrmm('L','L','C',diag,n1,n2,-cone,a(n2),n,a(0),n)

                    call la_ztrtri('U',diag,n2,a(n1),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ztrmm('R','U','N',diag,n1,n2,cone,a(n1),n,a(0),n)

                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0), t2 -> a(1), s -> a(0+n1*n1)
                    call la_ztrtri('U',diag,n1,a(0),n1,info)
                    if (info > 0) return
                    call la_ztrmm('L','U','N',diag,n1,n2,-cone,a(0),n1,a(n1*n1), &
                              n1)
                    call la_ztrtri('L',diag,n2,a(1),n1,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ztrmm('R','L','C',diag,n1,n2,cone,a(1),n1,a(n1*n1), &
                              n1)
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0+n2*n2), t2 -> a(0+n1*n2), s -> a(0)
                    call la_ztrtri('U',diag,n1,a(n2*n2),n2,info)
                    if (info > 0) return
                    call la_ztrmm('R','U','C',diag,n2,n1,-cone,a(n2*n2),n2,a(0), &
                              n2)
                    call la_ztrtri('L',diag,n2,a(n1*n2),n2,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_ztrmm('L','L','N',diag,n2,n1,cone,a(n1*n2),n2,a(0), &
                              n2)
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    call la_ztrtri('L',diag,k,a(1),n + 1,info)
                    if (info > 0) return
                    call la_ztrmm('R','L','N',diag,k,k,-cone,a(1),n + 1,a(k + 1),n + &
                              1)
                    call la_ztrtri('U',diag,k,a(0),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ztrmm('L','U','C',diag,k,k,cone,a(0),n + 1,a(k + 1),n + 1 &
                              )
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    call la_ztrtri('L',diag,k,a(k + 1),n + 1,info)
                    if (info > 0) return
                    call la_ztrmm('L','L','C',diag,k,k,-cone,a(k + 1),n + 1,a(0),n + &
                              1)
                    call la_ztrtri('U',diag,k,a(k),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ztrmm('R','U','N',diag,k,k,cone,a(k),n + 1,a(0),n + 1)

                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    call la_ztrtri('U',diag,k,a(k),k,info)
                    if (info > 0) return
                    call la_ztrmm('L','U','N',diag,k,k,-cone,a(k),k,a(k*(k + 1)), &
                               k)
                    call la_ztrtri('L',diag,k,a(0),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ztrmm('R','L','C',diag,k,k,cone,a(0),k,a(k*(k + 1)), &
                              k)
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    call la_ztrtri('U',diag,k,a(k*(k + 1)),k,info)
                    if (info > 0) return
                    call la_ztrmm('R','U','C',diag,k,k,-cone,a(k*(k + 1)),k,a(0), &
                               k)
                    call la_ztrtri('L',diag,k,a(k*k),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_ztrmm('L','L','N',diag,k,k,cone,a(k*k),k,a(0),k)

                 end if
              end if
           end if
           return
     end subroutine la_ztftri
     !> WTFTRI: computes the inverse of a triangular matrix A stored in RFP
     !> format.
     !> This is a Level 3 BLAS version of the algorithm.

     pure subroutine la_wtftri(transr,uplo,diag,n,a,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: transr,uplo,diag
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           complex(qp),intent(inout) :: a(0:*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lower,nisodd,normaltransr
           integer(ilp) :: n1,n2,k
           ! Intrinsic Functions
           intrinsic :: mod
           ! Executable Statements
           ! test the input parameters.
           info = 0
           normaltransr = la_lsame(transr,'N')
           lower = la_lsame(uplo,'L')
           if (.not. normaltransr .and. .not. la_lsame(transr,'C')) then
              info = -1
           else if (.not. lower .and. .not. la_lsame(uplo,'U')) then
              info = -2
           else if (.not. la_lsame(diag,'N') .and. .not. la_lsame(diag,'U')) &
                     then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WTFTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! if n is odd, set nisodd = .true.
           ! if n is even, set k = n/2 and nisodd = .false.
           if (mod(n,2) == 0) then
              k = n/2
              nisodd = .false.
           else
              nisodd = .true.
           end if
           ! set n1 and n2 depending on lower
           if (lower) then
              n2 = n/2
              n1 = n - n2
           else
              n1 = n/2
              n2 = n - n1
           end if
           ! start execution: there are eight cases
           if (nisodd) then
              ! n is odd
              if (normaltransr) then
                 ! n is odd and transr = 'n'
                 if (lower) then
                   ! srpa for lower, normal and n is odd ( a(0:n-1,0:n1-1) )
                   ! t1 -> a(0,0), t2 -> a(0,1), s -> a(n1,0)
                   ! t1 -> a(0), t2 -> a(n), s -> a(n1)
                    call la_wtrtri('L',diag,n1,a(0),n,info)
                    if (info > 0) return
                    call la_wtrmm('R','L','N',diag,n2,n1,-cone,a(0),n,a(n1),n)

                    call la_wtrtri('U',diag,n2,a(n),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_wtrmm('L','U','C',diag,n2,n1,cone,a(n),n,a(n1),n)

                 else
                   ! srpa for upper, normal and n is odd ( a(0:n-1,0:n2-1)
                   ! t1 -> a(n1+1,0), t2 -> a(n1,0), s -> a(0,0)
                   ! t1 -> a(n2), t2 -> a(n1), s -> a(0)
                    call la_wtrtri('L',diag,n1,a(n2),n,info)
                    if (info > 0) return
                    call la_wtrmm('L','L','C',diag,n1,n2,-cone,a(n2),n,a(0),n)

                    call la_wtrtri('U',diag,n2,a(n1),n,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_wtrmm('R','U','N',diag,n1,n2,cone,a(n1),n,a(0),n)

                 end if
              else
                 ! n is odd and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is odd
                    ! t1 -> a(0), t2 -> a(1), s -> a(0+n1*n1)
                    call la_wtrtri('U',diag,n1,a(0),n1,info)
                    if (info > 0) return
                    call la_wtrmm('L','U','N',diag,n1,n2,-cone,a(0),n1,a(n1*n1), &
                              n1)
                    call la_wtrtri('L',diag,n2,a(1),n1,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_wtrmm('R','L','C',diag,n1,n2,cone,a(1),n1,a(n1*n1), &
                              n1)
                 else
                    ! srpa for upper, transpose and n is odd
                    ! t1 -> a(0+n2*n2), t2 -> a(0+n1*n2), s -> a(0)
                    call la_wtrtri('U',diag,n1,a(n2*n2),n2,info)
                    if (info > 0) return
                    call la_wtrmm('R','U','C',diag,n2,n1,-cone,a(n2*n2),n2,a(0), &
                              n2)
                    call la_wtrtri('L',diag,n2,a(n1*n2),n2,info)
                    if (info > 0) info = info + n1
                    if (info > 0) return
                    call la_wtrmm('L','L','N',diag,n2,n1,cone,a(n1*n2),n2,a(0), &
                              n2)
                 end if
              end if
           else
              ! n is even
              if (normaltransr) then
                 ! n is even and transr = 'n'
                 if (lower) then
                    ! srpa for lower, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(1,0), t2 -> a(0,0), s -> a(k+1,0)
                    ! t1 -> a(1), t2 -> a(0), s -> a(k+1)
                    call la_wtrtri('L',diag,k,a(1),n + 1,info)
                    if (info > 0) return
                    call la_wtrmm('R','L','N',diag,k,k,-cone,a(1),n + 1,a(k + 1),n + &
                              1)
                    call la_wtrtri('U',diag,k,a(0),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_wtrmm('L','U','C',diag,k,k,cone,a(0),n + 1,a(k + 1),n + 1 &
                              )
                 else
                    ! srpa for upper, normal, and n is even ( a(0:n,0:k-1) )
                    ! t1 -> a(k+1,0) ,  t2 -> a(k,0),   s -> a(0,0)
                    ! t1 -> a(k+1), t2 -> a(k), s -> a(0)
                    call la_wtrtri('L',diag,k,a(k + 1),n + 1,info)
                    if (info > 0) return
                    call la_wtrmm('L','L','C',diag,k,k,-cone,a(k + 1),n + 1,a(0),n + &
                              1)
                    call la_wtrtri('U',diag,k,a(k),n + 1,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_wtrmm('R','U','N',diag,k,k,cone,a(k),n + 1,a(0),n + 1)

                 end if
              else
                 ! n is even and transr = 'c'
                 if (lower) then
                    ! srpa for lower, transpose and n is even (see paper)
                    ! t1 -> b(0,1), t2 -> b(0,0), s -> b(0,k+1)
                    ! t1 -> a(0+k), t2 -> a(0+0), s -> a(0+k*(k+1)); lda=k
                    call la_wtrtri('U',diag,k,a(k),k,info)
                    if (info > 0) return
                    call la_wtrmm('L','U','N',diag,k,k,-cone,a(k),k,a(k*(k + 1)), &
                               k)
                    call la_wtrtri('L',diag,k,a(0),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_wtrmm('R','L','C',diag,k,k,cone,a(0),k,a(k*(k + 1)), &
                              k)
                 else
                    ! srpa for upper, transpose and n is even (see paper)
                    ! t1 -> b(0,k+1),     t2 -> b(0,k),   s -> b(0,0)
                    ! t1 -> a(0+k*(k+1)), t2 -> a(0+k*k), s -> a(0+0)); lda=k
                    call la_wtrtri('U',diag,k,a(k*(k + 1)),k,info)
                    if (info > 0) return
                    call la_wtrmm('R','U','C',diag,k,k,-cone,a(k*(k + 1)),k,a(0), &
                               k)
                    call la_wtrtri('L',diag,k,a(k*k),k,info)
                    if (info > 0) info = info + k
                    if (info > 0) return
                    call la_wtrmm('L','L','N',diag,k,k,cone,a(k*k),k,a(0),k)

                 end if
              end if
           end if
           return
     end subroutine la_wtftri

     !> CTPCON: estimates the reciprocal of the condition number of a packed
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_ctpcon(norm,uplo,diag,n,ap,rcond,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: ap(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(sp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CTPCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_slamch('SAFE MINIMUM')*real(max(1,n),KIND=sp)
           ! compute the norm of the triangular matrix a.
           anorm = la_clantp(norm,uplo,diag,n,ap,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_clatps(uplo,'NO TRANSPOSE',diag,normin,n,ap,work,scale, &
                              rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_clatps(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,ap,work, &
                              scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_icamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_csrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_ctpcon
     !> ZTPCON: estimates the reciprocal of the condition number of a packed
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_ztpcon(norm,uplo,diag,n,ap,rcond,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: ap(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(dp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZTPCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_dlamch('SAFE MINIMUM')*real(max(1,n),KIND=dp)
           ! compute the norm of the triangular matrix a.
           anorm = la_zlantp(norm,uplo,diag,n,ap,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_zlatps(uplo,'NO TRANSPOSE',diag,normin,n,ap,work,scale, &
                              rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_zlatps(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,ap,work, &
                              scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_izamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_zdrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_ztpcon
     !> WTPCON: estimates the reciprocal of the condition number of a packed
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_wtpcon(norm,uplo,diag,n,ap,rcond,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: ap(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(qp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WTPCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_qlamch('SAFE MINIMUM')*real(max(1,n),KIND=qp)
           ! compute the norm of the triangular matrix a.
           anorm = la_wlantp(norm,uplo,diag,n,ap,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_wlatps(uplo,'NO TRANSPOSE',diag,normin,n,ap,work,scale, &
                              rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_wlatps(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,ap,work, &
                              scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_iwamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_wqrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_wtpcon

     !> CTRCON: estimates the reciprocal of the condition number of a
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_ctrcon(norm,uplo,diag,n,a,lda,rcond,work,rwork,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(out) :: rcond
           ! Array Arguments
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(sp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(sp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CTRCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_slamch('SAFE MINIMUM')*real(max(1,n),KIND=sp)
           ! compute the norm of the triangular matrix a.
           anorm = la_clantr(norm,uplo,diag,n,n,a,lda,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_clatrs(uplo,'NO TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                               rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_clatrs(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,a,lda,work, &
                               scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_icamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_csrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_ctrcon
     !> ZTRCON: estimates the reciprocal of the condition number of a
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_ztrcon(norm,uplo,diag,n,a,lda,rcond,work,rwork,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(out) :: rcond
           ! Array Arguments
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(dp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(dp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZTRCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_dlamch('SAFE MINIMUM')*real(max(1,n),KIND=dp)
           ! compute the norm of the triangular matrix a.
           anorm = la_zlantr(norm,uplo,diag,n,n,a,lda,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_zlatrs(uplo,'NO TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                               rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_zlatrs(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,a,lda,work, &
                               scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_izamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_zdrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_ztrcon
     !> WTRCON: estimates the reciprocal of the condition number of a
     !> triangular matrix A, in either the 1-norm or the infinity-norm.
     !> The norm of A is computed and an estimate is obtained for
     !> norm(inv(A)), then the reciprocal of the condition number is
     !> computed as
     !> RCOND = 1 / ( norm(A) * norm(inv(A)) ).

     subroutine la_wtrcon(norm,uplo,diag,n,a,lda,rcond,work,rwork,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: diag,norm,uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(out) :: rcond
           ! Array Arguments
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: nounit,onenrm,upper
           character :: normin
           integer(ilp) :: ix,kase,kase1
           real(qp) :: ainvnm,anorm,scale,smlnum,xnorm
           complex(qp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           onenrm = norm == '1' .or. la_lsame(norm,'O')
           nounit = la_lsame(diag,'N')
           if (.not. onenrm .and. .not. la_lsame(norm,'I')) then
              info = -1
           else if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -2
           else if (.not. nounit .and. .not. la_lsame(diag,'U')) then
              info = -3
           else if (n < 0) then
              info = -4
           else if (lda < max(1,n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WTRCON',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) then
              rcond = one
              return
           end if
           rcond = zero
           smlnum = la_qlamch('SAFE MINIMUM')*real(max(1,n),KIND=qp)
           ! compute the norm of the triangular matrix a.
           anorm = la_wlantr(norm,uplo,diag,n,n,a,lda,rwork)
           ! continue only if anorm > 0.
           if (anorm > zero) then
              ! estimate the norm of the inverse of a.
              ainvnm = zero
              normin = 'N'
              if (onenrm) then
                 kase1 = 1
              else
                 kase1 = 2
              end if
              kase = 0
              10 continue
              call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
              if (kase /= 0) then
                 if (kase == kase1) then
                    ! multiply by inv(a).
                    call la_wlatrs(uplo,'NO TRANSPOSE',diag,normin,n,a,lda,work,scale, &
                               rwork,info)
                 else
                    ! multiply by inv(a**h).
                    call la_wlatrs(uplo,'CONJUGATE TRANSPOSE',diag,normin,n,a,lda,work, &
                               scale,rwork,info)
                 end if
                 normin = 'Y'
                 ! multiply by 1/scale if doing so will not cause overflow.
                 if (scale /= one) then
                    ix = la_iwamax(n,work,1)
                    xnorm = cabs1(work(ix))
                    if (scale < xnorm*smlnum .or. scale == zero) go to 20
                    call la_wqrscl(n,scale,work,1)
                 end if
                 go to 10
              end if
              ! compute the estimate of the reciprocal condition number.
              if (ainvnm /= zero) rcond = (one/anorm)/ainvnm
           end if
           20 continue
           return
     end subroutine la_wtrcon

end module la_lapack_solve_tri_comp
