!> Hermitian indefinite components: rook, Aasen and rank-k variants
module la_lapack_solve_ldl_comp4
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_pac
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_scalar
     use la_lapack_solve_aux
     use la_lapack_solve_ldl_comp3
     use la_lapack_solve_lu
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slaqsy
     public :: la_dlaqsy
     public :: la_qlaqsy
     public :: la_chetf2_rk
     public :: la_chetf2_rook
     public :: la_chetri_rook
     public :: la_clahef_rk
     public :: la_clahef_rook
     public :: la_claqsy
     public :: la_chetrf_rk
     public :: la_chetrf_rook
     public :: la_chetrs_aa
     public :: la_chetrs_rook
     public :: la_checon_rook
     public :: la_chprfs
     public :: la_chetrf_aa
     public :: la_clahef_aa
     public :: la_zhetf2_rk
     public :: la_zhetf2_rook
     public :: la_zhetri_rook
     public :: la_zlahef_rk
     public :: la_zlahef_rook
     public :: la_zlaqsy
     public :: la_zhetrf_rk
     public :: la_zhetrf_rook
     public :: la_zhetrs_aa
     public :: la_zhetrs_rook
     public :: la_zhecon_rook
     public :: la_zhprfs
     public :: la_zhetrf_aa
     public :: la_zlahef_aa
     public :: la_whetf2_rk
     public :: la_whetf2_rook
     public :: la_whetri_rook
     public :: la_wlahef_rk
     public :: la_wlahef_rook
     public :: la_wlaqsy
     public :: la_whetrf_rk
     public :: la_whetrf_rook
     public :: la_whetrs_aa
     public :: la_whetrs_rook
     public :: la_whecon_rook
     public :: la_whprfs
     public :: la_whetrf_aa
     public :: la_wlahef_aa

     contains

     !> SLAQSY: equilibrates a symmetric matrix A using the scaling factors
     !> in the vector S.

     pure subroutine la_slaqsy(uplo,n,a,lda,s,scond,amax,equed)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(in) :: amax,scond
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: thresh = 0.1e+0_sp

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: cj,large,small
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = 1,j
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,n
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_slaqsy
     !> DLAQSY: equilibrates a symmetric matrix A using the scaling factors
     !> in the vector S.

     pure subroutine la_dlaqsy(uplo,n,a,lda,s,scond,amax,equed)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(in) :: amax,scond
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: thresh = 0.1e+0_dp

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: cj,large,small
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = 1,j
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,n
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_dlaqsy
     !> QLAQSY: equilibrates a symmetric matrix A using the scaling factors
     !> in the vector S.

     pure subroutine la_qlaqsy(uplo,n,a,lda,s,scond,amax,equed)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(in) :: amax,scond
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: s(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: thresh = 0.1e+0_qp

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: cj,large,small
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = 1,j
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,n
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_qlaqsy

     !> CHETF2_RK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
     !> A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.
     !> For more information see Further Details section.

     pure subroutine la_chetf2_rk(uplo,n,a,lda,e,ipiv,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: e(*)
        ! ======================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: done,upper
           integer(ilp) :: i,ii,imax,itemp,j,jmax,k,kk,kp,kstep,p
           real(sp) :: absakk,alpha,colmax,d,d11,d22,r1,stemp,rowmax,tt,sfmin
           complex(sp) :: d12,d21,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=sp)) + abs(aimag(z))
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
              call la_xerbla('CHETF2_RK',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_slamch('S')
           if (upper) then
              ! factorize a as u*d*u**h using the upper triangle of a
              ! initialize the first entry of array e, where superdiagonal
              ! elements of d are stored
              e(1) = czero
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 34
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_icamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero)) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=sp)
                 ! set e( k ) to zero
                 if (k > 1) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    12 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_icamax(k - imax,a(imax,imax + 1),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_icamax(imax - 1,a(1,imax),1)
                          stemp = cabs1(a(itemp,imax))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the leading submatrix a(1:k,1:k)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p > 1) call la_cswap(p - 1,a(1,k),1,a(1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = p + 1,k - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=sp)
                    a(k,k) = real(a(p,p),KIND=sp)
                    a(p,p) = r1
                    ! convert upper triangle of a into u form by applying
                    ! the interchanges in columns k+1:n.
                    if (k < n) call la_cswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the leading submatrix a(1:k,1:k)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp > 1) call la_cswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! (2) swap and conjugate middle parts
                    do j = kp + 1,kk - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=sp)
                    a(kk,kk) = real(a(kp,kp),KIND=sp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=sp)
                       ! (5) swap row elements
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                    ! convert upper triangle of a into u form by applying
                    ! the interchanges in columns k+1:n.
                    if (k < n) call la_cswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=sp)
                    if (kstep == 2) a(k - 1,k - 1) = real(a(k - 1,k - 1),KIND=sp)
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    if (k > 1) then
                       ! perform a rank-1 update of a(1:k-1,1:k-1) and
                       ! store u(k) in column k
                       if (abs(real(a(k,k),KIND=sp)) >= sfmin) then
                          ! perform a rank-1 update of a(1:k-1,1:k-1) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*1/d(k)*w(k)**t
                          d11 = one/real(a(k,k),KIND=sp)
                          call la_cher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                          ! store u(k) in column k
                          call la_csscal(k - 1,d11,a(1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=sp)
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_cher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                       end if
                       ! store the superdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( ( a(k-1)a(k) )*inv(d(k)) ) * ( a(k-1)a(k) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k > 2) then
                       ! d = |a12|
                       d = la_slapy2(real(a(k - 1,k),KIND=sp),aimag(a(k - 1,k)))

                       d11 = real(a(k,k)/d,KIND=sp)
                       d22 = real(a(k - 1,k - 1)/d,KIND=sp)
                       d12 = a(k - 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k - 2,1,-1
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wkm1 = tt*(d11*a(j,k - 1) - conjg(d12)*a(j,k))
                          wk = tt*(d22*a(j,k) - d12*a(j,k - 1))
                          ! perform a rank-2 update of a(1:k-2,1:k-2)
                          do i = j,1,-1
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k - 1) &
                                       /d)*conjg(wkm1)
                          end do
                          ! store u(k) and u(k-1) in cols k and k-1 for row j
                          a(j,k) = wk/d
                          a(j,k - 1) = wkm1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=sp),zero,KIND=sp)
                       end do
                    end if
                    ! copy superdiagonal elements of d(k) to e(k) and
                    ! zero out superdiagonal entry of a
                    e(k) = a(k - 1,k)
                    e(k - 1) = czero
                    a(k - 1,k) = czero
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              34 continue
           else
              ! factorize a as l*d*l**h using the lower triangle of a
              ! initialize the unused last entry of the subdiagonal array e.
              e(n) = czero
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 64
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_icamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=sp)
                 ! set e( k ) to zero
                 if (k < n) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    42 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_icamax(imax - k,a(imax,k),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_icamax(n - imax,a(imax + 1,imax),1)
                          stemp = cabs1(a(itemp,imax))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 42
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the trailing submatrix a(k:n,k:n)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p < n) call la_cswap(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = k + 1,p - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=sp)
                    a(k,k) = real(a(p,p),KIND=sp)
                    a(p,p) = r1
                    ! convert lower triangle of a into l form by applying
                    ! the interchanges in columns 1:k-1.
                    if (k > 1) call la_cswap(k - 1,a(k,1),lda,a(p,1),lda)
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the trailing submatrix a(k:n,k:n)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp < n) call la_cswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! (2) swap and conjugate middle parts
                    do j = kk + 1,kp - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=sp)
                    a(kk,kk) = real(a(kp,kp),KIND=sp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=sp)
                       ! (5) swap row elements
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                    ! convert lower triangle of a into l form by applying
                    ! the interchanges in columns 1:k-1.
                    if (k > 1) call la_cswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=sp)
                    if (kstep == 2) a(k + 1,k + 1) = real(a(k + 1,k + 1),KIND=sp)
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of a now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) and
                       ! store l(k) in column k
                       ! handle division by a small number
                       if (abs(real(a(k,k),KIND=sp)) >= sfmin) then
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                          d11 = one/real(a(k,k),KIND=sp)
                          call la_cher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                          ! store l(k) in column k
                          call la_csscal(n - k,d11,a(k + 1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=sp)
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_cher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       end if
                       ! store the subdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! perform a rank-2 update of a(k+2:n,k+2:n) as
                    ! a := a - ( l(k) l(k+1) ) * d(k) * ( l(k) l(k+1) )**t
                       ! = a - ( ( a(k)a(k+1) )*inv(d(k) ) * ( a(k)a(k+1) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k < n - 1) then
                       ! d = |a21|
                       d = la_slapy2(real(a(k + 1,k),KIND=sp),aimag(a(k + 1,k)))

                       d11 = real(a(k + 1,k + 1),KIND=sp)/d
                       d22 = real(a(k,k),KIND=sp)/d
                       d21 = a(k + 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k + 2,n
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wk = tt*(d11*a(j,k) - d21*a(j,k + 1))
                          wkp1 = tt*(d22*a(j,k + 1) - conjg(d21)*a(j,k))
                          ! perform a rank-2 update of a(k+2:n,k+2:n)
                          do i = j,n
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k + 1) &
                                       /d)*conjg(wkp1)
                          end do
                          ! store l(k) and l(k+1) in cols k and k+1 for row j
                          a(j,k) = wk/d
                          a(j,k + 1) = wkp1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=sp),zero,KIND=sp)
                       end do
                    end if
                    ! copy subdiagonal elements of d(k) to e(k) and
                    ! zero out subdiagonal entry of a
                    e(k) = a(k + 1,k)
                    e(k + 1) = czero
                    a(k + 1,k) = czero
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
              64 continue
           end if
           return
     end subroutine la_chetf2_rk
     !> ZHETF2_RK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
     !> A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.
     !> For more information see Further Details section.

     pure subroutine la_zhetf2_rk(uplo,n,a,lda,e,ipiv,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: e(*)
        ! ======================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: done,upper
           integer(ilp) :: i,ii,imax,itemp,j,jmax,k,kk,kp,kstep,p
           real(dp) :: absakk,alpha,colmax,d,d11,d22,r1,dtemp,rowmax,tt,sfmin
           complex(dp) :: d12,d21,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=dp)) + abs(aimag(z))
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
              call la_xerbla('ZHETF2_RK',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_dlamch('S')
           if (upper) then
              ! factorize a as u*d*u**h using the upper triangle of a
              ! initialize the first entry of array e, where superdiagonal
              ! elements of d are stored
              e(1) = czero
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 34
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_izamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero)) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=dp)
                 ! set e( k ) to zero
                 if (k > 1) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    12 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_izamax(k - imax,a(imax,imax + 1),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_izamax(imax - 1,a(1,imax),1)
                          dtemp = cabs1(a(itemp,imax))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the leading submatrix a(1:k,1:k)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p > 1) call la_zswap(p - 1,a(1,k),1,a(1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = p + 1,k - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=dp)
                    a(k,k) = real(a(p,p),KIND=dp)
                    a(p,p) = r1
                    ! convert upper triangle of a into u form by applying
                    ! the interchanges in columns k+1:n.
                    if (k < n) call la_zswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the leading submatrix a(1:k,1:k)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp > 1) call la_zswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! (2) swap and conjugate middle parts
                    do j = kp + 1,kk - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=dp)
                    a(kk,kk) = real(a(kp,kp),KIND=dp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=dp)
                       ! (5) swap row elements
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                    ! convert upper triangle of a into u form by applying
                    ! the interchanges in columns k+1:n.
                    if (k < n) call la_zswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=dp)
                    if (kstep == 2) a(k - 1,k - 1) = real(a(k - 1,k - 1),KIND=dp)
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    if (k > 1) then
                       ! perform a rank-1 update of a(1:k-1,1:k-1) and
                       ! store u(k) in column k
                       if (abs(real(a(k,k),KIND=dp)) >= sfmin) then
                          ! perform a rank-1 update of a(1:k-1,1:k-1) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*1/d(k)*w(k)**t
                          d11 = one/real(a(k,k),KIND=dp)
                          call la_zher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                          ! store u(k) in column k
                          call la_zdscal(k - 1,d11,a(1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=dp)
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_zher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                       end if
                       ! store the superdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( ( a(k-1)a(k) )*inv(d(k)) ) * ( a(k-1)a(k) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k > 2) then
                       ! d = |a12|
                       d = la_dlapy2(real(a(k - 1,k),KIND=dp),aimag(a(k - 1,k)))

                       d11 = real(a(k,k)/d,KIND=dp)
                       d22 = real(a(k - 1,k - 1)/d,KIND=dp)
                       d12 = a(k - 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k - 2,1,-1
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wkm1 = tt*(d11*a(j,k - 1) - conjg(d12)*a(j,k))
                          wk = tt*(d22*a(j,k) - d12*a(j,k - 1))
                          ! perform a rank-2 update of a(1:k-2,1:k-2)
                          do i = j,1,-1
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k - 1) &
                                       /d)*conjg(wkm1)
                          end do
                          ! store u(k) and u(k-1) in cols k and k-1 for row j
                          a(j,k) = wk/d
                          a(j,k - 1) = wkm1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=dp),zero,KIND=dp)
                       end do
                    end if
                    ! copy superdiagonal elements of d(k) to e(k) and
                    ! zero out superdiagonal entry of a
                    e(k) = a(k - 1,k)
                    e(k - 1) = czero
                    a(k - 1,k) = czero
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              34 continue
           else
              ! factorize a as l*d*l**h using the lower triangle of a
              ! initialize the unused last entry of the subdiagonal array e.
              e(n) = czero
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 64
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_izamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=dp)
                 ! set e( k ) to zero
                 if (k < n) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    42 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_izamax(imax - k,a(imax,k),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_izamax(n - imax,a(imax + 1,imax),1)
                          dtemp = cabs1(a(itemp,imax))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 42
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the trailing submatrix a(k:n,k:n)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p < n) call la_zswap(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = k + 1,p - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=dp)
                    a(k,k) = real(a(p,p),KIND=dp)
                    a(p,p) = r1
                    ! convert lower triangle of a into l form by applying
                    ! the interchanges in columns 1:k-1.
                    if (k > 1) call la_zswap(k - 1,a(k,1),lda,a(p,1),lda)
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the trailing submatrix a(k:n,k:n)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp < n) call la_zswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! (2) swap and conjugate middle parts
                    do j = kk + 1,kp - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=dp)
                    a(kk,kk) = real(a(kp,kp),KIND=dp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=dp)
                       ! (5) swap row elements
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                    ! convert lower triangle of a into l form by applying
                    ! the interchanges in columns 1:k-1.
                    if (k > 1) call la_zswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=dp)
                    if (kstep == 2) a(k + 1,k + 1) = real(a(k + 1,k + 1),KIND=dp)
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of a now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) and
                       ! store l(k) in column k
                       ! handle division by a small number
                       if (abs(real(a(k,k),KIND=dp)) >= sfmin) then
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                          d11 = one/real(a(k,k),KIND=dp)
                          call la_zher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                          ! store l(k) in column k
                          call la_zdscal(n - k,d11,a(k + 1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=dp)
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_zher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       end if
                       ! store the subdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! perform a rank-2 update of a(k+2:n,k+2:n) as
                    ! a := a - ( l(k) l(k+1) ) * d(k) * ( l(k) l(k+1) )**t
                       ! = a - ( ( a(k)a(k+1) )*inv(d(k) ) * ( a(k)a(k+1) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k < n - 1) then
                       ! d = |a21|
                       d = la_dlapy2(real(a(k + 1,k),KIND=dp),aimag(a(k + 1,k)))

                       d11 = real(a(k + 1,k + 1),KIND=dp)/d
                       d22 = real(a(k,k),KIND=dp)/d
                       d21 = a(k + 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k + 2,n
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wk = tt*(d11*a(j,k) - d21*a(j,k + 1))
                          wkp1 = tt*(d22*a(j,k + 1) - conjg(d21)*a(j,k))
                          ! perform a rank-2 update of a(k+2:n,k+2:n)
                          do i = j,n
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k + 1) &
                                       /d)*conjg(wkp1)
                          end do
                          ! store l(k) and l(k+1) in cols k and k+1 for row j
                          a(j,k) = wk/d
                          a(j,k + 1) = wkp1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=dp),zero,KIND=dp)
                       end do
                    end if
                    ! copy subdiagonal elements of d(k) to e(k) and
                    ! zero out subdiagonal entry of a
                    e(k) = a(k + 1,k)
                    e(k + 1) = czero
                    a(k + 1,k) = czero
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
              64 continue
           end if
           return
     end subroutine la_zhetf2_rk
     !> WHETF2_RK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
     !> A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.
     !> For more information see Further Details section.

     pure subroutine la_whetf2_rk(uplo,n,a,lda,e,ipiv,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: e(*)
        ! ======================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: done,upper
           integer(ilp) :: i,ii,imax,itemp,j,jmax,k,kk,kp,kstep,p
           real(qp) :: absakk,alpha,colmax,d,d11,d22,r1,qtemp,rowmax,tt,sfmin
           complex(qp) :: d12,d21,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=qp)) + abs(aimag(z))
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
              call la_xerbla('WHETF2_RK',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_qlamch('S')
           if (upper) then
              ! factorize a as u*d*u**h using the upper triangle of a
              ! initialize the first entry of array e, where superdiagonal
              ! elements of d are stored
              e(1) = czero
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 34
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iwamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero)) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=qp)
                 ! set e( k ) to zero
                 if (k > 1) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    12 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_iwamax(k - imax,a(imax,imax + 1),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_iwamax(imax - 1,a(1,imax),1)
                          qtemp = cabs1(a(itemp,imax))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the leading submatrix a(1:k,1:k)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p > 1) call la_wswap(p - 1,a(1,k),1,a(1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = p + 1,k - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=qp)
                    a(k,k) = real(a(p,p),KIND=qp)
                    a(p,p) = r1
                    ! convert upper triangle of a into u form by applying
                    ! the interchanges in columns k+1:n.
                    if (k < n) call la_wswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the leading submatrix a(1:k,1:k)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp > 1) call la_wswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! (2) swap and conjugate middle parts
                    do j = kp + 1,kk - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=qp)
                    a(kk,kk) = real(a(kp,kp),KIND=qp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=qp)
                       ! (5) swap row elements
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                    ! convert upper triangle of a into u form by applying
                    ! the interchanges in columns k+1:n.
                    if (k < n) call la_wswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=qp)
                    if (kstep == 2) a(k - 1,k - 1) = real(a(k - 1,k - 1),KIND=qp)
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    if (k > 1) then
                       ! perform a rank-1 update of a(1:k-1,1:k-1) and
                       ! store u(k) in column k
                       if (abs(real(a(k,k),KIND=qp)) >= sfmin) then
                          ! perform a rank-1 update of a(1:k-1,1:k-1) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*1/d(k)*w(k)**t
                          d11 = one/real(a(k,k),KIND=qp)
                          call la_wher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                          ! store u(k) in column k
                          call la_wqscal(k - 1,d11,a(1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=qp)
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_wher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                       end if
                       ! store the superdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( ( a(k-1)a(k) )*inv(d(k)) ) * ( a(k-1)a(k) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k > 2) then
                       ! d = |a12|
                       d = la_qlapy2(real(a(k - 1,k),KIND=qp),aimag(a(k - 1,k)))

                       d11 = real(a(k,k)/d,KIND=qp)
                       d22 = real(a(k - 1,k - 1)/d,KIND=qp)
                       d12 = a(k - 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k - 2,1,-1
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wkm1 = tt*(d11*a(j,k - 1) - conjg(d12)*a(j,k))
                          wk = tt*(d22*a(j,k) - d12*a(j,k - 1))
                          ! perform a rank-2 update of a(1:k-2,1:k-2)
                          do i = j,1,-1
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k - 1) &
                                       /d)*conjg(wkm1)
                          end do
                          ! store u(k) and u(k-1) in cols k and k-1 for row j
                          a(j,k) = wk/d
                          a(j,k - 1) = wkm1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=qp),zero,KIND=qp)
                       end do
                    end if
                    ! copy superdiagonal elements of d(k) to e(k) and
                    ! zero out superdiagonal entry of a
                    e(k) = a(k - 1,k)
                    e(k - 1) = czero
                    a(k - 1,k) = czero
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              34 continue
           else
              ! factorize a as l*d*l**h using the lower triangle of a
              ! initialize the unused last entry of the subdiagonal array e.
              e(n) = czero
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 64
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iwamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=qp)
                 ! set e( k ) to zero
                 if (k < n) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    42 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_iwamax(imax - k,a(imax,k),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_iwamax(n - imax,a(imax + 1,imax),1)
                          qtemp = cabs1(a(itemp,imax))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 42
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the trailing submatrix a(k:n,k:n)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p < n) call la_wswap(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = k + 1,p - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=qp)
                    a(k,k) = real(a(p,p),KIND=qp)
                    a(p,p) = r1
                    ! convert lower triangle of a into l form by applying
                    ! the interchanges in columns 1:k-1.
                    if (k > 1) call la_wswap(k - 1,a(k,1),lda,a(p,1),lda)
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the trailing submatrix a(k:n,k:n)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp < n) call la_wswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! (2) swap and conjugate middle parts
                    do j = kk + 1,kp - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=qp)
                    a(kk,kk) = real(a(kp,kp),KIND=qp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=qp)
                       ! (5) swap row elements
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                    ! convert lower triangle of a into l form by applying
                    ! the interchanges in columns 1:k-1.
                    if (k > 1) call la_wswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=qp)
                    if (kstep == 2) a(k + 1,k + 1) = real(a(k + 1,k + 1),KIND=qp)
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of a now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) and
                       ! store l(k) in column k
                       ! handle division by a small number
                       if (abs(real(a(k,k),KIND=qp)) >= sfmin) then
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                          d11 = one/real(a(k,k),KIND=qp)
                          call la_wher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                          ! store l(k) in column k
                          call la_wqscal(n - k,d11,a(k + 1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=qp)
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_wher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       end if
                       ! store the subdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! perform a rank-2 update of a(k+2:n,k+2:n) as
                    ! a := a - ( l(k) l(k+1) ) * d(k) * ( l(k) l(k+1) )**t
                       ! = a - ( ( a(k)a(k+1) )*inv(d(k) ) * ( a(k)a(k+1) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k < n - 1) then
                       ! d = |a21|
                       d = la_qlapy2(real(a(k + 1,k),KIND=qp),aimag(a(k + 1,k)))

                       d11 = real(a(k + 1,k + 1),KIND=qp)/d
                       d22 = real(a(k,k),KIND=qp)/d
                       d21 = a(k + 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k + 2,n
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wk = tt*(d11*a(j,k) - d21*a(j,k + 1))
                          wkp1 = tt*(d22*a(j,k + 1) - conjg(d21)*a(j,k))
                          ! perform a rank-2 update of a(k+2:n,k+2:n)
                          do i = j,n
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k + 1) &
                                       /d)*conjg(wkp1)
                          end do
                          ! store l(k) and l(k+1) in cols k and k+1 for row j
                          a(j,k) = wk/d
                          a(j,k + 1) = wkp1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=qp),zero,KIND=qp)
                       end do
                    end if
                    ! copy subdiagonal elements of d(k) to e(k) and
                    ! zero out subdiagonal entry of a
                    e(k) = a(k + 1,k)
                    e(k + 1) = czero
                    a(k + 1,k) = czero
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
              64 continue
           end if
           return
     end subroutine la_whetf2_rk

     !> CHETF2_ROOK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:
     !> A = U*D*U**H  or  A = L*D*L**H
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**H is the conjugate transpose of U, and D is
     !> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_chetf2_rook(uplo,n,a,lda,ipiv,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
        ! ======================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: done,upper
           integer(ilp) :: i,ii,imax,itemp,j,jmax,k,kk,kp,kstep,p
           real(sp) :: absakk,alpha,colmax,d,d11,d22,r1,stemp,rowmax,tt,sfmin
           complex(sp) :: d12,d21,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,aimag,cmplx,conjg,max,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=sp)) + abs(aimag(z))
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
              call la_xerbla('CHETF2_ROOK',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_slamch('S')
           if (upper) then
              ! factorize a as u*d*u**h using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_icamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero)) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=sp)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    12 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_icamax(k - imax,a(imax,imax + 1),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_icamax(imax - 1,a(1,imax),1)
                          stemp = cabs1(a(itemp,imax))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the leading submatrix a(1:k,1:k)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p > 1) call la_cswap(p - 1,a(1,k),1,a(1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = p + 1,k - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=sp)
                    a(k,k) = real(a(p,p),KIND=sp)
                    a(p,p) = r1
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the leading submatrix a(1:k,1:k)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp > 1) call la_cswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! (2) swap and conjugate middle parts
                    do j = kp + 1,kk - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=sp)
                    a(kk,kk) = real(a(kp,kp),KIND=sp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=sp)
                       ! (5) swap row elements
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=sp)
                    if (kstep == 2) a(k - 1,k - 1) = real(a(k - 1,k - 1),KIND=sp)
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    if (k > 1) then
                       ! perform a rank-1 update of a(1:k-1,1:k-1) and
                       ! store u(k) in column k
                       if (abs(real(a(k,k),KIND=sp)) >= sfmin) then
                          ! perform a rank-1 update of a(1:k-1,1:k-1) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*1/d(k)*w(k)**t
                          d11 = one/real(a(k,k),KIND=sp)
                          call la_cher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                          ! store u(k) in column k
                          call la_csscal(k - 1,d11,a(1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=sp)
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_cher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                       end if
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( ( a(k-1)a(k) )*inv(d(k)) ) * ( a(k-1)a(k) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k > 2) then
                       ! d = |a12|
                       d = la_slapy2(real(a(k - 1,k),KIND=sp),aimag(a(k - 1,k)))

                       d11 = real(a(k,k)/d,KIND=sp)
                       d22 = real(a(k - 1,k - 1)/d,KIND=sp)
                       d12 = a(k - 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k - 2,1,-1
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wkm1 = tt*(d11*a(j,k - 1) - conjg(d12)*a(j,k))
                          wk = tt*(d22*a(j,k) - d12*a(j,k - 1))
                          ! perform a rank-2 update of a(1:k-2,1:k-2)
                          do i = j,1,-1
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k - 1) &
                                       /d)*conjg(wkm1)
                          end do
                          ! store u(k) and u(k-1) in cols k and k-1 for row j
                          a(j,k) = wk/d
                          a(j,k - 1) = wkm1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=sp),zero,KIND=sp)
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**h using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_icamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=sp)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    42 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_icamax(imax - k,a(imax,k),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_icamax(n - imax,a(imax + 1,imax),1)
                          stemp = cabs1(a(itemp,imax))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 42
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the trailing submatrix a(k:n,k:n)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p < n) call la_cswap(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = k + 1,p - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=sp)
                    a(k,k) = real(a(p,p),KIND=sp)
                    a(p,p) = r1
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the trailing submatrix a(k:n,k:n)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp < n) call la_cswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! (2) swap and conjugate middle parts
                    do j = kk + 1,kp - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=sp)
                    a(kk,kk) = real(a(kp,kp),KIND=sp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=sp)
                       ! (5) swap row elements
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=sp)
                    if (kstep == 2) a(k + 1,k + 1) = real(a(k + 1,k + 1),KIND=sp)
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of a now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) and
                       ! store l(k) in column k
                       ! handle division by a small number
                       if (abs(real(a(k,k),KIND=sp)) >= sfmin) then
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                          d11 = one/real(a(k,k),KIND=sp)
                          call la_cher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                          ! store l(k) in column k
                          call la_csscal(n - k,d11,a(k + 1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=sp)
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_cher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       end if
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! perform a rank-2 update of a(k+2:n,k+2:n) as
                    ! a := a - ( l(k) l(k+1) ) * d(k) * ( l(k) l(k+1) )**t
                       ! = a - ( ( a(k)a(k+1) )*inv(d(k) ) * ( a(k)a(k+1) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k < n - 1) then
                       ! d = |a21|
                       d = la_slapy2(real(a(k + 1,k),KIND=sp),aimag(a(k + 1,k)))

                       d11 = real(a(k + 1,k + 1),KIND=sp)/d
                       d22 = real(a(k,k),KIND=sp)/d
                       d21 = a(k + 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k + 2,n
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wk = tt*(d11*a(j,k) - d21*a(j,k + 1))
                          wkp1 = tt*(d22*a(j,k + 1) - conjg(d21)*a(j,k))
                          ! perform a rank-2 update of a(k+2:n,k+2:n)
                          do i = j,n
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k + 1) &
                                       /d)*conjg(wkp1)
                          end do
                          ! store l(k) and l(k+1) in cols k and k+1 for row j
                          a(j,k) = wk/d
                          a(j,k + 1) = wkp1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=sp),zero,KIND=sp)
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_chetf2_rook
     !> ZHETF2_ROOK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:
     !> A = U*D*U**H  or  A = L*D*L**H
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**H is the conjugate transpose of U, and D is
     !> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_zhetf2_rook(uplo,n,a,lda,ipiv,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
        ! ======================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: done,upper
           integer(ilp) :: i,ii,imax,itemp,j,jmax,k,kk,kp,kstep,p
           real(dp) :: absakk,alpha,colmax,d,d11,d22,r1,dtemp,rowmax,tt,sfmin
           complex(dp) :: d12,d21,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=dp)) + abs(aimag(z))
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
              call la_xerbla('ZHETF2_ROOK',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_dlamch('S')
           if (upper) then
              ! factorize a as u*d*u**h using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_izamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero)) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=dp)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    12 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_izamax(k - imax,a(imax,imax + 1),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_izamax(imax - 1,a(1,imax),1)
                          dtemp = cabs1(a(itemp,imax))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the leading submatrix a(1:k,1:k)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p > 1) call la_zswap(p - 1,a(1,k),1,a(1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = p + 1,k - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=dp)
                    a(k,k) = real(a(p,p),KIND=dp)
                    a(p,p) = r1
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the leading submatrix a(1:k,1:k)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp > 1) call la_zswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! (2) swap and conjugate middle parts
                    do j = kp + 1,kk - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=dp)
                    a(kk,kk) = real(a(kp,kp),KIND=dp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=dp)
                       ! (5) swap row elements
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=dp)
                    if (kstep == 2) a(k - 1,k - 1) = real(a(k - 1,k - 1),KIND=dp)
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    if (k > 1) then
                       ! perform a rank-1 update of a(1:k-1,1:k-1) and
                       ! store u(k) in column k
                       if (abs(real(a(k,k),KIND=dp)) >= sfmin) then
                          ! perform a rank-1 update of a(1:k-1,1:k-1) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*1/d(k)*w(k)**t
                          d11 = one/real(a(k,k),KIND=dp)
                          call la_zher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                          ! store u(k) in column k
                          call la_zdscal(k - 1,d11,a(1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=dp)
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_zher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                       end if
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( ( a(k-1)a(k) )*inv(d(k)) ) * ( a(k-1)a(k) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k > 2) then
                       ! d = |a12|
                       d = la_dlapy2(real(a(k - 1,k),KIND=dp),aimag(a(k - 1,k)))

                       d11 = real(a(k,k)/d,KIND=dp)
                       d22 = real(a(k - 1,k - 1)/d,KIND=dp)
                       d12 = a(k - 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k - 2,1,-1
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wkm1 = tt*(d11*a(j,k - 1) - conjg(d12)*a(j,k))
                          wk = tt*(d22*a(j,k) - d12*a(j,k - 1))
                          ! perform a rank-2 update of a(1:k-2,1:k-2)
                          do i = j,1,-1
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k - 1) &
                                       /d)*conjg(wkm1)
                          end do
                          ! store u(k) and u(k-1) in cols k and k-1 for row j
                          a(j,k) = wk/d
                          a(j,k - 1) = wkm1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=dp),zero,KIND=dp)
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**h using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_izamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=dp)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    42 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_izamax(imax - k,a(imax,k),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_izamax(n - imax,a(imax + 1,imax),1)
                          dtemp = cabs1(a(itemp,imax))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 42
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the trailing submatrix a(k:n,k:n)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p < n) call la_zswap(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = k + 1,p - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=dp)
                    a(k,k) = real(a(p,p),KIND=dp)
                    a(p,p) = r1
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the trailing submatrix a(k:n,k:n)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp < n) call la_zswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! (2) swap and conjugate middle parts
                    do j = kk + 1,kp - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=dp)
                    a(kk,kk) = real(a(kp,kp),KIND=dp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=dp)
                       ! (5) swap row elements
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=dp)
                    if (kstep == 2) a(k + 1,k + 1) = real(a(k + 1,k + 1),KIND=dp)
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of a now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) and
                       ! store l(k) in column k
                       ! handle division by a small number
                       if (abs(real(a(k,k),KIND=dp)) >= sfmin) then
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                          d11 = one/real(a(k,k),KIND=dp)
                          call la_zher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                          ! store l(k) in column k
                          call la_zdscal(n - k,d11,a(k + 1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=dp)
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_zher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       end if
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! perform a rank-2 update of a(k+2:n,k+2:n) as
                    ! a := a - ( l(k) l(k+1) ) * d(k) * ( l(k) l(k+1) )**t
                       ! = a - ( ( a(k)a(k+1) )*inv(d(k) ) * ( a(k)a(k+1) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k < n - 1) then
                       ! d = |a21|
                       d = la_dlapy2(real(a(k + 1,k),KIND=dp),aimag(a(k + 1,k)))

                       d11 = real(a(k + 1,k + 1),KIND=dp)/d
                       d22 = real(a(k,k),KIND=dp)/d
                       d21 = a(k + 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k + 2,n
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wk = tt*(d11*a(j,k) - d21*a(j,k + 1))
                          wkp1 = tt*(d22*a(j,k + 1) - conjg(d21)*a(j,k))
                          ! perform a rank-2 update of a(k+2:n,k+2:n)
                          do i = j,n
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k + 1) &
                                       /d)*conjg(wkp1)
                          end do
                          ! store l(k) and l(k+1) in cols k and k+1 for row j
                          a(j,k) = wk/d
                          a(j,k + 1) = wkp1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=dp),zero,KIND=dp)
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_zhetf2_rook
     !> WHETF2_ROOK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method:
     !> A = U*D*U**H  or  A = L*D*L**H
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**H is the conjugate transpose of U, and D is
     !> Hermitian and block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_whetf2_rook(uplo,n,a,lda,ipiv,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
        ! ======================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: done,upper
           integer(ilp) :: i,ii,imax,itemp,j,jmax,k,kk,kp,kstep,p
           real(qp) :: absakk,alpha,colmax,d,d11,d22,r1,qtemp,rowmax,tt,sfmin
           complex(qp) :: d12,d21,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,conjg,aimag,max,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=qp)) + abs(aimag(z))
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
              call la_xerbla('WHETF2_ROOK',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_qlamch('S')
           if (upper) then
              ! factorize a as u*d*u**h using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iwamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero)) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=qp)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    12 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_iwamax(k - imax,a(imax,imax + 1),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_iwamax(imax - 1,a(1,imax),1)
                          qtemp = cabs1(a(itemp,imax))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the leading submatrix a(1:k,1:k)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p > 1) call la_wswap(p - 1,a(1,k),1,a(1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = p + 1,k - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=qp)
                    a(k,k) = real(a(p,p),KIND=qp)
                    a(p,p) = r1
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the leading submatrix a(1:k,1:k)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp > 1) call la_wswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! (2) swap and conjugate middle parts
                    do j = kp + 1,kk - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=qp)
                    a(kk,kk) = real(a(kp,kp),KIND=qp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=qp)
                       ! (5) swap row elements
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=qp)
                    if (kstep == 2) a(k - 1,k - 1) = real(a(k - 1,k - 1),KIND=qp)
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    if (k > 1) then
                       ! perform a rank-1 update of a(1:k-1,1:k-1) and
                       ! store u(k) in column k
                       if (abs(real(a(k,k),KIND=qp)) >= sfmin) then
                          ! perform a rank-1 update of a(1:k-1,1:k-1) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*1/d(k)*w(k)**t
                          d11 = one/real(a(k,k),KIND=qp)
                          call la_wher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                          ! store u(k) in column k
                          call la_wqscal(k - 1,d11,a(1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=qp)
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - u(k)*d(k)*u(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_wher(uplo,k - 1,-d11,a(1,k),1,a,lda)
                       end if
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( ( a(k-1)a(k) )*inv(d(k)) ) * ( a(k-1)a(k) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k > 2) then
                       ! d = |a12|
                       d = la_qlapy2(real(a(k - 1,k),KIND=qp),aimag(a(k - 1,k)))

                       d11 = real(a(k,k)/d,KIND=qp)
                       d22 = real(a(k - 1,k - 1)/d,KIND=qp)
                       d12 = a(k - 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k - 2,1,-1
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wkm1 = tt*(d11*a(j,k - 1) - conjg(d12)*a(j,k))
                          wk = tt*(d22*a(j,k) - d12*a(j,k - 1))
                          ! perform a rank-2 update of a(1:k-2,1:k-2)
                          do i = j,1,-1
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k - 1) &
                                       /d)*conjg(wkm1)
                          end do
                          ! store u(k) and u(k-1) in cols k and k-1 for row j
                          a(j,k) = wk/d
                          a(j,k - 1) = wkm1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=qp),zero,KIND=qp)
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**h using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              p = k
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(a(k,k),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iwamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(a(k,k),KIND=qp)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    42 continue
                       ! begin pivot search loop body
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_iwamax(imax - k,a(imax,k),lda)
                          rowmax = cabs1(a(imax,jmax))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_iwamax(n - imax,a(imax + 1,imax),1)
                          qtemp = cabs1(a(itemp,imax))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(a(imax,imax),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 42
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! for only a 2x2 pivot, interchange rows and columns k and p
                 ! in the trailing submatrix a(k:n,k:n)
                 if ((kstep == 2) .and. (p /= k)) then
                    ! (1) swap columnar parts
                    if (p < n) call la_wswap(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! (2) swap and conjugate middle parts
                    do j = k + 1,p - 1
                       t = conjg(a(j,k))
                       a(j,k) = conjg(a(p,j))
                       a(p,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(p,k) = conjg(a(p,k))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(k,k),KIND=qp)
                    a(k,k) = real(a(p,p),KIND=qp)
                    a(p,p) = r1
                 end if
                 ! for both 1x1 and 2x2 pivots, interchange rows and
                 ! columns kk and kp in the trailing submatrix a(k:n,k:n)
                 if (kp /= kk) then
                    ! (1) swap columnar parts
                    if (kp < n) call la_wswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! (2) swap and conjugate middle parts
                    do j = kk + 1,kp - 1
                       t = conjg(a(j,kk))
                       a(j,kk) = conjg(a(kp,j))
                       a(kp,j) = t
                    end do
                    ! (3) swap and conjugate corner elements at row-col interserction
                    a(kp,kk) = conjg(a(kp,kk))
                    ! (4) swap diagonal elements at row-col intersection
                    r1 = real(a(kk,kk),KIND=qp)
                    a(kk,kk) = real(a(kp,kp),KIND=qp)
                    a(kp,kp) = r1
                    if (kstep == 2) then
                       ! (*) make sure that diagonal element of pivot is real
                       a(k,k) = real(a(k,k),KIND=qp)
                       ! (5) swap row elements
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 else
                    ! (*) make sure that diagonal element of pivot is real
                    a(k,k) = real(a(k,k),KIND=qp)
                    if (kstep == 2) a(k + 1,k + 1) = real(a(k + 1,k + 1),KIND=qp)
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of a now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) and
                       ! store l(k) in column k
                       ! handle division by a small number
                       if (abs(real(a(k,k),KIND=qp)) >= sfmin) then
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                          d11 = one/real(a(k,k),KIND=qp)
                          call la_wher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                          ! store l(k) in column k
                          call la_wqscal(n - k,d11,a(k + 1,k),1)
                       else
                          ! store l(k) in column k
                          d11 = real(a(k,k),KIND=qp)
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/d11
                          end do
                          ! perform a rank-1 update of a(k+1:n,k+1:n) as
                          ! a := a - l(k)*d(k)*l(k)**t
                             ! = a - w(k)*(1/d(k))*w(k)**t
                             ! = a - (w(k)/d(k))*(d(k))*(w(k)/d(k))**t
                          call la_wher(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       end if
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! perform a rank-2 update of a(k+2:n,k+2:n) as
                    ! a := a - ( l(k) l(k+1) ) * d(k) * ( l(k) l(k+1) )**t
                       ! = a - ( ( a(k)a(k+1) )*inv(d(k) ) * ( a(k)a(k+1) )**t
                    ! and store l(k) and l(k+1) in columns k and k+1
                    if (k < n - 1) then
                       ! d = |a21|
                       d = la_qlapy2(real(a(k + 1,k),KIND=qp),aimag(a(k + 1,k)))

                       d11 = real(a(k + 1,k + 1),KIND=qp)/d
                       d22 = real(a(k,k),KIND=qp)/d
                       d21 = a(k + 1,k)/d
                       tt = one/(d11*d22 - one)
                       do j = k + 2,n
                          ! compute  d21 * ( w(k)w(k+1) ) * inv(d(k)) for row j
                          wk = tt*(d11*a(j,k) - d21*a(j,k + 1))
                          wkp1 = tt*(d22*a(j,k + 1) - conjg(d21)*a(j,k))
                          ! perform a rank-2 update of a(k+2:n,k+2:n)
                          do i = j,n
                             a(i,j) = a(i,j) - (a(i,k)/d)*conjg(wk) - (a(i,k + 1) &
                                       /d)*conjg(wkp1)
                          end do
                          ! store l(k) and l(k+1) in cols k and k+1 for row j
                          a(j,k) = wk/d
                          a(j,k + 1) = wkp1/d
                          ! (*) make sure that diagonal element of pivot is real
                          a(j,j) = cmplx(real(a(j,j),KIND=qp),zero,KIND=qp)
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_whetf2_rook

     !> CHETRI_ROOK: computes the inverse of a complex Hermitian indefinite matrix
     !> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
     !> CHETRF_ROOK.

     pure subroutine la_chetri_rook(uplo,n,a,lda,ipiv,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp,kstep
           real(sp) :: ak,akp1,d,t
           complex(sp) :: akkp1,temp
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real
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
              call la_xerbla('CHETRI_ROOK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == czero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == czero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**h.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 70
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/real(a(k,k),KIND=sp)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_ccopy(k - 1,a(1,k),1,work,1)
                    call la_chemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - real(la_cdotc(k - 1,work,1,a(1,k),1), &
                              KIND=sp)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = real(a(k,k),KIND=sp)/t
                 akp1 = real(a(k + 1,k + 1),KIND=sp)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_ccopy(k - 1,a(1,k),1,work,1)
                    call la_chemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - real(la_cdotc(k - 1,work,1,a(1,k),1), &
                              KIND=sp)
                    a(k,k + 1) = a(k,k + 1) - la_cdotc(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_ccopy(k - 1,a(1,k + 1),1,work,1)
                    call la_chemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - real(la_cdotc(k - 1,work,1,a(1,k + 1), &
                              1),KIND=sp)
                 end if
                 kstep = 2
              end if
              if (kstep == 1) then
                 ! interchange rows and columns k and ipiv(k) in the leading
                 ! submatrix a(1:k,1:k)
                 kp = ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_cswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              else
                 ! interchange rows and columns k and k+1 with -ipiv(k) and
                 ! -ipiv(k+1) in the leading submatrix a(k+1:n,k+1:n)
                 ! (1) interchange rows and columns k and -ipiv(k)
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_cswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
                 ! (2) interchange rows and columns k+1 and -ipiv(k+1)
                 k = k + 1
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_cswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              end if
              k = k + 1
              go to 30
              70 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**h.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              80 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 120
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/real(a(k,k),KIND=sp)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_ccopy(n - k,a(k + 1,k),1,work,1)
                    call la_chemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - real(la_cdotc(n - k,work,1,a(k + 1,k),1), &
                              KIND=sp)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = real(a(k - 1,k - 1),KIND=sp)/t
                 akp1 = real(a(k,k),KIND=sp)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_ccopy(n - k,a(k + 1,k),1,work,1)
                    call la_chemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - real(la_cdotc(n - k,work,1,a(k + 1,k),1), &
                              KIND=sp)
                    a(k,k - 1) = a(k,k - 1) - la_cdotc(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_ccopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_chemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - real(la_cdotc(n - k,work,1,a(k + 1,k - 1) &
                              ,1),KIND=sp)
                 end if
                 kstep = 2
              end if
              if (kstep == 1) then
                 ! interchange rows and columns k and ipiv(k) in the trailing
                 ! submatrix a(k:n,k:n)
                 kp = ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_cswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              else
                 ! interchange rows and columns k and k-1 with -ipiv(k) and
                 ! -ipiv(k-1) in the trailing submatrix a(k-1:n,k-1:n)
                 ! (1) interchange rows and columns k and -ipiv(k)
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_cswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
                 ! (2) interchange rows and columns k-1 and -ipiv(k-1)
                 k = k - 1
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_cswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              end if
              k = k - 1
              go to 80
              120 continue
           end if
           return
     end subroutine la_chetri_rook
     !> ZHETRI_ROOK: computes the inverse of a complex Hermitian indefinite matrix
     !> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
     !> ZHETRF_ROOK.

     pure subroutine la_zhetri_rook(uplo,n,a,lda,ipiv,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp,kstep
           real(dp) :: ak,akp1,d,t
           complex(dp) :: akkp1,temp
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real
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
              call la_xerbla('ZHETRI_ROOK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == czero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == czero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**h.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 70
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/real(a(k,k),KIND=dp)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_zcopy(k - 1,a(1,k),1,work,1)
                    call la_zhemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - real(la_zdotc(k - 1,work,1,a(1,k),1), &
                              KIND=dp)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = real(a(k,k),KIND=dp)/t
                 akp1 = real(a(k + 1,k + 1),KIND=dp)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_zcopy(k - 1,a(1,k),1,work,1)
                    call la_zhemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - real(la_zdotc(k - 1,work,1,a(1,k),1), &
                              KIND=dp)
                    a(k,k + 1) = a(k,k + 1) - la_zdotc(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_zcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_zhemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - real(la_zdotc(k - 1,work,1,a(1,k + 1), &
                              1),KIND=dp)
                 end if
                 kstep = 2
              end if
              if (kstep == 1) then
                 ! interchange rows and columns k and ipiv(k) in the leading
                 ! submatrix a(1:k,1:k)
                 kp = ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_zswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              else
                 ! interchange rows and columns k and k+1 with -ipiv(k) and
                 ! -ipiv(k+1) in the leading submatrix a(k+1:n,k+1:n)
                 ! (1) interchange rows and columns k and -ipiv(k)
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_zswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
                 ! (2) interchange rows and columns k+1 and -ipiv(k+1)
                 k = k + 1
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_zswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              end if
              k = k + 1
              go to 30
              70 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**h.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              80 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 120
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/real(a(k,k),KIND=dp)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_zcopy(n - k,a(k + 1,k),1,work,1)
                    call la_zhemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - real(la_zdotc(n - k,work,1,a(k + 1,k),1), &
                              KIND=dp)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = real(a(k - 1,k - 1),KIND=dp)/t
                 akp1 = real(a(k,k),KIND=dp)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_zcopy(n - k,a(k + 1,k),1,work,1)
                    call la_zhemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - real(la_zdotc(n - k,work,1,a(k + 1,k),1), &
                              KIND=dp)
                    a(k,k - 1) = a(k,k - 1) - la_zdotc(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_zcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_zhemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - real(la_zdotc(n - k,work,1,a(k + 1,k - 1) &
                              ,1),KIND=dp)
                 end if
                 kstep = 2
              end if
              if (kstep == 1) then
                 ! interchange rows and columns k and ipiv(k) in the trailing
                 ! submatrix a(k:n,k:n)
                 kp = ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_zswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              else
                 ! interchange rows and columns k and k-1 with -ipiv(k) and
                 ! -ipiv(k-1) in the trailing submatrix a(k-1:n,k-1:n)
                 ! (1) interchange rows and columns k and -ipiv(k)
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_zswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
                 ! (2) interchange rows and columns k-1 and -ipiv(k-1)
                 k = k - 1
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_zswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              end if
              k = k - 1
              go to 80
              120 continue
           end if
           return
     end subroutine la_zhetri_rook
     !> WHETRI_ROOK: computes the inverse of a complex Hermitian indefinite matrix
     !> A using the factorization A = U*D*U**H or A = L*D*L**H computed by
     !> WHETRF_ROOK.

     pure subroutine la_whetri_rook(uplo,n,a,lda,ipiv,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp,kstep
           real(qp) :: ak,akp1,d,t
           complex(qp) :: akkp1,temp
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real
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
              call la_xerbla('WHETRI_ROOK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == czero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == czero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**h.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 70
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/real(a(k,k),KIND=qp)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_wcopy(k - 1,a(1,k),1,work,1)
                    call la_whemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - real(la_wdotc(k - 1,work,1,a(1,k),1), &
                              KIND=qp)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = real(a(k,k),KIND=qp)/t
                 akp1 = real(a(k + 1,k + 1),KIND=qp)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_wcopy(k - 1,a(1,k),1,work,1)
                    call la_whemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - real(la_wdotc(k - 1,work,1,a(1,k),1), &
                              KIND=qp)
                    a(k,k + 1) = a(k,k + 1) - la_wdotc(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_wcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_whemv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - real(la_wdotc(k - 1,work,1,a(1,k + 1), &
                              1),KIND=qp)
                 end if
                 kstep = 2
              end if
              if (kstep == 1) then
                 ! interchange rows and columns k and ipiv(k) in the leading
                 ! submatrix a(1:k,1:k)
                 kp = ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_wswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              else
                 ! interchange rows and columns k and k+1 with -ipiv(k) and
                 ! -ipiv(k+1) in the leading submatrix a(k+1:n,k+1:n)
                 ! (1) interchange rows and columns k and -ipiv(k)
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_wswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
                 ! (2) interchange rows and columns k+1 and -ipiv(k+1)
                 k = k + 1
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp > 1) call la_wswap(kp - 1,a(1,k),1,a(1,kp),1)
                    do j = kp + 1,k - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              end if
              k = k + 1
              go to 30
              70 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**h.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              80 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 120
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/real(a(k,k),KIND=qp)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_wcopy(n - k,a(k + 1,k),1,work,1)
                    call la_whemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - real(la_wdotc(n - k,work,1,a(k + 1,k),1), &
                              KIND=qp)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = real(a(k - 1,k - 1),KIND=qp)/t
                 akp1 = real(a(k,k),KIND=qp)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_wcopy(n - k,a(k + 1,k),1,work,1)
                    call la_whemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - real(la_wdotc(n - k,work,1,a(k + 1,k),1), &
                              KIND=qp)
                    a(k,k - 1) = a(k,k - 1) - la_wdotc(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_wcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_whemv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - real(la_wdotc(n - k,work,1,a(k + 1,k - 1) &
                              ,1),KIND=qp)
                 end if
                 kstep = 2
              end if
              if (kstep == 1) then
                 ! interchange rows and columns k and ipiv(k) in the trailing
                 ! submatrix a(k:n,k:n)
                 kp = ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_wswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              else
                 ! interchange rows and columns k and k-1 with -ipiv(k) and
                 ! -ipiv(k-1) in the trailing submatrix a(k-1:n,k-1:n)
                 ! (1) interchange rows and columns k and -ipiv(k)
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_wswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
                 ! (2) interchange rows and columns k-1 and -ipiv(k-1)
                 k = k - 1
                 kp = -ipiv(k)
                 if (kp /= k) then
                    if (kp < n) call la_wswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                    do j = k + 1,kp - 1
                       temp = conjg(a(j,k))
                       a(j,k) = conjg(a(kp,j))
                       a(kp,j) = temp
                    end do
                    a(kp,k) = conjg(a(kp,k))
                    temp = a(k,k)
                    a(k,k) = a(kp,kp)
                    a(kp,kp) = temp
                 end if
              end if
              k = k - 1
              go to 80
              120 continue
           end if
           return
     end subroutine la_whetri_rook

     !> CLAHEF_RK: computes a partial factorization of a complex Hermitian
     !> matrix A using the bounded Bunch-Kaufman (rook) diagonal
     !> pivoting method. The partial factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L',
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> CLAHEF_RK is an auxiliary routine called by CHETRF_RK. It uses
     !> blocked code (calling Level 3 BLAS) to update the submatrix
     !> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').

     pure subroutine la_clahef_rk(uplo,n,nb,kb,a,lda,e,ipiv,w,ldw,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: w(ldw,*),e(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: imax,itemp,ii,j,jb,jj,jmax,k,kk,kkw,kp,kstep,kw,p
           real(sp) :: absakk,alpha,colmax,stemp,r1,rowmax,t,sfmin
           complex(sp) :: d11,d21,d22,z
           ! Intrinsic Functions
           intrinsic :: abs,conjg,aimag,max,min,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=sp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_slamch('S')
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11 (note that conjg(w) is actually stored)
              ! initialize the first entry of array e, where superdiagonal
              ! elements of d are stored
              e(1) = czero
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              k = n
              10 continue
              ! kw is the column of w which corresponds to column k of a
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              kstep = 1
              p = k
              ! copy column k of a to column kw of w and update it
              if (k > 1) call la_ccopy(k - 1,a(1,k),1,w(1,kw),1)
              w(k,kw) = real(a(k,k),KIND=sp)
              if (k < n) then
                 call la_cgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k,kw + 1), &
                           ldw,cone,w(1,kw),1)
                 w(k,kw) = real(w(k,kw),KIND=sp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,kw),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_icamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,kw),KIND=sp)
                 if (k > 1) call la_ccopy(k - 1,w(1,kw),1,a(1,k),1)
                 ! set e( k ) to zero
                 if (k > 1) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! lop until pivot found
                    done = .false.
                    12 continue
                       ! begin pivot search loop body
                       ! copy column imax to column kw-1 of w and update it
                       if (imax > 1) call la_ccopy(imax - 1,a(1,imax),1,w(1,kw - 1),1)

                       w(imax,kw - 1) = real(a(imax,imax),KIND=sp)
                       call la_ccopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                       call la_clacgv(k - imax,w(imax + 1,kw - 1),1)
                       if (k < n) then
                          call la_cgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                                    imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                          w(imax,kw - 1) = real(w(imax,kw - 1),KIND=sp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_icamax(k - imax,w(imax + 1,kw - 1),1)
                          rowmax = cabs1(w(jmax,kw - 1))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_icamax(imax - 1,w(1,kw - 1),1)
                          stemp = cabs1(w(itemp,kw - 1))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,kw - 1),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column kw-1 of w to column kw of w
                          call la_ccopy(k,w(1,kw - 1),1,w(1,kw),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_ccopy(k,w(1,kw - 1),1,w(1,kw),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns p and k.
                 ! updated column p is already stored in column kw of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column k to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k-1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=sp)
                    call la_ccopy(k - 1 - p,a(p + 1,k),1,a(p,p + 1),lda)
                    call la_clacgv(k - 1 - p,a(p,p + 1),lda)
                    if (p > 1) call la_ccopy(p - 1,a(1,k),1,a(1,p),1)
                    ! interchange rows k and p in the last k+1 to n columns of a
                    ! (columns k and k-1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_cswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                    call la_cswap(n - kk + 1,w(k,kkw),ldw,w(p,kkw),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=sp)
                    call la_ccopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    call la_clacgv(kk - 1 - kp,a(kp,kp + 1),lda)
                    if (kp > 1) call la_ccopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_cswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_cswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! (1) store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element u(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=sp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_ccopy(k,w(1,kw),1,a(1,k),1)
                    if (k > 1) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=sp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_csscal(k - 1,r1,a(1,k),1)
                       else
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(kw)
                       call la_clacgv(k - 1,w(1,kw),1)
                       ! store the superdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! (1) store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! (note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored)
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/conjg(d21)
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(real(d11*d22,KIND=sp) - one)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = t*((d11*w(j,kw - 1) - w(j,kw))/d21)
                          a(j,k) = t*((d22*w(j,kw) - w(j,kw - 1))/conjg(d21))
                       end do
                    end if
                    ! copy diagonal elements of d(k) to a,
                    ! copy superdiagonal element of d(k) to e(k) and
                    ! zero out superdiagonal entry of a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = czero
                    a(k,k) = w(k,kw)
                    e(k) = w(k - 1,kw)
                    e(k - 1) = czero
                    ! (2) conjugate columns w(kw) and w(kw-1)
                    call la_clacgv(k - 1,w(1,kw),1)
                    call la_clacgv(k - 2,w(1,kw - 1),1)
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**h = a11 - u12*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                    call la_cgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                 end do
                 ! update the rectangular superdiagonal block
                 if (j >= 2) call la_cgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a( &
                           1,k + 1),lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22 (note that conjg(w) is actually stored)
              ! initialize the unused last entry of the subdiagonal array e.
              e(n) = czero
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              kstep = 1
              p = k
              ! copy column k of a to column k of w and update column k of w
              w(k,k) = real(a(k,k),KIND=sp)
              if (k < n) call la_ccopy(n - k,a(k + 1,k),1,w(k + 1,k),1)
              if (k > 1) then
                 call la_cgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1), &
                           ldw,cone,w(k,k),1)
                 w(k,k) = real(w(k,k),KIND=sp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,k),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_icamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,k),KIND=sp)
                 if (k < n) call la_ccopy(n - k,w(k + 1,k),1,a(k + 1,k),1)
                 ! set e( k ) to zero
                 if (k < n) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    72 continue
                       ! begin pivot search loop body
                       ! copy column imax to column k+1 of w and update it
                       call la_ccopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                       call la_clacgv(imax - k,w(k,k + 1),1)
                       w(imax,k + 1) = real(a(imax,imax),KIND=sp)
                       if (imax < n) call la_ccopy(n - imax,a(imax + 1,imax),1,w(imax + 1,k + 1 &
                                 ),1)
                       if (k > 1) then
                          call la_cgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w( &
                                    imax,1),ldw,cone,w(k,k + 1),1)
                          w(imax,k + 1) = real(w(imax,k + 1),KIND=sp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_icamax(imax - k,w(k,k + 1),1)
                          rowmax = cabs1(w(jmax,k + 1))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_icamax(n - imax,w(imax + 1,k + 1),1)
                          stemp = cabs1(w(itemp,k + 1))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,k+1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,k + 1),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column k+1 of w to column k of w
                          call la_ccopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_ccopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 72
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns p and k (only for 2-by-2 pivot).
                 ! updated column p is already stored in column k of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column kk-1 to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k+1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=sp)
                    call la_ccopy(p - k - 1,a(k + 1,k),1,a(p,k + 1),lda)
                    call la_clacgv(p - k - 1,a(p,k + 1),lda)
                    if (p < n) call la_ccopy(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! interchange rows k and p in first k-1 columns of a
                    ! (columns k and k+1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in first kk columns of w.
                    if (k > 1) call la_cswap(k - 1,a(k,1),lda,a(p,1),lda)
                    call la_cswap(kk,w(k,1),ldw,w(p,1),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=sp)
                    call la_ccopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    call la_clacgv(kp - kk - 1,a(kp,kk + 1),lda)
                    if (kp < n) call la_ccopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (column k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_cswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_cswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! (1) store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=sp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_ccopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=sp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_csscal(n - k,r1,a(k + 1,k),1)
                       else
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(k)
                       call la_clacgv(n - k,w(k + 1,k),1)
                       ! store the subdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! (1) store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored.
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/conjg(d21)
                       t = one/(real(d11*d22,KIND=sp) - one)
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = t*((d11*w(j,k) - w(j,k + 1))/conjg(d21))
                          a(j,k + 1) = t*((d22*w(j,k + 1) - w(j,k))/d21)
                       end do
                    end if
                    ! copy diagonal elements of d(k) to a,
                    ! copy subdiagonal element of d(k) to e(k) and
                    ! zero out subdiagonal entry of a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = czero
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                    e(k) = w(k + 1,k)
                    e(k + 1) = czero
                    ! (2) conjugate columns w(k) and w(k+1)
                    call la_clacgv(n - k,w(k + 1,k),1)
                    call la_clacgv(n - k - 1,w(k + 2,k + 1),1)
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**h = a22 - l21*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                    call la_cgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_cgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_clahef_rk
     !> ZLAHEF_RK: computes a partial factorization of a complex Hermitian
     !> matrix A using the bounded Bunch-Kaufman (rook) diagonal
     !> pivoting method. The partial factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L',
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> ZLAHEF_RK is an auxiliary routine called by ZHETRF_RK. It uses
     !> blocked code (calling Level 3 BLAS) to update the submatrix
     !> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').

     pure subroutine la_zlahef_rk(uplo,n,nb,kb,a,lda,e,ipiv,w,ldw,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: w(ldw,*),e(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: imax,itemp,ii,j,jb,jj,jmax,k,kk,kkw,kp,kstep,kw,p
           real(dp) :: absakk,alpha,colmax,dtemp,r1,rowmax,t,sfmin
           complex(dp) :: d11,d21,d22,z
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,aimag,max,min,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=dp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_dlamch('S')
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11 (note that conjg(w) is actually stored)
              ! initialize the first entry of array e, where superdiagonal
              ! elements of d are stored
              e(1) = czero
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              k = n
              10 continue
              ! kw is the column of w which corresponds to column k of a
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              kstep = 1
              p = k
              ! copy column k of a to column kw of w and update it
              if (k > 1) call la_zcopy(k - 1,a(1,k),1,w(1,kw),1)
              w(k,kw) = real(a(k,k),KIND=dp)
              if (k < n) then
                 call la_zgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k,kw + 1), &
                           ldw,cone,w(1,kw),1)
                 w(k,kw) = real(w(k,kw),KIND=dp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,kw),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_izamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,kw),KIND=dp)
                 if (k > 1) call la_zcopy(k - 1,w(1,kw),1,a(1,k),1)
                 ! set e( k ) to zero
                 if (k > 1) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! lop until pivot found
                    done = .false.
                    12 continue
                       ! begin pivot search loop body
                       ! copy column imax to column kw-1 of w and update it
                       if (imax > 1) call la_zcopy(imax - 1,a(1,imax),1,w(1,kw - 1),1)

                       w(imax,kw - 1) = real(a(imax,imax),KIND=dp)
                       call la_zcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                       call la_zlacgv(k - imax,w(imax + 1,kw - 1),1)
                       if (k < n) then
                          call la_zgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                                    imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                          w(imax,kw - 1) = real(w(imax,kw - 1),KIND=dp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_izamax(k - imax,w(imax + 1,kw - 1),1)
                          rowmax = cabs1(w(jmax,kw - 1))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_izamax(imax - 1,w(1,kw - 1),1)
                          dtemp = cabs1(w(itemp,kw - 1))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,kw - 1),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column kw-1 of w to column kw of w
                          call la_zcopy(k,w(1,kw - 1),1,w(1,kw),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_zcopy(k,w(1,kw - 1),1,w(1,kw),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns p and k.
                 ! updated column p is already stored in column kw of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column k to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k-1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=dp)
                    call la_zcopy(k - 1 - p,a(p + 1,k),1,a(p,p + 1),lda)
                    call la_zlacgv(k - 1 - p,a(p,p + 1),lda)
                    if (p > 1) call la_zcopy(p - 1,a(1,k),1,a(1,p),1)
                    ! interchange rows k and p in the last k+1 to n columns of a
                    ! (columns k and k-1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_zswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                    call la_zswap(n - kk + 1,w(k,kkw),ldw,w(p,kkw),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=dp)
                    call la_zcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    call la_zlacgv(kk - 1 - kp,a(kp,kp + 1),lda)
                    if (kp > 1) call la_zcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_zswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_zswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! (1) store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element u(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=dp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_zcopy(k,w(1,kw),1,a(1,k),1)
                    if (k > 1) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=dp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_zdscal(k - 1,r1,a(1,k),1)
                       else
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(kw)
                       call la_zlacgv(k - 1,w(1,kw),1)
                       ! store the superdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! (1) store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! (note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored)
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/conjg(d21)
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(real(d11*d22,KIND=dp) - one)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = t*((d11*w(j,kw - 1) - w(j,kw))/d21)
                          a(j,k) = t*((d22*w(j,kw) - w(j,kw - 1))/conjg(d21))
                       end do
                    end if
                    ! copy diagonal elements of d(k) to a,
                    ! copy superdiagonal element of d(k) to e(k) and
                    ! zero out superdiagonal entry of a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = czero
                    a(k,k) = w(k,kw)
                    e(k) = w(k - 1,kw)
                    e(k - 1) = czero
                    ! (2) conjugate columns w(kw) and w(kw-1)
                    call la_zlacgv(k - 1,w(1,kw),1)
                    call la_zlacgv(k - 2,w(1,kw - 1),1)
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**h = a11 - u12*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                    call la_zgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                 end do
                 ! update the rectangular superdiagonal block
                 if (j >= 2) call la_zgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a( &
                           1,k + 1),lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22 (note that conjg(w) is actually stored)
              ! initialize the unused last entry of the subdiagonal array e.
              e(n) = czero
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              kstep = 1
              p = k
              ! copy column k of a to column k of w and update column k of w
              w(k,k) = real(a(k,k),KIND=dp)
              if (k < n) call la_zcopy(n - k,a(k + 1,k),1,w(k + 1,k),1)
              if (k > 1) then
                 call la_zgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1), &
                           ldw,cone,w(k,k),1)
                 w(k,k) = real(w(k,k),KIND=dp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,k),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_izamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,k),KIND=dp)
                 if (k < n) call la_zcopy(n - k,w(k + 1,k),1,a(k + 1,k),1)
                 ! set e( k ) to zero
                 if (k < n) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    72 continue
                       ! begin pivot search loop body
                       ! copy column imax to column k+1 of w and update it
                       call la_zcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                       call la_zlacgv(imax - k,w(k,k + 1),1)
                       w(imax,k + 1) = real(a(imax,imax),KIND=dp)
                       if (imax < n) call la_zcopy(n - imax,a(imax + 1,imax),1,w(imax + 1,k + 1 &
                                 ),1)
                       if (k > 1) then
                          call la_zgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w( &
                                    imax,1),ldw,cone,w(k,k + 1),1)
                          w(imax,k + 1) = real(w(imax,k + 1),KIND=dp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_izamax(imax - k,w(k,k + 1),1)
                          rowmax = cabs1(w(jmax,k + 1))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_izamax(n - imax,w(imax + 1,k + 1),1)
                          dtemp = cabs1(w(itemp,k + 1))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,k+1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,k + 1),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column k+1 of w to column k of w
                          call la_zcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_zcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 72
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns p and k (only for 2-by-2 pivot).
                 ! updated column p is already stored in column k of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column kk-1 to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k+1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=dp)
                    call la_zcopy(p - k - 1,a(k + 1,k),1,a(p,k + 1),lda)
                    call la_zlacgv(p - k - 1,a(p,k + 1),lda)
                    if (p < n) call la_zcopy(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! interchange rows k and p in first k-1 columns of a
                    ! (columns k and k+1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in first kk columns of w.
                    if (k > 1) call la_zswap(k - 1,a(k,1),lda,a(p,1),lda)
                    call la_zswap(kk,w(k,1),ldw,w(p,1),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=dp)
                    call la_zcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    call la_zlacgv(kp - kk - 1,a(kp,kk + 1),lda)
                    if (kp < n) call la_zcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (column k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_zswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_zswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! (1) store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=dp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_zcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=dp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_zdscal(n - k,r1,a(k + 1,k),1)
                       else
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(k)
                       call la_zlacgv(n - k,w(k + 1,k),1)
                       ! store the subdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! (1) store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored.
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/conjg(d21)
                       t = one/(real(d11*d22,KIND=dp) - one)
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = t*((d11*w(j,k) - w(j,k + 1))/conjg(d21))
                          a(j,k + 1) = t*((d22*w(j,k + 1) - w(j,k))/d21)
                       end do
                    end if
                    ! copy diagonal elements of d(k) to a,
                    ! copy subdiagonal element of d(k) to e(k) and
                    ! zero out subdiagonal entry of a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = czero
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                    e(k) = w(k + 1,k)
                    e(k + 1) = czero
                    ! (2) conjugate columns w(k) and w(k+1)
                    call la_zlacgv(n - k,w(k + 1,k),1)
                    call la_zlacgv(n - k - 1,w(k + 2,k + 1),1)
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**h = a22 - l21*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                    call la_zgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_zgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_zlahef_rk
     !> WLAHEF_RK: computes a partial factorization of a complex Hermitian
     !> matrix A using the bounded Bunch-Kaufman (rook) diagonal
     !> pivoting method. The partial factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L',
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> WLAHEF_RK is an auxiliary routine called by WHETRF_RK. It uses
     !> blocked code (calling Level 3 BLAS) to update the submatrix
     !> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').

     pure subroutine la_wlahef_rk(uplo,n,nb,kb,a,lda,e,ipiv,w,ldw,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: w(ldw,*),e(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: imax,itemp,ii,j,jb,jj,jmax,k,kk,kkw,kp,kstep,kw,p
           real(qp) :: absakk,alpha,colmax,qtemp,r1,rowmax,t,sfmin
           complex(qp) :: d11,d21,d22,z
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,aimag,max,min,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=qp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_qlamch('S')
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11 (note that conjg(w) is actually stored)
              ! initialize the first entry of array e, where superdiagonal
              ! elements of d are stored
              e(1) = czero
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              k = n
              10 continue
              ! kw is the column of w which corresponds to column k of a
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              kstep = 1
              p = k
              ! copy column k of a to column kw of w and update it
              if (k > 1) call la_wcopy(k - 1,a(1,k),1,w(1,kw),1)
              w(k,kw) = real(a(k,k),KIND=qp)
              if (k < n) then
                 call la_wgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k,kw + 1), &
                           ldw,cone,w(1,kw),1)
                 w(k,kw) = real(w(k,kw),KIND=qp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,kw),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iwamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,kw),KIND=qp)
                 if (k > 1) call la_wcopy(k - 1,w(1,kw),1,a(1,k),1)
                 ! set e( k ) to zero
                 if (k > 1) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! lop until pivot found
                    done = .false.
                    12 continue
                       ! begin pivot search loop body
                       ! copy column imax to column kw-1 of w and update it
                       if (imax > 1) call la_wcopy(imax - 1,a(1,imax),1,w(1,kw - 1),1)

                       w(imax,kw - 1) = real(a(imax,imax),KIND=qp)
                       call la_wcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                       call la_wlacgv(k - imax,w(imax + 1,kw - 1),1)
                       if (k < n) then
                          call la_wgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                                    imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                          w(imax,kw - 1) = real(w(imax,kw - 1),KIND=qp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_iwamax(k - imax,w(imax + 1,kw - 1),1)
                          rowmax = cabs1(w(jmax,kw - 1))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_iwamax(imax - 1,w(1,kw - 1),1)
                          qtemp = cabs1(w(itemp,kw - 1))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,kw - 1),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column kw-1 of w to column kw of w
                          call la_wcopy(k,w(1,kw - 1),1,w(1,kw),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_wcopy(k,w(1,kw - 1),1,w(1,kw),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns p and k.
                 ! updated column p is already stored in column kw of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column k to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k-1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=qp)
                    call la_wcopy(k - 1 - p,a(p + 1,k),1,a(p,p + 1),lda)
                    call la_wlacgv(k - 1 - p,a(p,p + 1),lda)
                    if (p > 1) call la_wcopy(p - 1,a(1,k),1,a(1,p),1)
                    ! interchange rows k and p in the last k+1 to n columns of a
                    ! (columns k and k-1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_wswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                    call la_wswap(n - kk + 1,w(k,kkw),ldw,w(p,kkw),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=qp)
                    call la_wcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    call la_wlacgv(kk - 1 - kp,a(kp,kp + 1),lda)
                    if (kp > 1) call la_wcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_wswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_wswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! (1) store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element u(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=qp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_wcopy(k,w(1,kw),1,a(1,k),1)
                    if (k > 1) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=qp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_wqscal(k - 1,r1,a(1,k),1)
                       else
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(kw)
                       call la_wlacgv(k - 1,w(1,kw),1)
                       ! store the superdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! (1) store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! (note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored)
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/conjg(d21)
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(real(d11*d22,KIND=qp) - one)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = t*((d11*w(j,kw - 1) - w(j,kw))/d21)
                          a(j,k) = t*((d22*w(j,kw) - w(j,kw - 1))/conjg(d21))
                       end do
                    end if
                    ! copy diagonal elements of d(k) to a,
                    ! copy superdiagonal element of d(k) to e(k) and
                    ! zero out superdiagonal entry of a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = czero
                    a(k,k) = w(k,kw)
                    e(k) = w(k - 1,kw)
                    e(k - 1) = czero
                    ! (2) conjugate columns w(kw) and w(kw-1)
                    call la_wlacgv(k - 1,w(1,kw),1)
                    call la_wlacgv(k - 2,w(1,kw - 1),1)
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**h = a11 - u12*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                    call la_wgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                 end do
                 ! update the rectangular superdiagonal block
                 if (j >= 2) call la_wgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a( &
                           1,k + 1),lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22 (note that conjg(w) is actually stored)
              ! initialize the unused last entry of the subdiagonal array e.
              e(n) = czero
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              kstep = 1
              p = k
              ! copy column k of a to column k of w and update column k of w
              w(k,k) = real(a(k,k),KIND=qp)
              if (k < n) call la_wcopy(n - k,a(k + 1,k),1,w(k + 1,k),1)
              if (k > 1) then
                 call la_wgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1), &
                           ldw,cone,w(k,k),1)
                 w(k,k) = real(w(k,k),KIND=qp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,k),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iwamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,k),KIND=qp)
                 if (k < n) call la_wcopy(n - k,w(k + 1,k),1,a(k + 1,k),1)
                 ! set e( k ) to zero
                 if (k < n) e(k) = czero
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    72 continue
                       ! begin pivot search loop body
                       ! copy column imax to column k+1 of w and update it
                       call la_wcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                       call la_wlacgv(imax - k,w(k,k + 1),1)
                       w(imax,k + 1) = real(a(imax,imax),KIND=qp)
                       if (imax < n) call la_wcopy(n - imax,a(imax + 1,imax),1,w(imax + 1,k + 1 &
                                 ),1)
                       if (k > 1) then
                          call la_wgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w( &
                                    imax,1),ldw,cone,w(k,k + 1),1)
                          w(imax,k + 1) = real(w(imax,k + 1),KIND=qp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_iwamax(imax - k,w(k,k + 1),1)
                          rowmax = cabs1(w(jmax,k + 1))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_iwamax(n - imax,w(imax + 1,k + 1),1)
                          qtemp = cabs1(w(itemp,k + 1))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,k+1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,k + 1),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column k+1 of w to column k of w
                          call la_wcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_wcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 72
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns p and k (only for 2-by-2 pivot).
                 ! updated column p is already stored in column k of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column kk-1 to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k+1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=qp)
                    call la_wcopy(p - k - 1,a(k + 1,k),1,a(p,k + 1),lda)
                    call la_wlacgv(p - k - 1,a(p,k + 1),lda)
                    if (p < n) call la_wcopy(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! interchange rows k and p in first k-1 columns of a
                    ! (columns k and k+1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in first kk columns of w.
                    if (k > 1) call la_wswap(k - 1,a(k,1),lda,a(p,1),lda)
                    call la_wswap(kk,w(k,1),ldw,w(p,1),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=qp)
                    call la_wcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    call la_wlacgv(kp - kk - 1,a(kp,kk + 1),lda)
                    if (kp < n) call la_wcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (column k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_wswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_wswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! (1) store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=qp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_wcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=qp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_wqscal(n - k,r1,a(k + 1,k),1)
                       else
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(k)
                       call la_wlacgv(n - k,w(k + 1,k),1)
                       ! store the subdiagonal element of d in array e
                       e(k) = czero
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! (1) store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored.
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/conjg(d21)
                       t = one/(real(d11*d22,KIND=qp) - one)
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = t*((d11*w(j,k) - w(j,k + 1))/conjg(d21))
                          a(j,k + 1) = t*((d22*w(j,k + 1) - w(j,k))/d21)
                       end do
                    end if
                    ! copy diagonal elements of d(k) to a,
                    ! copy subdiagonal element of d(k) to e(k) and
                    ! zero out subdiagonal entry of a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = czero
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                    e(k) = w(k + 1,k)
                    e(k + 1) = czero
                    ! (2) conjugate columns w(k) and w(k+1)
                    call la_wlacgv(n - k,w(k + 1,k),1)
                    call la_wlacgv(n - k - 1,w(k + 2,k + 1),1)
                 end if
                 ! end column k is nonsingular
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**h = a22 - l21*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                    call la_wgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_wgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_wlahef_rk

     !> CLAHEF_ROOK: computes a partial factorization of a complex Hermitian
     !> matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting
     !> method. The partial factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0      I     )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**H denotes the conjugate transpose of U.
     !> CLAHEF_ROOK is an auxiliary routine called by CHETRF_ROOK. It uses
     !> blocked code (calling Level 3 BLAS) to update the submatrix
     !> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').

     pure subroutine la_clahef_rook(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: imax,itemp,ii,j,jb,jj,jmax,jp1,jp2,k,kk,kkw,kp,kstep,kw, &
                     p
           real(sp) :: absakk,alpha,colmax,stemp,r1,rowmax,t,sfmin
           complex(sp) :: d11,d21,d22,z
           ! Intrinsic Functions
           intrinsic :: abs,conjg,aimag,max,min,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=sp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_slamch('S')
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11 (note that conjg(w) is actually stored)
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              k = n
              10 continue
              ! kw is the column of w which corresponds to column k of a
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              kstep = 1
              p = k
              ! copy column k of a to column kw of w and update it
              if (k > 1) call la_ccopy(k - 1,a(1,k),1,w(1,kw),1)
              w(k,kw) = real(a(k,k),KIND=sp)
              if (k < n) then
                 call la_cgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k,kw + 1), &
                           ldw,cone,w(1,kw),1)
                 w(k,kw) = real(w(k,kw),KIND=sp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,kw),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_icamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,kw),KIND=sp)
                 if (k > 1) call la_ccopy(k - 1,w(1,kw),1,a(1,k),1)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! lop until pivot found
                    done = .false.
                    12 continue
                       ! begin pivot search loop body
                       ! copy column imax to column kw-1 of w and update it
                       if (imax > 1) call la_ccopy(imax - 1,a(1,imax),1,w(1,kw - 1),1)

                       w(imax,kw - 1) = real(a(imax,imax),KIND=sp)
                       call la_ccopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                       call la_clacgv(k - imax,w(imax + 1,kw - 1),1)
                       if (k < n) then
                          call la_cgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                                    imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                          w(imax,kw - 1) = real(w(imax,kw - 1),KIND=sp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_icamax(k - imax,w(imax + 1,kw - 1),1)
                          rowmax = cabs1(w(jmax,kw - 1))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_icamax(imax - 1,w(1,kw - 1),1)
                          stemp = cabs1(w(itemp,kw - 1))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,kw - 1),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column kw-1 of w to column kw of w
                          call la_ccopy(k,w(1,kw - 1),1,w(1,kw),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_ccopy(k,w(1,kw - 1),1,w(1,kw),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns p and k.
                 ! updated column p is already stored in column kw of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column k to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k-1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=sp)
                    call la_ccopy(k - 1 - p,a(p + 1,k),1,a(p,p + 1),lda)
                    call la_clacgv(k - 1 - p,a(p,p + 1),lda)
                    if (p > 1) call la_ccopy(p - 1,a(1,k),1,a(1,p),1)
                    ! interchange rows k and p in the last k+1 to n columns of a
                    ! (columns k and k-1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_cswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                    call la_cswap(n - kk + 1,w(k,kkw),ldw,w(p,kkw),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=sp)
                    call la_ccopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    call la_clacgv(kk - 1 - kp,a(kp,kp + 1),lda)
                    if (kp > 1) call la_ccopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_cswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_cswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! (1) store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element u(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=sp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_ccopy(k,w(1,kw),1,a(1,k),1)
                    if (k > 1) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=sp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_csscal(k - 1,r1,a(1,k),1)
                       else
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(kw)
                       call la_clacgv(k - 1,w(1,kw),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! (1) store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! (note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored)
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/conjg(d21)
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(real(d11*d22,KIND=sp) - one)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = t*((d11*w(j,kw - 1) - w(j,kw))/d21)
                          a(j,k) = t*((d22*w(j,kw) - w(j,kw - 1))/conjg(d21))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                    ! (2) conjugate columns w(kw) and w(kw-1)
                    call la_clacgv(k - 1,w(1,kw),1)
                    call la_clacgv(k - 2,w(1,kw - 1),1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**h = a11 - u12*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                    call la_cgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                 end do
                 ! update the rectangular superdiagonal block
                 if (j >= 2) call la_cgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a( &
                           1,k + 1),lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in of rows in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows j and jp2
                 ! (or j and jp2, and j+1 and jp1) at each step j
                 kstep = 1
                 jp1 = 1
                 ! (here, j is a diagonal index)
                 jj = j
                 jp2 = ipiv(j)
                 if (jp2 < 0) then
                    jp2 = -jp2
                    ! (here, j is a diagonal index)
                    j = j + 1
                    jp1 = -ipiv(j)
                    kstep = 2
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp2 /= jj .and. j <= n) call la_cswap(n - j + 1,a(jp2,j),lda,a(jj,j), &
                           lda)
                 jj = jj + 1
                 if (kstep == 2 .and. jp1 /= jj .and. j <= n) call la_cswap(n - j + 1,a(jp1,j), &
                           lda,a(jj,j),lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22 (note that conjg(w) is actually stored)
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              kstep = 1
              p = k
              ! copy column k of a to column k of w and update column k of w
              w(k,k) = real(a(k,k),KIND=sp)
              if (k < n) call la_ccopy(n - k,a(k + 1,k),1,w(k + 1,k),1)
              if (k > 1) then
                 call la_cgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1), &
                           ldw,cone,w(k,k),1)
                 w(k,k) = real(w(k,k),KIND=sp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,k),KIND=sp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_icamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,k),KIND=sp)
                 if (k < n) call la_ccopy(n - k,w(k + 1,k),1,a(k + 1,k),1)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    72 continue
                       ! begin pivot search loop body
                       ! copy column imax to column k+1 of w and update it
                       call la_ccopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                       call la_clacgv(imax - k,w(k,k + 1),1)
                       w(imax,k + 1) = real(a(imax,imax),KIND=sp)
                       if (imax < n) call la_ccopy(n - imax,a(imax + 1,imax),1,w(imax + 1,k + 1 &
                                 ),1)
                       if (k > 1) then
                          call la_cgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w( &
                                    imax,1),ldw,cone,w(k,k + 1),1)
                          w(imax,k + 1) = real(w(imax,k + 1),KIND=sp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_icamax(imax - k,w(k,k + 1),1)
                          rowmax = cabs1(w(jmax,k + 1))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_icamax(n - imax,w(imax + 1,k + 1),1)
                          stemp = cabs1(w(itemp,k + 1))
                          if (stemp > rowmax) then
                             rowmax = stemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,k+1 ),KIND=sp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,k + 1),KIND=sp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column k+1 of w to column k of w
                          call la_ccopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_ccopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 72
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns p and k (only for 2-by-2 pivot).
                 ! updated column p is already stored in column k of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column kk-1 to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k+1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=sp)
                    call la_ccopy(p - k - 1,a(k + 1,k),1,a(p,k + 1),lda)
                    call la_clacgv(p - k - 1,a(p,k + 1),lda)
                    if (p < n) call la_ccopy(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! interchange rows k and p in first k-1 columns of a
                    ! (columns k and k+1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in first kk columns of w.
                    if (k > 1) call la_cswap(k - 1,a(k,1),lda,a(p,1),lda)
                    call la_cswap(kk,w(k,1),ldw,w(p,1),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=sp)
                    call la_ccopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    call la_clacgv(kp - kk - 1,a(kp,kk + 1),lda)
                    if (kp < n) call la_ccopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (column k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_cswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_cswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! (1) store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=sp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_ccopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=sp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_csscal(n - k,r1,a(k + 1,k),1)
                       else
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(k)
                       call la_clacgv(n - k,w(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! (1) store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored.
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/conjg(d21)
                       t = one/(real(d11*d22,KIND=sp) - one)
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = t*((d11*w(j,k) - w(j,k + 1))/conjg(d21))
                          a(j,k + 1) = t*((d22*w(j,k + 1) - w(j,k))/d21)
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                    ! (2) conjugate columns w(k) and w(k+1)
                    call la_clacgv(n - k,w(k + 1,k),1)
                    call la_clacgv(n - k - 1,w(k + 2,k + 1),1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**h = a22 - l21*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                    call la_cgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=sp)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_cgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows j and jp2
                 ! (or j and jp2, and j-1 and jp1) at each step j
                 kstep = 1
                 jp1 = 1
                 ! (here, j is a diagonal index)
                 jj = j
                 jp2 = ipiv(j)
                 if (jp2 < 0) then
                    jp2 = -jp2
                    ! (here, j is a diagonal index)
                    j = j - 1
                    jp1 = -ipiv(j)
                    kstep = 2
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp2 /= jj .and. j >= 1) call la_cswap(j,a(jp2,1),lda,a(jj,1),lda)

                 jj = jj - 1
                 if (kstep == 2 .and. jp1 /= jj .and. j >= 1) call la_cswap(j,a(jp1,1),lda,a( &
                            jj,1),lda)
              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_clahef_rook
     !> ZLAHEF_ROOK: computes a partial factorization of a complex Hermitian
     !> matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting
     !> method. The partial factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0      I     )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**H denotes the conjugate transpose of U.
     !> ZLAHEF_ROOK is an auxiliary routine called by ZHETRF_ROOK. It uses
     !> blocked code (calling Level 3 BLAS) to update the submatrix
     !> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').

     pure subroutine la_zlahef_rook(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: imax,itemp,ii,j,jb,jj,jmax,jp1,jp2,k,kk,kkw,kp,kstep,kw, &
                     p
           real(dp) :: absakk,alpha,colmax,dtemp,r1,rowmax,t,sfmin
           complex(dp) :: d11,d21,d22,z
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,aimag,max,min,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=dp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_dlamch('S')
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11 (note that conjg(w) is actually stored)
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              k = n
              10 continue
              ! kw is the column of w which corresponds to column k of a
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              kstep = 1
              p = k
              ! copy column k of a to column kw of w and update it
              if (k > 1) call la_zcopy(k - 1,a(1,k),1,w(1,kw),1)
              w(k,kw) = real(a(k,k),KIND=dp)
              if (k < n) then
                 call la_zgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k,kw + 1), &
                           ldw,cone,w(1,kw),1)
                 w(k,kw) = real(w(k,kw),KIND=dp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,kw),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_izamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,kw),KIND=dp)
                 if (k > 1) call la_zcopy(k - 1,w(1,kw),1,a(1,k),1)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! lop until pivot found
                    done = .false.
                    12 continue
                       ! begin pivot search loop body
                       ! copy column imax to column kw-1 of w and update it
                       if (imax > 1) call la_zcopy(imax - 1,a(1,imax),1,w(1,kw - 1),1)

                       w(imax,kw - 1) = real(a(imax,imax),KIND=dp)
                       call la_zcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                       call la_zlacgv(k - imax,w(imax + 1,kw - 1),1)
                       if (k < n) then
                          call la_zgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                                    imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                          w(imax,kw - 1) = real(w(imax,kw - 1),KIND=dp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_izamax(k - imax,w(imax + 1,kw - 1),1)
                          rowmax = cabs1(w(jmax,kw - 1))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_izamax(imax - 1,w(1,kw - 1),1)
                          dtemp = cabs1(w(itemp,kw - 1))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,kw - 1),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column kw-1 of w to column kw of w
                          call la_zcopy(k,w(1,kw - 1),1,w(1,kw),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_zcopy(k,w(1,kw - 1),1,w(1,kw),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns p and k.
                 ! updated column p is already stored in column kw of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column k to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k-1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=dp)
                    call la_zcopy(k - 1 - p,a(p + 1,k),1,a(p,p + 1),lda)
                    call la_zlacgv(k - 1 - p,a(p,p + 1),lda)
                    if (p > 1) call la_zcopy(p - 1,a(1,k),1,a(1,p),1)
                    ! interchange rows k and p in the last k+1 to n columns of a
                    ! (columns k and k-1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_zswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                    call la_zswap(n - kk + 1,w(k,kkw),ldw,w(p,kkw),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=dp)
                    call la_zcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    call la_zlacgv(kk - 1 - kp,a(kp,kp + 1),lda)
                    if (kp > 1) call la_zcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_zswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_zswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! (1) store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element u(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=dp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_zcopy(k,w(1,kw),1,a(1,k),1)
                    if (k > 1) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=dp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_zdscal(k - 1,r1,a(1,k),1)
                       else
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(kw)
                       call la_zlacgv(k - 1,w(1,kw),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! (1) store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! (note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored)
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/conjg(d21)
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(real(d11*d22,KIND=dp) - one)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = t*((d11*w(j,kw - 1) - w(j,kw))/d21)
                          a(j,k) = t*((d22*w(j,kw) - w(j,kw - 1))/conjg(d21))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                    ! (2) conjugate columns w(kw) and w(kw-1)
                    call la_zlacgv(k - 1,w(1,kw),1)
                    call la_zlacgv(k - 2,w(1,kw - 1),1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**h = a11 - u12*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                    call la_zgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                 end do
                 ! update the rectangular superdiagonal block
                 if (j >= 2) call la_zgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a( &
                           1,k + 1),lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in of rows in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows j and jp2
                 ! (or j and jp2, and j+1 and jp1) at each step j
                 kstep = 1
                 jp1 = 1
                 ! (here, j is a diagonal index)
                 jj = j
                 jp2 = ipiv(j)
                 if (jp2 < 0) then
                    jp2 = -jp2
                    ! (here, j is a diagonal index)
                    j = j + 1
                    jp1 = -ipiv(j)
                    kstep = 2
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp2 /= jj .and. j <= n) call la_zswap(n - j + 1,a(jp2,j),lda,a(jj,j), &
                           lda)
                 jj = jj + 1
                 if (kstep == 2 .and. jp1 /= jj .and. j <= n) call la_zswap(n - j + 1,a(jp1,j), &
                           lda,a(jj,j),lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22 (note that conjg(w) is actually stored)
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              kstep = 1
              p = k
              ! copy column k of a to column k of w and update column k of w
              w(k,k) = real(a(k,k),KIND=dp)
              if (k < n) call la_zcopy(n - k,a(k + 1,k),1,w(k + 1,k),1)
              if (k > 1) then
                 call la_zgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1), &
                           ldw,cone,w(k,k),1)
                 w(k,k) = real(w(k,k),KIND=dp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,k),KIND=dp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_izamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,k),KIND=dp)
                 if (k < n) call la_zcopy(n - k,w(k + 1,k),1,a(k + 1,k),1)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    72 continue
                       ! begin pivot search loop body
                       ! copy column imax to column k+1 of w and update it
                       call la_zcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                       call la_zlacgv(imax - k,w(k,k + 1),1)
                       w(imax,k + 1) = real(a(imax,imax),KIND=dp)
                       if (imax < n) call la_zcopy(n - imax,a(imax + 1,imax),1,w(imax + 1,k + 1 &
                                 ),1)
                       if (k > 1) then
                          call la_zgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w( &
                                    imax,1),ldw,cone,w(k,k + 1),1)
                          w(imax,k + 1) = real(w(imax,k + 1),KIND=dp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_izamax(imax - k,w(k,k + 1),1)
                          rowmax = cabs1(w(jmax,k + 1))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_izamax(n - imax,w(imax + 1,k + 1),1)
                          dtemp = cabs1(w(itemp,k + 1))
                          if (dtemp > rowmax) then
                             rowmax = dtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,k+1 ),KIND=dp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,k + 1),KIND=dp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column k+1 of w to column k of w
                          call la_zcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_zcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 72
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns p and k (only for 2-by-2 pivot).
                 ! updated column p is already stored in column k of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column kk-1 to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k+1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=dp)
                    call la_zcopy(p - k - 1,a(k + 1,k),1,a(p,k + 1),lda)
                    call la_zlacgv(p - k - 1,a(p,k + 1),lda)
                    if (p < n) call la_zcopy(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! interchange rows k and p in first k-1 columns of a
                    ! (columns k and k+1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in first kk columns of w.
                    if (k > 1) call la_zswap(k - 1,a(k,1),lda,a(p,1),lda)
                    call la_zswap(kk,w(k,1),ldw,w(p,1),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=dp)
                    call la_zcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    call la_zlacgv(kp - kk - 1,a(kp,kk + 1),lda)
                    if (kp < n) call la_zcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (column k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_zswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_zswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! (1) store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=dp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_zcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=dp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_zdscal(n - k,r1,a(k + 1,k),1)
                       else
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(k)
                       call la_zlacgv(n - k,w(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! (1) store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored.
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/conjg(d21)
                       t = one/(real(d11*d22,KIND=dp) - one)
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = t*((d11*w(j,k) - w(j,k + 1))/conjg(d21))
                          a(j,k + 1) = t*((d22*w(j,k + 1) - w(j,k))/d21)
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                    ! (2) conjugate columns w(k) and w(k+1)
                    call la_zlacgv(n - k,w(k + 1,k),1)
                    call la_zlacgv(n - k - 1,w(k + 2,k + 1),1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**h = a22 - l21*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                    call la_zgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=dp)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_zgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows j and jp2
                 ! (or j and jp2, and j-1 and jp1) at each step j
                 kstep = 1
                 jp1 = 1
                 ! (here, j is a diagonal index)
                 jj = j
                 jp2 = ipiv(j)
                 if (jp2 < 0) then
                    jp2 = -jp2
                    ! (here, j is a diagonal index)
                    j = j - 1
                    jp1 = -ipiv(j)
                    kstep = 2
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp2 /= jj .and. j >= 1) call la_zswap(j,a(jp2,1),lda,a(jj,1),lda)

                 jj = jj - 1
                 if (kstep == 2 .and. jp1 /= jj .and. j >= 1) call la_zswap(j,a(jp1,1),lda,a( &
                            jj,1),lda)
              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_zlahef_rook
     !> WLAHEF_ROOK: computes a partial factorization of a complex Hermitian
     !> matrix A using the bounded Bunch-Kaufman ("rook") diagonal pivoting
     !> method. The partial factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I      0     )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**H U22**H )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**H L21**H )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0      I     )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**H denotes the conjugate transpose of U.
     !> WLAHEF_ROOK is an auxiliary routine called by WHETRF_ROOK. It uses
     !> blocked code (calling Level 3 BLAS) to update the submatrix
     !> A11 (if UPLO = 'U') or A22 (if UPLO = 'L').

     pure subroutine la_wlahef_rook(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: done
           integer(ilp) :: imax,itemp,ii,j,jb,jj,jmax,jp1,jp2,k,kk,kkw,kp,kstep,kw, &
                     p
           real(qp) :: absakk,alpha,colmax,qtemp,r1,rowmax,t,sfmin
           complex(qp) :: d11,d21,d22,z
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,aimag,max,min,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=qp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           ! compute machine safe minimum
           sfmin = la_qlamch('S')
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11 (note that conjg(w) is actually stored)
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              k = n
              10 continue
              ! kw is the column of w which corresponds to column k of a
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              kstep = 1
              p = k
              ! copy column k of a to column kw of w and update it
              if (k > 1) call la_wcopy(k - 1,a(1,k),1,w(1,kw),1)
              w(k,kw) = real(a(k,k),KIND=qp)
              if (k < n) then
                 call la_wgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k,kw + 1), &
                           ldw,cone,w(1,kw),1)
                 w(k,kw) = real(w(k,kw),KIND=qp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,kw),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iwamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,kw),KIND=qp)
                 if (k > 1) call la_wcopy(k - 1,w(1,kw),1,a(1,k),1)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! lop until pivot found
                    done = .false.
                    12 continue
                       ! begin pivot search loop body
                       ! copy column imax to column kw-1 of w and update it
                       if (imax > 1) call la_wcopy(imax - 1,a(1,imax),1,w(1,kw - 1),1)

                       w(imax,kw - 1) = real(a(imax,imax),KIND=qp)
                       call la_wcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                       call la_wlacgv(k - imax,w(imax + 1,kw - 1),1)
                       if (k < n) then
                          call la_wgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                                    imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                          w(imax,kw - 1) = real(w(imax,kw - 1),KIND=qp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = imax + la_iwamax(k - imax,w(imax + 1,kw - 1),1)
                          rowmax = cabs1(w(jmax,kw - 1))
                       else
                          rowmax = zero
                       end if
                       if (imax > 1) then
                          itemp = la_iwamax(imax - 1,w(1,kw - 1),1)
                          qtemp = cabs1(w(itemp,kw - 1))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,kw-1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,kw - 1),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column kw-1 of w to column kw of w
                          call la_wcopy(k,w(1,kw - 1),1,w(1,kw),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k-1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_wcopy(k,w(1,kw - 1),1,w(1,kw),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 12
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns p and k.
                 ! updated column p is already stored in column kw of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column k to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k-1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=qp)
                    call la_wcopy(k - 1 - p,a(p + 1,k),1,a(p,p + 1),lda)
                    call la_wlacgv(k - 1 - p,a(p,p + 1),lda)
                    if (p > 1) call la_wcopy(p - 1,a(1,k),1,a(1,p),1)
                    ! interchange rows k and p in the last k+1 to n columns of a
                    ! (columns k and k-1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_wswap(n - k,a(k,k + 1),lda,a(p,k + 1),lda)
                    call la_wswap(n - kk + 1,w(k,kkw),ldw,w(p,kkw),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=qp)
                    call la_wcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    call la_wlacgv(kk - 1 - kp,a(kp,kp + 1),lda)
                    if (kp > 1) call la_wcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_wswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_wswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! (1) store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element u(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=qp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_wcopy(k,w(1,kw),1,a(1,k),1)
                    if (k > 1) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=qp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_wqscal(k - 1,r1,a(1,k),1)
                       else
                          do ii = 1,k - 1
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(kw)
                       call la_wlacgv(k - 1,w(1,kw),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! (1) store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! (note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored)
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/conjg(d21)
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(real(d11*d22,KIND=qp) - one)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = t*((d11*w(j,kw - 1) - w(j,kw))/d21)
                          a(j,k) = t*((d22*w(j,kw) - w(j,kw - 1))/conjg(d21))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                    ! (2) conjugate columns w(kw) and w(kw-1)
                    call la_wlacgv(k - 1,w(1,kw),1)
                    call la_wlacgv(k - 2,w(1,kw - 1),1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**h = a11 - u12*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                    call la_wgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                 end do
                 ! update the rectangular superdiagonal block
                 if (j >= 2) call la_wgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a( &
                           1,k + 1),lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in of rows in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows j and jp2
                 ! (or j and jp2, and j+1 and jp1) at each step j
                 kstep = 1
                 jp1 = 1
                 ! (here, j is a diagonal index)
                 jj = j
                 jp2 = ipiv(j)
                 if (jp2 < 0) then
                    jp2 = -jp2
                    ! (here, j is a diagonal index)
                    j = j + 1
                    jp1 = -ipiv(j)
                    kstep = 2
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp2 /= jj .and. j <= n) call la_wswap(n - j + 1,a(jp2,j),lda,a(jj,j), &
                           lda)
                 jj = jj + 1
                 if (kstep == 2 .and. jp1 /= jj .and. j <= n) call la_wswap(n - j + 1,a(jp1,j), &
                           lda,a(jj,j),lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22 (note that conjg(w) is actually stored)
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              kstep = 1
              p = k
              ! copy column k of a to column k of w and update column k of w
              w(k,k) = real(a(k,k),KIND=qp)
              if (k < n) call la_wcopy(n - k,a(k + 1,k),1,w(k + 1,k),1)
              if (k > 1) then
                 call la_wgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1), &
                           ldw,cone,w(k,k),1)
                 w(k,k) = real(w(k,k),KIND=qp)
              end if
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(real(w(k,k),KIND=qp))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iwamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
                 a(k,k) = real(w(k,k),KIND=qp)
                 if (k < n) call la_wcopy(n - k,w(k + 1,k),1,a(k + 1,k),1)
              else
                 ! ============================================================
                 ! begin pivot search
                 ! case(1)
                 ! equivalent to testing for absakk>=alpha*colmax
                 ! (used to handle nan and inf)
                 if (.not. (absakk < alpha*colmax)) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    done = .false.
                    ! loop until pivot found
                    72 continue
                       ! begin pivot search loop body
                       ! copy column imax to column k+1 of w and update it
                       call la_wcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                       call la_wlacgv(imax - k,w(k,k + 1),1)
                       w(imax,k + 1) = real(a(imax,imax),KIND=qp)
                       if (imax < n) call la_wcopy(n - imax,a(imax + 1,imax),1,w(imax + 1,k + 1 &
                                 ),1)
                       if (k > 1) then
                          call la_wgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w( &
                                    imax,1),ldw,cone,w(k,k + 1),1)
                          w(imax,k + 1) = real(w(imax,k + 1),KIND=qp)
                       end if
                       ! jmax is the column-index of the largest off-diagonal
                       ! element in row imax, and rowmax is its absolute value.
                       ! determine both rowmax and jmax.
                       if (imax /= k) then
                          jmax = k - 1 + la_iwamax(imax - k,w(k,k + 1),1)
                          rowmax = cabs1(w(jmax,k + 1))
                       else
                          rowmax = zero
                       end if
                       if (imax < n) then
                          itemp = imax + la_iwamax(n - imax,w(imax + 1,k + 1),1)
                          qtemp = cabs1(w(itemp,k + 1))
                          if (qtemp > rowmax) then
                             rowmax = qtemp
                             jmax = itemp
                          end if
                       end if
                       ! case(2)
                       ! equivalent to testing for
                       ! abs( real( w( imax,k+1 ),KIND=qp) )>=alpha*rowmax
                       ! (used to handle nan and inf)
                       if (.not. (abs(real(w(imax,k + 1),KIND=qp)) < alpha*rowmax)) &
                                 then
                          ! interchange rows and columns k and imax,
                          ! use 1-by-1 pivot block
                          kp = imax
                          ! copy column k+1 of w to column k of w
                          call la_wcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                          done = .true.
                       ! case(3)
                       ! equivalent to testing for rowmax==colmax,
                       ! (used to handle nan and inf)
                       else if ((p == jmax) .or. (rowmax <= colmax)) then
                          ! interchange rows and columns k+1 and imax,
                          ! use 2-by-2 pivot block
                          kp = imax
                          kstep = 2
                          done = .true.
                       ! case(4)
                       else
                          ! pivot not found: set params and repeat
                          p = imax
                          colmax = rowmax
                          imax = jmax
                          ! copy updated jmaxth (next imaxth) column to kth of w
                          call la_wcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                       end if
                       ! end pivot search loop body
                    if (.not. done) goto 72
                 end if
                 ! end pivot search
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns p and k (only for 2-by-2 pivot).
                 ! updated column p is already stored in column k of w.
                 if ((kstep == 2) .and. (p /= k)) then
                    ! copy non-updated column kk-1 to column p of submatrix a
                    ! at step k. no need to copy element into columns
                    ! k and k+1 of a for 2-by-2 pivot, since these columns
                    ! will be later overwritten.
                    a(p,p) = real(a(k,k),KIND=qp)
                    call la_wcopy(p - k - 1,a(k + 1,k),1,a(p,k + 1),lda)
                    call la_wlacgv(p - k - 1,a(p,k + 1),lda)
                    if (p < n) call la_wcopy(n - p,a(p + 1,k),1,a(p + 1,p),1)
                    ! interchange rows k and p in first k-1 columns of a
                    ! (columns k and k+1 of a for 2-by-2 pivot will be
                    ! later overwritten). interchange rows k and p
                    ! in first kk columns of w.
                    if (k > 1) call la_wswap(k - 1,a(k,1),lda,a(p,1),lda)
                    call la_wswap(kk,w(k,1),ldw,w(p,1),ldw)
                 end if
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = real(a(kk,kk),KIND=qp)
                    call la_wcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    call la_wlacgv(kp - kk - 1,a(kp,kk + 1),lda)
                    if (kp < n) call la_wcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (column k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_wswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_wswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! (1) store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    ! (note: no need to use for hermitian matrix
                    ! a( k, k ) = real( w( k, k),KIND=qp) to separately copy diagonal
                    ! element d(k,k) from w (potentially saves only one load))
                    call la_wcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       ! (note: no need to check if a(k,k) is not zero,
                        ! since that was ensured earlier in pivot search:
                        ! case a(k,k) = 0 falls into 2x2 pivot case(3))
                       ! handle division by a small number
                       t = real(a(k,k),KIND=qp)
                       if (abs(t) >= sfmin) then
                          r1 = one/t
                          call la_wqscal(n - k,r1,a(k + 1,k),1)
                       else
                          do ii = k + 1,n
                             a(ii,k) = a(ii,k)/t
                          end do
                       end if
                       ! (2) conjugate column w(k)
                       call la_wlacgv(n - k,w(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! (1) store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored.
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! factor out the columns of the inverse of 2-by-2 pivot
                       ! block d, so that each column contains 1, to reduce the
                       ! number of flops when we multiply panel
                       ! ( w(kw-1) w(kw) ) by this inverse, i.e. by d**(-1).
                       ! d**(-1) = ( d11 cj(d21) )**(-1) =
                                 ! ( d21    d22 )
                       ! = 1/(d11*d22-|d21|**2) * ( ( d22) (-cj(d21) ) ) =
                                                ! ( (-d21) (     d11 ) )
                       ! = 1/(|d21|**2) * 1/((d11/cj(d21))*(d22/d21)-1) *
                         ! * ( d21*( d22/d21 ) conj(d21)*(           - 1 ) ) =
                           ! (     (      -1 )           ( d11/conj(d21) ) )
                       ! = 1/(|d21|**2) * 1/(d22*d11-1) *
                         ! * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                           ! (     (  -1 )           ( d22 ) )
                       ! = (1/|d21|**2) * t * ( d21*( d11 ) conj(d21)*(  -1 ) ) =
                                            ! (     (  -1 )           ( d22 ) )
                       ! = ( (t/conj(d21))*( d11 ) (t/d21)*(  -1 ) ) =
                         ! (               (  -1 )         ( d22 ) )
                       ! handle division by a small number. (note: order of
                       ! operations is important)
                       ! = ( t*(( d11 )/conj(d21)) t*((  -1 )/d21 ) )
                         ! (   ((  -1 )          )   (( d22 )     ) ),
                       ! where d11 = d22/d21,
                             ! d22 = d11/conj(d21),
                             ! d21 = d21,
                             ! t = 1/(d22*d11-1).
                       ! (note: no need to check for division by zero,
                        ! since that was ensured earlier in pivot search:
                        ! (a) d21 != 0 in 2x2 pivot case(4),
                            ! since |d21| should be larger than |d11| and |d22|;
                        ! (b) (d22*d11 - 1) != 0, since from (a),
                            ! both |d11| < 1, |d22| < 1, hence |d22*d11| << 1.)
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/conjg(d21)
                       t = one/(real(d11*d22,KIND=qp) - one)
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = t*((d11*w(j,k) - w(j,k + 1))/conjg(d21))
                          a(j,k + 1) = t*((d22*w(j,k + 1) - w(j,k))/d21)
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                    ! (2) conjugate columns w(k) and w(k+1)
                    call la_wlacgv(n - k,w(k + 1,k),1)
                    call la_wlacgv(n - k - 1,w(k + 2,k + 1),1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -p
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**h = a22 - l21*w**h
              ! computing blocks of nb columns at a time (note that conjg(w) is
              ! actually stored)
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                    call la_wgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                    a(jj,jj) = real(a(jj,jj),KIND=qp)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_wgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows j and jp2
                 ! (or j and jp2, and j-1 and jp1) at each step j
                 kstep = 1
                 jp1 = 1
                 ! (here, j is a diagonal index)
                 jj = j
                 jp2 = ipiv(j)
                 if (jp2 < 0) then
                    jp2 = -jp2
                    ! (here, j is a diagonal index)
                    j = j - 1
                    jp1 = -ipiv(j)
                    kstep = 2
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp2 /= jj .and. j >= 1) call la_wswap(j,a(jp2,1),lda,a(jj,1),lda)

                 jj = jj - 1
                 if (kstep == 2 .and. jp1 /= jj .and. j >= 1) call la_wswap(j,a(jp1,1),lda,a( &
                            jj,1),lda)
              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_wlahef_rook

     !> CLAQSY: equilibrates a symmetric matrix A using the scaling factors
     !> in the vector S.

     pure subroutine la_claqsy(uplo,n,a,lda,s,scond,amax,equed)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(in) :: amax,scond
           ! Array Arguments
           real(sp),intent(in) :: s(*)
           complex(sp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: thresh = 0.1e+0_sp

           ! Local Scalars
           integer(ilp) :: i,j
           real(sp) :: cj,large,small
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_slamch('SAFE MINIMUM')/la_slamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = 1,j
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,n
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_claqsy
     !> ZLAQSY: equilibrates a symmetric matrix A using the scaling factors
     !> in the vector S.

     pure subroutine la_zlaqsy(uplo,n,a,lda,s,scond,amax,equed)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(in) :: amax,scond
           ! Array Arguments
           real(dp),intent(in) :: s(*)
           complex(dp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: thresh = 0.1e+0_dp

           ! Local Scalars
           integer(ilp) :: i,j
           real(dp) :: cj,large,small
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_dlamch('SAFE MINIMUM')/la_dlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = 1,j
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,n
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_zlaqsy
     !> WLAQSY: equilibrates a symmetric matrix A using the scaling factors
     !> in the vector S.

     pure subroutine la_wlaqsy(uplo,n,a,lda,s,scond,amax,equed)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(out) :: equed
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(in) :: amax,scond
           ! Array Arguments
           real(qp),intent(in) :: s(*)
           complex(qp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: thresh = 0.1e+0_qp

           ! Local Scalars
           integer(ilp) :: i,j
           real(qp) :: cj,large,small
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              equed = 'N'
              return
           end if
           ! initialize large and small.
           small = la_qlamch('SAFE MINIMUM')/la_qlamch('PRECISION')
           large = one/small
           if (scond >= thresh .and. amax >= small .and. amax <= large) then
              ! no equilibration
              equed = 'N'
           else
              ! replace a by diag(s) * a * diag(s).
              if (la_lsame(uplo,'U')) then
                 ! upper triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = 1,j
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              else
                 ! lower triangle of a is stored.
                 do j = 1,n
                    cj = s(j)
                    do i = j,n
                       a(i,j) = cj*s(i)*a(i,j)
                    end do
                 end do
              end if
              equed = 'Y'
           end if
           return
     end subroutine la_wlaqsy

     !> CHETRF_RK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
     !> A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.
     !> For more information see Further Details section.

     pure subroutine la_chetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: i,iinfo,ip,iws,k,kb,ldwork,lwkopt,nb,nbmin
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < 1 .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! determine the block size
              nb = la_ilaenv(1,'CHETRF_RK',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHETRF_RK',-info)
              return
           else if (lquery) then
              return
           end if
           nbmin = 2
           ldwork = n
           if (nb > 1 .and. nb < n) then
              iws = ldwork*nb
              if (lwork < iws) then
                 nb = max(lwork/ldwork,1)
                 nbmin = max(2,la_ilaenv(2,'CHETRF_RK',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_clahef_rk;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 15
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_clahef_rk(uplo,k,nb,kb,a,lda,e,ipiv,work,ldwork,iinfo)

              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_chetf2_rk(uplo,k,a,lda,e,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! no need to adjust ipiv
              ! apply permutations to the leading panel 1:k-1
              ! read ipiv from the last block factored, i.e.
              ! indices  k-kb+1:k and apply row permutations to the
              ! last k+1 colunms k+1:n after that block
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              if (k < n) then
                 do i = k, (k - kb + 1),-1
                    ip = abs(ipiv(i))
                    if (ip /= i) then
                       call la_cswap(n - k,a(i,k + 1),lda,a(ip,k + 1),lda)
                    end if
                 end do
              end if
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
              ! this label is the exit from main loop over k decreasing
              ! from n to 1 in steps of kb
              15 continue
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_clahef_rk;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 35
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_clahef_rk(uplo,n - k + 1,nb,kb,a(k,k),lda,e(k),ipiv(k), &
                           work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_chetf2_rk(uplo,n - k + 1,a(k,k),lda,e(k),ipiv(k),iinfo)

                 kb = n - k + 1
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
              ! adjust ipiv
              do i = k,k + kb - 1
                 if (ipiv(i) > 0) then
                    ipiv(i) = ipiv(i) + k - 1
                 else
                    ipiv(i) = ipiv(i) - k + 1
                 end if
              end do
              ! apply permutations to the leading panel 1:k-1
              ! read ipiv from the last block factored, i.e.
              ! indices  k:k+kb-1 and apply row permutations to the
              ! first k-1 colunms 1:k-1 before that block
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              if (k > 1) then
                 do i = k, (k + kb - 1),1
                    ip = abs(ipiv(i))
                    if (ip /= i) then
                       call la_cswap(k - 1,a(i,1),lda,a(ip,1),lda)
                    end if
                 end do
              end if
              ! increase k and return to the start of the main loop
              k = k + kb
              go to 20
              ! this label is the exit from main loop over k increasing
              ! from 1 to n in steps of kb
              35 continue
           ! end lower
           end if
           work(1) = lwkopt
           return
     end subroutine la_chetrf_rk
     !> ZHETRF_RK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
     !> A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.
     !> For more information see Further Details section.

     pure subroutine la_zhetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: i,iinfo,ip,iws,k,kb,ldwork,lwkopt,nb,nbmin
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < 1 .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! determine the block size
              nb = la_ilaenv(1,'ZHETRF_RK',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHETRF_RK',-info)
              return
           else if (lquery) then
              return
           end if
           nbmin = 2
           ldwork = n
           if (nb > 1 .and. nb < n) then
              iws = ldwork*nb
              if (lwork < iws) then
                 nb = max(lwork/ldwork,1)
                 nbmin = max(2,la_ilaenv(2,'ZHETRF_RK',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_zlahef_rk;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 15
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_zlahef_rk(uplo,k,nb,kb,a,lda,e,ipiv,work,ldwork,iinfo)

              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_zhetf2_rk(uplo,k,a,lda,e,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! no need to adjust ipiv
              ! apply permutations to the leading panel 1:k-1
              ! read ipiv from the last block factored, i.e.
              ! indices  k-kb+1:k and apply row permutations to the
              ! last k+1 colunms k+1:n after that block
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              if (k < n) then
                 do i = k, (k - kb + 1),-1
                    ip = abs(ipiv(i))
                    if (ip /= i) then
                       call la_zswap(n - k,a(i,k + 1),lda,a(ip,k + 1),lda)
                    end if
                 end do
              end if
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
              ! this label is the exit from main loop over k decreasing
              ! from n to 1 in steps of kb
              15 continue
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_zlahef_rk;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 35
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_zlahef_rk(uplo,n - k + 1,nb,kb,a(k,k),lda,e(k),ipiv(k), &
                           work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_zhetf2_rk(uplo,n - k + 1,a(k,k),lda,e(k),ipiv(k),iinfo)

                 kb = n - k + 1
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
              ! adjust ipiv
              do i = k,k + kb - 1
                 if (ipiv(i) > 0) then
                    ipiv(i) = ipiv(i) + k - 1
                 else
                    ipiv(i) = ipiv(i) - k + 1
                 end if
              end do
              ! apply permutations to the leading panel 1:k-1
              ! read ipiv from the last block factored, i.e.
              ! indices  k:k+kb-1 and apply row permutations to the
              ! first k-1 colunms 1:k-1 before that block
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              if (k > 1) then
                 do i = k, (k + kb - 1),1
                    ip = abs(ipiv(i))
                    if (ip /= i) then
                       call la_zswap(k - 1,a(i,1),lda,a(ip,1),lda)
                    end if
                 end do
              end if
              ! increase k and return to the start of the main loop
              k = k + kb
              go to 20
              ! this label is the exit from main loop over k increasing
              ! from 1 to n in steps of kb
              35 continue
           ! end lower
           end if
           work(1) = lwkopt
           return
     end subroutine la_zhetrf_rk
     !> WHETRF_RK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman (rook) diagonal pivoting method:
     !> A = P*U*D*(U**H)*(P**T) or A = P*L*D*(L**H)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**H (or L**H) is the conjugate of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is Hermitian and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.
     !> For more information see Further Details section.

     pure subroutine la_whetrf_rk(uplo,n,a,lda,e,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: e(*),work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: i,iinfo,ip,iws,k,kb,ldwork,lwkopt,nb,nbmin
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < 1 .and. .not. lquery) then
              info = -8
           end if
           if (info == 0) then
              ! determine the block size
              nb = la_ilaenv(1,'WHETRF_RK',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHETRF_RK',-info)
              return
           else if (lquery) then
              return
           end if
           nbmin = 2
           ldwork = n
           if (nb > 1 .and. nb < n) then
              iws = ldwork*nb
              if (lwork < iws) then
                 nb = max(lwork/ldwork,1)
                 nbmin = max(2,la_ilaenv(2,'WHETRF_RK',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_wlahef_rk;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 15
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_wlahef_rk(uplo,k,nb,kb,a,lda,e,ipiv,work,ldwork,iinfo)

              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_whetf2_rk(uplo,k,a,lda,e,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! no need to adjust ipiv
              ! apply permutations to the leading panel 1:k-1
              ! read ipiv from the last block factored, i.e.
              ! indices  k-kb+1:k and apply row permutations to the
              ! last k+1 colunms k+1:n after that block
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              if (k < n) then
                 do i = k, (k - kb + 1),-1
                    ip = abs(ipiv(i))
                    if (ip /= i) then
                       call la_wswap(n - k,a(i,k + 1),lda,a(ip,k + 1),lda)
                    end if
                 end do
              end if
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
              ! this label is the exit from main loop over k decreasing
              ! from n to 1 in steps of kb
              15 continue
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_wlahef_rk;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 35
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_wlahef_rk(uplo,n - k + 1,nb,kb,a(k,k),lda,e(k),ipiv(k), &
                           work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_whetf2_rk(uplo,n - k + 1,a(k,k),lda,e(k),ipiv(k),iinfo)

                 kb = n - k + 1
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
              ! adjust ipiv
              do i = k,k + kb - 1
                 if (ipiv(i) > 0) then
                    ipiv(i) = ipiv(i) + k - 1
                 else
                    ipiv(i) = ipiv(i) - k + 1
                 end if
              end do
              ! apply permutations to the leading panel 1:k-1
              ! read ipiv from the last block factored, i.e.
              ! indices  k:k+kb-1 and apply row permutations to the
              ! first k-1 colunms 1:k-1 before that block
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              if (k > 1) then
                 do i = k, (k + kb - 1),1
                    ip = abs(ipiv(i))
                    if (ip /= i) then
                       call la_wswap(k - 1,a(i,1),lda,a(ip,1),lda)
                    end if
                 end do
              end if
              ! increase k and return to the start of the main loop
              k = k + kb
              go to 20
              ! this label is the exit from main loop over k increasing
              ! from 1 to n in steps of kb
              35 continue
           ! end lower
           end if
           work(1) = lwkopt
           return
     end subroutine la_whetrf_rk

     !> CHETRF_ROOK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
     !> The form of the factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_chetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: iinfo,iws,j,k,kb,ldwork,lwkopt,nb,nbmin
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < 1 .and. .not. lquery) then
              info = -7
           end if
           if (info == 0) then
              ! determine the block size
              nb = la_ilaenv(1,'CHETRF_ROOK',uplo,n,-1,-1,-1)
              lwkopt = max(1,n*nb)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHETRF_ROOK',-info)
              return
           else if (lquery) then
              return
           end if
           nbmin = 2
           ldwork = n
           if (nb > 1 .and. nb < n) then
              iws = ldwork*nb
              if (lwork < iws) then
                 nb = max(lwork/ldwork,1)
                 nbmin = max(2,la_ilaenv(2,'CHETRF_ROOK',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_clahef_rook;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_clahef_rook(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)

              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_chetf2_rook(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! no need to adjust ipiv
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_clahef_rook;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_clahef_rook(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work, &
                           ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_chetf2_rook(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
                 kb = n - k + 1
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
              ! adjust ipiv
              do j = k,k + kb - 1
                 if (ipiv(j) > 0) then
                    ipiv(j) = ipiv(j) + k - 1
                 else
                    ipiv(j) = ipiv(j) - k + 1
                 end if
              end do
              ! increase k and return to the start of the main loop
              k = k + kb
              go to 20
           end if
           40 continue
           work(1) = lwkopt
           return
     end subroutine la_chetrf_rook
     !> ZHETRF_ROOK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
     !> The form of the factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_zhetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: iinfo,iws,j,k,kb,ldwork,lwkopt,nb,nbmin
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < 1 .and. .not. lquery) then
              info = -7
           end if
           if (info == 0) then
              ! determine the block size
              nb = la_ilaenv(1,'ZHETRF_ROOK',uplo,n,-1,-1,-1)
              lwkopt = max(1,n*nb)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHETRF_ROOK',-info)
              return
           else if (lquery) then
              return
           end if
           nbmin = 2
           ldwork = n
           if (nb > 1 .and. nb < n) then
              iws = ldwork*nb
              if (lwork < iws) then
                 nb = max(lwork/ldwork,1)
                 nbmin = max(2,la_ilaenv(2,'ZHETRF_ROOK',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_zlahef_rook;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_zlahef_rook(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)

              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_zhetf2_rook(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! no need to adjust ipiv
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_zlahef_rook;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_zlahef_rook(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work, &
                           ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_zhetf2_rook(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
                 kb = n - k + 1
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
              ! adjust ipiv
              do j = k,k + kb - 1
                 if (ipiv(j) > 0) then
                    ipiv(j) = ipiv(j) + k - 1
                 else
                    ipiv(j) = ipiv(j) - k + 1
                 end if
              end do
              ! increase k and return to the start of the main loop
              k = k + kb
              go to 20
           end if
           40 continue
           work(1) = lwkopt
           return
     end subroutine la_zhetrf_rook
     !> WHETRF_ROOK: computes the factorization of a complex Hermitian matrix A
     !> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
     !> The form of the factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is Hermitian and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_whetrf_rook(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: iinfo,iws,j,k,kb,ldwork,lwkopt,nb,nbmin
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < 1 .and. .not. lquery) then
              info = -7
           end if
           if (info == 0) then
              ! determine the block size
              nb = la_ilaenv(1,'WHETRF_ROOK',uplo,n,-1,-1,-1)
              lwkopt = max(1,n*nb)
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHETRF_ROOK',-info)
              return
           else if (lquery) then
              return
           end if
           nbmin = 2
           ldwork = n
           if (nb > 1 .and. nb < n) then
              iws = ldwork*nb
              if (lwork < iws) then
                 nb = max(lwork/ldwork,1)
                 nbmin = max(2,la_ilaenv(2,'WHETRF_ROOK',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_wlahef_rook;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_wlahef_rook(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)

              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_whetf2_rook(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! no need to adjust ipiv
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_wlahef_rook;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_wlahef_rook(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work, &
                           ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_whetf2_rook(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
                 kb = n - k + 1
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
              ! adjust ipiv
              do j = k,k + kb - 1
                 if (ipiv(j) > 0) then
                    ipiv(j) = ipiv(j) + k - 1
                 else
                    ipiv(j) = ipiv(j) - k + 1
                 end if
              end do
              ! increase k and return to the start of the main loop
              k = k + kb
              go to 20
           end if
           40 continue
           work(1) = lwkopt
           return
     end subroutine la_whetrf_rook

     !> CHETRS_AA: solves a system of linear equations A*X = B with a complex
     !> hermitian matrix A using the factorization A = U**H*T*U or
     !> A = L*T*L**H computed by CHETRF_AA.

     pure subroutine la_chetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        use la_constants_sp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,nrhs,lda,ldb,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           logical(lk) :: lquery,upper
           integer(ilp) :: k,kp,lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(1,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CHETRS_AA',-info)
              return
           else if (lquery) then
              lwkopt = (3*n - 2)
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u**h*t*u.
              ! 1) forward substitution with u**h
              if (n > 1) then
                 ! pivot, p**t * b -> b
                 k = 1
                 do while (k <= n)
                    kp = ipiv(k)
                    if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                    k = k + 1
                 end do
                 ! compute u**h \ b -> b    [ (u**h \p**t * b) ]
                 call la_ctrsm('L','U','C','U',n - 1,nrhs,cone,a(1,2),lda,b(2,1), &
                            ldb)
              end if
              ! 2) solve with triangular matrix t
              ! compute t \ b -> b   [ t \ (u**h \p**t * b) ]
              call la_clacpy('F',1,n,a(1,1),lda + 1,work(n),1)
              if (n > 1) then
                  call la_clacpy('F',1,n - 1,a(1,2),lda + 1,work(2*n),1)
                  call la_clacpy('F',1,n - 1,a(1,2),lda + 1,work(1),1)
                  call la_clacgv(n - 1,work(1),1)
              end if
              call la_cgtsv(n,nrhs,work(1),work(n),work(2*n),b,ldb,info)
              ! 3) backward substitution with u
              if (n > 1) then
                 ! compute u \ b -> b   [ u \ (t \ (u**h \p**t * b) ) ]
                 call la_ctrsm('L','U','N','U',n - 1,nrhs,cone,a(1,2),lda,b(2,1), &
                           ldb)
                 ! pivot, p * b  -> b [ p * (u \ (t \ (u**h \p**t * b) )) ]
                 k = n
                 do while (k >= 1)
                    kp = ipiv(k)
                    if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                    k = k - 1
                 end do
              end if
           else
              ! solve a*x = b, where a = l*t*l**h.
              ! 1) forward substitution with l
              if (n > 1) then
                 ! pivot, p**t * b -> b
                 k = 1
                 do while (k <= n)
                    kp = ipiv(k)
                    if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                    k = k + 1
                 end do
                 ! compute l \ b -> b    [ (l \p**t * b) ]
                 call la_ctrsm('L','L','N','U',n - 1,nrhs,cone,a(2,1),lda,b(2,1), &
                           ldb)
              end if
              ! 2) solve with triangular matrix t
              ! compute t \ b -> b   [ t \ (l \p**t * b) ]
              call la_clacpy('F',1,n,a(1,1),lda + 1,work(n),1)
              if (n > 1) then
                  call la_clacpy('F',1,n - 1,a(2,1),lda + 1,work(1),1)
                  call la_clacpy('F',1,n - 1,a(2,1),lda + 1,work(2*n),1)
                  call la_clacgv(n - 1,work(2*n),1)
              end if
              call la_cgtsv(n,nrhs,work(1),work(n),work(2*n),b,ldb,info)
              ! 3) backward substitution with l**h
              if (n > 1) then
                 ! compute (l**h \ b) -> b   [ l**h \ (t \ (l \p**t * b) ) ]
                 call la_ctrsm('L','L','C','U',n - 1,nrhs,cone,a(2,1),lda,b(2,1), &
                            ldb)
                 ! pivot, p * b -> b  [ p * (l**h \ (t \ (l \p**t * b) )) ]
                 k = n
                 do while (k >= 1)
                    kp = ipiv(k)
                    if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                    k = k - 1
                 end do
              end if
           end if
           return
     end subroutine la_chetrs_aa
     !> ZHETRS_AA: solves a system of linear equations A*X = B with a complex
     !> hermitian matrix A using the factorization A = U**H*T*U or
     !> A = L*T*L**H computed by ZHETRF_AA.

     pure subroutine la_zhetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        use la_constants_dp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,nrhs,lda,ldb,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           logical(lk) :: lquery,upper
           integer(ilp) :: k,kp,lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(1,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZHETRS_AA',-info)
              return
           else if (lquery) then
              lwkopt = (3*n - 2)
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u**h*t*u.
              ! 1) forward substitution with u**h
              if (n > 1) then
                 ! pivot, p**t * b -> b
                 do k = 1,n
                    kp = ipiv(k)
                    if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
                 ! compute u**h \ b -> b    [ (u**h \p**t * b) ]
                 call la_ztrsm('L','U','C','U',n - 1,nrhs,cone,a(1,2),lda,b(2,1), &
                            ldb)
              end if
              ! 2) solve with triangular matrix t
              ! compute t \ b -> b   [ t \ (u**h \p**t * b) ]
              call la_zlacpy('F',1,n,a(1,1),lda + 1,work(n),1)
              if (n > 1) then
                  call la_zlacpy('F',1,n - 1,a(1,2),lda + 1,work(2*n),1)
                  call la_zlacpy('F',1,n - 1,a(1,2),lda + 1,work(1),1)
                  call la_zlacgv(n - 1,work(1),1)
              end if
              call la_zgtsv(n,nrhs,work(1),work(n),work(2*n),b,ldb,info)
              ! 3) backward substitution with u
              if (n > 1) then
                 ! compute u \ b -> b   [ u \ (t \ (u**h \p**t * b) ) ]
                 call la_ztrsm('L','U','N','U',n - 1,nrhs,cone,a(1,2),lda,b(2,1), &
                           ldb)
                 ! pivot, p * b  [ p * (u**h \ (t \ (u \p**t * b) )) ]
                 do k = n,1,-1
                    kp = ipiv(k)
                    if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
              end if
           else
              ! solve a*x = b, where a = l*t*l**h.
              ! 1) forward substitution with l
              if (n > 1) then
                 ! pivot, p**t * b -> b
                 do k = 1,n
                    kp = ipiv(k)
                    if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
                 ! compute l \ b -> b    [ (l \p**t * b) ]
                 call la_ztrsm('L','L','N','U',n - 1,nrhs,cone,a(2,1),lda,b(2,1), &
                           ldb)
              end if
              ! 2) solve with triangular matrix t
              ! compute t \ b -> b   [ t \ (l \p**t * b) ]
              call la_zlacpy('F',1,n,a(1,1),lda + 1,work(n),1)
              if (n > 1) then
                  call la_zlacpy('F',1,n - 1,a(2,1),lda + 1,work(1),1)
                  call la_zlacpy('F',1,n - 1,a(2,1),lda + 1,work(2*n),1)
                  call la_zlacgv(n - 1,work(2*n),1)
              end if
              call la_zgtsv(n,nrhs,work(1),work(n),work(2*n),b,ldb,info)
              ! 3) backward substitution with l**h
              if (n > 1) then
                 ! compute l**h \ b -> b   [ l**h \ (t \ (l \p**t * b) ) ]
                 call la_ztrsm('L','L','C','U',n - 1,nrhs,cone,a(2,1),lda,b(2,1), &
                            ldb)
                 ! pivot, p * b  [ p * (l**h \ (t \ (l \p**t * b) )) ]
                 do k = n,1,-1
                    kp = ipiv(k)
                    if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
              end if
           end if
           return
     end subroutine la_zhetrs_aa
     !> WHETRS_AA: solves a system of linear equations A*X = B with a complex
     !> hermitian matrix A using the factorization A = U**H*T*U or
     !> A = L*T*L**H computed by WHETRF_AA.

     pure subroutine la_whetrs_aa(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,lwork,info)
        use la_constants_qp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,nrhs,lda,ldb,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           logical(lk) :: lquery,upper
           integer(ilp) :: k,kp,lwkopt
           ! Intrinsic Functions
           intrinsic :: max
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           else if (lwork < max(1,3*n - 2) .and. .not. lquery) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WHETRS_AA',-info)
              return
           else if (lquery) then
              lwkopt = (3*n - 2)
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u**h*t*u.
              ! 1) forward substitution with u**h
              if (n > 1) then
                 ! pivot, p**t * b -> b
                 do k = 1,n
                    kp = ipiv(k)
                    if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
                 ! compute u**h \ b -> b    [ (u**h \p**t * b) ]
                 call la_wtrsm('L','U','C','U',n - 1,nrhs,cone,a(1,2),lda,b(2,1), &
                            ldb)
              end if
              ! 2) solve with triangular matrix t
              ! compute t \ b -> b   [ t \ (u**h \p**t * b) ]
              call la_wlacpy('F',1,n,a(1,1),lda + 1,work(n),1)
              if (n > 1) then
                  call la_wlacpy('F',1,n - 1,a(1,2),lda + 1,work(2*n),1)
                  call la_wlacpy('F',1,n - 1,a(1,2),lda + 1,work(1),1)
                  call la_wlacgv(n - 1,work(1),1)
              end if
              call la_wgtsv(n,nrhs,work(1),work(n),work(2*n),b,ldb,info)
              ! 3) backward substitution with u
              if (n > 1) then
                 ! compute u \ b -> b   [ u \ (t \ (u**h \p**t * b) ) ]
                 call la_wtrsm('L','U','N','U',n - 1,nrhs,cone,a(1,2),lda,b(2,1), &
                           ldb)
                 ! pivot, p * b  [ p * (u**h \ (t \ (u \p**t * b) )) ]
                 do k = n,1,-1
                    kp = ipiv(k)
                    if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
              end if
           else
              ! solve a*x = b, where a = l*t*l**h.
              ! 1) forward substitution with l
              if (n > 1) then
                 ! pivot, p**t * b -> b
                 do k = 1,n
                    kp = ipiv(k)
                    if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
                 ! compute l \ b -> b    [ (l \p**t * b) ]
                 call la_wtrsm('L','L','N','U',n - 1,nrhs,cone,a(2,1),lda,b(2,1), &
                           ldb)
              end if
              ! 2) solve with triangular matrix t
              ! compute t \ b -> b   [ t \ (l \p**t * b) ]
              call la_wlacpy('F',1,n,a(1,1),lda + 1,work(n),1)
              if (n > 1) then
                  call la_wlacpy('F',1,n - 1,a(2,1),lda + 1,work(1),1)
                  call la_wlacpy('F',1,n - 1,a(2,1),lda + 1,work(2*n),1)
                  call la_wlacgv(n - 1,work(2*n),1)
              end if
              call la_wgtsv(n,nrhs,work(1),work(n),work(2*n),b,ldb,info)
              ! 3) backward substitution with l**h
              if (n > 1) then
                 ! compute l**h \ b -> b   [ l**h \ (t \ (l \p**t * b) ) ]
                 call la_wtrsm('L','L','C','U',n - 1,nrhs,cone,a(2,1),lda,b(2,1), &
                            ldb)
                 ! pivot, p * b  [ p * (l**h \ (t \ (l \p**t * b) )) ]
                 do k = n,1,-1
                    kp = ipiv(k)
                    if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end do
              end if
           end if
           return
     end subroutine la_whetrs_aa

     !> CHETRS_ROOK: solves a system of linear equations A*X = B with a complex
     !> Hermitian matrix A using the factorization A = U*D*U**H or
     !> A = L*D*L**H computed by CHETRF_ROOK.

     pure subroutine la_chetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(sp) :: s
           complex(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: conjg,max,real
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('CHETRS_ROOK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**h.
              ! first solve u*d*x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              10 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 30
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_cgeru(k - 1,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 ! multiply by the inverse of the diagonal block.
                 s = real(cone,KIND=sp)/real(a(k,k),KIND=sp)
                 call la_csscal(nrhs,s,b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k), then k-1 and -ipiv(k-1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k - 1)
                 if (kp /= k - 1) call la_cswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_cgeru(k - 2,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 call la_cgeru(k - 2,nrhs,-cone,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/conjg(akm1k)
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/conjg(akm1k)
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**h *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**h(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 if (k > 1) then
                    call la_clacgv(nrhs,b(k,1),ldb)
                    call la_cgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k), &
                              1,cone,b(k,1),ldb)
                    call la_clacgv(nrhs,b(k,1),ldb)
                 end if
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**h(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k > 1) then
                    call la_clacgv(nrhs,b(k,1),ldb)
                    call la_cgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k), &
                              1,cone,b(k,1),ldb)
                    call la_clacgv(nrhs,b(k,1),ldb)
                    call la_clacgv(nrhs,b(k + 1,1),ldb)
                    call la_cgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1) &
                              ,1,cone,b(k + 1,1),ldb)
                    call la_clacgv(nrhs,b(k + 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k), then k+1 and -ipiv(k+1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k + 1)
                 if (kp /= k + 1) call la_cswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**h.
              ! first solve l*d*x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              60 continue
              ! if k > n, exit from loop.
              if (k > n) go to 80
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_cgeru(n - k,nrhs,-cone,a(k + 1,k),1,b(k,1),ldb,b( &
                           k + 1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 s = real(cone,KIND=sp)/real(a(k,k),KIND=sp)
                 call la_csscal(nrhs,s,b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k), then k+1 and -ipiv(k+1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k + 1)
                 if (kp /= k + 1) call la_cswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_cgeru(n - k - 1,nrhs,-cone,a(k + 2,k),1,b(k,1),ldb,b(k + 2, &
                              1),ldb)
                    call la_cgeru(n - k - 1,nrhs,-cone,a(k + 2,k + 1),1,b(k + 1,1),ldb,b( &
                              k + 2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/conjg(akm1k)
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/conjg(akm1k)
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**h *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**h(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) then
                    call la_clacgv(nrhs,b(k,1),ldb)
                    call la_cgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k),1,cone,b(k,1),ldb)
                    call la_clacgv(nrhs,b(k,1),ldb)
                 end if
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**h(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_clacgv(nrhs,b(k,1),ldb)
                    call la_cgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k),1,cone,b(k,1),ldb)
                    call la_clacgv(nrhs,b(k,1),ldb)
                    call la_clacgv(nrhs,b(k - 1,1),ldb)
                    call la_cgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k - 1),1,cone,b(k - 1,1),ldb)
                    call la_clacgv(nrhs,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k), then k-1 and -ipiv(k-1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k - 1)
                 if (kp /= k - 1) call la_cswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_chetrs_rook
     !> ZHETRS_ROOK: solves a system of linear equations A*X = B with a complex
     !> Hermitian matrix A using the factorization A = U*D*U**H or
     !> A = L*D*L**H computed by ZHETRF_ROOK.

     pure subroutine la_zhetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(dp) :: s
           complex(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: conjg,max,real
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('ZHETRS_ROOK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**h.
              ! first solve u*d*x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              10 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 30
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_zgeru(k - 1,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 ! multiply by the inverse of the diagonal block.
                 s = real(cone,KIND=dp)/real(a(k,k),KIND=dp)
                 call la_zdscal(nrhs,s,b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k), then k-1 and -ipiv(k-1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k - 1)
                 if (kp /= k - 1) call la_zswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_zgeru(k - 2,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 call la_zgeru(k - 2,nrhs,-cone,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/conjg(akm1k)
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/conjg(akm1k)
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**h *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**h(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 if (k > 1) then
                    call la_zlacgv(nrhs,b(k,1),ldb)
                    call la_zgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k), &
                              1,cone,b(k,1),ldb)
                    call la_zlacgv(nrhs,b(k,1),ldb)
                 end if
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**h(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k > 1) then
                    call la_zlacgv(nrhs,b(k,1),ldb)
                    call la_zgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k), &
                              1,cone,b(k,1),ldb)
                    call la_zlacgv(nrhs,b(k,1),ldb)
                    call la_zlacgv(nrhs,b(k + 1,1),ldb)
                    call la_zgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1) &
                              ,1,cone,b(k + 1,1),ldb)
                    call la_zlacgv(nrhs,b(k + 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k), then k+1 and -ipiv(k+1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k + 1)
                 if (kp /= k + 1) call la_zswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**h.
              ! first solve l*d*x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              60 continue
              ! if k > n, exit from loop.
              if (k > n) go to 80
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_zgeru(n - k,nrhs,-cone,a(k + 1,k),1,b(k,1),ldb,b( &
                           k + 1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 s = real(cone,KIND=dp)/real(a(k,k),KIND=dp)
                 call la_zdscal(nrhs,s,b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k), then k+1 and -ipiv(k+1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k + 1)
                 if (kp /= k + 1) call la_zswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_zgeru(n - k - 1,nrhs,-cone,a(k + 2,k),1,b(k,1),ldb,b(k + 2, &
                              1),ldb)
                    call la_zgeru(n - k - 1,nrhs,-cone,a(k + 2,k + 1),1,b(k + 1,1),ldb,b( &
                              k + 2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/conjg(akm1k)
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/conjg(akm1k)
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**h *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**h(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) then
                    call la_zlacgv(nrhs,b(k,1),ldb)
                    call la_zgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k),1,cone,b(k,1),ldb)
                    call la_zlacgv(nrhs,b(k,1),ldb)
                 end if
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**h(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_zlacgv(nrhs,b(k,1),ldb)
                    call la_zgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k),1,cone,b(k,1),ldb)
                    call la_zlacgv(nrhs,b(k,1),ldb)
                    call la_zlacgv(nrhs,b(k - 1,1),ldb)
                    call la_zgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k - 1),1,cone,b(k - 1,1),ldb)
                    call la_zlacgv(nrhs,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k), then k-1 and -ipiv(k-1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k - 1)
                 if (kp /= k - 1) call la_zswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_zhetrs_rook
     !> WHETRS_ROOK: solves a system of linear equations A*X = B with a complex
     !> Hermitian matrix A using the factorization A = U*D*U**H or
     !> A = L*D*L**H computed by WHETRF_ROOK.

     pure subroutine la_whetrs_rook(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(qp) :: s
           complex(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: conjg,max,real
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           else if (ldb < max(1,n)) then
              info = -8
           end if
           if (info /= 0) then
              call la_xerbla('WHETRS_ROOK',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**h.
              ! first solve u*d*x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              10 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 30
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_wgeru(k - 1,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 ! multiply by the inverse of the diagonal block.
                 s = real(cone,KIND=qp)/real(a(k,k),KIND=qp)
                 call la_wqscal(nrhs,s,b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k), then k-1 and -ipiv(k-1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k - 1)
                 if (kp /= k - 1) call la_wswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_wgeru(k - 2,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 call la_wgeru(k - 2,nrhs,-cone,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/conjg(akm1k)
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/conjg(akm1k)
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**h *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**h(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 if (k > 1) then
                    call la_wlacgv(nrhs,b(k,1),ldb)
                    call la_wgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k), &
                              1,cone,b(k,1),ldb)
                    call la_wlacgv(nrhs,b(k,1),ldb)
                 end if
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**h(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k > 1) then
                    call la_wlacgv(nrhs,b(k,1),ldb)
                    call la_wgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k), &
                              1,cone,b(k,1),ldb)
                    call la_wlacgv(nrhs,b(k,1),ldb)
                    call la_wlacgv(nrhs,b(k + 1,1),ldb)
                    call la_wgemv('CONJUGATE TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1) &
                              ,1,cone,b(k + 1,1),ldb)
                    call la_wlacgv(nrhs,b(k + 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k), then k+1 and -ipiv(k+1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k + 1)
                 if (kp /= k + 1) call la_wswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**h.
              ! first solve l*d*x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              60 continue
              ! if k > n, exit from loop.
              if (k > n) go to 80
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_wgeru(n - k,nrhs,-cone,a(k + 1,k),1,b(k,1),ldb,b( &
                           k + 1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 s = real(cone,KIND=qp)/real(a(k,k),KIND=qp)
                 call la_wqscal(nrhs,s,b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k), then k+1 and -ipiv(k+1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k + 1)
                 if (kp /= k + 1) call la_wswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_wgeru(n - k - 1,nrhs,-cone,a(k + 2,k),1,b(k,1),ldb,b(k + 2, &
                              1),ldb)
                    call la_wgeru(n - k - 1,nrhs,-cone,a(k + 2,k + 1),1,b(k + 1,1),ldb,b( &
                              k + 2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/conjg(akm1k)
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/conjg(akm1k)
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**h *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**h(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) then
                    call la_wlacgv(nrhs,b(k,1),ldb)
                    call la_wgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k),1,cone,b(k,1),ldb)
                    call la_wlacgv(nrhs,b(k,1),ldb)
                 end if
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**h(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_wlacgv(nrhs,b(k,1),ldb)
                    call la_wgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k),1,cone,b(k,1),ldb)
                    call la_wlacgv(nrhs,b(k,1),ldb)
                    call la_wlacgv(nrhs,b(k - 1,1),ldb)
                    call la_wgemv('CONJUGATE TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb, &
                              a(k + 1,k - 1),1,cone,b(k - 1,1),ldb)
                    call la_wlacgv(nrhs,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k), then k-1 and -ipiv(k-1)
                 kp = -ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 kp = -ipiv(k - 1)
                 if (kp /= k - 1) call la_wswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_whetrs_rook

     !> CHECON_ROOK: estimates the reciprocal of the condition number of a complex
     !> Hermitian matrix A using the factorization A = U*D*U**H or
     !> A = L*D*L**H computed by CHETRF_ROOK.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_checon_rook(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(in) :: anorm
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,kase
           real(sp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
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
           else if (anorm < zero) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('CHECON_ROOK',-info)
              return
           end if
           ! quick return if possible
           rcond = zero
           if (n == 0) then
              rcond = one
              return
           else if (anorm <= zero) then
              return
           end if
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. a(i,i) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do i = 1,n
                 if (ipiv(i) > 0 .and. a(i,i) == zero) return
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**h) or inv(u*d*u**h).
              call la_chetrs_rook(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_checon_rook
     !> ZHECON_ROOK: estimates the reciprocal of the condition number of a complex
     !> Hermitian matrix A using the factorization A = U*D*U**H or
     !> A = L*D*L**H computed by CHETRF_ROOK.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_zhecon_rook(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(in) :: anorm
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,kase
           real(dp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
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
           else if (anorm < zero) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('ZHECON_ROOK',-info)
              return
           end if
           ! quick return if possible
           rcond = zero
           if (n == 0) then
              rcond = one
              return
           else if (anorm <= zero) then
              return
           end if
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. a(i,i) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do i = 1,n
                 if (ipiv(i) > 0 .and. a(i,i) == zero) return
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**h) or inv(u*d*u**h).
              call la_zhetrs_rook(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_zhecon_rook
     !> WHECON_ROOK: estimates the reciprocal of the condition number of a complex
     !> Hermitian matrix A using the factorization A = U*D*U**H or
     !> A = L*D*L**H computed by ZHETRF_ROOK.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_whecon_rook(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(in) :: anorm
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,kase
           real(qp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
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
           else if (anorm < zero) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('WHECON_ROOK',-info)
              return
           end if
           ! quick return if possible
           rcond = zero
           if (n == 0) then
              rcond = one
              return
           else if (anorm <= zero) then
              return
           end if
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. a(i,i) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do i = 1,n
                 if (ipiv(i) > 0 .and. a(i,i) == zero) return
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**h) or inv(u*d*u**h).
              call la_whetrs_rook(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_whecon_rook

     !> CHPRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is Hermitian indefinite
     !> and packed, and provides error bounds and backward error estimates
     !> for the solution.

     pure subroutine la_chprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                rwork,info)
        use la_constants_sp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: afp(*),ap(*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,ik,j,k,kase,kk,nz
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
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('CHPRFS',-info)
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
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_slamch('EPSILON')
           safmin = la_slamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_140: do j = 1,nrhs
              count = 1
              lstres = three
              20 continue
              ! loop until stopping criterion is satisfied.
              ! compute residual r = b - a * x
              call la_ccopy(n,b(1,j),1,work,1)
              call la_chpmv(uplo,n,-cone,ap,x(1,j),1,cone,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              kk = 1
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    ik = kk
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(ap(ik))*xk
                       s = s + cabs1(ap(ik))*cabs1(x(i,j))
                       ik = ik + 1
                    end do
                    rwork(k) = rwork(k) + abs(real(ap(kk + k - 1),KIND=sp))*xk + s
                    kk = kk + k
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + abs(real(ap(kk),KIND=sp))*xk
                    ik = kk + 1
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(ap(ik))*xk
                       s = s + cabs1(ap(ik))*cabs1(x(i,j))
                       ik = ik + 1
                    end do
                    rwork(k) = rwork(k) + s
                    kk = kk + (n - k + 1)
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_chptrs(uplo,n,1,afp,ipiv,work,n,info)
                 call la_caxpy(n,cone,work,1,x(1,j),1)
                 lstres = berr(j)
                 count = count + 1
                 go to 20
              end if
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(a))*
                 ! ( abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(a) is the inverse of a
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(a)*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(a)*abs(x) + abs(b) is less than safe2.
              ! use la_clacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_clacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**h).
                    call la_chptrs(uplo,n,1,afp,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_chptrs(uplo,n,1,afp,ipiv,work,n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_chprfs
     !> ZHPRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is Hermitian indefinite
     !> and packed, and provides error bounds and backward error estimates
     !> for the solution.

     pure subroutine la_zhprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                rwork,info)
        use la_constants_dp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: afp(*),ap(*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,ik,j,k,kase,kk,nz
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
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('ZHPRFS',-info)
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
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_dlamch('EPSILON')
           safmin = la_dlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_140: do j = 1,nrhs
              count = 1
              lstres = three
              20 continue
              ! loop until stopping criterion is satisfied.
              ! compute residual r = b - a * x
              call la_zcopy(n,b(1,j),1,work,1)
              call la_zhpmv(uplo,n,-cone,ap,x(1,j),1,cone,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              kk = 1
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    ik = kk
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(ap(ik))*xk
                       s = s + cabs1(ap(ik))*cabs1(x(i,j))
                       ik = ik + 1
                    end do
                    rwork(k) = rwork(k) + abs(real(ap(kk + k - 1),KIND=dp))*xk + s
                    kk = kk + k
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + abs(real(ap(kk),KIND=dp))*xk
                    ik = kk + 1
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(ap(ik))*xk
                       s = s + cabs1(ap(ik))*cabs1(x(i,j))
                       ik = ik + 1
                    end do
                    rwork(k) = rwork(k) + s
                    kk = kk + (n - k + 1)
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_zhptrs(uplo,n,1,afp,ipiv,work,n,info)
                 call la_zaxpy(n,cone,work,1,x(1,j),1)
                 lstres = berr(j)
                 count = count + 1
                 go to 20
              end if
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(a))*
                 ! ( abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(a) is the inverse of a
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(a)*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(a)*abs(x) + abs(b) is less than safe2.
              ! use la_zlacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_zlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**h).
                    call la_zhptrs(uplo,n,1,afp,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_zhptrs(uplo,n,1,afp,ipiv,work,n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_zhprfs
     !> WHPRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is Hermitian indefinite
     !> and packed, and provides error bounds and backward error estimates
     !> for the solution.

     pure subroutine la_whprfs(uplo,n,nrhs,ap,afp,ipiv,b,ldb,x,ldx,ferr,berr,work, &
                rwork,info)
        use la_constants_qp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: afp(*),ap(*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,ik,j,k,kase,kk,nz
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
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (nrhs < 0) then
              info = -3
           else if (ldb < max(1,n)) then
              info = -8
           else if (ldx < max(1,n)) then
              info = -10
           end if
           if (info /= 0) then
              call la_xerbla('WHPRFS',-info)
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
           ! nz = maximum number of nonzero elements in each row of a, plus 1
           nz = n + 1
           eps = la_qlamch('EPSILON')
           safmin = la_qlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_140: do j = 1,nrhs
              count = 1
              lstres = three
              20 continue
              ! loop until stopping criterion is satisfied.
              ! compute residual r = b - a * x
              call la_wcopy(n,b(1,j),1,work,1)
              call la_whpmv(uplo,n,-cone,ap,x(1,j),1,cone,work,1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 rwork(i) = cabs1(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              kk = 1
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    ik = kk
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(ap(ik))*xk
                       s = s + cabs1(ap(ik))*cabs1(x(i,j))
                       ik = ik + 1
                    end do
                    rwork(k) = rwork(k) + abs(real(ap(kk + k - 1),KIND=qp))*xk + s
                    kk = kk + k
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + abs(real(ap(kk),KIND=qp))*xk
                    ik = kk + 1
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(ap(ik))*xk
                       s = s + cabs1(ap(ik))*cabs1(x(i,j))
                       ik = ik + 1
                    end do
                    rwork(k) = rwork(k) + s
                    kk = kk + (n - k + 1)
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_whptrs(uplo,n,1,afp,ipiv,work,n,info)
                 call la_waxpy(n,cone,work,1,x(1,j),1)
                 lstres = berr(j)
                 count = count + 1
                 go to 20
              end if
              ! bound error from formula
              ! norm(x - xtrue) / norm(x) .le. ferr =
              ! norm( abs(inv(a))*
                 ! ( abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) ))) / norm(x)
              ! where
                ! norm(z) is the magnitude of the largest component of z
                ! inv(a) is the inverse of a
                ! abs(z) is the componentwise absolute value of the matrix or
                   ! vector z
                ! nz is the maximum number of nonzeros in any row of a, plus 1
                ! eps is machine epsilon
              ! the i-th component of abs(r)+nz*eps*(abs(a)*abs(x)+abs(b))
              ! is incremented by safe1 if the i-th component of
              ! abs(a)*abs(x) + abs(b) is less than safe2.
              ! use la_wlacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (rwork(i) > safe2) then
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i)
                 else
                    rwork(i) = cabs1(work(i)) + nz*eps*rwork(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_wlacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**h).
                    call la_whptrs(uplo,n,1,afp,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_whptrs(uplo,n,1,afp,ipiv,work,n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,cabs1(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_whprfs

     !> CHETRF_AA: computes the factorization of a complex hermitian matrix A
     !> using the Aasen's algorithm.  The form of the factorization is
     !> A = U**H*T*U  or  A = L*T*L**H
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is a hermitian tridiagonal matrix.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_chetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,lda,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: j,lwkopt
           integer(ilp) :: nb,mj,nj,k1,k2,j1,j2,j3,jb
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
           ! Executable Statements
           ! determine the block size
           nb = la_ilaenv(1,'CHETRF_AA',uplo,n,-1,-1,-1)
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < (2*n) .and. .not. lquery) then
              info = -7
           end if
           if (info == 0) then
              lwkopt = (nb + 1)*n
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CHETRF_AA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return
           if (n == 0) then
               return
           end if
           ipiv(1) = 1
           if (n == 1) then
              a(1,1) = real(a(1,1),KIND=sp)
              return
           end if
           ! adjust block size based on the workspace size
           if (lwork < ((1 + nb)*n)) then
              nb = (lwork - n)/n
           end if
           if (upper) then
              ! .....................................................
              ! factorize a as u**h*d*u using the upper triangle of a
              ! .....................................................
              ! copy first row a(1, 1:n) into h(1:n) (stored in work(1:n))
              call la_ccopy(n,a(1,1),lda,work(1),1)
              ! j is the main loop index, increasing from 1 to n in steps of
              ! jb, where jb is the number of columns factorized by la_clahef;
              ! jb is either nb, or n-j+1 for the last block
              j = 0
              10 continue
              if (j >= n) go to 20
              ! each step of the main loop
               ! j is the last column of the previous panel
               ! j1 is the first column of the current panel
               ! k1 identifies if the previous column of the panel has been
                ! explicitly stored, e.g., k1=1 for the first panel, and
                ! k1=0 for the rest
              j1 = j + 1
              jb = min(n - j1 + 1,nb)
              k1 = max(1,j) - j
              ! panel factorization
              call la_clahef_aa(uplo,2 - k1,n - j,jb,a(max(1,j),j + 1),lda,ipiv(j + 1), &
                        work,n,work(n*nb + 1))
              ! adjust ipiv and apply it back (j-th step picks (j+1)-th pivot)
              do j2 = j + 2,min(n,j + jb + 1)
                 ipiv(j2) = ipiv(j2) + j
                 if ((j2 /= ipiv(j2)) .and. ((j1 - k1) > 2)) then
                    call la_cswap(j1 - k1 - 2,a(1,j2),1,a(1,ipiv(j2)),1)
                 end if
              end do
              j = j + jb
              ! trailing submatrix update, where
               ! the row a(j1-1, j2-1:n) stores u(j1, j2+1:n) and
               ! work stores the current block of the auxiriarly matrix h
              if (j < n) then
                ! if the first panel and jb=1 (nb=1), then nothing to do
                 if (j1 > 1 .or. jb > 1) then
                    ! merge rank-1 update with blas-3 update
                    alpha = conjg(a(j,j + 1))
                    a(j,j + 1) = cone
                    call la_ccopy(n - j,a(j - 1,j + 1),lda,work((j + 1 - j1 + 1) + jb*n),1)

                    call la_cscal(n - j,alpha,work((j + 1 - j1 + 1) + jb*n),1)
                    ! k1 identifies if the previous column of the panel has been
                     ! explicitly stored, e.g., k1=0 and k2=1 for the first panel,
                     ! and k1=1 and k2=0 for the rest
                    if (j1 > 1) then
                       ! not first panel
                       k2 = 1
                    else
                       ! first panel
                       k2 = 0
                       ! first update skips the first column
                       jb = jb - 1
                    end if
                    do j2 = j + 1,n,nb
                       nj = min(nb,n - j2 + 1)
                       ! update (j2, j2) diagonal block with la_cgemv
                       j3 = j2
                       do mj = nj - 1,1,-1
                          call la_cgemm('CONJUGATE TRANSPOSE','TRANSPOSE',1,mj,jb + 1,-cone, &
                           a(j1 - k2,j3),lda,work((j3 - j1 + 1) + k1*n),n,cone,a(j3,j3),lda)

                          j3 = j3 + 1
                       end do
                       ! update off-diagonal block of j2-th block row with la_cgemm
                       call la_cgemm('CONJUGATE TRANSPOSE','TRANSPOSE',nj,n - j3 + 1,jb + 1,- &
                       cone,a(j1 - k2,j2),lda,work((j3 - j1 + 1) + k1*n),n,cone,a(j2,j3),lda &
                                 )
                    end do
                    ! recover t( j, j+1 )
                    a(j,j + 1) = conjg(alpha)
                 end if
                 ! work(j+1, 1) stores h(j+1, 1)
                 call la_ccopy(n - j,a(j + 1,j + 1),lda,work(1),1)
              end if
              go to 10
           else
              ! .....................................................
              ! factorize a as l*d*l**h using the lower triangle of a
              ! .....................................................
              ! copy first column a(1:n, 1) into h(1:n, 1)
               ! (stored in work(1:n))
              call la_ccopy(n,a(1,1),1,work(1),1)
              ! j is the main loop index, increasing from 1 to n in steps of
              ! jb, where jb is the number of columns factorized by la_clahef;
              ! jb is either nb, or n-j+1 for the last block
              j = 0
              11 continue
              if (j >= n) go to 20
              ! each step of the main loop
               ! j is the last column of the previous panel
               ! j1 is the first column of the current panel
               ! k1 identifies if the previous column of the panel has been
                ! explicitly stored, e.g., k1=1 for the first panel, and
                ! k1=0 for the rest
              j1 = j + 1
              jb = min(n - j1 + 1,nb)
              k1 = max(1,j) - j
              ! panel factorization
              call la_clahef_aa(uplo,2 - k1,n - j,jb,a(j + 1,max(1,j)),lda,ipiv(j + 1), &
                        work,n,work(n*nb + 1))
              ! adjust ipiv and apply it back (j-th step picks (j+1)-th pivot)
              do j2 = j + 2,min(n,j + jb + 1)
                 ipiv(j2) = ipiv(j2) + j
                 if ((j2 /= ipiv(j2)) .and. ((j1 - k1) > 2)) then
                    call la_cswap(j1 - k1 - 2,a(j2,1),lda,a(ipiv(j2),1),lda)
                 end if
              end do
              j = j + jb
              ! trailing submatrix update, where
                ! a(j2+1, j1-1) stores l(j2+1, j1) and
                ! work(j2+1, 1) stores h(j2+1, 1)
              if (j < n) then
                ! if the first panel and jb=1 (nb=1), then nothing to do
                 if (j1 > 1 .or. jb > 1) then
                    ! merge rank-1 update with blas-3 update
                    alpha = conjg(a(j + 1,j))
                    a(j + 1,j) = cone
                    call la_ccopy(n - j,a(j + 1,j - 1),1,work((j + 1 - j1 + 1) + jb*n),1)
                    call la_cscal(n - j,alpha,work((j + 1 - j1 + 1) + jb*n),1)
                    ! k1 identifies if the previous column of the panel has been
                     ! explicitly stored, e.g., k1=0 and k2=1 for the first panel,
                     ! and k1=1 and k2=0 for the rest
                    if (j1 > 1) then
                       ! not first panel
                       k2 = 1
                    else
                       ! first panel
                       k2 = 0
                       ! first update skips the first column
                       jb = jb - 1
                    end if
                    do j2 = j + 1,n,nb
                       nj = min(nb,n - j2 + 1)
                       ! update (j2, j2) diagonal block with la_cgemv
                       j3 = j2
                       do mj = nj - 1,1,-1
                          call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',mj,1,jb + 1,- &
                          cone,work((j3 - j1 + 1) + k1*n),n,a(j3,j1 - k2),lda,cone,a(j3,j3), &
                                    lda)
                          j3 = j3 + 1
                       end do
                       ! update off-diagonal block of j2-th block column with la_cgemm
                       call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',n - j3 + 1,nj,jb + 1,- &
                       cone,work((j3 - j1 + 1) + k1*n),n,a(j2,j1 - k2),lda,cone,a(j3,j2),lda &
                                 )
                    end do
                    ! recover t( j+1, j )
                    a(j + 1,j) = conjg(alpha)
                 end if
                 ! work(j+1, 1) stores h(j+1, 1)
                 call la_ccopy(n - j,a(j + 1,j + 1),1,work(1),1)
              end if
              go to 11
           end if
           20 continue
           work(1) = lwkopt
           return
     end subroutine la_chetrf_aa
     !> ZHETRF_AA: computes the factorization of a complex hermitian matrix A
     !> using the Aasen's algorithm.  The form of the factorization is
     !> A = U**H*T*U  or  A = L*T*L**H
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is a hermitian tridiagonal matrix.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_zhetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,lda,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: j,lwkopt
           integer(ilp) :: nb,mj,nj,k1,k2,j1,j2,j3,jb
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
           ! Executable Statements
           ! determine the block size
           nb = la_ilaenv(1,'ZHETRF_AA',uplo,n,-1,-1,-1)
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -7
           end if
           if (info == 0) then
              lwkopt = (nb + 1)*n
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZHETRF_AA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return
           if (n == 0) then
               return
           end if
           ipiv(1) = 1
           if (n == 1) then
              a(1,1) = real(a(1,1),KIND=dp)
              return
           end if
           ! adjust block size based on the workspace size
           if (lwork < ((1 + nb)*n)) then
              nb = (lwork - n)/n
           end if
           if (upper) then
              ! .....................................................
              ! factorize a as u**h*d*u using the upper triangle of a
              ! .....................................................
              ! copy first row a(1, 1:n) into h(1:n) (stored in work(1:n))
              call la_zcopy(n,a(1,1),lda,work(1),1)
              ! j is the main loop index, increasing from 1 to n in steps of
              ! jb, where jb is the number of columns factorized by la_zlahef;
              ! jb is either nb, or n-j+1 for the last block
              j = 0
              10 continue
              if (j >= n) go to 20
              ! each step of the main loop
               ! j is the last column of the previous panel
               ! j1 is the first column of the current panel
               ! k1 identifies if the previous column of the panel has been
                ! explicitly stored, e.g., k1=1 for the first panel, and
                ! k1=0 for the rest
              j1 = j + 1
              jb = min(n - j1 + 1,nb)
              k1 = max(1,j) - j
              ! panel factorization
              call la_zlahef_aa(uplo,2 - k1,n - j,jb,a(max(1,j),j + 1),lda,ipiv(j + 1), &
                        work,n,work(n*nb + 1))
              ! adjust ipiv and apply it back (j-th step picks (j+1)-th pivot)
              do j2 = j + 2,min(n,j + jb + 1)
                 ipiv(j2) = ipiv(j2) + j
                 if ((j2 /= ipiv(j2)) .and. ((j1 - k1) > 2)) then
                    call la_zswap(j1 - k1 - 2,a(1,j2),1,a(1,ipiv(j2)),1)
                 end if
              end do
              j = j + jb
              ! trailing submatrix update, where
               ! the row a(j1-1, j2-1:n) stores u(j1, j2+1:n) and
               ! work stores the current block of the auxiriarly matrix h
              if (j < n) then
                ! if the first panel and jb=1 (nb=1), then nothing to do
                 if (j1 > 1 .or. jb > 1) then
                    ! merge rank-1 update with blas-3 update
                    alpha = conjg(a(j,j + 1))
                    a(j,j + 1) = cone
                    call la_zcopy(n - j,a(j - 1,j + 1),lda,work((j + 1 - j1 + 1) + jb*n),1)

                    call la_zscal(n - j,alpha,work((j + 1 - j1 + 1) + jb*n),1)
                    ! k1 identifies if the previous column of the panel has been
                     ! explicitly stored, e.g., k1=0 and k2=1 for the first panel,
                     ! and k1=1 and k2=0 for the rest
                    if (j1 > 1) then
                       ! not first panel
                       k2 = 1
                    else
                       ! first panel
                       k2 = 0
                       ! first update skips the first column
                       jb = jb - 1
                    end if
                    do j2 = j + 1,n,nb
                       nj = min(nb,n - j2 + 1)
                       ! update (j2, j2) diagonal block with la_zgemv
                       j3 = j2
                       do mj = nj - 1,1,-1
                          call la_zgemm('CONJUGATE TRANSPOSE','TRANSPOSE',1,mj,jb + 1,-cone, &
                           a(j1 - k2,j3),lda,work((j3 - j1 + 1) + k1*n),n,cone,a(j3,j3),lda)

                          j3 = j3 + 1
                       end do
                       ! update off-diagonal block of j2-th block row with la_zgemm
                       call la_zgemm('CONJUGATE TRANSPOSE','TRANSPOSE',nj,n - j3 + 1,jb + 1,- &
                       cone,a(j1 - k2,j2),lda,work((j3 - j1 + 1) + k1*n),n,cone,a(j2,j3),lda &
                                 )
                    end do
                    ! recover t( j, j+1 )
                    a(j,j + 1) = conjg(alpha)
                 end if
                 ! work(j+1, 1) stores h(j+1, 1)
                 call la_zcopy(n - j,a(j + 1,j + 1),lda,work(1),1)
              end if
              go to 10
           else
              ! .....................................................
              ! factorize a as l*d*l**h using the lower triangle of a
              ! .....................................................
              ! copy first column a(1:n, 1) into h(1:n, 1)
               ! (stored in work(1:n))
              call la_zcopy(n,a(1,1),1,work(1),1)
              ! j is the main loop index, increasing from 1 to n in steps of
              ! jb, where jb is the number of columns factorized by la_zlahef;
              ! jb is either nb, or n-j+1 for the last block
              j = 0
              11 continue
              if (j >= n) go to 20
              ! each step of the main loop
               ! j is the last column of the previous panel
               ! j1 is the first column of the current panel
               ! k1 identifies if the previous column of the panel has been
                ! explicitly stored, e.g., k1=1 for the first panel, and
                ! k1=0 for the rest
              j1 = j + 1
              jb = min(n - j1 + 1,nb)
              k1 = max(1,j) - j
              ! panel factorization
              call la_zlahef_aa(uplo,2 - k1,n - j,jb,a(j + 1,max(1,j)),lda,ipiv(j + 1), &
                        work,n,work(n*nb + 1))
              ! adjust ipiv and apply it back (j-th step picks (j+1)-th pivot)
              do j2 = j + 2,min(n,j + jb + 1)
                 ipiv(j2) = ipiv(j2) + j
                 if ((j2 /= ipiv(j2)) .and. ((j1 - k1) > 2)) then
                    call la_zswap(j1 - k1 - 2,a(j2,1),lda,a(ipiv(j2),1),lda)
                 end if
              end do
              j = j + jb
              ! trailing submatrix update, where
                ! a(j2+1, j1-1) stores l(j2+1, j1) and
                ! work(j2+1, 1) stores h(j2+1, 1)
              if (j < n) then
                ! if the first panel and jb=1 (nb=1), then nothing to do
                 if (j1 > 1 .or. jb > 1) then
                    ! merge rank-1 update with blas-3 update
                    alpha = conjg(a(j + 1,j))
                    a(j + 1,j) = cone
                    call la_zcopy(n - j,a(j + 1,j - 1),1,work((j + 1 - j1 + 1) + jb*n),1)
                    call la_zscal(n - j,alpha,work((j + 1 - j1 + 1) + jb*n),1)
                    ! k1 identifies if the previous column of the panel has been
                     ! explicitly stored, e.g., k1=0 and k2=1 for the first panel,
                     ! and k1=1 and k2=0 for the rest
                    if (j1 > 1) then
                       ! not first panel
                       k2 = 1
                    else
                       ! first panel
                       k2 = 0
                       ! first update skips the first column
                       jb = jb - 1
                    end if
                    do j2 = j + 1,n,nb
                       nj = min(nb,n - j2 + 1)
                       ! update (j2, j2) diagonal block with la_zgemv
                       j3 = j2
                       do mj = nj - 1,1,-1
                          call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',mj,1,jb + 1,- &
                          cone,work((j3 - j1 + 1) + k1*n),n,a(j3,j1 - k2),lda,cone,a(j3,j3), &
                                    lda)
                          j3 = j3 + 1
                       end do
                       ! update off-diagonal block of j2-th block column with la_zgemm
                       call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',n - j3 + 1,nj,jb + 1,- &
                       cone,work((j3 - j1 + 1) + k1*n),n,a(j2,j1 - k2),lda,cone,a(j3,j2),lda &
                                 )
                    end do
                    ! recover t( j+1, j )
                    a(j + 1,j) = conjg(alpha)
                 end if
                 ! work(j+1, 1) stores h(j+1, 1)
                 call la_zcopy(n - j,a(j + 1,j + 1),1,work(1),1)
              end if
              go to 11
           end if
           20 continue
           work(1) = lwkopt
           return
     end subroutine la_zhetrf_aa
     !> WHETRF_AA: computes the factorization of a complex hermitian matrix A
     !> using the Aasen's algorithm.  The form of the factorization is
     !> A = U**H*T*U  or  A = L*T*L**H
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and T is a hermitian tridiagonal matrix.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_whetrf_aa(uplo,n,a,lda,ipiv,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,lda,lwork
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,upper
           integer(ilp) :: j,lwkopt
           integer(ilp) :: nb,mj,nj,k1,k2,j1,j2,j3,jb
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
           ! Executable Statements
           ! determine the block size
           nb = la_ilaenv(1,'WHETRF_AA',uplo,n,-1,-1,-1)
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           lquery = (lwork == -1)
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           else if (lwork < max(1,2*n) .and. .not. lquery) then
              info = -7
           end if
           if (info == 0) then
              lwkopt = (nb + 1)*n
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WHETRF_AA',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return
           if (n == 0) then
               return
           end if
           ipiv(1) = 1
           if (n == 1) then
              a(1,1) = real(a(1,1),KIND=qp)
              return
           end if
           ! adjust block size based on the workspace size
           if (lwork < ((1 + nb)*n)) then
              nb = (lwork - n)/n
           end if
           if (upper) then
              ! .....................................................
              ! factorize a as u**h*d*u using the upper triangle of a
              ! .....................................................
              ! copy first row a(1, 1:n) into h(1:n) (stored in work(1:n))
              call la_wcopy(n,a(1,1),lda,work(1),1)
              ! j is the main loop index, increasing from 1 to n in steps of
              ! jb, where jb is the number of columns factorized by la_wlahef;
              ! jb is either nb, or n-j+1 for the last block
              j = 0
              10 continue
              if (j >= n) go to 20
              ! each step of the main loop
               ! j is the last column of the previous panel
               ! j1 is the first column of the current panel
               ! k1 identifies if the previous column of the panel has been
                ! explicitly stored, e.g., k1=1 for the first panel, and
                ! k1=0 for the rest
              j1 = j + 1
              jb = min(n - j1 + 1,nb)
              k1 = max(1,j) - j
              ! panel factorization
              call la_wlahef_aa(uplo,2 - k1,n - j,jb,a(max(1,j),j + 1),lda,ipiv(j + 1), &
                        work,n,work(n*nb + 1))
              ! adjust ipiv and apply it back (j-th step picks (j+1)-th pivot)
              do j2 = j + 2,min(n,j + jb + 1)
                 ipiv(j2) = ipiv(j2) + j
                 if ((j2 /= ipiv(j2)) .and. ((j1 - k1) > 2)) then
                    call la_wswap(j1 - k1 - 2,a(1,j2),1,a(1,ipiv(j2)),1)
                 end if
              end do
              j = j + jb
              ! trailing submatrix update, where
               ! the row a(j1-1, j2-1:n) stores u(j1, j2+1:n) and
               ! work stores the current block of the auxiriarly matrix h
              if (j < n) then
                ! if the first panel and jb=1 (nb=1), then nothing to do
                 if (j1 > 1 .or. jb > 1) then
                    ! merge rank-1 update with blas-3 update
                    alpha = conjg(a(j,j + 1))
                    a(j,j + 1) = cone
                    call la_wcopy(n - j,a(j - 1,j + 1),lda,work((j + 1 - j1 + 1) + jb*n),1)

                    call la_wscal(n - j,alpha,work((j + 1 - j1 + 1) + jb*n),1)
                    ! k1 identifies if the previous column of the panel has been
                     ! explicitly stored, e.g., k1=0 and k2=1 for the first panel,
                     ! and k1=1 and k2=0 for the rest
                    if (j1 > 1) then
                       ! not first panel
                       k2 = 1
                    else
                       ! first panel
                       k2 = 0
                       ! first update skips the first column
                       jb = jb - 1
                    end if
                    do j2 = j + 1,n,nb
                       nj = min(nb,n - j2 + 1)
                       ! update (j2, j2) diagonal block with la_wgemv
                       j3 = j2
                       do mj = nj - 1,1,-1
                          call la_wgemm('CONJUGATE TRANSPOSE','TRANSPOSE',1,mj,jb + 1,-cone, &
                           a(j1 - k2,j3),lda,work((j3 - j1 + 1) + k1*n),n,cone,a(j3,j3),lda)

                          j3 = j3 + 1
                       end do
                       ! update off-diagonal block of j2-th block row with la_wgemm
                       call la_wgemm('CONJUGATE TRANSPOSE','TRANSPOSE',nj,n - j3 + 1,jb + 1,- &
                       cone,a(j1 - k2,j2),lda,work((j3 - j1 + 1) + k1*n),n,cone,a(j2,j3),lda &
                                 )
                    end do
                    ! recover t( j, j+1 )
                    a(j,j + 1) = conjg(alpha)
                 end if
                 ! work(j+1, 1) stores h(j+1, 1)
                 call la_wcopy(n - j,a(j + 1,j + 1),lda,work(1),1)
              end if
              go to 10
           else
              ! .....................................................
              ! factorize a as l*d*l**h using the lower triangle of a
              ! .....................................................
              ! copy first column a(1:n, 1) into h(1:n, 1)
               ! (stored in work(1:n))
              call la_wcopy(n,a(1,1),1,work(1),1)
              ! j is the main loop index, increasing from 1 to n in steps of
              ! jb, where jb is the number of columns factorized by la_wlahef;
              ! jb is either nb, or n-j+1 for the last block
              j = 0
              11 continue
              if (j >= n) go to 20
              ! each step of the main loop
               ! j is the last column of the previous panel
               ! j1 is the first column of the current panel
               ! k1 identifies if the previous column of the panel has been
                ! explicitly stored, e.g., k1=1 for the first panel, and
                ! k1=0 for the rest
              j1 = j + 1
              jb = min(n - j1 + 1,nb)
              k1 = max(1,j) - j
              ! panel factorization
              call la_wlahef_aa(uplo,2 - k1,n - j,jb,a(j + 1,max(1,j)),lda,ipiv(j + 1), &
                        work,n,work(n*nb + 1))
              ! adjust ipiv and apply it back (j-th step picks (j+1)-th pivot)
              do j2 = j + 2,min(n,j + jb + 1)
                 ipiv(j2) = ipiv(j2) + j
                 if ((j2 /= ipiv(j2)) .and. ((j1 - k1) > 2)) then
                    call la_wswap(j1 - k1 - 2,a(j2,1),lda,a(ipiv(j2),1),lda)
                 end if
              end do
              j = j + jb
              ! trailing submatrix update, where
                ! a(j2+1, j1-1) stores l(j2+1, j1) and
                ! work(j2+1, 1) stores h(j2+1, 1)
              if (j < n) then
                ! if the first panel and jb=1 (nb=1), then nothing to do
                 if (j1 > 1 .or. jb > 1) then
                    ! merge rank-1 update with blas-3 update
                    alpha = conjg(a(j + 1,j))
                    a(j + 1,j) = cone
                    call la_wcopy(n - j,a(j + 1,j - 1),1,work((j + 1 - j1 + 1) + jb*n),1)
                    call la_wscal(n - j,alpha,work((j + 1 - j1 + 1) + jb*n),1)
                    ! k1 identifies if the previous column of the panel has been
                     ! explicitly stored, e.g., k1=0 and k2=1 for the first panel,
                     ! and k1=1 and k2=0 for the rest
                    if (j1 > 1) then
                       ! not first panel
                       k2 = 1
                    else
                       ! first panel
                       k2 = 0
                       ! first update skips the first column
                       jb = jb - 1
                    end if
                    do j2 = j + 1,n,nb
                       nj = min(nb,n - j2 + 1)
                       ! update (j2, j2) diagonal block with la_wgemv
                       j3 = j2
                       do mj = nj - 1,1,-1
                          call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',mj,1,jb + 1,- &
                          cone,work((j3 - j1 + 1) + k1*n),n,a(j3,j1 - k2),lda,cone,a(j3,j3), &
                                    lda)
                          j3 = j3 + 1
                       end do
                       ! update off-diagonal block of j2-th block column with la_wgemm
                       call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',n - j3 + 1,nj,jb + 1,- &
                       cone,work((j3 - j1 + 1) + k1*n),n,a(j2,j1 - k2),lda,cone,a(j3,j2),lda &
                                 )
                    end do
                    ! recover t( j+1, j )
                    a(j + 1,j) = conjg(alpha)
                 end if
                 ! work(j+1, 1) stores h(j+1, 1)
                 call la_wcopy(n - j,a(j + 1,j + 1),1,work(1),1)
              end if
              go to 11
           end if
           20 continue
           work(1) = lwkopt
           return
     end subroutine la_whetrf_aa

     !> CLAHEF_AA: factorizes a panel of a complex hermitian matrix A using
     !> the Aasen's algorithm. The panel consists of a set of NB rows of A
     !> when UPLO is U, or a set of NB columns when UPLO is L.
     !> In order to factorize the panel, the Aasen's algorithm requires the
     !> last row, or column, of the previous panel. The first row, or column,
     !> of A is set to be the first row, or column, of an identity matrix,
     !> which is used to factorize the first panel.
     !> The resulting J-th row of U, or J-th column of L, is stored in the
     !> (J-1)-th row, or column, of A (without the unit diagonals), while
     !> the diagonal and subdiagonal of A are overwritten by those of T.

     pure subroutine la_clahef_aa(uplo,j1,m,nb,a,lda,ipiv,h,ldh,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: m,nb,j1,lda,ldh
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*),h(ldh,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k,k1,i1,i2,mj
           complex(sp) :: piv,alpha
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
           ! Executable Statements
           j = 1
           ! k1 is the first column of the panel to be factorized
           ! i.e.,  k1 is 2 for the first block column, and 1 for the rest of the blocks
           k1 = (2 - j1) + 1
           if (la_lsame(uplo,'U')) then
              ! .....................................................
              ! factorize a as u**t*d*u using the upper triangle of a
              ! .....................................................
              10 continue
              if (j > min(m,nb)) go to 20
              ! k is the column to be factorized
               ! when being called from la_chetrf_aa,
               ! > for the first block column, j1 is 1, hence j1+j-1 is j,
               ! > for the rest of the columns, j1 is 2, and j1+j-1 is j+1,
              k = j1 + j - 1
              if (j == m) then
                  ! only need to compute t(j, j)
                  mj = 1
              else
                  mj = m - j + 1
              end if
              ! h(j:n, j) := a(j, j:n) - h(j:n, 1:(j-1)) * l(j1:(j-1), j),
               ! where h(j:n, j) has been initialized to be a(j, j:n)
              if (k > 2) then
              ! k is the column to be factorized
               ! > for the first block column, k is j, skipping the first two
                 ! columns
               ! > for the rest of the columns, k is j+1, skipping only the
                 ! first column
                 call la_clacgv(j - k1,a(1,j),1)
                 call la_cgemv('NO TRANSPOSE',mj,j - k1,-cone,h(j,k1),ldh,a(1,j),1, &
                           cone,h(j,j),1)
                 call la_clacgv(j - k1,a(1,j),1)
              end if
              ! copy h(i:n, i) into work
              call la_ccopy(mj,h(j,j),1,work(1),1)
              if (j > k1) then
                 ! compute work := work - l(j-1, j:n) * t(j-1,j),
                  ! where a(j-1, j) stores t(j-1, j) and a(j-2, j:n) stores u(j-1, j:n)
                 alpha = -conjg(a(k - 1,j))
                 call la_caxpy(mj,alpha,a(k - 2,j),lda,work(1),1)
              end if
              ! set a(j, j) = t(j, j)
              a(k,j) = real(work(1),KIND=sp)
              if (j < m) then
                 ! compute work(2:n) = t(j, j) l(j, (j+1):n)
                  ! where a(j, j) stores t(j, j) and a(j-1, (j+1):n) stores u(j, (j+1):n)
                 if (k > 1) then
                    alpha = -a(k,j)
                    call la_caxpy(m - j,alpha,a(k - 1,j + 1),lda,work(2),1)
                 end if
                 ! find max(|work(2:n)|)
                 i2 = la_icamax(m - j,work(2),1) + 1
                 piv = work(i2)
                 ! apply hermitian pivot
                 if ((i2 /= 2) .and. (piv /= 0)) then
                    ! swap work(i1) and work(i2)
                    i1 = 2
                    work(i2) = work(i1)
                    work(i1) = piv
                    ! swap a(i1, i1+1:n) with a(i1+1:n, i2)
                    i1 = i1 + j - 1
                    i2 = i2 + j - 1
                    call la_cswap(i2 - i1 - 1,a(j1 + i1 - 1,i1 + 1),lda,a(j1 + i1,i2),1)

                    call la_clacgv(i2 - i1,a(j1 + i1 - 1,i1 + 1),lda)
                    call la_clacgv(i2 - i1 - 1,a(j1 + i1,i2),1)
                    ! swap a(i1, i2+1:n) with a(i2, i2+1:n)
                    if (i2 < m) call la_cswap(m - i2,a(j1 + i1 - 1,i2 + 1),lda,a(j1 + i2 - 1,i2 + 1), &
                               lda)
                    ! swap a(i1, i1) with a(i2,i2)
                    piv = a(i1 + j1 - 1,i1)
                    a(j1 + i1 - 1,i1) = a(j1 + i2 - 1,i2)
                    a(j1 + i2 - 1,i2) = piv
                    ! swap h(i1, 1:j1) with h(i2, 1:j1)
                    call la_cswap(i1 - 1,h(i1,1),ldh,h(i2,1),ldh)
                    ipiv(i1) = i2
                    if (i1 > (k1 - 1)) then
                       ! swap l(1:i1-1, i1) with l(1:i1-1, i2),
                        ! skipping the first column
                       call la_cswap(i1 - k1 + 1,a(1,i1),1,a(1,i2),1)
                    end if
                 else
                    ipiv(j + 1) = j + 1
                 end if
                 ! set a(j, j+1) = t(j, j+1)
                 a(k,j + 1) = work(2)
                 if (j < nb) then
                    ! copy a(j+1:n, j+1) into h(j:n, j),
                    call la_ccopy(m - j,a(k + 1,j + 1),lda,h(j + 1,j + 1),1)
                 end if
                 ! compute l(j+2, j+1) = work( 3:n ) / t(j, j+1),
                  ! where a(j, j+1) = t(j, j+1) and a(j+2:n, j) = l(j+2:n, j+1)
                 if (j < (m - 1)) then
                    if (a(k,j + 1) /= czero) then
                       alpha = cone/a(k,j + 1)
                       call la_ccopy(m - j - 1,work(3),1,a(k,j + 2),lda)
                       call la_cscal(m - j - 1,alpha,a(k,j + 2),lda)
                    else
                       call la_claset('FULL',1,m - j - 1,czero,czero,a(k,j + 2),lda)

                    end if
                 end if
              end if
              j = j + 1
              go to 10
              20 continue
           else
              ! .....................................................
              ! factorize a as l*d*l**t using the lower triangle of a
              ! .....................................................
              30 continue
              if (j > min(m,nb)) go to 40
              ! k is the column to be factorized
               ! when being called from la_chetrf_aa,
               ! > for the first block column, j1 is 1, hence j1+j-1 is j,
               ! > for the rest of the columns, j1 is 2, and j1+j-1 is j+1,
              k = j1 + j - 1
              if (j == m) then
                  ! only need to compute t(j, j)
                  mj = 1
              else
                  mj = m - j + 1
              end if
              ! h(j:n, j) := a(j:n, j) - h(j:n, 1:(j-1)) * l(j, j1:(j-1))^t,
               ! where h(j:n, j) has been initialized to be a(j:n, j)
              if (k > 2) then
              ! k is the column to be factorized
               ! > for the first block column, k is j, skipping the first two
                 ! columns
               ! > for the rest of the columns, k is j+1, skipping only the
                 ! first column
                 call la_clacgv(j - k1,a(j,1),lda)
                 call la_cgemv('NO TRANSPOSE',mj,j - k1,-cone,h(j,k1),ldh,a(j,1), &
                           lda,cone,h(j,j),1)
                 call la_clacgv(j - k1,a(j,1),lda)
              end if
              ! copy h(j:n, j) into work
              call la_ccopy(mj,h(j,j),1,work(1),1)
              if (j > k1) then
                 ! compute work := work - l(j:n, j-1) * t(j-1,j),
                  ! where a(j-1, j) = t(j-1, j) and a(j, j-2) = l(j, j-1)
                 alpha = -conjg(a(j,k - 1))
                 call la_caxpy(mj,alpha,a(j,k - 2),1,work(1),1)
              end if
              ! set a(j, j) = t(j, j)
              a(j,k) = real(work(1),KIND=sp)
              if (j < m) then
                 ! compute work(2:n) = t(j, j) l((j+1):n, j)
                  ! where a(j, j) = t(j, j) and a((j+1):n, j-1) = l((j+1):n, j)
                 if (k > 1) then
                    alpha = -a(j,k)
                    call la_caxpy(m - j,alpha,a(j + 1,k - 1),1,work(2),1)
                 end if
                 ! find max(|work(2:n)|)
                 i2 = la_icamax(m - j,work(2),1) + 1
                 piv = work(i2)
                 ! apply hermitian pivot
                 if ((i2 /= 2) .and. (piv /= 0)) then
                    ! swap work(i1) and work(i2)
                    i1 = 2
                    work(i2) = work(i1)
                    work(i1) = piv
                    ! swap a(i1+1:n, i1) with a(i2, i1+1:n)
                    i1 = i1 + j - 1
                    i2 = i2 + j - 1
                    call la_cswap(i2 - i1 - 1,a(i1 + 1,j1 + i1 - 1),1,a(i2,j1 + i1),lda)

                    call la_clacgv(i2 - i1,a(i1 + 1,j1 + i1 - 1),1)
                    call la_clacgv(i2 - i1 - 1,a(i2,j1 + i1),lda)
                    ! swap a(i2+1:n, i1) with a(i2+1:n, i2)
                    if (i2 < m) call la_cswap(m - i2,a(i2 + 1,j1 + i1 - 1),1,a(i2 + 1,j1 + i2 - 1), &
                              1)
                    ! swap a(i1, i1) with a(i2, i2)
                    piv = a(i1,j1 + i1 - 1)
                    a(i1,j1 + i1 - 1) = a(i2,j1 + i2 - 1)
                    a(i2,j1 + i2 - 1) = piv
                    ! swap h(i1, i1:j1) with h(i2, i2:j1)
                    call la_cswap(i1 - 1,h(i1,1),ldh,h(i2,1),ldh)
                    ipiv(i1) = i2
                    if (i1 > (k1 - 1)) then
                       ! swap l(1:i1-1, i1) with l(1:i1-1, i2),
                        ! skipping the first column
                       call la_cswap(i1 - k1 + 1,a(i1,1),lda,a(i2,1),lda)
                    end if
                 else
                    ipiv(j + 1) = j + 1
                 end if
                 ! set a(j+1, j) = t(j+1, j)
                 a(j + 1,k) = work(2)
                 if (j < nb) then
                    ! copy a(j+1:n, j+1) into h(j+1:n, j),
                    call la_ccopy(m - j,a(j + 1,k + 1),1,h(j + 1,j + 1),1)
                 end if
                 ! compute l(j+2, j+1) = work( 3:n ) / t(j, j+1),
                  ! where a(j, j+1) = t(j, j+1) and a(j+2:n, j) = l(j+2:n, j+1)
                 if (j < (m - 1)) then
                    if (a(j + 1,k) /= czero) then
                       alpha = cone/a(j + 1,k)
                       call la_ccopy(m - j - 1,work(3),1,a(j + 2,k),1)
                       call la_cscal(m - j - 1,alpha,a(j + 2,k),1)
                    else
                       call la_claset('FULL',m - j - 1,1,czero,czero,a(j + 2,k),lda)

                    end if
                 end if
              end if
              j = j + 1
              go to 30
              40 continue
           end if
           return
     end subroutine la_clahef_aa
     !> DLAHEF_AA factorizes a panel of a complex hermitian matrix A using
     !> the Aasen's algorithm. The panel consists of a set of NB rows of A
     !> when UPLO is U, or a set of NB columns when UPLO is L.
     !> In order to factorize the panel, the Aasen's algorithm requires the
     !> last row, or column, of the previous panel. The first row, or column,
     !> of A is set to be the first row, or column, of an identity matrix,
     !> which is used to factorize the first panel.
     !> The resulting J-th row of U, or J-th column of L, is stored in the
     !> (J-1)-th row, or column, of A (without the unit diagonals), while
     !> the diagonal and subdiagonal of A are overwritten by those of T.

     pure subroutine la_zlahef_aa(uplo,j1,m,nb,a,lda,ipiv,h,ldh,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: m,nb,j1,lda,ldh
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*),h(ldh,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k,k1,i1,i2,mj
           complex(dp) :: piv,alpha
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
           ! Executable Statements
           j = 1
           ! k1 is the first column of the panel to be factorized
           ! i.e.,  k1 is 2 for the first block column, and 1 for the rest of the blocks
           k1 = (2 - j1) + 1
           if (la_lsame(uplo,'U')) then
              ! .....................................................
              ! factorize a as u**t*d*u using the upper triangle of a
              ! .....................................................
              10 continue
              if (j > min(m,nb)) go to 20
              ! k is the column to be factorized
               ! when being called from la_zhetrf_aa,
               ! > for the first block column, j1 is 1, hence j1+j-1 is j,
               ! > for the rest of the columns, j1 is 2, and j1+j-1 is j+1,
              k = j1 + j - 1
              if (j == m) then
                  ! only need to compute t(j, j)
                  mj = 1
              else
                  mj = m - j + 1
              end if
              ! h(j:n, j) := a(j, j:n) - h(j:n, 1:(j-1)) * l(j1:(j-1), j),
               ! where h(j:n, j) has been initialized to be a(j, j:n)
              if (k > 2) then
              ! k is the column to be factorized
               ! > for the first block column, k is j, skipping the first two
                 ! columns
               ! > for the rest of the columns, k is j+1, skipping only the
                 ! first column
                 call la_zlacgv(j - k1,a(1,j),1)
                 call la_zgemv('NO TRANSPOSE',mj,j - k1,-cone,h(j,k1),ldh,a(1,j),1, &
                           cone,h(j,j),1)
                 call la_zlacgv(j - k1,a(1,j),1)
              end if
              ! copy h(i:n, i) into work
              call la_zcopy(mj,h(j,j),1,work(1),1)
              if (j > k1) then
                 ! compute work := work - l(j-1, j:n) * t(j-1,j),
                  ! where a(j-1, j) stores t(j-1, j) and a(j-2, j:n) stores u(j-1, j:n)
                 alpha = -conjg(a(k - 1,j))
                 call la_zaxpy(mj,alpha,a(k - 2,j),lda,work(1),1)
              end if
              ! set a(j, j) = t(j, j)
              a(k,j) = real(work(1),KIND=dp)
              if (j < m) then
                 ! compute work(2:n) = t(j, j) l(j, (j+1):n)
                  ! where a(j, j) stores t(j, j) and a(j-1, (j+1):n) stores u(j, (j+1):n)
                 if (k > 1) then
                    alpha = -a(k,j)
                    call la_zaxpy(m - j,alpha,a(k - 1,j + 1),lda,work(2),1)
                 end if
                 ! find max(|work(2:n)|)
                 i2 = la_izamax(m - j,work(2),1) + 1
                 piv = work(i2)
                 ! apply hermitian pivot
                 if ((i2 /= 2) .and. (piv /= 0)) then
                    ! swap work(i1) and work(i2)
                    i1 = 2
                    work(i2) = work(i1)
                    work(i1) = piv
                    ! swap a(i1, i1+1:n) with a(i1+1:n, i2)
                    i1 = i1 + j - 1
                    i2 = i2 + j - 1
                    call la_zswap(i2 - i1 - 1,a(j1 + i1 - 1,i1 + 1),lda,a(j1 + i1,i2),1)

                    call la_zlacgv(i2 - i1,a(j1 + i1 - 1,i1 + 1),lda)
                    call la_zlacgv(i2 - i1 - 1,a(j1 + i1,i2),1)
                    ! swap a(i1, i2+1:n) with a(i2, i2+1:n)
                    if (i2 < m) call la_zswap(m - i2,a(j1 + i1 - 1,i2 + 1),lda,a(j1 + i2 - 1,i2 + 1), &
                               lda)
                    ! swap a(i1, i1) with a(i2,i2)
                    piv = a(i1 + j1 - 1,i1)
                    a(j1 + i1 - 1,i1) = a(j1 + i2 - 1,i2)
                    a(j1 + i2 - 1,i2) = piv
                    ! swap h(i1, 1:j1) with h(i2, 1:j1)
                    call la_zswap(i1 - 1,h(i1,1),ldh,h(i2,1),ldh)
                    ipiv(i1) = i2
                    if (i1 > (k1 - 1)) then
                       ! swap l(1:i1-1, i1) with l(1:i1-1, i2),
                        ! skipping the first column
                       call la_zswap(i1 - k1 + 1,a(1,i1),1,a(1,i2),1)
                    end if
                 else
                    ipiv(j + 1) = j + 1
                 end if
                 ! set a(j, j+1) = t(j, j+1)
                 a(k,j + 1) = work(2)
                 if (j < nb) then
                    ! copy a(j+1:n, j+1) into h(j:n, j),
                    call la_zcopy(m - j,a(k + 1,j + 1),lda,h(j + 1,j + 1),1)
                 end if
                 ! compute l(j+2, j+1) = work( 3:n ) / t(j, j+1),
                  ! where a(j, j+1) = t(j, j+1) and a(j+2:n, j) = l(j+2:n, j+1)
                 if (j < (m - 1)) then
                    if (a(k,j + 1) /= czero) then
                       alpha = cone/a(k,j + 1)
                       call la_zcopy(m - j - 1,work(3),1,a(k,j + 2),lda)
                       call la_zscal(m - j - 1,alpha,a(k,j + 2),lda)
                    else
                       call la_zlaset('FULL',1,m - j - 1,czero,czero,a(k,j + 2),lda)

                    end if
                 end if
              end if
              j = j + 1
              go to 10
              20 continue
           else
              ! .....................................................
              ! factorize a as l*d*l**t using the lower triangle of a
              ! .....................................................
              30 continue
              if (j > min(m,nb)) go to 40
              ! k is the column to be factorized
               ! when being called from la_zhetrf_aa,
               ! > for the first block column, j1 is 1, hence j1+j-1 is j,
               ! > for the rest of the columns, j1 is 2, and j1+j-1 is j+1,
              k = j1 + j - 1
              if (j == m) then
                  ! only need to compute t(j, j)
                  mj = 1
              else
                  mj = m - j + 1
              end if
              ! h(j:n, j) := a(j:n, j) - h(j:n, 1:(j-1)) * l(j, j1:(j-1))^t,
               ! where h(j:n, j) has been initialized to be a(j:n, j)
              if (k > 2) then
              ! k is the column to be factorized
               ! > for the first block column, k is j, skipping the first two
                 ! columns
               ! > for the rest of the columns, k is j+1, skipping only the
                 ! first column
                 call la_zlacgv(j - k1,a(j,1),lda)
                 call la_zgemv('NO TRANSPOSE',mj,j - k1,-cone,h(j,k1),ldh,a(j,1), &
                           lda,cone,h(j,j),1)
                 call la_zlacgv(j - k1,a(j,1),lda)
              end if
              ! copy h(j:n, j) into work
              call la_zcopy(mj,h(j,j),1,work(1),1)
              if (j > k1) then
                 ! compute work := work - l(j:n, j-1) * t(j-1,j),
                  ! where a(j-1, j) = t(j-1, j) and a(j, j-2) = l(j, j-1)
                 alpha = -conjg(a(j,k - 1))
                 call la_zaxpy(mj,alpha,a(j,k - 2),1,work(1),1)
              end if
              ! set a(j, j) = t(j, j)
              a(j,k) = real(work(1),KIND=dp)
              if (j < m) then
                 ! compute work(2:n) = t(j, j) l((j+1):n, j)
                  ! where a(j, j) = t(j, j) and a((j+1):n, j-1) = l((j+1):n, j)
                 if (k > 1) then
                    alpha = -a(j,k)
                    call la_zaxpy(m - j,alpha,a(j + 1,k - 1),1,work(2),1)
                 end if
                 ! find max(|work(2:n)|)
                 i2 = la_izamax(m - j,work(2),1) + 1
                 piv = work(i2)
                 ! apply hermitian pivot
                 if ((i2 /= 2) .and. (piv /= 0)) then
                    ! swap work(i1) and work(i2)
                    i1 = 2
                    work(i2) = work(i1)
                    work(i1) = piv
                    ! swap a(i1+1:n, i1) with a(i2, i1+1:n)
                    i1 = i1 + j - 1
                    i2 = i2 + j - 1
                    call la_zswap(i2 - i1 - 1,a(i1 + 1,j1 + i1 - 1),1,a(i2,j1 + i1),lda)

                    call la_zlacgv(i2 - i1,a(i1 + 1,j1 + i1 - 1),1)
                    call la_zlacgv(i2 - i1 - 1,a(i2,j1 + i1),lda)
                    ! swap a(i2+1:n, i1) with a(i2+1:n, i2)
                    if (i2 < m) call la_zswap(m - i2,a(i2 + 1,j1 + i1 - 1),1,a(i2 + 1,j1 + i2 - 1), &
                              1)
                    ! swap a(i1, i1) with a(i2, i2)
                    piv = a(i1,j1 + i1 - 1)
                    a(i1,j1 + i1 - 1) = a(i2,j1 + i2 - 1)
                    a(i2,j1 + i2 - 1) = piv
                    ! swap h(i1, i1:j1) with h(i2, i2:j1)
                    call la_zswap(i1 - 1,h(i1,1),ldh,h(i2,1),ldh)
                    ipiv(i1) = i2
                    if (i1 > (k1 - 1)) then
                       ! swap l(1:i1-1, i1) with l(1:i1-1, i2),
                        ! skipping the first column
                       call la_zswap(i1 - k1 + 1,a(i1,1),lda,a(i2,1),lda)
                    end if
                 else
                    ipiv(j + 1) = j + 1
                 end if
                 ! set a(j+1, j) = t(j+1, j)
                 a(j + 1,k) = work(2)
                 if (j < nb) then
                    ! copy a(j+1:n, j+1) into h(j+1:n, j),
                    call la_zcopy(m - j,a(j + 1,k + 1),1,h(j + 1,j + 1),1)
                 end if
                 ! compute l(j+2, j+1) = work( 3:n ) / t(j, j+1),
                  ! where a(j, j+1) = t(j, j+1) and a(j+2:n, j) = l(j+2:n, j+1)
                 if (j < (m - 1)) then
                    if (a(j + 1,k) /= czero) then
                       alpha = cone/a(j + 1,k)
                       call la_zcopy(m - j - 1,work(3),1,a(j + 2,k),1)
                       call la_zscal(m - j - 1,alpha,a(j + 2,k),1)
                    else
                       call la_zlaset('FULL',m - j - 1,1,czero,czero,a(j + 2,k),lda)

                    end if
                 end if
              end if
              j = j + 1
              go to 30
              40 continue
           end if
           return
     end subroutine la_zlahef_aa
     !> DLAHEF_AA factorizes a panel of a complex hermitian matrix A using
     !> the Aasen's algorithm. The panel consists of a set of NB rows of A
     !> when UPLO is U, or a set of NB columns when UPLO is L.
     !> In order to factorize the panel, the Aasen's algorithm requires the
     !> last row, or column, of the previous panel. The first row, or column,
     !> of A is set to be the first row, or column, of an identity matrix,
     !> which is used to factorize the first panel.
     !> The resulting J-th row of U, or J-th column of L, is stored in the
     !> (J-1)-th row, or column, of A (without the unit diagonals), while
     !> the diagonal and subdiagonal of A are overwritten by those of T.

     pure subroutine la_wlahef_aa(uplo,j1,m,nb,a,lda,ipiv,h,ldh,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: m,nb,j1,lda,ldh
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*),h(ldh,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: j,k,k1,i1,i2,mj
           complex(qp) :: piv,alpha
           ! Intrinsic Functions
           intrinsic :: real,conjg,max
           ! Executable Statements
           j = 1
           ! k1 is the first column of the panel to be factorized
           ! i.e.,  k1 is 2 for the first block column, and 1 for the rest of the blocks
           k1 = (2 - j1) + 1
           if (la_lsame(uplo,'U')) then
              ! .....................................................
              ! factorize a as u**t*d*u using the upper triangle of a
              ! .....................................................
              10 continue
              if (j > min(m,nb)) go to 20
              ! k is the column to be factorized
               ! when being called from la_whetrf_aa,
               ! > for the first block column, j1 is 1, hence j1+j-1 is j,
               ! > for the rest of the columns, j1 is 2, and j1+j-1 is j+1,
              k = j1 + j - 1
              if (j == m) then
                  ! only need to compute t(j, j)
                  mj = 1
              else
                  mj = m - j + 1
              end if
              ! h(j:n, j) := a(j, j:n) - h(j:n, 1:(j-1)) * l(j1:(j-1), j),
               ! where h(j:n, j) has been initialized to be a(j, j:n)
              if (k > 2) then
              ! k is the column to be factorized
               ! > for the first block column, k is j, skipping the first two
                 ! columns
               ! > for the rest of the columns, k is j+1, skipping only the
                 ! first column
                 call la_wlacgv(j - k1,a(1,j),1)
                 call la_wgemv('NO TRANSPOSE',mj,j - k1,-cone,h(j,k1),ldh,a(1,j),1, &
                           cone,h(j,j),1)
                 call la_wlacgv(j - k1,a(1,j),1)
              end if
              ! copy h(i:n, i) into work
              call la_wcopy(mj,h(j,j),1,work(1),1)
              if (j > k1) then
                 ! compute work := work - l(j-1, j:n) * t(j-1,j),
                  ! where a(j-1, j) stores t(j-1, j) and a(j-2, j:n) stores u(j-1, j:n)
                 alpha = -conjg(a(k - 1,j))
                 call la_waxpy(mj,alpha,a(k - 2,j),lda,work(1),1)
              end if
              ! set a(j, j) = t(j, j)
              a(k,j) = real(work(1),KIND=qp)
              if (j < m) then
                 ! compute work(2:n) = t(j, j) l(j, (j+1):n)
                  ! where a(j, j) stores t(j, j) and a(j-1, (j+1):n) stores u(j, (j+1):n)
                 if (k > 1) then
                    alpha = -a(k,j)
                    call la_waxpy(m - j,alpha,a(k - 1,j + 1),lda,work(2),1)
                 end if
                 ! find max(|work(2:n)|)
                 i2 = la_iwamax(m - j,work(2),1) + 1
                 piv = work(i2)
                 ! apply hermitian pivot
                 if ((i2 /= 2) .and. (piv /= 0)) then
                    ! swap work(i1) and work(i2)
                    i1 = 2
                    work(i2) = work(i1)
                    work(i1) = piv
                    ! swap a(i1, i1+1:n) with a(i1+1:n, i2)
                    i1 = i1 + j - 1
                    i2 = i2 + j - 1
                    call la_wswap(i2 - i1 - 1,a(j1 + i1 - 1,i1 + 1),lda,a(j1 + i1,i2),1)

                    call la_wlacgv(i2 - i1,a(j1 + i1 - 1,i1 + 1),lda)
                    call la_wlacgv(i2 - i1 - 1,a(j1 + i1,i2),1)
                    ! swap a(i1, i2+1:n) with a(i2, i2+1:n)
                    if (i2 < m) call la_wswap(m - i2,a(j1 + i1 - 1,i2 + 1),lda,a(j1 + i2 - 1,i2 + 1), &
                               lda)
                    ! swap a(i1, i1) with a(i2,i2)
                    piv = a(i1 + j1 - 1,i1)
                    a(j1 + i1 - 1,i1) = a(j1 + i2 - 1,i2)
                    a(j1 + i2 - 1,i2) = piv
                    ! swap h(i1, 1:j1) with h(i2, 1:j1)
                    call la_wswap(i1 - 1,h(i1,1),ldh,h(i2,1),ldh)
                    ipiv(i1) = i2
                    if (i1 > (k1 - 1)) then
                       ! swap l(1:i1-1, i1) with l(1:i1-1, i2),
                        ! skipping the first column
                       call la_wswap(i1 - k1 + 1,a(1,i1),1,a(1,i2),1)
                    end if
                 else
                    ipiv(j + 1) = j + 1
                 end if
                 ! set a(j, j+1) = t(j, j+1)
                 a(k,j + 1) = work(2)
                 if (j < nb) then
                    ! copy a(j+1:n, j+1) into h(j:n, j),
                    call la_wcopy(m - j,a(k + 1,j + 1),lda,h(j + 1,j + 1),1)
                 end if
                 ! compute l(j+2, j+1) = work( 3:n ) / t(j, j+1),
                  ! where a(j, j+1) = t(j, j+1) and a(j+2:n, j) = l(j+2:n, j+1)
                 if (j < (m - 1)) then
                    if (a(k,j + 1) /= czero) then
                       alpha = cone/a(k,j + 1)
                       call la_wcopy(m - j - 1,work(3),1,a(k,j + 2),lda)
                       call la_wscal(m - j - 1,alpha,a(k,j + 2),lda)
                    else
                       call la_wlaset('FULL',1,m - j - 1,czero,czero,a(k,j + 2),lda)

                    end if
                 end if
              end if
              j = j + 1
              go to 10
              20 continue
           else
              ! .....................................................
              ! factorize a as l*d*l**t using the lower triangle of a
              ! .....................................................
              30 continue
              if (j > min(m,nb)) go to 40
              ! k is the column to be factorized
               ! when being called from la_whetrf_aa,
               ! > for the first block column, j1 is 1, hence j1+j-1 is j,
               ! > for the rest of the columns, j1 is 2, and j1+j-1 is j+1,
              k = j1 + j - 1
              if (j == m) then
                  ! only need to compute t(j, j)
                  mj = 1
              else
                  mj = m - j + 1
              end if
              ! h(j:n, j) := a(j:n, j) - h(j:n, 1:(j-1)) * l(j, j1:(j-1))^t,
               ! where h(j:n, j) has been initialized to be a(j:n, j)
              if (k > 2) then
              ! k is the column to be factorized
               ! > for the first block column, k is j, skipping the first two
                 ! columns
               ! > for the rest of the columns, k is j+1, skipping only the
                 ! first column
                 call la_wlacgv(j - k1,a(j,1),lda)
                 call la_wgemv('NO TRANSPOSE',mj,j - k1,-cone,h(j,k1),ldh,a(j,1), &
                           lda,cone,h(j,j),1)
                 call la_wlacgv(j - k1,a(j,1),lda)
              end if
              ! copy h(j:n, j) into work
              call la_wcopy(mj,h(j,j),1,work(1),1)
              if (j > k1) then
                 ! compute work := work - l(j:n, j-1) * t(j-1,j),
                  ! where a(j-1, j) = t(j-1, j) and a(j, j-2) = l(j, j-1)
                 alpha = -conjg(a(j,k - 1))
                 call la_waxpy(mj,alpha,a(j,k - 2),1,work(1),1)
              end if
              ! set a(j, j) = t(j, j)
              a(j,k) = real(work(1),KIND=qp)
              if (j < m) then
                 ! compute work(2:n) = t(j, j) l((j+1):n, j)
                  ! where a(j, j) = t(j, j) and a((j+1):n, j-1) = l((j+1):n, j)
                 if (k > 1) then
                    alpha = -a(j,k)
                    call la_waxpy(m - j,alpha,a(j + 1,k - 1),1,work(2),1)
                 end if
                 ! find max(|work(2:n)|)
                 i2 = la_iwamax(m - j,work(2),1) + 1
                 piv = work(i2)
                 ! apply hermitian pivot
                 if ((i2 /= 2) .and. (piv /= 0)) then
                    ! swap work(i1) and work(i2)
                    i1 = 2
                    work(i2) = work(i1)
                    work(i1) = piv
                    ! swap a(i1+1:n, i1) with a(i2, i1+1:n)
                    i1 = i1 + j - 1
                    i2 = i2 + j - 1
                    call la_wswap(i2 - i1 - 1,a(i1 + 1,j1 + i1 - 1),1,a(i2,j1 + i1),lda)

                    call la_wlacgv(i2 - i1,a(i1 + 1,j1 + i1 - 1),1)
                    call la_wlacgv(i2 - i1 - 1,a(i2,j1 + i1),lda)
                    ! swap a(i2+1:n, i1) with a(i2+1:n, i2)
                    if (i2 < m) call la_wswap(m - i2,a(i2 + 1,j1 + i1 - 1),1,a(i2 + 1,j1 + i2 - 1), &
                              1)
                    ! swap a(i1, i1) with a(i2, i2)
                    piv = a(i1,j1 + i1 - 1)
                    a(i1,j1 + i1 - 1) = a(i2,j1 + i2 - 1)
                    a(i2,j1 + i2 - 1) = piv
                    ! swap h(i1, i1:j1) with h(i2, i2:j1)
                    call la_wswap(i1 - 1,h(i1,1),ldh,h(i2,1),ldh)
                    ipiv(i1) = i2
                    if (i1 > (k1 - 1)) then
                       ! swap l(1:i1-1, i1) with l(1:i1-1, i2),
                        ! skipping the first column
                       call la_wswap(i1 - k1 + 1,a(i1,1),lda,a(i2,1),lda)
                    end if
                 else
                    ipiv(j + 1) = j + 1
                 end if
                 ! set a(j+1, j) = t(j+1, j)
                 a(j + 1,k) = work(2)
                 if (j < nb) then
                    ! copy a(j+1:n, j+1) into h(j+1:n, j),
                    call la_wcopy(m - j,a(j + 1,k + 1),1,h(j + 1,j + 1),1)
                 end if
                 ! compute l(j+2, j+1) = work( 3:n ) / t(j, j+1),
                  ! where a(j, j+1) = t(j, j+1) and a(j+2:n, j) = l(j+2:n, j+1)
                 if (j < (m - 1)) then
                    if (a(j + 1,k) /= czero) then
                       alpha = cone/a(j + 1,k)
                       call la_wcopy(m - j - 1,work(3),1,a(j + 2,k),1)
                       call la_wscal(m - j - 1,alpha,a(j + 2,k),1)
                    else
                       call la_wlaset('FULL',m - j - 1,1,czero,czero,a(j + 2,k),lda)

                    end if
                 end if
              end if
              j = j + 1
              go to 30
              40 continue
           end if
           return
     end subroutine la_wlahef_aa

end module la_lapack_solve_ldl_comp4
