!> Symmetric indefinite components: Bunch-Kaufman factorization, solve, inverse
module la_lapack_solve_ldl_comp
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_blas_level2_pac
     use la_blas_level2_sym
     use la_blas_level3_gen
     use la_blas_level3_tri
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_l2
     use la_lapack_blas_like_scalar
     use la_lapack_solve_aux
     use la_lapack_solve_ldl_comp2
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slasyf
     public :: la_ssptrf
     public :: la_ssyconv
     public :: la_ssyequb
     public :: la_ssyswapr
     public :: la_ssytri
     public :: la_ssytrs
     public :: la_ssytrs2
     public :: la_ssytrs_3
     public :: la_sspcon
     public :: la_ssycon
     public :: la_ssyrfs
     public :: la_ssytf2
     public :: la_ssytrf
     public :: la_dlasyf
     public :: la_dsptrf
     public :: la_dsyconv
     public :: la_dsyequb
     public :: la_dsyswapr
     public :: la_dsytri
     public :: la_dsytrs
     public :: la_dsytrs2
     public :: la_dsytrs_3
     public :: la_dspcon
     public :: la_dsycon
     public :: la_dsyrfs
     public :: la_dsytf2
     public :: la_dsytrf
#ifdef LA_WITH_XDP
     public :: la_xlasyf
     public :: la_xsptrf
     public :: la_xsyconv
     public :: la_xsyequb
     public :: la_xsyswapr
     public :: la_xsytri
     public :: la_xsytrs
     public :: la_xsytrs2
     public :: la_xsytrs_3
     public :: la_xspcon
     public :: la_xsycon
     public :: la_xsyrfs
     public :: la_xsytf2
     public :: la_xsytrf
#endif
#ifdef LA_WITH_QP
     public :: la_qlasyf
     public :: la_qsptrf
     public :: la_qsyconv
     public :: la_qsyequb
     public :: la_qsyswapr
     public :: la_qsytri
     public :: la_qsytrs
     public :: la_qsytrs2
     public :: la_qsytrs_3
     public :: la_qspcon
     public :: la_qsycon
     public :: la_qsyrfs
     public :: la_qsytf2
     public :: la_qsytrf
#endif
     public :: la_clasyf
     public :: la_csptrf
     public :: la_csyconv
     public :: la_csyequb
     public :: la_csyswapr
     public :: la_csytf2
     public :: la_csytrf
     public :: la_csytri
     public :: la_csytrs
     public :: la_csytrs2
     public :: la_csytrs_3
     public :: la_cla_herpvgrw
     public :: la_cspcon
     public :: la_csycon
     public :: la_csyrfs
     public :: la_zlasyf
     public :: la_zsptrf
     public :: la_zsyconv
     public :: la_zsyequb
     public :: la_zsyswapr
     public :: la_zsytf2
     public :: la_zsytrf
     public :: la_zsytri
     public :: la_zsytrs
     public :: la_zsytrs2
     public :: la_zsytrs_3
     public :: la_zla_herpvgrw
     public :: la_zspcon
     public :: la_zsycon
     public :: la_zsyrfs
#ifdef LA_WITH_XDP
     public :: la_ylasyf
     public :: la_ysptrf
     public :: la_ysyconv
     public :: la_ysyequb
     public :: la_ysyswapr
     public :: la_ysytf2
     public :: la_ysytrf
     public :: la_ysytri
     public :: la_ysytrs
     public :: la_ysytrs2
     public :: la_ysytrs_3
     public :: la_yla_herpvgrw
     public :: la_yspcon
     public :: la_ysycon
     public :: la_ysyrfs
#endif
#ifdef LA_WITH_QP
     public :: la_wlasyf
     public :: la_wsptrf
     public :: la_wsyconv
     public :: la_wsyequb
     public :: la_wsyswapr
     public :: la_wsytf2
     public :: la_wsytrf
     public :: la_wsytri
     public :: la_wsytrs
     public :: la_wsytrs2
     public :: la_wsytrs_3
     public :: la_wla_herpvgrw
     public :: la_wspcon
     public :: la_wsycon
     public :: la_wsyrfs
#endif

     contains

     !> SLASYF: computes a partial factorization of a real symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> SLASYF is an auxiliary routine called by SSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_slasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
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
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(sp) :: absakk,alpha,colmax,d11,d21,d22,r1,rowmax,t
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_scopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_sgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w(k,kw + &
                        1),ldw,one,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_isamax(k - 1,w(1,kw),1)
                 colmax = abs(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_scopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_scopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_sgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w( &
                              imax,kw + 1),ldw,one,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_isamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = abs(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_isamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,abs(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_scopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_scopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    if (kp > 1) call la_scopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_sswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_sswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_scopy(k,w(1,kw),1,a(1,k),1)
                    r1 = one/a(k,k)
                    call la_sscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_sgemv('NO TRANSPOSE',jj - j + 1,n - k,-one,a(j,k + 1),lda,w(jj, &
                              kw + 1),ldw,one,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_sgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-one,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,one,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_sswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_scopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_sgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(k,1),ldw, &
                        one,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_isamax(n - k,w(k + 1,k),1)
                 colmax = abs(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_scopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_scopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_sgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(imax, &
                              1),ldw,one,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_isamax(imax - k,w(k,k + 1),1)
                    rowmax = abs(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_isamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,abs(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_scopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_scopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_scopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_sswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_sswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_scopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = one/a(k,k)
                       call la_sscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_sgemv('NO TRANSPOSE',j + jb - jj,k - 1,-one,a(jj,1),lda,w(jj, &
                              1),ldw,one,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_sgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           one,a(j + jb,1),lda,w(j,1),ldw,one,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_sswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_slasyf
     !> DLASYF: computes a partial factorization of a real symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_dlasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
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
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(dp) :: absakk,alpha,colmax,d11,d21,d22,r1,rowmax,t
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_dcopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_dgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w(k,kw + &
                        1),ldw,one,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_idamax(k - 1,w(1,kw),1)
                 colmax = abs(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_dcopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_dcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_dgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w( &
                              imax,kw + 1),ldw,one,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_idamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = abs(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_idamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,abs(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_dcopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_dcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    if (kp > 1) call la_dcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_dswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_dswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_dcopy(k,w(1,kw),1,a(1,k),1)
                    r1 = one/a(k,k)
                    call la_dscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_dgemv('NO TRANSPOSE',jj - j + 1,n - k,-one,a(j,k + 1),lda,w(jj, &
                              kw + 1),ldw,one,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_dgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-one,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,one,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_dswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_dcopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_dgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(k,1),ldw, &
                        one,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_idamax(n - k,w(k + 1,k),1)
                 colmax = abs(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_dcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_dcopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_dgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(imax, &
                              1),ldw,one,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_idamax(imax - k,w(k,k + 1),1)
                    rowmax = abs(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_idamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,abs(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_dcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_dcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_dcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_dswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_dswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_dcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = one/a(k,k)
                       call la_dscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_dgemv('NO TRANSPOSE',j + jb - jj,k - 1,-one,a(jj,1),lda,w(jj, &
                              1),ldw,one,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_dgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           one,a(j + jb,1),lda,w(j,1),ldw,one,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_dswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_dlasyf
#ifdef LA_WITH_XDP
     !> XLASYF: computes a partial factorization of a real symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> XLASYF is an auxiliary routine called by XSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_xlasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sevten = 17.0e+0_xdp

           ! Local Scalars
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(xdp) :: absakk,alpha,colmax,d11,d21,d22,r1,rowmax,t
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_xcopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_xgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w(k,kw + &
                        1),ldw,one,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_ixamax(k - 1,w(1,kw),1)
                 colmax = abs(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_xcopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_xcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_xgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w( &
                              imax,kw + 1),ldw,one,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_ixamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = abs(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_ixamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,abs(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_xcopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_xcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    if (kp > 1) call la_xcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_xswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_xswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_xcopy(k,w(1,kw),1,a(1,k),1)
                    r1 = one/a(k,k)
                    call la_xscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_xgemv('NO TRANSPOSE',jj - j + 1,n - k,-one,a(j,k + 1),lda,w(jj, &
                              kw + 1),ldw,one,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_xgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-one,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,one,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_xswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_xcopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_xgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(k,1),ldw, &
                        one,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_ixamax(n - k,w(k + 1,k),1)
                 colmax = abs(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_xcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_xcopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_xgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(imax, &
                              1),ldw,one,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_ixamax(imax - k,w(k,k + 1),1)
                    rowmax = abs(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_ixamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,abs(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_xcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_xcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_xcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_xswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_xswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_xcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = one/a(k,k)
                       call la_xscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_xgemv('NO TRANSPOSE',j + jb - jj,k - 1,-one,a(jj,1),lda,w(jj, &
                              1),ldw,one,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_xgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           one,a(j + jb,1),lda,w(j,1),ldw,one,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_xswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_xlasyf
#endif
#ifdef LA_WITH_QP
     !> QLASYF: computes a partial factorization of a real symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) (  D   0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) (  0  A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> QLASYF is an auxiliary routine called by QSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_qlasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
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
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(qp) :: absakk,alpha,colmax,d11,d21,d22,r1,rowmax,t
           ! Intrinsic Functions
           intrinsic :: abs,max,min,sqrt
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_qcopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_qgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w(k,kw + &
                        1),ldw,one,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iqamax(k - 1,w(1,kw),1)
                 colmax = abs(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_qcopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_qcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_qgemv('NO TRANSPOSE',k,n - k,-one,a(1,k + 1),lda,w( &
                              imax,kw + 1),ldw,one,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_iqamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = abs(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_iqamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,abs(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_qcopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_qcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    if (kp > 1) call la_qcopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_qswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_qswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_qcopy(k,w(1,kw),1,a(1,k),1)
                    r1 = one/a(k,k)
                    call la_qscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_qgemv('NO TRANSPOSE',jj - j + 1,n - k,-one,a(j,k + 1),lda,w(jj, &
                              kw + 1),ldw,one,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_qgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-one,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,one,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_qswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_qcopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_qgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(k,1),ldw, &
                        one,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iqamax(n - k,w(k + 1,k),1)
                 colmax = abs(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_qcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_qcopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_qgemv('NO TRANSPOSE',n - k + 1,k - 1,-one,a(k,1),lda,w(imax, &
                              1),ldw,one,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_iqamax(imax - k,w(k,k + 1),1)
                    rowmax = abs(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_iqamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,abs(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_qcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_qcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_qcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_qswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_qswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_qcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = one/a(k,k)
                       call la_qscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_qgemv('NO TRANSPOSE',j + jb - jj,k - 1,-one,a(jj,1),lda,w(jj, &
                              1),ldw,one,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_qgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           one,a(j + jb,1),lda,w(j,1),ldw,one,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_qswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_qlasyf
#endif

     !> SSPTRF: computes the factorization of a real symmetric matrix A stored
     !> in packed format using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_ssptrf(uplo,n,ap,ipiv,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(sp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('SSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_isamax(k - 1,ap(kc),1)
                 colmax = abs(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_isamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_sswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/ap(kc + k - 1)
                    call la_sspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_sscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_isamax(n - k,ap(kc + 1),1)
                 colmax = abs(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_isamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_sswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = one/ap(kc)
                       call la_sspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_sscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_ssptrf
     !> DSPTRF: computes the factorization of a real symmetric matrix A stored
     !> in packed format using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_dsptrf(uplo,n,ap,ipiv,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(dp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('DSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_idamax(k - 1,ap(kc),1)
                 colmax = abs(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_idamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_dswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/ap(kc + k - 1)
                    call la_dspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_dscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_idamax(n - k,ap(kc + 1),1)
                 colmax = abs(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_idamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_dswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = one/ap(kc)
                       call la_dspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_dscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_dsptrf
#ifdef LA_WITH_XDP
     !> XSPTRF: computes the factorization of a real symmetric matrix A stored
     !> in packed format using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_xsptrf(uplo,n,ap,ipiv,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sevten = 17.0e+0_xdp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(xdp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('XSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_ixamax(k - 1,ap(kc),1)
                 colmax = abs(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_ixamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_xswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/ap(kc + k - 1)
                    call la_xspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_xscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_ixamax(n - k,ap(kc + 1),1)
                 colmax = abs(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_ixamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_xswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = one/ap(kc)
                       call la_xspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_xscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_xsptrf
#endif
#ifdef LA_WITH_QP
     !> QSPTRF: computes the factorization of a real symmetric matrix A stored
     !> in packed format using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_qsptrf(uplo,n,ap,ipiv,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(qp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('QSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_iqamax(k - 1,ap(kc),1)
                 colmax = abs(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_iqamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_qswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/ap(kc + k - 1)
                    call la_qspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_qscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_iqamax(n - k,ap(kc + 1),1)
                 colmax = abs(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (abs(ap(kx)) > rowmax) then
                          rowmax = abs(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_iqamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,abs(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_qswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = one/ap(kc)
                       call la_qspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_qscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_qsptrf
#endif

     !> SSYCONV: convert A given by TRF into L and D and vice-versa.
     !> Get Non-diag elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_ssyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           real(sp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('SSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
            ! a is upper
            ! convert a (a is upper)
              ! convert value
              if (convert) then
                 i = n
                 e(1) = zero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = zero
                       a(i - 1,i) = zero
                       i = i - 1
                    else
                       e(i) = zero
                    end if
                    i = i - 1
                 end do
              ! convert permutations
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                    end if
                 else
                   ip = -ipiv(i)
                    if (i < n) then
                  do j = i + 1,n
                      temp = a(ip,j)
                      a(ip,j) = a(i - 1,j)
                      a(i - 1,j) = temp
                  end do
                     end if
                     i = i - 1
                end if
                i = i - 1
             end do
              else
            ! revert a (a is upper)
              ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
              ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
            ! a is lower
              if (convert) then
            ! convert a (a is lower)
              ! convert value
                 i = 1
                 e(n) = zero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = zero
                       a(i + 1,i) = zero
                       i = i + 1
                    else
                       e(i) = zero
                    end if
                    i = i + 1
                 end do
              ! convert permutations
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i > 1) then
                    do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i,j)
                      a(i,j) = temp
                    end do
                    end if
                 else
                   ip = -ipiv(i)
                   if (i > 1) then
                   do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i + 1,j)
                      a(i + 1,j) = temp
                   end do
                   end if
                   i = i + 1
                end if
                i = i + 1
             end do
              else
            ! revert a (a is lower)
              ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
              ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_ssyconv
     !> DSYCONV: convert A given by TRF into L and D and vice-versa.
     !> Get Non-diag elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_dsyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           real(dp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
            ! a is upper
            ! convert a (a is upper)
              ! convert value
              if (convert) then
                 i = n
                 e(1) = zero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = zero
                       a(i - 1,i) = zero
                       i = i - 1
                    else
                       e(i) = zero
                    end if
                    i = i - 1
                 end do
              ! convert permutations
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                    end if
                 else
                   ip = -ipiv(i)
                    if (i < n) then
                  do j = i + 1,n
                      temp = a(ip,j)
                      a(ip,j) = a(i - 1,j)
                      a(i - 1,j) = temp
                  end do
                     end if
                     i = i - 1
                end if
                i = i - 1
             end do
              else
            ! revert a (a is upper)
              ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
              ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
            ! a is lower
              if (convert) then
            ! convert a (a is lower)
              ! convert value
                 i = 1
                 e(n) = zero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = zero
                       a(i + 1,i) = zero
                       i = i + 1
                    else
                       e(i) = zero
                    end if
                    i = i + 1
                 end do
              ! convert permutations
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i > 1) then
                    do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i,j)
                      a(i,j) = temp
                    end do
                    end if
                 else
                   ip = -ipiv(i)
                   if (i > 1) then
                   do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i + 1,j)
                      a(i + 1,j) = temp
                   end do
                   end if
                   i = i + 1
                end if
                i = i + 1
             end do
              else
            ! revert a (a is lower)
              ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
              ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_dsyconv
#ifdef LA_WITH_XDP
     !> XSYCONV: convert A given by TRF into L and D and vice-versa.
     !> Get Non-diag elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_xsyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           real(xdp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('XSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
            ! a is upper
            ! convert a (a is upper)
              ! convert value
              if (convert) then
                 i = n
                 e(1) = zero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = zero
                       a(i - 1,i) = zero
                       i = i - 1
                    else
                       e(i) = zero
                    end if
                    i = i - 1
                 end do
              ! convert permutations
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                    end if
                 else
                   ip = -ipiv(i)
                    if (i < n) then
                  do j = i + 1,n
                      temp = a(ip,j)
                      a(ip,j) = a(i - 1,j)
                      a(i - 1,j) = temp
                  end do
                     end if
                     i = i - 1
                end if
                i = i - 1
             end do
              else
            ! revert a (a is upper)
              ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
              ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
            ! a is lower
              if (convert) then
            ! convert a (a is lower)
              ! convert value
                 i = 1
                 e(n) = zero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = zero
                       a(i + 1,i) = zero
                       i = i + 1
                    else
                       e(i) = zero
                    end if
                    i = i + 1
                 end do
              ! convert permutations
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i > 1) then
                    do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i,j)
                      a(i,j) = temp
                    end do
                    end if
                 else
                   ip = -ipiv(i)
                   if (i > 1) then
                   do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i + 1,j)
                      a(i + 1,j) = temp
                   end do
                   end if
                   i = i + 1
                end if
                i = i + 1
             end do
              else
            ! revert a (a is lower)
              ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
              ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_xsyconv
#endif
#ifdef LA_WITH_QP
     !> QSYCONV: convert A given by TRF into L and D and vice-versa.
     !> Get Non-diag elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_qsyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           real(qp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
            ! a is upper
            ! convert a (a is upper)
              ! convert value
              if (convert) then
                 i = n
                 e(1) = zero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = zero
                       a(i - 1,i) = zero
                       i = i - 1
                    else
                       e(i) = zero
                    end if
                    i = i - 1
                 end do
              ! convert permutations
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                    end if
                 else
                   ip = -ipiv(i)
                    if (i < n) then
                  do j = i + 1,n
                      temp = a(ip,j)
                      a(ip,j) = a(i - 1,j)
                      a(i - 1,j) = temp
                  end do
                     end if
                     i = i - 1
                end if
                i = i - 1
             end do
              else
            ! revert a (a is upper)
              ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
              ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
            ! a is lower
              if (convert) then
            ! convert a (a is lower)
              ! convert value
                 i = 1
                 e(n) = zero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = zero
                       a(i + 1,i) = zero
                       i = i + 1
                    else
                       e(i) = zero
                    end if
                    i = i + 1
                 end do
              ! convert permutations
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    ip = ipiv(i)
                    if (i > 1) then
                    do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i,j)
                      a(i,j) = temp
                    end do
                    end if
                 else
                   ip = -ipiv(i)
                   if (i > 1) then
                   do j = 1,i - 1
                      temp = a(ip,j)
                      a(ip,j) = a(i + 1,j)
                      a(i + 1,j) = temp
                   end do
                   end if
                   i = i + 1
                end if
                i = i + 1
             end do
              else
            ! revert a (a is lower)
              ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
              ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_qsyconv
#endif

     !> SSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_ssyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(sp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('SSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_sp/s(j)
           end do
           tol = one/sqrt(2.0_sp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + abs(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + abs(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*work(i)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_slassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = abs(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(work(i) - t*si)
                 c0 = -(t*si)*si + 2*work(i)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + work(i))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_slamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_slamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_ssyequb
     !> DSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_dsyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(dp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('DSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_dp/s(j)
           end do
           tol = one/sqrt(2.0_dp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + abs(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + abs(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*work(i)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_dlassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = abs(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(work(i) - t*si)
                 c0 = -(t*si)*si + 2*work(i)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + work(i))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_dlamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_dlamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_dsyequb
#ifdef LA_WITH_XDP
     !> XSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_xsyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(xdp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(xdp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('XSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_xdp/s(j)
           end do
           tol = one/sqrt(2.0_xdp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + abs(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + abs(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*work(i)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_xlassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = abs(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(work(i) - t*si)
                 c0 = -(t*si)*si + 2*work(i)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + work(i))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_xlamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_xlamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_xsyequb
#endif
#ifdef LA_WITH_QP
     !> QSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_qsyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: s(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(qp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('QSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),abs(a(j,j)))
                 amax = max(amax,abs(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),abs(a(i,j)))
                    s(j) = max(s(j),abs(a(i,j)))
                    amax = max(amax,abs(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_qp/s(j)
           end do
           tol = one/sqrt(2.0_qp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + abs(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + abs(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + abs(a(i,j))*s(j)
                       work(j) = work(j) + abs(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*work(i)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_qlassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = abs(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(work(i) - t*si)
                 c0 = -(t*si)*si + 2*work(i)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = abs(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = abs(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + work(i))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_qlamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_qlamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_qsyequb
#endif

     !> SSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_ssyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(sp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_sswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_sswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_ssyswapr
     !> DSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_dsyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(dp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_dswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_dswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_dsyswapr
#ifdef LA_WITH_XDP
     !> XSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_xsyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(xdp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_xswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_xswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_xsyswapr
#endif
#ifdef LA_WITH_QP
     !> QSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_qsyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           real(qp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_qswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_qswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_qsyswapr
#endif

     !> SSYTRI: computes the inverse of a real symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> SSYTRF.

     pure subroutine la_ssytri(uplo,n,a,lda,ipiv,work,info)
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
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: k,kp,kstep
           real(sp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('SSYTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_scopy(k - 1,a(1,k),1,work,1)
                    call la_ssymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_sdot(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_scopy(k - 1,a(1,k),1,work,1)
                    call la_ssymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_sdot(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_sdot(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_scopy(k - 1,a(1,k + 1),1,work,1)
                    call la_ssymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_sdot(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_sswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_sswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_scopy(n - k,a(k + 1,k),1,work,1)
                    call la_ssymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_sdot(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_scopy(n - k,a(k + 1,k),1,work,1)
                    call la_ssymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_sdot(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_sdot(n - k,a(k + 1,k),1,a(k + 1,k - 1),1)

                    call la_scopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_ssymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_sdot(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_sswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_sswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_ssytri
     !> DSYTRI: computes the inverse of a real symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> DSYTRF.

     pure subroutine la_dsytri(uplo,n,a,lda,ipiv,work,info)
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
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: k,kp,kstep
           real(dp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('DSYTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_dcopy(k - 1,a(1,k),1,work,1)
                    call la_dsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_ddot(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_dcopy(k - 1,a(1,k),1,work,1)
                    call la_dsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_ddot(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_ddot(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_dcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_dsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_ddot(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_dswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_dswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_dcopy(n - k,a(k + 1,k),1,work,1)
                    call la_dsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_ddot(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_dcopy(n - k,a(k + 1,k),1,work,1)
                    call la_dsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_ddot(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_ddot(n - k,a(k + 1,k),1,a(k + 1,k - 1),1)

                    call la_dcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_dsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_ddot(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_dswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_dswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_dsytri
#ifdef LA_WITH_XDP
     !> XSYTRI: computes the inverse of a real symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> XSYTRF.

     pure subroutine la_xsytri(uplo,n,a,lda,ipiv,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: k,kp,kstep
           real(xdp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('XSYTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_xcopy(k - 1,a(1,k),1,work,1)
                    call la_xsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_xdot(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_xcopy(k - 1,a(1,k),1,work,1)
                    call la_xsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_xdot(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_xdot(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_xcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_xsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_xdot(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_xswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_xswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_xcopy(n - k,a(k + 1,k),1,work,1)
                    call la_xsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_xdot(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_xcopy(n - k,a(k + 1,k),1,work,1)
                    call la_xsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_xdot(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_xdot(n - k,a(k + 1,k),1,a(k + 1,k - 1),1)

                    call la_xcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_xsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_xdot(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_xswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_xswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_xsytri
#endif
#ifdef LA_WITH_QP
     !> QSYTRI: computes the inverse of a real symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> QSYTRF.

     pure subroutine la_qsytri(uplo,n,a,lda,ipiv,work,info)
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
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: k,kp,kstep
           real(qp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('QSYTRI',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! check that the diagonal matrix d is nonsingular.
           if (upper) then
              ! upper triangular storage: examine d from bottom to top
              do info = n,1,-1
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              do info = 1,n
                 if (ipiv(info) > 0 .and. a(info,info) == zero) return
              end do
           end if
           info = 0
           if (upper) then
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_qcopy(k - 1,a(1,k),1,work,1)
                    call la_qsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_qdot(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k + 1))
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - one)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_qcopy(k - 1,a(1,k),1,work,1)
                    call la_qsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k),1)

                    a(k,k) = a(k,k) - la_qdot(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_qdot(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_qcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_qsymv(uplo,k - 1,-one,a,lda,work,1,zero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_qdot(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_qswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_qswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = one/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_qcopy(n - k,a(k + 1,k),1,work,1)
                    call la_qsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_qdot(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = abs(a(k,k - 1))
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - one)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_qcopy(n - k,a(k + 1,k),1,work,1)
                    call la_qsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k),1)
                    a(k,k) = a(k,k) - la_qdot(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_qdot(n - k,a(k + 1,k),1,a(k + 1,k - 1),1)

                    call la_qcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_qsymv(uplo,n - k,-one,a(k + 1,k + 1),lda,work,1,zero,a(k + 1, &
                              k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_qdot(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_qswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_qswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_qsytri
#endif

     !> SSYTRS: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by SSYTRF.

     pure subroutine la_ssytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
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
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('SSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_sger(k - 1,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 ! multiply by the inverse of the diagonal block.
                 call la_sscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k - 1) call la_sswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_sger(k - 2,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 call la_sger(k - 2,nrhs,-one,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_sgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_sgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 call la_sgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k + 1),1,one,b( &
                           k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_sger(n - k,nrhs,-one,a(k + 1,k),1,b(k,1),ldb,b(k + &
                           1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 call la_sscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k + 1) call la_sswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_sger(n - k - 1,nrhs,-one,a(k + 2,k),1,b(k,1),ldb,b(k + 2,1 &
                              ),ldb)
                    call la_sger(n - k - 1,nrhs,-one,a(k + 2,k + 1),1,b(k + 1,1),ldb,b(k + &
                              2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_sgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + &
                           1,k),1,one,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_sgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k), &
                               1,one,b(k,1),ldb)
                    call la_sgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k - 1 &
                              ),1,one,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_ssytrs
     !> DSYTRS: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by DSYTRF.

     pure subroutine la_dsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
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
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('DSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_dger(k - 1,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 ! multiply by the inverse of the diagonal block.
                 call la_dscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k - 1) call la_dswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_dger(k - 2,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 call la_dger(k - 2,nrhs,-one,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_dgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_dgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 call la_dgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k + 1),1,one,b( &
                           k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_dger(n - k,nrhs,-one,a(k + 1,k),1,b(k,1),ldb,b(k + &
                           1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 call la_dscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k + 1) call la_dswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_dger(n - k - 1,nrhs,-one,a(k + 2,k),1,b(k,1),ldb,b(k + 2,1 &
                              ),ldb)
                    call la_dger(n - k - 1,nrhs,-one,a(k + 2,k + 1),1,b(k + 1,1),ldb,b(k + &
                              2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_dgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + &
                           1,k),1,one,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_dgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k), &
                               1,one,b(k,1),ldb)
                    call la_dgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k - 1 &
                              ),1,one,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_dsytrs
#ifdef LA_WITH_XDP
     !> XSYTRS: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by XSYTRF.

     pure subroutine la_xsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(xdp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('XSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_xger(k - 1,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 ! multiply by the inverse of the diagonal block.
                 call la_xscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k - 1) call la_xswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_xger(k - 2,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 call la_xger(k - 2,nrhs,-one,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_xgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_xgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 call la_xgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k + 1),1,one,b( &
                           k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_xger(n - k,nrhs,-one,a(k + 1,k),1,b(k,1),ldb,b(k + &
                           1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 call la_xscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k + 1) call la_xswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_xger(n - k - 1,nrhs,-one,a(k + 2,k),1,b(k,1),ldb,b(k + 2,1 &
                              ),ldb)
                    call la_xger(n - k - 1,nrhs,-one,a(k + 2,k + 1),1,b(k + 1,1),ldb,b(k + &
                              2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_xgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + &
                           1,k),1,one,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_xgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k), &
                               1,one,b(k,1),ldb)
                    call la_xgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k - 1 &
                              ),1,one,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_xsytrs
#endif
#ifdef LA_WITH_QP
     !> QSYTRS: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by QSYTRF.

     pure subroutine la_qsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
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
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           real(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('QSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_qger(k - 1,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 ! multiply by the inverse of the diagonal block.
                 call la_qscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k - 1) call la_qswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_qger(k - 2,nrhs,-one,a(1,k),1,b(k,1),ldb,b(1,1),ldb)

                 call la_qger(k - 2,nrhs,-one,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_qgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_qgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k),1,one,b(k, &
                           1),ldb)
                 call la_qgemv('TRANSPOSE',k - 1,nrhs,-one,b,ldb,a(1,k + 1),1,one,b( &
                           k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_qger(n - k,nrhs,-one,a(k + 1,k),1,b(k,1),ldb,b(k + &
                           1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 call la_qscal(nrhs,one/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k + 1) call la_qswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_qger(n - k - 1,nrhs,-one,a(k + 2,k),1,b(k,1),ldb,b(k + 2,1 &
                              ),ldb)
                    call la_qger(n - k - 1,nrhs,-one,a(k + 2,k + 1),1,b(k + 1,1),ldb,b(k + &
                              2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - one
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_qgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + &
                           1,k),1,one,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_qgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k), &
                               1,one,b(k,1),ldb)
                    call la_qgemv('TRANSPOSE',n - k,nrhs,-one,b(k + 1,1),ldb,a(k + 1,k - 1 &
                              ),1,one,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_qsytrs
#endif

     !> SSYTRS2: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by SSYTRF and converted by SSYCONV.

     pure subroutine la_ssytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
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
           real(sp),intent(inout) :: a(lda,*),b(ldb,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           real(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('SSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_ssyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_sswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_strsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_sscal(nrhs,one/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_strsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_sswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_sswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_strsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_sscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_strsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_sswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_ssyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_ssytrs2
     !> DSYTRS2: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by DSYTRF and converted by DSYCONV.

     pure subroutine la_dsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
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
           real(dp),intent(inout) :: a(lda,*),b(ldb,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           real(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('DSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_dsyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_dswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_dtrsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_dscal(nrhs,one/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_dtrsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_dswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_dswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_dtrsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_dscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_dtrsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_dswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_dsyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_dsytrs2
#ifdef LA_WITH_XDP
     !> XSYTRS2: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by XSYTRF and converted by XSYCONV.

     pure subroutine la_xsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           real(xdp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('XSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_xsyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_xswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_xtrsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_xscal(nrhs,one/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_xtrsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_xswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_xswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_xtrsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_xscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_xtrsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_xswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_xsyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_xsytrs2
#endif
#ifdef LA_WITH_QP
     !> QSYTRS2: solves a system of linear equations A*X = B with a real
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by QSYTRF and converted by QSYCONV.

     pure subroutine la_qsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
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
           real(qp),intent(inout) :: a(lda,*),b(ldb,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           real(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('QSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_qsyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_qswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_qtrsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_qscal(nrhs,one/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_qtrsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_qswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_qswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_qtrsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_qscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - one
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_qtrsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_qswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_qsyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_qsytrs2
#endif

     !> SSYTRS_3: solves a system of linear equations A * X = B with a real
     !> symmetric matrix A using the factorization computed
     !> by SSYTRF_RK or SSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_ssytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
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
           real(sp),intent(in) :: a(lda,*),e(*)
           real(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           real(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('SSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_strsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_sscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_strsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_strsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_sscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_strsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_sswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_ssytrs_3
     !> DSYTRS_3: solves a system of linear equations A * X = B with a real
     !> symmetric matrix A using the factorization computed
     !> by DSYTRF_RK or DSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_dsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
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
           real(dp),intent(in) :: a(lda,*),e(*)
           real(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           real(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('DSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_dtrsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_dscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_dtrsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_dtrsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_dscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_dtrsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_dswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_dsytrs_3
#ifdef LA_WITH_XDP
     !> XSYTRS_3: solves a system of linear equations A * X = B with a real
     !> symmetric matrix A using the factorization computed
     !> by XSYTRF_RK or DSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_xsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(xdp),intent(in) :: a(lda,*),e(*)
           real(xdp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           real(xdp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('XSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_xtrsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_xscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_xtrsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_xtrsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_xscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_xtrsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_xswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_xsytrs_3
#endif
#ifdef LA_WITH_QP
     !> QSYTRS_3: solves a system of linear equations A * X = B with a real
     !> symmetric matrix A using the factorization computed
     !> by QSYTRF_RK or DSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_qsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
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
           real(qp),intent(in) :: a(lda,*),e(*)
           real(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           real(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('QSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_qtrsm('L','U','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_qscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_qtrsm('L','U','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_qtrsm('L','L','N','U',n,nrhs,one,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_qscal(nrhs,one/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - one
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_qtrsm('L','L','T','U',n,nrhs,one,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_qswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_qsytrs_3
#endif

     !> SSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric packed matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by SSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_sspcon(uplo,n,ap,ipiv,anorm,rcond,work,iwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(sp),intent(in) :: anorm
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: ap(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(sp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('SSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_slacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_ssptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_sspcon
     !> DSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric packed matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by DSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_dspcon(uplo,n,ap,ipiv,anorm,rcond,work,iwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(dp),intent(in) :: anorm
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: ap(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(dp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('DSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_dlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_dsptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_dspcon
#ifdef LA_WITH_XDP
     !> XSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric packed matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by XSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_xspcon(uplo,n,ap,ipiv,anorm,rcond,work,iwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(xdp),intent(in) :: anorm
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: ap(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(xdp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('XSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_xlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_xsptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_xspcon
#endif
#ifdef LA_WITH_QP
     !> QSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric packed matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by QSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_qspcon(uplo,n,ap,ipiv,anorm,rcond,work,iwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(qp),intent(in) :: anorm
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: ap(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(qp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('QSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_qlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_qsptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_qspcon
#endif

     !> SSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by SSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_ssycon(uplo,n,a,lda,ipiv,anorm,rcond,work,iwork,info)
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
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*)
           real(sp),intent(out) :: work(*)
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
              call la_xerbla('SSYCON',-info)
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
           call la_slacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_ssytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_ssycon
     !> DSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by DSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_dsycon(uplo,n,a,lda,ipiv,anorm,rcond,work,iwork,info)
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
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*)
           real(dp),intent(out) :: work(*)
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
              call la_xerbla('DSYCON',-info)
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
           call la_dlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_dsytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_dsycon
#ifdef LA_WITH_XDP
     !> XSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by XSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_xsycon(uplo,n,a,lda,ipiv,anorm,rcond,work,iwork,info)
        use la_constants_xdp

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(xdp),intent(in) :: anorm
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,kase
           real(xdp) :: ainvnm
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
              call la_xerbla('XSYCON',-info)
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
           call la_xlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_xsytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_xsycon
#endif
#ifdef LA_WITH_QP
     !> QSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a real symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by QSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_qsycon(uplo,n,a,lda,ipiv,anorm,rcond,work,iwork,info)
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
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*)
           real(qp),intent(out) :: work(*)
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
              call la_xerbla('QSYCON',-info)
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
           call la_qlacn2(n,work(n + 1),work,iwork,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_qsytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_qsycon
#endif

     !> SSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_ssyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,iwork,info)
        use la_constants_sp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           real(sp),intent(out) :: berr(*),ferr(*),work(*)
           real(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
           real(sp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('SSYRFS',-info)
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
              call la_scopy(n,b(1,j),1,work(n + 1),1)
              call la_ssymv(uplo,n,-one,a,lda,x(1,j),1,one,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    do i = 1,k - 1
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + abs(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    work(k) = work(k) + abs(a(k,k))*xk
                    do i = k + 1,n
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + s
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_ssytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 call la_saxpy(n,one,work(n + 1),1,x(1,j),1)
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
              ! use la_slacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_slacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**t).
                    call la_ssytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_ssytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_ssyrfs
     !> DSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_dsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,iwork,info)
        use la_constants_dp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           real(dp),intent(out) :: berr(*),ferr(*),work(*)
           real(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
           real(dp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('DSYRFS',-info)
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
              call la_dcopy(n,b(1,j),1,work(n + 1),1)
              call la_dsymv(uplo,n,-one,a,lda,x(1,j),1,one,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    do i = 1,k - 1
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + abs(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    work(k) = work(k) + abs(a(k,k))*xk
                    do i = k + 1,n
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + s
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_dsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 call la_daxpy(n,one,work(n + 1),1,x(1,j),1)
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
              ! use la_dlacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_dlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**t).
                    call la_dsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_dsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_dsyrfs
#ifdef LA_WITH_XDP
     !> XSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_xsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,iwork,info)
        use la_constants_xdp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(xdp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           real(xdp),intent(out) :: berr(*),ferr(*),work(*)
           real(xdp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
           real(xdp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('XSYRFS',-info)
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
           eps = la_xlamch('EPSILON')
           safmin = la_xlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_140: do j = 1,nrhs
              count = 1
              lstres = three
              20 continue
              ! loop until stopping criterion is satisfied.
              ! compute residual r = b - a * x
              call la_xcopy(n,b(1,j),1,work(n + 1),1)
              call la_xsymv(uplo,n,-one,a,lda,x(1,j),1,one,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    do i = 1,k - 1
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + abs(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    work(k) = work(k) + abs(a(k,k))*xk
                    do i = k + 1,n
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + s
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_xsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 call la_xaxpy(n,one,work(n + 1),1,x(1,j),1)
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
              ! use la_xlacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_xlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**t).
                    call la_xsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_xsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_xsyrfs
#endif
#ifdef LA_WITH_QP
     !> QSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_qsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,iwork,info)
        use la_constants_qp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           real(qp),intent(out) :: berr(*),ferr(*),work(*)
           real(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
           real(qp) :: eps,lstres,s,safe1,safe2,safmin,xk
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,max
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('QSYRFS',-info)
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
              call la_qcopy(n,b(1,j),1,work(n + 1),1)
              call la_qsymv(uplo,n,-one,a,lda,x(1,j),1,one,work(n + 1),1)
              ! compute componentwise relative backward error from formula
              ! max(i) ( abs(r(i)) / ( abs(a)*abs(x) + abs(b) )(i) )
              ! where abs(z) is the componentwise absolute value of the matrix
              ! or vector z.  if the i-th component of the denominator is less
              ! than safe2, then safe1 is added to the i-th components of the
              ! numerator and denominator before dividing.
              do i = 1,n
                 work(i) = abs(b(i,j))
              end do
              ! compute abs(a)*abs(x) + abs(b).
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    do i = 1,k - 1
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + abs(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = abs(x(k,j))
                    work(k) = work(k) + abs(a(k,k))*xk
                    do i = k + 1,n
                       work(i) = work(i) + abs(a(i,k))*xk
                       s = s + abs(a(i,k))*abs(x(i,j))
                    end do
                    work(k) = work(k) + s
                 end do
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
              ! test stopping criterion. continue iterating if
                 ! 1) the residual berr(j) is larger than machine epsilon, and
                 ! 2) berr(j) decreased by at least a factor of 2 during the
                    ! last iteration, and
                 ! 3) at most itmax iterations tried.
              if (berr(j) > eps .and. two*berr(j) <= lstres .and. count <= itmax) then
                 ! update solution and try again.
                 call la_qsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 call la_qaxpy(n,one,work(n + 1),1,x(1,j),1)
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
              ! use la_qlacn2 to estimate the infinity-norm of the matrix
                 ! inv(a) * diag(w),
              ! where w = abs(r) + nz*eps*( abs(a)*abs(x)+abs(b) )))
              do i = 1,n
                 if (work(i) > safe2) then
                    work(i) = abs(work(n + i)) + nz*eps*work(i)
                 else
                    work(i) = abs(work(n + i)) + nz*eps*work(i) + safe1
                 end if
              end do
              kase = 0
              100 continue
              call la_qlacn2(n,work(2*n + 1),work(n + 1),iwork,ferr(j),kase,isave)

              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**t).
                    call la_qsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(n + i) = work(i)*work(n + i)
                    end do
                    call la_qsytrs(uplo,n,1,af,ldaf,ipiv,work(n + 1),n,info)
                 end if
                 go to 100
              end if
              ! normalize error.
              lstres = zero
              do i = 1,n
                 lstres = max(lstres,abs(x(i,j)))
              end do
              if (lstres /= zero) ferr(j) = ferr(j)/lstres
           end do loop_140
           return
     end subroutine la_qsyrfs
#endif

     !> SSYTF2: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_ssytf2(uplo,n,a,lda,ipiv,info)
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
           real(sp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(sp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
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
              call la_xerbla('SSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_isamax(k - 1,a(1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_sisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_isamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_isamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_sswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_sswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/a(k,k)
                    call la_ssyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_sscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_isamax(n - k,a(k + 1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_sisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_isamax(imax - k,a(imax,k),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_isamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_sswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_sswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       d11 = one/a(k,k)
                       call la_ssyr(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_sscal(n - k,d11,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( (a(k) a(k+1))*d(k)**(-1) ) * (a(k) a(k+1))**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_ssytf2
     !> DSYTF2: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_dsytf2(uplo,n,a,lda,ipiv,info)
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
           real(dp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(dp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
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
              call la_xerbla('DSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_idamax(k - 1,a(1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_disnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_idamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_idamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_dswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_dswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/a(k,k)
                    call la_dsyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_dscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_idamax(n - k,a(k + 1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_disnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_idamax(imax - k,a(imax,k),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_idamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_dswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_dswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       d11 = one/a(k,k)
                       call la_dsyr(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_dscal(n - k,d11,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( (a(k) a(k+1))*d(k)**(-1) ) * (a(k) a(k+1))**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_dsytf2
#ifdef LA_WITH_XDP
     !> XSYTF2: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_xsytf2(uplo,n,a,lda,ipiv,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sevten = 17.0e+0_xdp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(xdp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
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
              call la_xerbla('XSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_ixamax(k - 1,a(1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_xisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_ixamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_ixamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_xswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_xswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/a(k,k)
                    call la_xsyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_xscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_ixamax(n - k,a(k + 1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_xisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_ixamax(imax - k,a(imax,k),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_ixamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_xswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_xswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       d11 = one/a(k,k)
                       call la_xsyr(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_xscal(n - k,d11,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( (a(k) a(k+1))*d(k)**(-1) ) * (a(k) a(k+1))**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_xsytf2
#endif
#ifdef LA_WITH_QP
     !> QSYTF2: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_qsytf2(uplo,n,a,lda,ipiv,info)
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
           real(qp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(qp) :: absakk,alpha,colmax,d11,d12,d21,d22,r1,rowmax,t,wk,wkm1, &
                     wkp1
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
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
              call la_xerbla('QSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iqamax(k - 1,a(1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_qisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_iqamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_iqamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_qswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_qswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = one/a(k,k)
                    call la_qsyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_qscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = one/(d11*d22 - one)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = abs(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iqamax(n - k,a(k + 1,k),1)
                 colmax = abs(a(imax,k))
              else
                 colmax = zero
              end if
              if ((max(absakk,colmax) == zero) .or. la_qisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_iqamax(imax - k,a(imax,k),lda)
                    rowmax = abs(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_iqamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,abs(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (abs(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_qswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_qswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       d11 = one/a(k,k)
                       call la_qsyr(uplo,n - k,-d11,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_qscal(n - k,d11,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( (a(k) a(k+1))*d(k)**(-1) ) * (a(k) a(k+1))**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = one/(d11*d22 - one)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_qsytf2
#endif

     !> SSYTRF: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U**T*D*U  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_ssytrf(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: work(*)
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
              nb = la_ilaenv(1,'SSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'SSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u**t*d*u using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_slasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_slasyf(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_ssytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_slasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_slasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,ldwork, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_ssytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_ssytrf
     !> DSYTRF: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U**T*D*U  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_dsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: work(*)
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
              nb = la_ilaenv(1,'DSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'DSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u**t*d*u using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_dlasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_dlasyf(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_dsytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_dlasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_dlasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,ldwork, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_dsytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_dsytrf
#ifdef LA_WITH_XDP
     !> XSYTRF: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U**T*D*U  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_xsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: work(*)
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
              nb = la_ilaenv(1,'XSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'XSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u**t*d*u using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_xlasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_xlasyf(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_xsytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_xlasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_xlasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,ldwork, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_xsytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_xsytrf
#endif
#ifdef LA_WITH_QP
     !> QSYTRF: computes the factorization of a real symmetric matrix A using
     !> the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U**T*D*U  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_qsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: work(*)
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
              nb = la_ilaenv(1,'QSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'QSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u**t*d*u using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_qlasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_qlasyf(uplo,k,nb,kb,a,lda,ipiv,work,ldwork,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_qsytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_qlasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_qlasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,ldwork, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_qsytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_qsytrf
#endif

     !> CLASYF: computes a partial factorization of a complex symmetric matrix
     !> A using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) ( 0   A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**T denotes the transpose of U.
     !> CLASYF is an auxiliary routine called by CSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_clasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
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
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(sp) :: absakk,alpha,colmax,rowmax
           complex(sp) :: d11,d21,d22,r1,t,z
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,min,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=sp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_ccopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_cgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k, &
                        kw + 1),ldw,cone,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,kw))
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
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_ccopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_ccopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_cgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                               imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_icamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = cabs1(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_icamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_ccopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_ccopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
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
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_ccopy(k,w(1,kw),1,a(1,k),1)
                    r1 = cone/a(k,k)
                    call la_cscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = cone/(d11*d22 - cone)
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       d21 = t/d21
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_cgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_cgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_cswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_ccopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_cgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1),ldw, &
                         cone,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,k))
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
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_ccopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_ccopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_cgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(imax, &
                              1),ldw,cone,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_icamax(imax - k,w(k,k + 1),1)
                    rowmax = cabs1(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_icamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_ccopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_ccopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_ccopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_cswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_cswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_ccopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = cone/a(k,k)
                       call la_cscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_cgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_cgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_cswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_clasyf
     !> ZLASYF: computes a partial factorization of a complex symmetric matrix
     !> A using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) ( 0   A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**T denotes the transpose of U.
     !> ZLASYF is an auxiliary routine called by ZSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_zlasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
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
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(dp) :: absakk,alpha,colmax,rowmax
           complex(dp) :: d11,d21,d22,r1,t,z
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=dp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_zcopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_zgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k, &
                        kw + 1),ldw,cone,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
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
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_zcopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_zcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_zgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                               imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_izamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = cabs1(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_izamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_zcopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_zcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
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
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_zcopy(k,w(1,kw),1,a(1,k),1)
                    r1 = cone/a(k,k)
                    call la_zscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_zgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_zgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_zswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_zcopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_zgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1),ldw, &
                         cone,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
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
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_zcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_zcopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_zgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(imax, &
                              1),ldw,cone,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_izamax(imax - k,w(k,k + 1),1)
                    rowmax = cabs1(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_izamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_zcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_zcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_zcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_zswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_zswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_zcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = cone/a(k,k)
                       call la_zscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_zgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_zgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_zswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_zlasyf
#ifdef LA_WITH_XDP
     !> YLASYF: computes a partial factorization of a complex symmetric matrix
     !> A using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) ( 0   A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**T denotes the transpose of U.
     !> YLASYF is an auxiliary routine called by YSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_ylasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info,kb
           integer(ilp),intent(in) :: lda,ldw,n,nb
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: w(ldw,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sevten = 17.0e+0_xdp

           ! Local Scalars
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(xdp) :: absakk,alpha,colmax,rowmax
           complex(xdp) :: d11,d21,d22,r1,t,z
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min,sqrt
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=xdp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_ycopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_ygemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k, &
                        kw + 1),ldw,cone,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
              if (k > 1) then
                 imax = la_iyamax(k - 1,w(1,kw),1)
                 colmax = cabs1(w(imax,kw))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_ycopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_ycopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_ygemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                               imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_iyamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = cabs1(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_iyamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_ycopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_ycopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    if (kp > 1) call la_ycopy(kp - 1,a(1,kk),1,a(1,kp),1)
                    ! interchange rows kk and kp in last k+1 to n columns of a
                    ! (columns k (or k and k-1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in last kkw to nb columns of w.
                    if (k < n) call la_yswap(n - k,a(kk,k + 1),lda,a(kp,k + 1),lda)
                    call la_yswap(n - kk + 1,w(kk,kkw),ldw,w(kp,kkw),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column kw of w now holds
                    ! w(kw) = u(k)*d(k),
                    ! where u(k) is the k-th column of u
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_ycopy(k,w(1,kw),1,a(1,k),1)
                    r1 = cone/a(k,k)
                    call la_yscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_ygemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_ygemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_yswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_ycopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_ygemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1),ldw, &
                         cone,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              if (k < n) then
                 imax = k + la_iyamax(n - k,w(k + 1,k),1)
                 colmax = cabs1(w(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero or underflow: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_ycopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_ycopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_ygemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(imax, &
                              1),ldw,cone,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_iyamax(imax - k,w(k,k + 1),1)
                    rowmax = cabs1(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_iyamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_ycopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_ycopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_ycopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_yswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_yswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_ycopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = cone/a(k,k)
                       call la_yscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_ygemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_ygemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_yswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_ylasyf
#endif
#ifdef LA_WITH_QP
     !> WLASYF: computes a partial factorization of a complex symmetric matrix
     !> A using the Bunch-Kaufman diagonal pivoting method. The partial
     !> factorization has the form:
     !> A  =  ( I  U12 ) ( A11  0  ) (  I       0    )  if UPLO = 'U', or:
     !> ( 0  U22 ) (  0   D  ) ( U12**T U22**T )
     !> A  =  ( L11  0 ) ( D    0  ) ( L11**T L21**T )  if UPLO = 'L'
     !> ( L21  I ) ( 0   A22 ) (  0       I    )
     !> where the order of D is at most NB. The actual order is returned in
     !> the argument KB, and is either NB or NB-1, or N if N <= NB.
     !> Note that U**T denotes the transpose of U.
     !> WLASYF is an auxiliary routine called by WSYTRF. It uses blocked code
     !> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
     !> A22 (if UPLO = 'L').

     pure subroutine la_wlasyf(uplo,n,nb,kb,a,lda,ipiv,w,ldw,info)
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
           integer(ilp) :: imax,j,jb,jj,jmax,jp,k,kk,kkw,kp,kstep,kw
           real(qp) :: absakk,alpha,colmax,rowmax
           complex(qp) :: d11,d21,d22,r1,t,z
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=qp)) + abs(aimag(z))
           ! Executable Statements
           info = 0
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (la_lsame(uplo,'U')) then
              ! factorize the trailing columns of a using the upper triangle
              ! of a and working backwards, and compute the matrix w = u12*d
              ! for use in updating a11
              ! k is the main loop index, decreasing from n in steps of 1 or 2
              ! kw is the column of w which corresponds to column k of a
              k = n
              10 continue
              kw = nb + k - n
              ! exit from loop
              if ((k <= n - nb + 1 .and. nb < n) .or. k < 1) go to 30
              ! copy column k of a to column kw of w and update it
              call la_wcopy(k,a(1,k),1,w(1,kw),1)
              if (k < n) call la_wgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w(k, &
                        kw + 1),ldw,cone,w(1,kw),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,kw))
              ! imax is the row-index of the largest off-diagonal element in
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
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column kw-1 of w and update it
                    call la_wcopy(imax,a(1,imax),1,w(1,kw - 1),1)
                    call la_wcopy(k - imax,a(imax,imax + 1),lda,w(imax + 1,kw - 1),1)

                    if (k < n) call la_wgemv('NO TRANSPOSE',k,n - k,-cone,a(1,k + 1),lda,w( &
                               imax,kw + 1),ldw,cone,w(1,kw - 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_iwamax(k - imax,w(imax + 1,kw - 1),1)
                    rowmax = cabs1(w(jmax,kw - 1))
                    if (imax > 1) then
                       jmax = la_iwamax(imax - 1,w(1,kw - 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,kw - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,kw - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column kw-1 of w to column kw of w
                       call la_wcopy(k,w(1,kw - 1),1,w(1,kw),1)
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k - kstep + 1
                 ! kkw is the column of w which corresponds to column kk of a
                 kkw = nb + kk - n
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kkw of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k-1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_wcopy(kk - 1 - kp,a(kp + 1,kk),1,a(kp,kp + 1),lda)
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
                    ! store subdiag. elements of column u(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! note: diagonal element u(k,k) is a unit element
                    ! and not stored.
                       ! a(k,k) := d(k,k) = w(k,kw)
                       ! a(1:k-1,k) := u(1:k-1,k) = w(1:k-1,kw)/d(k,k)
                    call la_wcopy(k,w(1,kw),1,a(1,k),1)
                    r1 = cone/a(k,k)
                    call la_wscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns kw and kw-1 of w now hold
                    ! ( w(kw-1) w(kw) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! store u(1:k-2,k-1) and u(1:k-2,k) and 2-by-2
                    ! block d(k-1:k,k-1:k) in columns k-1 and k of a.
                    ! note: 2-by-2 diagonal block u(k-1:k,k-1:k) is a unit
                    ! block and not stored.
                       ! a(k-1:k,k-1:k) := d(k-1:k,k-1:k) = w(k-1:k,kw-1:kw)
                       ! a(1:k-2,k-1:k) := u(1:k-2,k:k-1:k) =
                       ! = w(1:k-2,kw-1:kw) * ( d(k-1:k,k-1:k)**(-1) )
                    if (k > 2) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(kw-1) w(kw) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k - 1,kw)
                       d11 = w(k,kw)/d21
                       d22 = w(k - 1,kw - 1)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k-1) and a(k) as
                       ! dot products of rows of ( w(kw-1) w(kw) ) and columns
                       ! of d**(-1)
                       do j = 1,k - 2
                          a(j,k - 1) = d21*(d11*w(j,kw - 1) - w(j,kw))
                          a(j,k) = d21*(d22*w(j,kw) - w(j,kw - 1))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k - 1,k - 1) = w(k - 1,kw - 1)
                    a(k - 1,k) = w(k - 1,kw)
                    a(k,k) = w(k,kw)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
              30 continue
              ! update the upper triangle of a11 (= a(1:k,1:k)) as
              ! a11 := a11 - u12*d*u12**t = a11 - u12*w**t
              ! computing blocks of nb columns at a time
              do j = ((k - 1)/nb)*nb + 1,1,-nb
                 jb = min(nb,k - j + 1)
                 ! update the upper triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_wgemv('NO TRANSPOSE',jj - j + 1,n - k,-cone,a(j,k + 1),lda,w(jj, &
                               kw + 1),ldw,cone,a(j,jj),1)
                 end do
                 ! update the rectangular superdiagonal block
                 call la_wgemm('NO TRANSPOSE','TRANSPOSE',j - 1,jb,n - k,-cone,a(1,k + 1), &
                           lda,w(j,kw + 1),ldw,cone,a(1,j),lda)
              end do
              ! put u12 in standard form by partially undoing the interchanges
              ! in columns k+1:n looping backwards from k+1 to n
              j = k + 1
              60 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j + 1
                 end if
                 ! (note: here, j is used to determine row length. length n-j+1
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j + 1
                 if (jp /= jj .and. j <= n) call la_wswap(n - j + 1,a(jp,j),lda,a(jj,j), &
                           lda)
              if (j < n) go to 60
              ! set kb to the number of columns factorized
              kb = n - k
           else
              ! factorize the leading columns of a using the lower triangle
              ! of a and working forwards, and compute the matrix w = l21*d
              ! for use in updating a22
              ! k is the main loop index, increasing from 1 in steps of 1 or 2
              k = 1
              70 continue
              ! exit from loop
              if ((k >= nb .and. nb < n) .or. k > n) go to 90
              ! copy column k of a to column k of w and update it
              call la_wcopy(n - k + 1,a(k,k),1,w(k,k),1)
              call la_wgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(k,1),ldw, &
                         cone,w(k,k),1)
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(w(k,k))
              ! imax is the row-index of the largest off-diagonal element in
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
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! copy column imax to column k+1 of w and update it
                    call la_wcopy(imax - k,a(imax,k),lda,w(k,k + 1),1)
                    call la_wcopy(n - imax + 1,a(imax,imax),1,w(imax,k + 1),1)
                    call la_wgemv('NO TRANSPOSE',n - k + 1,k - 1,-cone,a(k,1),lda,w(imax, &
                              1),ldw,cone,w(k,k + 1),1)
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_iwamax(imax - k,w(k,k + 1),1)
                    rowmax = cabs1(w(jmax,k + 1))
                    if (imax < n) then
                       jmax = imax + la_iwamax(n - imax,w(imax + 1,k + 1),1)
                       rowmax = max(rowmax,cabs1(w(jmax,k + 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(w(imax,k + 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                       ! copy column k+1 of w to column k of w
                       call la_wcopy(n - k + 1,w(k,k + 1),1,w(k,k),1)
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 ! ============================================================
                 ! kk is the column of a where pivoting step stopped
                 kk = k + kstep - 1
                 ! interchange rows and columns kp and kk.
                 ! updated column kp is already stored in column kk of w.
                 if (kp /= kk) then
                    ! copy non-updated column kk to column kp of submatrix a
                    ! at step k. no need to copy element into column k
                    ! (or k and k+1 for 2-by-2 pivot) of a, since these columns
                    ! will be later overwritten.
                    a(kp,kp) = a(kk,kk)
                    call la_wcopy(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    if (kp < n) call la_wcopy(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    ! interchange rows kk and kp in first k-1 columns of a
                    ! (columns k (or k and k+1 for 2-by-2 pivot) of a will be
                    ! later overwritten). interchange rows kk and kp
                    ! in first kk columns of w.
                    if (k > 1) call la_wswap(k - 1,a(kk,1),lda,a(kp,1),lda)
                    call la_wswap(kk,w(kk,1),ldw,w(kp,1),ldw)
                 end if
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k of w now holds
                    ! w(k) = l(k)*d(k),
                    ! where l(k) is the k-th column of l
                    ! store subdiag. elements of column l(k)
                    ! and 1-by-1 block d(k) in column k of a.
                    ! (note: diagonal element l(k,k) is a unit element
                    ! and not stored)
                       ! a(k,k) := d(k,k) = w(k,k)
                       ! a(k+1:n,k) := l(k+1:n,k) = w(k+1:n,k)/d(k,k)
                    call la_wcopy(n - k + 1,w(k,k),1,a(k,k),1)
                    if (k < n) then
                       r1 = cone/a(k,k)
                       call la_wscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 of w now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    ! store l(k+2:n,k) and l(k+2:n,k+1) and 2-by-2
                    ! block d(k:k+1,k:k+1) in columns k and k+1 of a.
                    ! (note: 2-by-2 diagonal block l(k:k+1,k:k+1) is a unit
                    ! block and not stored)
                       ! a(k:k+1,k:k+1) := d(k:k+1,k:k+1) = w(k:k+1,k:k+1)
                       ! a(k+2:n,k:k+1) := l(k+2:n,k:k+1) =
                       ! = w(k+2:n,k:k+1) * ( d(k:k+1,k:k+1)**(-1) )
                    if (k < n - 1) then
                       ! compose the columns of the inverse of 2-by-2 pivot
                       ! block d in the following way to reduce the number
                       ! of flops when we myltiply panel ( w(k) w(k+1) ) by
                       ! this inverse
                       ! d**(-1) = ( d11 d21 )**(-1) =
                                 ! ( d21 d22 )
                       ! = 1/(d11*d22-d21**2) * ( ( d22 ) (-d21 ) ) =
                                              ! ( (-d21 ) ( d11 ) )
                       ! = 1/d21 * 1/((d11/d21)*(d22/d21)-1) *
                         ! * ( ( d22/d21 ) (      -1 ) ) =
                           ! ( (      -1 ) ( d11/d21 ) )
                       ! = 1/d21 * 1/(d22*d11-1) * ( ( d11 ) (  -1 ) ) =
                                                 ! ( ( -1  ) ( d22 ) )
                       ! = 1/d21 * t * ( ( d11 ) (  -1 ) )
                                     ! ( (  -1 ) ( d22 ) )
                       ! = d21 * ( ( d11 ) (  -1 ) )
                               ! ( (  -1 ) ( d22 ) )
                       d21 = w(k + 1,k)
                       d11 = w(k + 1,k + 1)/d21
                       d22 = w(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       ! update elements in columns a(k) and a(k+1) as
                       ! dot products of rows of ( w(k) w(k+1) ) and columns
                       ! of d**(-1)
                       do j = k + 2,n
                          a(j,k) = d21*(d11*w(j,k) - w(j,k + 1))
                          a(j,k + 1) = d21*(d22*w(j,k + 1) - w(j,k))
                       end do
                    end if
                    ! copy d(k) to a
                    a(k,k) = w(k,k)
                    a(k + 1,k) = w(k + 1,k)
                    a(k + 1,k + 1) = w(k + 1,k + 1)
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 70
              90 continue
              ! update the lower triangle of a22 (= a(k:n,k:n)) as
              ! a22 := a22 - l21*d*l21**t = a22 - l21*w**t
              ! computing blocks of nb columns at a time
              do j = k,n,nb
                 jb = min(nb,n - j + 1)
                 ! update the lower triangle of the diagonal block
                 do jj = j,j + jb - 1
                    call la_wgemv('NO TRANSPOSE',j + jb - jj,k - 1,-cone,a(jj,1),lda,w(jj, &
                               1),ldw,cone,a(jj,jj),1)
                 end do
                 ! update the rectangular subdiagonal block
                 if (j + jb <= n) call la_wgemm('NO TRANSPOSE','TRANSPOSE',n - j - jb + 1,jb,k - 1,- &
                           cone,a(j + jb,1),lda,w(j,1),ldw,cone,a(j + jb,j),lda)
              end do
              ! put l21 in standard form by partially undoing the interchanges
              ! of rows in columns 1:k-1 looping backwards from k-1 to 1
              j = k - 1
              120 continue
                 ! undo the interchanges (if any) of rows jj and jp at each
                 ! step j
                 ! (here, j is a diagonal index)
                 jj = j
                 jp = ipiv(j)
                 if (jp < 0) then
                    jp = -jp
                    ! (here, j is a diagonal index)
                    j = j - 1
                 end if
                 ! (note: here, j is used to determine row length. length j
                 ! of the rows to swap back doesn't include diagonal element)
                 j = j - 1
                 if (jp /= jj .and. j >= 1) call la_wswap(j,a(jp,1),lda,a(jj,1),lda)

              if (j > 1) go to 120
              ! set kb to the number of columns factorized
              kb = k - 1
           end if
           return
     end subroutine la_wlasyf
#endif

     !> CSPTRF: computes the factorization of a complex symmetric matrix A
     !> stored in packed format using the Bunch-Kaufman diagonal pivoting
     !> method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_csptrf(uplo,n,ap,ipiv,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(sp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(sp) :: absakk,alpha,colmax,rowmax
           complex(sp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,zdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real,sqrt
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
           end if
           if (info /= 0) then
              call la_xerbla('CSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_icamax(k - 1,ap(kc),1)
                 colmax = cabs1(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_icamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_cswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/ap(kc + k - 1)
                    call la_cspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_cscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_icamax(n - k,ap(kc + 1),1)
                 colmax = cabs1(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_icamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_cswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/ap(kc)
                       call la_cspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_cscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_csptrf
     !> ZSPTRF: computes the factorization of a complex symmetric matrix A
     !> stored in packed format using the Bunch-Kaufman diagonal pivoting
     !> method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_zsptrf(uplo,n,ap,ipiv,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(dp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(dp) :: absakk,alpha,colmax,rowmax
           complex(dp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
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
           end if
           if (info /= 0) then
              call la_xerbla('ZSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_izamax(k - 1,ap(kc),1)
                 colmax = cabs1(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_izamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_zswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/ap(kc + k - 1)
                    call la_zspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_zscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_izamax(n - k,ap(kc + 1),1)
                 colmax = cabs1(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_izamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_zswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/ap(kc)
                       call la_zspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_zscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_zsptrf
#ifdef LA_WITH_XDP
     !> YSPTRF: computes the factorization of a complex symmetric matrix A
     !> stored in packed format using the Bunch-Kaufman diagonal pivoting
     !> method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_ysptrf(uplo,n,ap,ipiv,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sevten = 17.0e+0_xdp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(xdp) :: absakk,alpha,colmax,rowmax
           complex(xdp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('YSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_iyamax(k - 1,ap(kc),1)
                 colmax = cabs1(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_iyamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_yswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/ap(kc + k - 1)
                    call la_yspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_yscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_iyamax(n - k,ap(kc + 1),1)
                 colmax = cabs1(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_iyamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_yswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/ap(kc)
                       call la_yspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_yscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_ysptrf
#endif
#ifdef LA_WITH_QP
     !> WSPTRF: computes the factorization of a complex symmetric matrix A
     !> stored in packed format using the Bunch-Kaufman diagonal pivoting
     !> method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.

     pure subroutine la_wsptrf(uplo,n,ap,ipiv,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(qp),intent(inout) :: ap(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kc,kk,knc,kp,kpc,kstep,kx,npp
           real(qp) :: absakk,alpha,colmax,rowmax
           complex(qp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
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
           end if
           if (info /= 0) then
              call la_xerbla('WSPTRF',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              kc = (n - 1)*n/2 + 1
              10 continue
              knc = kc
              ! if k < 1, exit from loop
              if (k < 1) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc + k - 1))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k > 1) then
                 imax = la_iwamax(k - 1,ap(kc),1)
                 colmax = cabs1(ap(kc + imax - 1))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    rowmax = zero
                    jmax = imax
                    kx = imax*(imax + 1)/2 + imax
                    do j = imax + 1,k
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + j
                    end do
                    kpc = (imax - 1)*imax/2 + 1
                    if (imax > 1) then
                       jmax = la_iwamax(imax - 1,ap(kpc),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - 1)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc + imax - 1)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kstep == 2) knc = knc - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_wswap(kp - 1,ap(knc),1,ap(kpc),1)
                    kx = kpc + kp - 1
                    do j = kp + 1,kk - 1
                       kx = kx + j - 1
                       t = ap(knc + j - 1)
                       ap(knc + j - 1) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc + kk - 1)
                    ap(knc + kk - 1) = ap(kpc + kp - 1)
                    ap(kpc + kp - 1) = t
                    if (kstep == 2) then
                       t = ap(kc + k - 2)
                       ap(kc + k - 2) = ap(kc + kp - 1)
                       ap(kc + kp - 1) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/ap(kc + k - 1)
                    call la_wspr(uplo,k - 1,-r1,ap(kc),1,ap)
                    ! store u(k) in column k
                    call la_wscal(k - 1,r1,ap(kc),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = ap(k - 1 + (k - 1)*k/2)
                       d22 = ap(k - 1 + (k - 2)*(k - 1)/2)/d12
                       d11 = ap(k + (k - 1)*k/2)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*ap(j + (k - 2)*(k - 1)/2) - ap(j + (k - 1)*k/2))

                          wk = d12*(d22*ap(j + (k - 1)*k/2) - ap(j + (k - 2)*(k - 1)/2))

                          do i = j,1,-1
                             ap(i + (j - 1)*j/2) = ap(i + (j - 1)*j/2) - ap(i + (k - 1)*k/2) &
                                       *wk - ap(i + (k - 2)*(k - 1)/2)*wkm1
                          end do
                          ap(j + (k - 1)*k/2) = wk
                          ap(j + (k - 2)*(k - 1)/2) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              kc = knc - k
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              kc = 1
              npp = n*(n + 1)/2
              60 continue
              knc = kc
              ! if k > n, exit from loop
              if (k > n) go to 110
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(ap(kc))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value
              if (k < n) then
                 imax = k + la_iwamax(n - k,ap(kc + 1),1)
                 colmax = cabs1(ap(kc + imax - k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero) then
                 ! column k is zero: set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    rowmax = zero
                    kx = kc + imax - k
                    do j = k,imax - 1
                       if (cabs1(ap(kx)) > rowmax) then
                          rowmax = cabs1(ap(kx))
                          jmax = j
                       end if
                       kx = kx + n - j
                    end do
                    kpc = npp - (n - imax + 1)*(n - imax + 2)/2 + 1
                    if (imax < n) then
                       jmax = imax + la_iwamax(n - imax,ap(kpc + 1),1)
                       rowmax = max(rowmax,cabs1(ap(kpc + jmax - imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(ap(kpc)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kstep == 2) knc = knc + n - k + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_wswap(n - kp,ap(knc + kp - kk + 1),1,ap(kpc + 1),1)

                    kx = knc + kp - kk
                    do j = kk + 1,kp - 1
                       kx = kx + n - j + 1
                       t = ap(knc + j - kk)
                       ap(knc + j - kk) = ap(kx)
                       ap(kx) = t
                    end do
                    t = ap(knc)
                    ap(knc) = ap(kpc)
                    ap(kpc) = t
                    if (kstep == 2) then
                       t = ap(kc + 1)
                       ap(kc + 1) = ap(kc + kp - k)
                       ap(kc + kp - k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/ap(kc)
                       call la_wspr(uplo,n - k,-r1,ap(kc + 1),1,ap(kc + n - k + 1))
                       ! store l(k) in column k
                       call la_wscal(n - k,r1,ap(kc + 1),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k): columns k and k+1 now hold
                    ! ( w(k) w(k+1) ) = ( l(k) l(k+1) )*d(k)
                    ! where l(k) and l(k+1) are the k-th and (k+1)-th columns
                    ! of l
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = ap(k + 1 + (k - 1)*(2*n - k)/2)
                       d11 = ap(k + 1 + k*(2*n - k - 1)/2)/d21
                       d22 = ap(k + (k - 1)*(2*n - k)/2)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*ap(j + (k - 1)*(2*n - k)/2) - ap(j + k*(2*n - k - 1)/2))

                          wkp1 = d21*(d22*ap(j + k*(2*n - k - 1)/2) - ap(j + (k - 1)*(2*n - k)/2) &
                                     )
                          do i = j,n
                             ap(i + (j - 1)*(2*n - j)/2) = ap(i + (j - 1)*(2*n - j)/2) - ap( &
                                       i + (k - 1)*(2*n - k)/2)*wk - ap(i + k*(2*n - k - 1)/2)*wkp1
                          end do
                          ap(j + (k - 1)*(2*n - k)/2) = wk
                          ap(j + k*(2*n - k - 1)/2) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              kc = knc + n - k + 2
              go to 60
           end if
           110 continue
           return
     end subroutine la_wsptrf
#endif

     !> CSYCONV: convert A given by TRF into L and D and vice-versa.
     !> Get Non-diag elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_csyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           complex(sp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! a is upper
                 ! convert a (a is upper)
                 ! convert value
              if (convert) then
                 i = n
                 e(1) = czero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = czero
                       a(i - 1,i) = czero
                       i = i - 1
                    else
                       e(i) = czero
                    end if
                    i = i - 1
                 end do
                 ! convert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i,j)
                            a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                             temp = a(ip,j)
                             a(ip,j) = a(i - 1,j)
                             a(i - 1,j) = temp
                          end do
                       end if
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              else
                 ! revert a (a is upper)
                 ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
                 ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
              ! a is lower
              if (convert) then
                 ! convert a (a is lower)
                 ! convert value
                 i = 1
                 e(n) = czero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = czero
                       a(i + 1,i) = czero
                       i = i + 1
                    else
                       e(i) = czero
                    end if
                    i = i + 1
                 end do
                 ! convert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i,j)
                             a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i + 1,j)
                             a(i + 1,j) = temp
                          end do
                       end if
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              else
                 ! revert a (a is lower)
                 ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
                 ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_csyconv
     !> ZSYCONV: converts A given by ZHETRF into L and D or vice-versa.
     !> Get nondiagonal elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_zsyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           complex(dp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! a is upper
              if (convert) then
                 ! convert a (a is upper)
                 ! convert value
                 i = n
                 e(1) = czero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = czero
                       a(i - 1,i) = czero
                       i = i - 1
                    else
                       e(i) = czero
                    end if
                    i = i - 1
                 end do
                 ! convert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i,j)
                            a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                             temp = a(ip,j)
                             a(ip,j) = a(i - 1,j)
                             a(i - 1,j) = temp
                          end do
                       end if
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              else
                 ! revert a (a is upper)
                 ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
                 ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
              ! a is lower
              if (convert) then
                 ! convert a (a is lower)
                 ! convert value
                 i = 1
                 e(n) = czero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = czero
                       a(i + 1,i) = czero
                       i = i + 1
                    else
                       e(i) = czero
                    end if
                    i = i + 1
                 end do
                 ! convert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i,j)
                             a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i + 1,j)
                             a(i + 1,j) = temp
                          end do
                       end if
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              else
                 ! revert a (a is lower)
                 ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
                 ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_zsyconv
#ifdef LA_WITH_XDP
     !> YSYCONV: converts A given by YHETRF into L and D or vice-versa.
     !> Get nondiagonal elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_ysyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           complex(xdp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('YSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! a is upper
              if (convert) then
                 ! convert a (a is upper)
                 ! convert value
                 i = n
                 e(1) = czero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = czero
                       a(i - 1,i) = czero
                       i = i - 1
                    else
                       e(i) = czero
                    end if
                    i = i - 1
                 end do
                 ! convert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i,j)
                            a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                             temp = a(ip,j)
                             a(ip,j) = a(i - 1,j)
                             a(i - 1,j) = temp
                          end do
                       end if
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              else
                 ! revert a (a is upper)
                 ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
                 ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
              ! a is lower
              if (convert) then
                 ! convert a (a is lower)
                 ! convert value
                 i = 1
                 e(n) = czero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = czero
                       a(i + 1,i) = czero
                       i = i + 1
                    else
                       e(i) = czero
                    end if
                    i = i + 1
                 end do
                 ! convert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i,j)
                             a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i + 1,j)
                             a(i + 1,j) = temp
                          end do
                       end if
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              else
                 ! revert a (a is lower)
                 ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
                 ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_ysyconv
#endif
#ifdef LA_WITH_QP
     !> WSYCONV: converts A given by WHETRF into L and D or vice-versa.
     !> Get nondiagonal elements of D (returned in workspace) and
     !> apply or reverse permutation done in TRF.

     pure subroutine la_wsyconv(uplo,way,n,a,lda,ipiv,e,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo,way
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: e(*)
        ! =====================================================================

           ! External Subroutines
           logical(lk) :: upper,convert
           integer(ilp) :: i,ip,j
           complex(qp) :: temp
           ! Executable Statements
           info = 0
           upper = la_lsame(uplo,'U')
           convert = la_lsame(way,'C')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (.not. convert .and. .not. la_lsame(way,'R')) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (lda < max(1,n)) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WSYCONV',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (upper) then
              ! a is upper
              if (convert) then
                 ! convert a (a is upper)
                 ! convert value
                 i = n
                 e(1) = czero
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       e(i) = a(i - 1,i)
                       e(i - 1) = czero
                       a(i - 1,i) = czero
                       i = i - 1
                    else
                       e(i) = czero
                    end if
                    i = i - 1
                 end do
                 ! convert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i,j)
                            a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i < n) then
                          do j = i + 1,n
                             temp = a(ip,j)
                             a(ip,j) = a(i - 1,j)
                             a(i - 1,j) = temp
                          end do
                       end if
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              else
                 ! revert a (a is upper)
                 ! revert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i < n) then
                       do j = i + 1,n
                         temp = a(ip,j)
                         a(ip,j) = a(i,j)
                         a(i,j) = temp
                       end do
                       end if
                    else
                      ip = -ipiv(i)
                      i = i + 1
                      if (i < n) then
                         do j = i + 1,n
                            temp = a(ip,j)
                            a(ip,j) = a(i - 1,j)
                            a(i - 1,j) = temp
                         end do
                      end if
                    end if
                    i = i + 1
                 end do
                 ! revert value
                 i = n
                 do while (i > 1)
                    if (ipiv(i) < 0) then
                       a(i - 1,i) = e(i)
                       i = i - 1
                    end if
                    i = i - 1
                 end do
              end if
           else
              ! a is lower
              if (convert) then
                 ! convert a (a is lower)
                 ! convert value
                 i = 1
                 e(n) = czero
                 do while (i <= n)
                    if (i < n .and. ipiv(i) < 0) then
                       e(i) = a(i + 1,i)
                       e(i + 1) = czero
                       a(i + 1,i) = czero
                       i = i + 1
                    else
                       e(i) = czero
                    end if
                    i = i + 1
                 end do
                 ! convert permutations
                 i = 1
                 do while (i <= n)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i,j)
                             a(i,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(ip,j)
                             a(ip,j) = a(i + 1,j)
                             a(i + 1,j) = temp
                          end do
                       end if
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              else
                 ! revert a (a is lower)
                 ! revert permutations
                 i = n
                 do while (i >= 1)
                    if (ipiv(i) > 0) then
                       ip = ipiv(i)
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i,j)
                             a(i,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    else
                       ip = -ipiv(i)
                       i = i - 1
                       if (i > 1) then
                          do j = 1,i - 1
                             temp = a(i + 1,j)
                             a(i + 1,j) = a(ip,j)
                             a(ip,j) = temp
                          end do
                       end if
                    end if
                    i = i - 1
                 end do
                 ! revert value
                 i = 1
                 do while (i <= n - 1)
                    if (ipiv(i) < 0) then
                       a(i + 1,i) = e(i)
                       i = i + 1
                    end if
                    i = i + 1
                 end do
              end if
           end if
           return
     end subroutine la_wsyconv
#endif

     !> CSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_csyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_sp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(sp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           complex(sp),intent(in) :: a(lda,*)
           complex(sp),intent(out) :: work(*)
           real(sp),intent(out) :: s(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(sp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           complex(sp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,aimag,int,log,max,min,real,sqrt
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('CSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_sp/s(j)
           end do
           tol = one/sqrt(2.0_sp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + real(s(i)*work(i),KIND=sp)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_classq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = cabs1(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = real(n - 2,KIND=sp)*(real(work(i),KIND=sp) - t*si)
                 c0 = -(t*si)*si + 2*real(work(i),KIND=sp)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + real(work(i),KIND=sp))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_slamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_slamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_csyequb
     !> ZSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_zsyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_dp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(dp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           complex(dp),intent(in) :: a(lda,*)
           complex(dp),intent(out) :: work(*)
           real(dp),intent(out) :: s(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(dp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,int,log,max,min,sqrt
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('ZSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_dp/s(j)
           end do
           tol = one/sqrt(2.0_dp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*real(work(i),KIND=dp)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_zlassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = cabs1(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(real(work(i),KIND=dp) - t*si)
                 c0 = -(t*si)*si + 2*real(work(i),KIND=dp)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + real(work(i),KIND=dp))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_dlamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_dlamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_zsyequb
#ifdef LA_WITH_XDP
     !> YSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_ysyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_xdp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(xdp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(out) :: work(*)
           real(xdp),intent(out) :: s(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(xdp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           complex(xdp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,int,log,max,min,sqrt
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('YSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_xdp/s(j)
           end do
           tol = one/sqrt(2.0_xdp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*real(work(i),KIND=xdp)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_ylassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = cabs1(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(real(work(i),KIND=xdp) - t*si)
                 c0 = -(t*si)*si + 2*real(work(i),KIND=xdp)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + real(work(i),KIND=xdp))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_xlamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_xlamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_ysyequb
#endif
#ifdef LA_WITH_QP
     !> WSYEQUB: computes row and column scalings intended to equilibrate a
     !> symmetric matrix A (with respect to the Euclidean norm) and reduce
     !> its condition number. The scale factors S are computed by the BIN
     !> algorithm (see references) so that the scaled matrix B with elements
     !> B(i,j) = S(i)*A(i,j)*S(j) has a condition number within a factor N of
     !> the smallest possible condition number over all possible diagonal
     !> scalings.

     pure subroutine la_wsyequb(uplo,n,a,lda,s,scond,amax,work,info)
        use la_constants_qp,only:zero,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(qp),intent(out) :: amax,scond
           character,intent(in) :: uplo
           ! Array Arguments
           complex(qp),intent(in) :: a(lda,*)
           complex(qp),intent(out) :: work(*)
           real(qp),intent(out) :: s(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: max_iter = 100

           ! Local Scalars
           integer(ilp) :: i,j,iter
           real(qp) :: avg,std,tol,c0,c1,c2,t,u,si,d,base,smin,smax,smlnum,bignum, &
                     scale,sumsq
           logical(lk) :: up
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,int,log,max,min,sqrt
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (.not. (la_lsame(uplo,'U') .or. la_lsame(uplo,'L'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,n)) then
              info = -4
           end if
           if (info /= 0) then
              call la_xerbla('WSYEQUB',-info)
              return
           end if
           up = la_lsame(uplo,'U')
           amax = zero
           ! quick return if possible.
           if (n == 0) then
              scond = one
              return
           end if
           do i = 1,n
              s(i) = zero
           end do
           amax = zero
           if (up) then
              do j = 1,n
                 do i = 1,j - 1
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
              end do
           else
              do j = 1,n
                 s(j) = max(s(j),cabs1(a(j,j)))
                 amax = max(amax,cabs1(a(j,j)))
                 do i = j + 1,n
                    s(i) = max(s(i),cabs1(a(i,j)))
                    s(j) = max(s(j),cabs1(a(i,j)))
                    amax = max(amax,cabs1(a(i,j)))
                 end do
              end do
           end if
           do j = 1,n
              s(j) = 1.0_qp/s(j)
           end do
           tol = one/sqrt(2.0_qp*n)
           do iter = 1,max_iter
              scale = zero
              sumsq = zero
              ! beta = |a|s
              do i = 1,n
                 work(i) = zero
              end do
              if (up) then
                 do j = 1,n
                    do i = 1,j - 1
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                 end do
              else
                 do j = 1,n
                    work(j) = work(j) + cabs1(a(j,j))*s(j)
                    do i = j + 1,n
                       work(i) = work(i) + cabs1(a(i,j))*s(j)
                       work(j) = work(j) + cabs1(a(i,j))*s(i)
                    end do
                 end do
              end if
              ! avg = s^t beta / n
              avg = zero
              do i = 1,n
                 avg = avg + s(i)*real(work(i),KIND=qp)
              end do
              avg = avg/n
              std = zero
              do i = n + 1,2*n
                 work(i) = s(i - n)*work(i - n) - avg
              end do
              call la_wlassq(n,work(n + 1),1,scale,sumsq)
              std = scale*sqrt(sumsq/n)
              if (std < tol*avg) goto 999
              do i = 1,n
                 t = cabs1(a(i,i))
                 si = s(i)
                 c2 = (n - 1)*t
                 c1 = (n - 2)*(real(work(i),KIND=qp) - t*si)
                 c0 = -(t*si)*si + 2*real(work(i),KIND=qp)*si - n*avg
                 d = c1*c1 - 4*c0*c2
                 if (d <= 0) then
                    info = -1
                    return
                 end if
                 si = -2*c0/(c1 + sqrt(d))
                 d = si - s(i)
                 u = zero
                 if (up) then
                    do j = 1,i
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 else
                    do j = 1,i
                       t = cabs1(a(i,j))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                    do j = i + 1,n
                       t = cabs1(a(j,i))
                       u = u + s(j)*t
                       work(j) = work(j) + d*t
                    end do
                 end if
                 avg = avg + (u + real(work(i),KIND=qp))*d/n
                 s(i) = si
              end do
           end do
           999 continue
           smlnum = la_qlamch('SAFEMIN')
           bignum = one/smlnum
           smin = bignum
           smax = zero
           t = one/sqrt(avg)
           base = la_qlamch('B')
           u = one/log(base)
           do i = 1,n
              s(i) = base**int(u*log(s(i)*t),KIND=ilp)
              smin = min(smin,s(i))
              smax = max(smax,s(i))
           end do
           scond = max(smin,smlnum)/min(smax,bignum)
     end subroutine la_wsyequb
#endif

     !> CSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_csyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           complex(sp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_cswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_cswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_csyswapr
     !> ZSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_zsyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           complex(dp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_zswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_zswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_zsyswapr
#ifdef LA_WITH_XDP
     !> YSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_ysyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           complex(xdp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_yswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_yswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_ysyswapr
#endif
#ifdef LA_WITH_QP
     !> WSYSWAPR: applies an elementary permutation on the rows and the columns of
     !> a symmetric matrix.

     pure subroutine la_wsyswapr(uplo,n,a,lda,i1,i2)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: i1,i2,lda,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,n)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i
           complex(qp) :: tmp
           ! Executable Statements
           upper = la_lsame(uplo,'U')
           if (upper) then
               ! upper
               ! first swap
                ! - swap column i1 and i2 from i1 to i1-1
              call la_wswap(i1 - 1,a(1,i1),1,a(1,i2),1)
                ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap row i1 from i1+1 to i2-1 with col i2 from i1+1 to i2-1
              tmp = a(i1,i1)
              a(i1,i1) = a(i2,i2)
              a(i2,i2) = tmp
              do i = 1,i2 - i1 - 1
                 tmp = a(i1,i1 + i)
                 a(i1,i1 + i) = a(i1 + i,i2)
                 a(i1 + i,i2) = tmp
              end do
                ! third swap
                ! - swap row i1 and i2 from i2+1 to n
              do i = i2 + 1,n
                 tmp = a(i1,i)
                 a(i1,i) = a(i2,i)
                 a(i2,i) = tmp
              end do
             else
               ! lower
               ! first swap
                ! - swap row i1 and i2 from i1 to i1-1
              call la_wswap(i1 - 1,a(i1,1),lda,a(i2,1),lda)
               ! second swap :
                ! - swap a(i1,i1) and a(i2,i2)
                ! - swap col i1 from i1+1 to i2-1 with row i2 from i1+1 to i2-1
               tmp = a(i1,i1)
               a(i1,i1) = a(i2,i2)
               a(i2,i2) = tmp
               do i = 1,i2 - i1 - 1
                  tmp = a(i1 + i,i1)
                  a(i1 + i,i1) = a(i2,i1 + i)
                  a(i2,i1 + i) = tmp
               end do
               ! third swap
                ! - swap col i1 and i2 from i2+1 to n
               do i = i2 + 1,n
                  tmp = a(i,i1)
                  a(i,i1) = a(i,i2)
                  a(i,i2) = tmp
               end do
           end if
     end subroutine la_wsyswapr
#endif

     !> CSYTF2: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_csytf2(uplo,n,a,lda,ipiv,info)
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
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: sevten = 17.0e+0_sp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(sp) :: absakk,alpha,colmax,rowmax
           complex(sp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,aimag,max,real,sqrt
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
              call la_xerbla('CSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_icamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_sisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_icamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_icamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_cswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_cswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/a(k,k)
                    call la_csyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_cscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_icamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_sisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_icamax(imax - k,a(imax,k),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_icamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_cswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_cswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/a(k,k)
                       call la_csyr(uplo,n - k,-r1,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_cscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_csytf2
     !> ZSYTF2: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_zsytf2(uplo,n,a,lda,ipiv,info)
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
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: sevten = 17.0e+0_dp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(dp) :: absakk,alpha,colmax,rowmax
           complex(dp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
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
              call la_xerbla('ZSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_izamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_disnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_izamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_izamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_zswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_zswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/a(k,k)
                    call la_zsyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_zscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_izamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_disnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_izamax(imax - k,a(imax,k),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_izamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_zswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_zswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/a(k,k)
                       call la_zsyr(uplo,n - k,-r1,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_zscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_zsytf2
#ifdef LA_WITH_XDP
     !> YSYTF2: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_ysytf2(uplo,n,a,lda,ipiv,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*)
        ! =====================================================================
           ! Parameters
           real(xdp),parameter :: sevten = 17.0e+0_xdp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(xdp) :: absakk,alpha,colmax,rowmax
           complex(xdp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(z) = abs(real(z,KIND=xdp)) + abs(aimag(z))
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
              call la_xerbla('YSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iyamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_xisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_iyamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_iyamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_yswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_yswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/a(k,k)
                    call la_ysyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_yscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iyamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_xisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_iyamax(imax - k,a(imax,k),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_iyamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_yswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_yswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/a(k,k)
                       call la_ysyr(uplo,n - k,-r1,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_yscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_ysytf2
#endif
#ifdef LA_WITH_QP
     !> WSYTF2: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method:
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, U**T is the transpose of U, and D is symmetric and
     !> block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the unblocked version of the algorithm, calling Level 2 BLAS.

     pure subroutine la_wsytf2(uplo,n,a,lda,ipiv,info)
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
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: sevten = 17.0e+0_qp

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,imax,j,jmax,k,kk,kp,kstep
           real(qp) :: absakk,alpha,colmax,rowmax
           complex(qp) :: d11,d12,d21,d22,r1,t,wk,wkm1,wkp1,z
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,sqrt
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
              call la_xerbla('WSYTF2',-info)
              return
           end if
           ! initialize alpha for use in choosing pivot block size.
           alpha = (one + sqrt(sevten))/eight
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k > 1) then
                 imax = la_iwamax(k - 1,a(1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_qisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = imax + la_iwamax(k - imax,a(imax,imax + 1),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax > 1) then
                       jmax = la_iwamax(imax - 1,a(1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k-1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k - kstep + 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the leading
                    ! submatrix a(1:k,1:k)
                    call la_wswap(kp - 1,a(1,kk),1,a(1,kp),1)
                    call la_wswap(kk - kp - 1,a(kp + 1,kk),1,a(kp,kp + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k - 1,k)
                       a(k - 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the leading submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = u(k)*d(k)
                    ! where u(k) is the k-th column of u
                    ! perform a rank-1 update of a(1:k-1,1:k-1) as
                    ! a := a - u(k)*d(k)*u(k)**t = a - w(k)*1/d(k)*w(k)**t
                    r1 = cone/a(k,k)
                    call la_wsyr(uplo,k - 1,-r1,a(1,k),1,a,lda)
                    ! store u(k) in column k
                    call la_wscal(k - 1,r1,a(1,k),1)
                 else
                    ! 2-by-2 pivot block d(k): columns k and k-1 now hold
                    ! ( w(k-1) w(k) ) = ( u(k-1) u(k) )*d(k)
                    ! where u(k) and u(k-1) are the k-th and (k-1)-th columns
                    ! of u
                    ! perform a rank-2 update of a(1:k-2,1:k-2) as
                    ! a := a - ( u(k-1) u(k) )*d(k)*( u(k-1) u(k) )**t
                       ! = a - ( w(k-1) w(k) )*inv(d(k))*( w(k-1) w(k) )**t
                    if (k > 2) then
                       d12 = a(k - 1,k)
                       d22 = a(k - 1,k - 1)/d12
                       d11 = a(k,k)/d12
                       t = cone/(d11*d22 - cone)
                       d12 = t/d12
                       do j = k - 2,1,-1
                          wkm1 = d12*(d11*a(j,k - 1) - a(j,k))
                          wk = d12*(d22*a(j,k) - a(j,k - 1))
                          do i = j,1,-1
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k - 1)*wkm1
                          end do
                          a(j,k) = wk
                          a(j,k - 1) = wkm1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k - 1) = -kp
              end if
              ! decrease k and return to the start of the main loop
              k = k - kstep
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2
              k = 1
              40 continue
              ! if k > n, exit from loop
              if (k > n) go to 70
              kstep = 1
              ! determine rows and columns to be interchanged and whether
              ! a 1-by-1 or 2-by-2 pivot block will be used
              absakk = cabs1(a(k,k))
              ! imax is the row-index of the largest off-diagonal element in
              ! column k, and colmax is its absolute value.
              ! determine both colmax and imax.
              if (k < n) then
                 imax = k + la_iwamax(n - k,a(k + 1,k),1)
                 colmax = cabs1(a(imax,k))
              else
                 colmax = zero
              end if
              if (max(absakk,colmax) == zero .or. la_qisnan(absakk)) then
                 ! column k is zero or underflow, or contains a nan:
                 ! set info and continue
                 if (info == 0) info = k
                 kp = k
              else
                 if (absakk >= alpha*colmax) then
                    ! no interchange, use 1-by-1 pivot block
                    kp = k
                 else
                    ! jmax is the column-index of the largest off-diagonal
                    ! element in row imax, and rowmax is its absolute value
                    jmax = k - 1 + la_iwamax(imax - k,a(imax,k),lda)
                    rowmax = cabs1(a(imax,jmax))
                    if (imax < n) then
                       jmax = imax + la_iwamax(n - imax,a(imax + 1,imax),1)
                       rowmax = max(rowmax,cabs1(a(jmax,imax)))
                    end if
                    if (absakk >= alpha*colmax*(colmax/rowmax)) then
                       ! no interchange, use 1-by-1 pivot block
                       kp = k
                    else if (cabs1(a(imax,imax)) >= alpha*rowmax) then
                       ! interchange rows and columns k and imax, use 1-by-1
                       ! pivot block
                       kp = imax
                    else
                       ! interchange rows and columns k+1 and imax, use 2-by-2
                       ! pivot block
                       kp = imax
                       kstep = 2
                    end if
                 end if
                 kk = k + kstep - 1
                 if (kp /= kk) then
                    ! interchange rows and columns kk and kp in the trailing
                    ! submatrix a(k:n,k:n)
                    if (kp < n) call la_wswap(n - kp,a(kp + 1,kk),1,a(kp + 1,kp),1)

                    call la_wswap(kp - kk - 1,a(kk + 1,kk),1,a(kp,kk + 1),lda)
                    t = a(kk,kk)
                    a(kk,kk) = a(kp,kp)
                    a(kp,kp) = t
                    if (kstep == 2) then
                       t = a(k + 1,k)
                       a(k + 1,k) = a(kp,k)
                       a(kp,k) = t
                    end if
                 end if
                 ! update the trailing submatrix
                 if (kstep == 1) then
                    ! 1-by-1 pivot block d(k): column k now holds
                    ! w(k) = l(k)*d(k)
                    ! where l(k) is the k-th column of l
                    if (k < n) then
                       ! perform a rank-1 update of a(k+1:n,k+1:n) as
                       ! a := a - l(k)*d(k)*l(k)**t = a - w(k)*(1/d(k))*w(k)**t
                       r1 = cone/a(k,k)
                       call la_wsyr(uplo,n - k,-r1,a(k + 1,k),1,a(k + 1,k + 1),lda)

                       ! store l(k) in column k
                       call la_wscal(n - k,r1,a(k + 1,k),1)
                    end if
                 else
                    ! 2-by-2 pivot block d(k)
                    if (k < n - 1) then
                       ! perform a rank-2 update of a(k+2:n,k+2:n) as
                       ! a := a - ( l(k) l(k+1) )*d(k)*( l(k) l(k+1) )**t
                          ! = a - ( w(k) w(k+1) )*inv(d(k))*( w(k) w(k+1) )**t
                       ! where l(k) and l(k+1) are the k-th and (k+1)-th
                       ! columns of l
                       d21 = a(k + 1,k)
                       d11 = a(k + 1,k + 1)/d21
                       d22 = a(k,k)/d21
                       t = cone/(d11*d22 - cone)
                       d21 = t/d21
                       do j = k + 2,n
                          wk = d21*(d11*a(j,k) - a(j,k + 1))
                          wkp1 = d21*(d22*a(j,k + 1) - a(j,k))
                          do i = j,n
                             a(i,j) = a(i,j) - a(i,k)*wk - a(i,k + 1)*wkp1
                          end do
                          a(j,k) = wk
                          a(j,k + 1) = wkp1
                       end do
                    end if
                 end if
              end if
              ! store details of the interchanges in ipiv
              if (kstep == 1) then
                 ipiv(k) = kp
              else
                 ipiv(k) = -kp
                 ipiv(k + 1) = -kp
              end if
              ! increase k and return to the start of the main loop
              k = k + kstep
              go to 40
           end if
           70 continue
           return
     end subroutine la_wsytf2
#endif

     !> CSYTRF: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_csytrf(uplo,n,a,lda,ipiv,work,lwork,info)
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
              nb = la_ilaenv(1,'CSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'CSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_clasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_clasyf(uplo,k,nb,kb,a,lda,ipiv,work,n,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_csytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_clasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_clasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,n, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_csytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_csytrf
     !> ZSYTRF: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_zsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
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
              nb = la_ilaenv(1,'ZSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'ZSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_zlasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_zlasyf(uplo,k,nb,kb,a,lda,ipiv,work,n,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_zsytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_zlasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_zlasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,n, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_zsytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_zsytrf
#ifdef LA_WITH_XDP
     !> YSYTRF: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_ysytrf(uplo,n,a,lda,ipiv,work,lwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: work(*)
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
              nb = la_ilaenv(1,'YSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'YSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_ylasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_ylasyf(uplo,k,nb,kb,a,lda,ipiv,work,n,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_ysytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_ylasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_ylasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,n, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_ysytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_ysytrf
#endif
#ifdef LA_WITH_QP
     !> WSYTRF: computes the factorization of a complex symmetric matrix A
     !> using the Bunch-Kaufman diagonal pivoting method.  The form of the
     !> factorization is
     !> A = U*D*U**T  or  A = L*D*L**T
     !> where U (or L) is a product of permutation and unit upper (lower)
     !> triangular matrices, and D is symmetric and block diagonal with
     !> 1-by-1 and 2-by-2 diagonal blocks.
     !> This is the blocked version of the algorithm, calling Level 3 BLAS.

     pure subroutine la_wsytrf(uplo,n,a,lda,ipiv,work,lwork,info)
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
              nb = la_ilaenv(1,'WSYTRF',uplo,n,-1,-1,-1)
              lwkopt = n*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WSYTRF',-info)
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
                 nbmin = max(2,la_ilaenv(2,'WSYTRF',uplo,n,-1,-1,-1))
              end if
           else
              iws = 1
           end if
           if (nb < nbmin) nb = n
           if (upper) then
              ! factorize a as u*d*u**t using the upper triangle of a
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! kb, where kb is the number of columns factorized by la_wlasyf;
              ! kb is either nb or nb-1, or k for the last block
              k = n
              10 continue
              ! if k < 1, exit from loop
              if (k < 1) go to 40
              if (k > nb) then
                 ! factorize columns k-kb+1:k of a and use blocked code to
                 ! update columns 1:k-kb
                 call la_wlasyf(uplo,k,nb,kb,a,lda,ipiv,work,n,iinfo)
              else
                 ! use unblocked code to factorize columns 1:k of a
                 call la_wsytf2(uplo,k,a,lda,ipiv,iinfo)
                 kb = k
              end if
              ! set info on the first occurrence of a zero pivot
              if (info == 0 .and. iinfo > 0) info = iinfo
              ! decrease k and return to the start of the main loop
              k = k - kb
              go to 10
           else
              ! factorize a as l*d*l**t using the lower triangle of a
              ! k is the main loop index, increasing from 1 to n in steps of
              ! kb, where kb is the number of columns factorized by la_wlasyf;
              ! kb is either nb or nb-1, or n-k+1 for the last block
              k = 1
              20 continue
              ! if k > n, exit from loop
              if (k > n) go to 40
              if (k <= n - nb) then
                 ! factorize columns k:k+kb-1 of a and use blocked code to
                 ! update columns k+kb:n
                 call la_wlasyf(uplo,n - k + 1,nb,kb,a(k,k),lda,ipiv(k),work,n, &
                           iinfo)
              else
                 ! use unblocked code to factorize columns k:n of a
                 call la_wsytf2(uplo,n - k + 1,a(k,k),lda,ipiv(k),iinfo)
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
     end subroutine la_wsytrf
#endif

     !> CSYTRI: computes the inverse of a complex symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> CSYTRF.

     pure subroutine la_csytri(uplo,n,a,lda,ipiv,work,info)
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
           integer(ilp) :: k,kp,kstep
           complex(sp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('CSYTRI',-info)
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
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_ccopy(k - 1,a(1,k),1,work,1)
                    call la_csymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_cdotu(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k + 1)
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_ccopy(k - 1,a(1,k),1,work,1)
                    call la_csymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_cdotu(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_cdotu(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_ccopy(k - 1,a(1,k + 1),1,work,1)
                    call la_csymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_cdotu(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_cswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_cswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_ccopy(n - k,a(k + 1,k),1,work,1)
                    call la_csymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_cdotu(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k - 1)
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_ccopy(n - k,a(k + 1,k),1,work,1)
                    call la_csymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_cdotu(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_cdotu(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_ccopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_csymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_cdotu(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_cswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_cswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_csytri
     !> ZSYTRI: computes the inverse of a complex symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> ZSYTRF.

     pure subroutine la_zsytri(uplo,n,a,lda,ipiv,work,info)
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
           integer(ilp) :: k,kp,kstep
           complex(dp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('ZSYTRI',-info)
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
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_zcopy(k - 1,a(1,k),1,work,1)
                    call la_zsymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_zdotu(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k + 1)
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_zcopy(k - 1,a(1,k),1,work,1)
                    call la_zsymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_zdotu(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_zdotu(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_zcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_zsymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_zdotu(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_zswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_zswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_zcopy(n - k,a(k + 1,k),1,work,1)
                    call la_zsymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_zdotu(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k - 1)
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_zcopy(n - k,a(k + 1,k),1,work,1)
                    call la_zsymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_zdotu(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_zdotu(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_zcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_zsymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_zdotu(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_zswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_zswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_zsytri
#ifdef LA_WITH_XDP
     !> YSYTRI: computes the inverse of a complex symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> YSYTRF.

     pure subroutine la_ysytri(uplo,n,a,lda,ipiv,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: k,kp,kstep
           complex(xdp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('YSYTRI',-info)
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
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_ycopy(k - 1,a(1,k),1,work,1)
                    call la_ysymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_ydotu(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k + 1)
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_ycopy(k - 1,a(1,k),1,work,1)
                    call la_ysymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_ydotu(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_ydotu(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_ycopy(k - 1,a(1,k + 1),1,work,1)
                    call la_ysymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_ydotu(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_yswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_yswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_ycopy(n - k,a(k + 1,k),1,work,1)
                    call la_ysymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_ydotu(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k - 1)
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_ycopy(n - k,a(k + 1,k),1,work,1)
                    call la_ysymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_ydotu(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_ydotu(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_ycopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_ysymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_ydotu(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_yswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_yswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_ysytri
#endif
#ifdef LA_WITH_QP
     !> WSYTRI: computes the inverse of a complex symmetric indefinite matrix
     !> A using the factorization A = U*D*U**T or A = L*D*L**T computed by
     !> WSYTRF.

     pure subroutine la_wsytri(uplo,n,a,lda,ipiv,work,info)
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
           integer(ilp) :: k,kp,kstep
           complex(qp) :: ak,akkp1,akp1,d,t,temp
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              call la_xerbla('WSYTRI',-info)
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
              ! compute inv(a) from the factorization a = u*d*u**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              30 continue
              ! if k > n, exit from loop.
              if (k > n) go to 40
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k > 1) then
                    call la_wcopy(k - 1,a(1,k),1,work,1)
                    call la_wsymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_wdotu(k - 1,work,1,a(1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k + 1)
                 ak = a(k,k)/t
                 akp1 = a(k + 1,k + 1)/t
                 akkp1 = a(k,k + 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k,k) = akp1/d
                 a(k + 1,k + 1) = ak/d
                 a(k,k + 1) = -akkp1/d
                 ! compute columns k and k+1 of the inverse.
                 if (k > 1) then
                    call la_wcopy(k - 1,a(1,k),1,work,1)
                    call la_wsymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k),1)

                    a(k,k) = a(k,k) - la_wdotu(k - 1,work,1,a(1,k),1)
                    a(k,k + 1) = a(k,k + 1) - la_wdotu(k - 1,a(1,k),1,a(1,k + 1),1)

                    call la_wcopy(k - 1,a(1,k + 1),1,work,1)
                    call la_wsymv(uplo,k - 1,-cone,a,lda,work,1,czero,a(1,k + 1),1)

                    a(k + 1,k + 1) = a(k + 1,k + 1) - la_wdotu(k - 1,work,1,a(1,k + 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the leading
                 ! submatrix a(1:k+1,1:k+1)
                 call la_wswap(kp - 1,a(1,k),1,a(1,kp),1)
                 call la_wswap(k - kp - 1,a(kp + 1,k),1,a(kp,kp + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k + 1)
                    a(k,k + 1) = a(kp,k + 1)
                    a(kp,k + 1) = temp
                 end if
              end if
              k = k + kstep
              go to 30
              40 continue
           else
              ! compute inv(a) from the factorization a = l*d*l**t.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              50 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 60
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! invert the diagonal block.
                 a(k,k) = cone/a(k,k)
                 ! compute column k of the inverse.
                 if (k < n) then
                    call la_wcopy(n - k,a(k + 1,k),1,work,1)
                    call la_wsymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_wdotu(n - k,work,1,a(k + 1,k),1)
                 end if
                 kstep = 1
              else
                 ! 2 x 2 diagonal block
                 ! invert the diagonal block.
                 t = a(k,k - 1)
                 ak = a(k - 1,k - 1)/t
                 akp1 = a(k,k)/t
                 akkp1 = a(k,k - 1)/t
                 d = t*(ak*akp1 - cone)
                 a(k - 1,k - 1) = akp1/d
                 a(k,k) = ak/d
                 a(k,k - 1) = -akkp1/d
                 ! compute columns k-1 and k of the inverse.
                 if (k < n) then
                    call la_wcopy(n - k,a(k + 1,k),1,work,1)
                    call la_wsymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k),1)
                    a(k,k) = a(k,k) - la_wdotu(n - k,work,1,a(k + 1,k),1)
                    a(k,k - 1) = a(k,k - 1) - la_wdotu(n - k,a(k + 1,k),1,a(k + 1,k - 1),1 &
                              )
                    call la_wcopy(n - k,a(k + 1,k - 1),1,work,1)
                    call la_wsymv(uplo,n - k,-cone,a(k + 1,k + 1),lda,work,1,czero,a(k + &
                              1,k - 1),1)
                    a(k - 1,k - 1) = a(k - 1,k - 1) - la_wdotu(n - k,work,1,a(k + 1,k - 1),1)

                 end if
                 kstep = 2
              end if
              kp = abs(ipiv(k))
              if (kp /= k) then
                 ! interchange rows and columns k and kp in the trailing
                 ! submatrix a(k-1:n,k-1:n)
                 if (kp < n) call la_wswap(n - kp,a(kp + 1,k),1,a(kp + 1,kp),1)
                 call la_wswap(kp - k - 1,a(k + 1,k),1,a(kp,k + 1),lda)
                 temp = a(k,k)
                 a(k,k) = a(kp,kp)
                 a(kp,kp) = temp
                 if (kstep == 2) then
                    temp = a(k,k - 1)
                    a(k,k - 1) = a(kp,k - 1)
                    a(kp,k - 1) = temp
                 end if
              end if
              k = k - kstep
              go to 50
              60 continue
           end if
           return
     end subroutine la_wsytri
#endif

     !> CSYTRS: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by CSYTRF.

     pure subroutine la_csytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
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
           complex(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('CSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 call la_cscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
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
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_cgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_cgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 call la_cgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1),1,cone,b( &
                            k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 call la_cscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
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
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_cgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + &
                           1,k),1,cone,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_cgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k) &
                              ,1,cone,b(k,1),ldb)
                    call la_cgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k - &
                              1),1,cone,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_csytrs
     !> ZSYTRS: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by ZSYTRF.

     pure subroutine la_zsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
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
           complex(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('ZSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 call la_zscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
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
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_zgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_zgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 call la_zgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1),1,cone,b( &
                            k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 call la_zscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
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
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_zgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + &
                           1,k),1,cone,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_zgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k) &
                              ,1,cone,b(k,1),ldb)
                    call la_zgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k - &
                              1),1,cone,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_zsytrs
#ifdef LA_WITH_XDP
     !> YSYTRS: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by YSYTRF.

     pure subroutine la_ysytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: j,k,kp
           complex(xdp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('YSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_ygeru(k - 1,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 ! multiply by the inverse of the diagonal block.
                 call la_yscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k - 1) call la_yswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(u(k)), where u(k) is the transformation
                 ! stored in columns k-1 and k of a.
                 call la_ygeru(k - 2,nrhs,-cone,a(1,k),1,b(k,1),ldb,b(1,1),ldb &
                           )
                 call la_ygeru(k - 2,nrhs,-cone,a(1,k - 1),1,b(k - 1,1),ldb,b(1,1), &
                           ldb)
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k - 1,k)
                 akm1 = a(k - 1,k - 1)/akm1k
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_ygemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_ygemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 call la_ygemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1),1,cone,b( &
                            k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_ygeru(n - k,nrhs,-cone,a(k + 1,k),1,b(k,1),ldb,b( &
                           k + 1,1),ldb)
                 ! multiply by the inverse of the diagonal block.
                 call la_yscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k + 1) call la_yswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)
                 ! multiply by inv(l(k)), where l(k) is the transformation
                 ! stored in columns k and k+1 of a.
                 if (k < n - 1) then
                    call la_ygeru(n - k - 1,nrhs,-cone,a(k + 2,k),1,b(k,1),ldb,b(k + 2, &
                              1),ldb)
                    call la_ygeru(n - k - 1,nrhs,-cone,a(k + 2,k + 1),1,b(k + 1,1),ldb,b( &
                              k + 2,1),ldb)
                 end if
                 ! multiply by the inverse of the diagonal block.
                 akm1k = a(k + 1,k)
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_ygemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + &
                           1,k),1,cone,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_ygemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k) &
                              ,1,cone,b(k,1),ldb)
                    call la_ygemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k - &
                              1),1,cone,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_ysytrs
#endif
#ifdef LA_WITH_QP
     !> WSYTRS: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by WSYTRF.

     pure subroutine la_wsytrs(uplo,n,nrhs,a,lda,ipiv,b,ldb,info)
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
           complex(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('WSYTRS',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
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
                 call la_wscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
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
                 ak = a(k,k)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k - 1,j)/akm1k
                    bk = b(k,j)/akm1k
                    b(k - 1,j) = (ak*bkm1 - bk)/denom
                    b(k,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k - 2
              end if
              go to 10
              30 continue
              ! next solve u**t *x = b, overwriting b with x.
              ! k is the main loop index, increasing from 1 to n in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = 1
              40 continue
              ! if k > n, exit from loop.
              if (k > n) go to 50
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(u**t(k)), where u(k) is the transformation
                 ! stored in column k of a.
                 call la_wgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(u**t(k+1)), where u(k+1) is the transformation
                 ! stored in columns k and k+1 of a.
                 call la_wgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k),1,cone,b( &
                           k,1),ldb)
                 call la_wgemv('TRANSPOSE',k - 1,nrhs,-cone,b,ldb,a(1,k + 1),1,cone,b( &
                            k + 1,1),ldb)
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 2
              end if
              go to 40
              50 continue
           else
              ! solve a*x = b, where a = l*d*l**t.
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
                 call la_wscal(nrhs,cone/a(k,k),b(k,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k+1 and -ipiv(k).
                 kp = -ipiv(k)
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
                 akm1 = a(k,k)/akm1k
                 ak = a(k + 1,k + 1)/akm1k
                 denom = akm1*ak - cone
                 do j = 1,nrhs
                    bkm1 = b(k,j)/akm1k
                    bk = b(k + 1,j)/akm1k
                    b(k,j) = (ak*bkm1 - bk)/denom
                    b(k + 1,j) = (akm1*bk - bkm1)/denom
                 end do
                 k = k + 2
              end if
              go to 60
              80 continue
              ! next solve l**t *x = b, overwriting b with x.
              ! k is the main loop index, decreasing from n to 1 in steps of
              ! 1 or 2, depending on the size of the diagonal blocks.
              k = n
              90 continue
              ! if k < 1, exit from loop.
              if (k < 1) go to 100
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! multiply by inv(l**t(k)), where l(k) is the transformation
                 ! stored in column k of a.
                 if (k < n) call la_wgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + &
                           1,k),1,cone,b(k,1),ldb)
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! multiply by inv(l**t(k-1)), where l(k-1) is the transformation
                 ! stored in columns k-1 and k of a.
                 if (k < n) then
                    call la_wgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k) &
                              ,1,cone,b(k,1),ldb)
                    call la_wgemv('TRANSPOSE',n - k,nrhs,-cone,b(k + 1,1),ldb,a(k + 1,k - &
                              1),1,cone,b(k - 1,1),ldb)
                 end if
                 ! interchange rows k and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 2
              end if
              go to 90
              100 continue
           end if
           return
     end subroutine la_wsytrs
#endif

     !> CSYTRS2: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by CSYTRF and converted by CSYCONV.

     pure subroutine la_csytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
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
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           complex(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('CSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_csyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_cswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_ctrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_cscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_ctrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_cswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_cswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_ctrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_cscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_ctrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_cswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_csyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_csytrs2
     !> ZSYTRS2: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by ZSYTRF and converted by ZSYCONV.

     pure subroutine la_zsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
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
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           complex(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('ZSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_zsyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_zswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_ztrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_zscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_ztrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_zswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_zswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_ztrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_zscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_ztrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_zswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_zsyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_zsytrs2
#ifdef LA_WITH_XDP
     !> YSYTRS2: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by YSYTRF and converted by YSYCONV.

     pure subroutine la_ysytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           complex(xdp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('YSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_ysyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_yswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_ytrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_yscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_ytrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_yswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_yswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_ytrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_yscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_ytrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_yswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_ysyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_ysytrs2
#endif
#ifdef LA_WITH_QP
     !> WSYTRS2: solves a system of linear equations A*X = B with a complex
     !> symmetric matrix A using the factorization A = U*D*U**T or
     !> A = L*D*L**T computed by WSYTRF and converted by WSYCONV.

     pure subroutine la_wsytrs2(uplo,n,nrhs,a,lda,ipiv,b,ldb,work,info)
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
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,iinfo,j,k,kp
           complex(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: max
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
              call la_xerbla('WSYTRS2',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           ! convert a
           call la_wsyconv(uplo,'C',n,a,lda,ipiv,work,iinfo)
           if (upper) then
              ! solve a*x = b, where a = u*d*u**t.
             ! p**t * b
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (kp == -ipiv(k - 1)) call la_wswap(nrhs,b(k - 1,1),ldb,b(kp,1),ldb &
                           )
                 k = k - 2
              end if
             end do
        ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
             call la_wtrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                   call la_wscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 elseif (i > 1) then
                    if (ipiv(i - 1) == ipiv(i)) then
                       akm1k = work(i)
                       akm1 = a(i - 1,i - 1)/akm1k
                       ak = a(i,i)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i - 1,j)/akm1k
                          bk = b(i,j)/akm1k
                          b(i - 1,j) = (ak*bkm1 - bk)/denom
                          b(i,j) = (akm1*bk - bkm1)/denom
                       end do
                    i = i - 1
                    end if
                 end if
                 i = i - 1
              end do
            ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_wtrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k < n .and. kp == -ipiv(k + 1)) call la_wswap(nrhs,b(k,1),ldb,b(kp, &
                            1),ldb)
                 k = k + 2
              end if
             end do
           else
              ! solve a*x = b, where a = l*d*l**t.
             ! p**t * b
             k = 1
             do while (k <= n)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k + 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k and -ipiv(k+1).
                 kp = -ipiv(k + 1)
                 if (kp == -ipiv(k)) call la_wswap(nrhs,b(k + 1,1),ldb,b(kp,1),ldb)

                 k = k + 2
              end if
             end do
        ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
             call la_wtrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
        ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                   call la_wscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else
                       akm1k = work(i)
                       akm1 = a(i,i)/akm1k
                       ak = a(i + 1,i + 1)/akm1k
                       denom = akm1*ak - cone
                       do j = 1,nrhs
                          bkm1 = b(i,j)/akm1k
                          bk = b(i + 1,j)/akm1k
                          b(i,j) = (ak*bkm1 - bk)/denom
                          b(i + 1,j) = (akm1*bk - bkm1)/denom
                       end do
                       i = i + 1
                 end if
                 i = i + 1
              end do
        ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
             call la_wtrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
             ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
             k = n
             do while (k >= 1)
              if (ipiv(k) > 0) then
                 ! 1 x 1 diagonal block
                 ! interchange rows k and ipiv(k).
                 kp = ipiv(k)
                 if (kp /= k) call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 k = k - 1
              else
                 ! 2 x 2 diagonal block
                 ! interchange rows k-1 and -ipiv(k).
                 kp = -ipiv(k)
                 if (k > 1 .and. kp == -ipiv(k - 1)) call la_wswap(nrhs,b(k,1),ldb,b(kp, &
                           1),ldb)
                 k = k - 2
              end if
             end do
           end if
           ! revert a
           call la_wsyconv(uplo,'R',n,a,lda,ipiv,work,iinfo)
           return
     end subroutine la_wsytrs2
#endif

     !> CSYTRS_3: solves a system of linear equations A * X = B with a complex
     !> symmetric matrix A using the factorization computed
     !> by CSYTRF_RK or CSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_csytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
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
           complex(sp),intent(in) :: a(lda,*),e(*)
           complex(sp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           complex(sp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('CSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_ctrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_cscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_ctrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv( i ) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_ctrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_cscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_ctrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_cswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_csytrs_3
     !> ZSYTRS_3: solves a system of linear equations A * X = B with a complex
     !> symmetric matrix A using the factorization computed
     !> by ZSYTRF_RK or ZSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_zsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
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
           complex(dp),intent(in) :: a(lda,*),e(*)
           complex(dp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           complex(dp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('ZSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_ztrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_zscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_ztrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_ztrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_zscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_ztrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_zswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_zsytrs_3
#ifdef LA_WITH_XDP
     !> YSYTRS_3: solves a system of linear equations A * X = B with a complex
     !> symmetric matrix A using the factorization computed
     !> by YSYTRF_RK or ZSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_ysytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldb,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(in) :: a(lda,*),e(*)
           complex(xdp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           complex(xdp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('YSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_ytrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_yscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_ytrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_ytrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_yscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_ytrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_yswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_ysytrs_3
#endif
#ifdef LA_WITH_QP
     !> WSYTRS_3: solves a system of linear equations A * X = B with a complex
     !> symmetric matrix A using the factorization computed
     !> by WSYTRF_RK or ZSYTRF_BK:
     !> A = P*U*D*(U**T)*(P**T) or A = P*L*D*(L**T)*(P**T),
     !> where U (or L) is unit upper (or lower) triangular matrix,
     !> U**T (or L**T) is the transpose of U (or L), P is a permutation
     !> matrix, P**T is the transpose of P, and D is symmetric and block
     !> diagonal with 1-by-1 and 2-by-2 diagonal blocks.
     !> This algorithm is using Level 3 BLAS.

     pure subroutine la_wsytrs_3(uplo,n,nrhs,a,lda,e,ipiv,b,ldb,info)
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
           complex(qp),intent(in) :: a(lda,*),e(*)
           complex(qp),intent(inout) :: b(ldb,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,j,k,kp
           complex(qp) :: ak,akm1,akm1k,bk,bkm1,denom
           ! Intrinsic Functions
           intrinsic :: abs,max
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
              info = -9
           end if
           if (info /= 0) then
              call la_xerbla('WSYTRS_3',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. nrhs == 0) return
           if (upper) then
              ! begin upper
              ! solve a*x = b, where a = u*d*u**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (u \p**t * b) -> b    [ (u \p**t * b) ]
              call la_wtrsm('L','U','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (u \p**t * b) ]
              i = n
              do while (i >= 1)
                 if (ipiv(i) > 0) then
                    call la_wscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i > 1) then
                    akm1k = e(i)
                    akm1 = a(i - 1,i - 1)/akm1k
                    ak = a(i,i)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i - 1,j)/akm1k
                       bk = b(i,j)/akm1k
                       b(i - 1,j) = (ak*bkm1 - bk)/denom
                       b(i,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i - 1
                 end if
                 i = i - 1
              end do
              ! compute (u**t \ b) -> b   [ u**t \ (d \ (u \p**t * b) ) ]
              call la_wtrsm('L','U','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (u**t \ (d \ (u \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for upper case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
           else
              ! begin lower
              ! solve a*x = b, where a = l*d*l**t.
              ! p**t * b
              ! interchange rows k and ipiv(k) of matrix b in the same order
              ! that the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with increment 1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = 1,n,1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! compute (l \p**t * b) -> b    [ (l \p**t * b) ]
              call la_wtrsm('L','L','N','U',n,nrhs,cone,a,lda,b,ldb)
              ! compute d \ b -> b   [ d \ (l \p**t * b) ]
              i = 1
              do while (i <= n)
                 if (ipiv(i) > 0) then
                    call la_wscal(nrhs,cone/a(i,i),b(i,1),ldb)
                 else if (i < n) then
                    akm1k = e(i)
                    akm1 = a(i,i)/akm1k
                    ak = a(i + 1,i + 1)/akm1k
                    denom = akm1*ak - cone
                    do j = 1,nrhs
                       bkm1 = b(i,j)/akm1k
                       bk = b(i + 1,j)/akm1k
                       b(i,j) = (ak*bkm1 - bk)/denom
                       b(i + 1,j) = (akm1*bk - bkm1)/denom
                    end do
                    i = i + 1
                 end if
                 i = i + 1
              end do
              ! compute (l**t \ b) -> b   [ l**t \ (d \ (l \p**t * b) ) ]
              call la_wtrsm('L','L','T','U',n,nrhs,cone,a,lda,b,ldb)
              ! p * b  [ p * (l**t \ (d \ (l \p**t * b) )) ]
              ! interchange rows k and ipiv(k) of matrix b in reverse order
              ! from the formation order of ipiv(i) vector for lower case.
              ! (we can do the simple loop over ipiv with decrement -1,
              ! since the abs value of ipiv(i) represents the row index
              ! of the interchange with row i in both 1x1 and 2x2 pivot cases)
              do k = n,1,-1
                 kp = abs(ipiv(k))
                 if (kp /= k) then
                    call la_wswap(nrhs,b(k,1),ldb,b(kp,1),ldb)
                 end if
              end do
              ! end lower
           end if
           return
     end subroutine la_wsytrs_3
#endif

     !> CLA_HERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(sp) function la_cla_herpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(sp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper,lsame
           complex(sp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(sp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=sp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_csytrs.
           ! calls to la_sswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= 0.0_sp) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= 0.0_sp) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_cla_herpvgrw = rpvgrw
     end function la_cla_herpvgrw
     !> ZLA_HERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(dp) function la_zla_herpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(dp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper,lsame
           complex(dp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(dp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=dp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_zsytrs.
           ! calls to la_sswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_zla_herpvgrw = rpvgrw
     end function la_zla_herpvgrw
#ifdef LA_WITH_XDP
     !> YLA_HERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(xdp) function la_yla_herpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(in) :: a(lda,*),af(ldaf,*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(xdp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper,lsame
           complex(xdp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_ysytrs.
           ! calls to la_dswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_yla_herpvgrw = rpvgrw
     end function la_yla_herpvgrw
#endif
#ifdef LA_WITH_QP
     !> WLA_HERPVGRW: computes the reciprocal pivot growth factor
     !> norm(A)/norm(U). The "max absolute element" norm is used. If this is
     !> much less than 1, the stability of the LU factorization of the
     !> (equilibrated) matrix A could be poor. This also means that the
     !> solution X, estimated condition numbers, and error bounds could be
     !> unreliable.

     real(qp) function la_wla_herpvgrw(uplo,n,info,a,lda,af,ldaf,ipiv,work)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(in) :: n,info,lda,ldaf
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: ncols,i,j,k,kp
           real(qp) :: amax,umax,rpvgrw,tmp
           logical(lk) :: upper,lsame
           complex(qp) :: zdum
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max,min
           ! Statement Functions
           real(qp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=qp)) + abs(aimag(zdum))
           ! Executable Statements
           upper = la_lsame('UPPER',uplo)
           if (info == 0) then
              if (upper) then
                 ncols = 1
              else
                 ncols = n
              end if
           else
              ncols = info
           end if
           rpvgrw = one
           do i = 1,2*n
              work(i) = zero
           end do
           ! find the max magnitude entry of each column of a.  compute the max
           ! for all n columns so we can apply the pivot permutation while
           ! looping below.  assume a full factorization is the common case.
           if (upper) then
              do j = 1,n
                 do i = 1,j
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           else
              do j = 1,n
                 do i = j,n
                    work(n + i) = max(cabs1(a(i,j)),work(n + i))
                    work(n + j) = max(cabs1(a(i,j)),work(n + j))
                 end do
              end do
           end if
           ! now find the max magnitude entry of each column of u or l.  also
           ! permute the magnitudes of a above so they're in the same order as
           ! the factor.
           ! the iteration orders and permutations were copied from la_wsytrs.
           ! calls to la_dswap would be severe overkill.
           if (upper) then
              k = n
              do while (k < ncols .and. k > 0)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = 1,k
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k - 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k - 1)
                    work(n + k - 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = 1,k - 1
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k - 1) = max(cabs1(af(i,k - 1)),work(k - 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k - 2
                 end if
              end do
              k = ncols
              do while (k <= n)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k + 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k + 2
                 end if
              end do
           else
              k = 1
              do while (k <= ncols)
                 if (ipiv(k) > 0) then
                    ! 1x1 pivot
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    do i = k,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                    end do
                    k = k + 1
                 else
                    ! 2x2 pivot
                    kp = -ipiv(k)
                    tmp = work(n + k + 1)
                    work(n + k + 1) = work(n + kp)
                    work(n + kp) = tmp
                    do i = k + 1,n
                       work(k) = max(cabs1(af(i,k)),work(k))
                       work(k + 1) = max(cabs1(af(i,k + 1)),work(k + 1))
                    end do
                    work(k) = max(cabs1(af(k,k)),work(k))
                    k = k + 2
                 end if
              end do
              k = ncols
              do while (k >= 1)
                 if (ipiv(k) > 0) then
                    kp = ipiv(k)
                    if (kp /= k) then
                       tmp = work(n + k)
                       work(n + k) = work(n + kp)
                       work(n + kp) = tmp
                    end if
                    k = k - 1
                 else
                    kp = -ipiv(k)
                    tmp = work(n + k)
                    work(n + k) = work(n + kp)
                    work(n + kp) = tmp
                    k = k - 2
                 end if
              end do
           end if
           ! compute the *inverse* of the max element growth factor.  dividing
           ! by zero would imply the largest entry of the factor's column is
           ! zero.  than can happen when either the column of a is zero or
           ! massive pivots made the factor underflow to zero.  neither counts
           ! as growth in itself, so simply ignore terms with zero
           ! denominators.
           if (upper) then
              do i = ncols,n
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           else
              do i = 1,ncols
                 umax = work(i)
                 amax = work(n + i)
                 if (umax /= zero) then
                    rpvgrw = min(amax/umax,rpvgrw)
                 end if
              end do
           end if
           la_wla_herpvgrw = rpvgrw
     end function la_wla_herpvgrw
#endif

     !> CSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric packed matrix A using the
     !> factorization A = U*D*U**T or A = L*D*L**T computed by CSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_cspcon(uplo,n,ap,ipiv,anorm,rcond,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(sp),intent(in) :: anorm
           real(sp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(sp),intent(in) :: ap(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(sp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('CSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_clacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_csptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_cspcon
     !> ZSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric packed matrix A using the
     !> factorization A = U*D*U**T or A = L*D*L**T computed by ZSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_zspcon(uplo,n,ap,ipiv,anorm,rcond,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(dp),intent(in) :: anorm
           real(dp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(dp),intent(in) :: ap(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(dp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('ZSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_zlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_zsptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_zspcon
#ifdef LA_WITH_XDP
     !> YSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric packed matrix A using the
     !> factorization A = U*D*U**T or A = L*D*L**T computed by YSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_yspcon(uplo,n,ap,ipiv,anorm,rcond,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(xdp),intent(in) :: anorm
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(in) :: ap(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(xdp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('YSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_ylacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_ysptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_yspcon
#endif
#ifdef LA_WITH_QP
     !> WSPCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric packed matrix A using the
     !> factorization A = U*D*U**T or A = L*D*L**T computed by WSPTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_wspcon(uplo,n,ap,ipiv,anorm,rcond,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           real(qp),intent(in) :: anorm
           real(qp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(qp),intent(in) :: ap(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,ip,kase
           real(qp) :: ainvnm
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           upper = la_lsame(uplo,'U')
           if (.not. upper .and. .not. la_lsame(uplo,'L')) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (anorm < zero) then
              info = -5
           end if
           if (info /= 0) then
              call la_xerbla('WSPCON',-info)
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
              ip = n*(n + 1)/2
              do i = n,1,-1
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip - i
              end do
           else
              ! lower triangular storage: examine d from top to bottom.
              ip = 1
              do i = 1,n
                 if (ipiv(i) > 0 .and. ap(ip) == zero) return
                 ip = ip + n - i + 1
              end do
           end if
           ! estimate the 1-norm of the inverse.
           kase = 0
           30 continue
           call la_wlacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_wsptrs(uplo,n,1,ap,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_wspcon
#endif

     !> CSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by CSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_csycon(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
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
              call la_xerbla('CSYCON',-info)
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
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_csytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_csycon
     !> ZSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by ZSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_zsycon(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
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
              call la_xerbla('ZSYCON',-info)
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
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_zsytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_zsycon
#ifdef LA_WITH_XDP
     !> YSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by YSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_ysycon(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,n
           real(xdp),intent(in) :: anorm
           real(xdp),intent(out) :: rcond
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           complex(xdp),intent(in) :: a(lda,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: i,kase
           real(xdp) :: ainvnm
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
              call la_xerbla('YSYCON',-info)
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
           call la_ylacn2(n,work(n + 1),work,ainvnm,kase,isave)
           if (kase /= 0) then
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_ysytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_ysycon
#endif
#ifdef LA_WITH_QP
     !> WSYCON: estimates the reciprocal of the condition number (in the
     !> 1-norm) of a complex symmetric matrix A using the factorization
     !> A = U*D*U**T or A = L*D*L**T computed by WSYTRF.
     !> An estimate is obtained for norm(inv(A)), and the reciprocal of the
     !> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).

     pure subroutine la_wsycon(uplo,n,a,lda,ipiv,anorm,rcond,work,info)
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
              call la_xerbla('WSYCON',-info)
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
              ! multiply by inv(l*d*l**t) or inv(u*d*u**t).
              call la_wsytrs(uplo,n,1,a,lda,ipiv,work,n,info)
              go to 30
           end if
           ! compute the estimate of the reciprocal condition number.
           if (ainvnm /= zero) rcond = (one/ainvnm)/anorm
           return
     end subroutine la_wsycon
#endif

     !> CSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_csyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,rwork,info)
        use la_constants_sp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(sp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(sp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('CSYRFS',-info)
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
              call la_csymv(uplo,n,-cone,a,lda,x(1,j),1,cone,work,1)
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
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + s
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
                 call la_csytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
                    ! multiply by diag(w)*inv(a**t).
                    call la_csytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_csytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
     end subroutine la_csyrfs
     !> ZSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_zsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,rwork,info)
        use la_constants_dp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(dp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(dp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('ZSYRFS',-info)
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
              call la_zsymv(uplo,n,-cone,a,lda,x(1,j),1,cone,work,1)
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
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + s
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
                 call la_zsytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
                    ! multiply by diag(w)*inv(a**t).
                    call la_zsytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_zsytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
     end subroutine la_zsyrfs
#ifdef LA_WITH_XDP
     !> YSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_ysyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,rwork,info)
        use la_constants_xdp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(xdp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(xdp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(xdp),intent(out) :: work(*)
           complex(xdp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
           real(xdp) :: eps,lstres,s,safe1,safe2,safmin,xk
           complex(xdp) :: zdum
           ! Local Arrays
           integer(ilp) :: isave(3)
           ! Intrinsic Functions
           intrinsic :: abs,real,aimag,max
           ! Statement Functions
           real(xdp) :: cabs1
           ! Statement Function Definitions
           cabs1(zdum) = abs(real(zdum,KIND=xdp)) + abs(aimag(zdum))
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('YSYRFS',-info)
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
           eps = la_xlamch('EPSILON')
           safmin = la_xlamch('SAFE MINIMUM')
           safe1 = nz*safmin
           safe2 = safe1/eps
           ! do for each right hand side
           loop_140: do j = 1,nrhs
              count = 1
              lstres = three
              20 continue
              ! loop until stopping criterion is satisfied.
              ! compute residual r = b - a * x
              call la_ycopy(n,b(1,j),1,work,1)
              call la_ysymv(uplo,n,-cone,a,lda,x(1,j),1,cone,work,1)
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
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + s
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
                 call la_ysytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
                 call la_yaxpy(n,cone,work,1,x(1,j),1)
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
              ! use la_ylacn2 to estimate the infinity-norm of the matrix
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
              call la_ylacn2(n,work(n + 1),work,ferr(j),kase,isave)
              if (kase /= 0) then
                 if (kase == 1) then
                    ! multiply by diag(w)*inv(a**t).
                    call la_ysytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_ysytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
     end subroutine la_ysyrfs
#endif
#ifdef LA_WITH_QP
     !> WSYRFS: improves the computed solution to a system of linear
     !> equations when the coefficient matrix is symmetric indefinite, and
     !> provides error bounds and backward error estimates for the solution.

     pure subroutine la_wsyrfs(uplo,n,nrhs,a,lda,af,ldaf,ipiv,b,ldb,x,ldx,ferr, &
               berr,work,rwork,info)
        use la_constants_qp,only:zero,two,three,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: uplo
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldaf,ldb,ldx,n,nrhs
           ! Array Arguments
           integer(ilp),intent(in) :: ipiv(*)
           real(qp),intent(out) :: berr(*),ferr(*),rwork(*)
           complex(qp),intent(in) :: a(lda,*),af(ldaf,*),b(ldb,*)
           complex(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: x(ldx,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: itmax = 5

           ! Local Scalars
           logical(lk) :: upper
           integer(ilp) :: count,i,j,k,kase,nz
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
           else if (lda < max(1,n)) then
              info = -5
           else if (ldaf < max(1,n)) then
              info = -7
           else if (ldb < max(1,n)) then
              info = -10
           else if (ldx < max(1,n)) then
              info = -12
           end if
           if (info /= 0) then
              call la_xerbla('WSYRFS',-info)
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
              call la_wsymv(uplo,n,-cone,a,lda,x(1,j),1,cone,work,1)
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
              if (upper) then
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    do i = 1,k - 1
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk + s
                 end do
              else
                 do k = 1,n
                    s = zero
                    xk = cabs1(x(k,j))
                    rwork(k) = rwork(k) + cabs1(a(k,k))*xk
                    do i = k + 1,n
                       rwork(i) = rwork(i) + cabs1(a(i,k))*xk
                       s = s + cabs1(a(i,k))*cabs1(x(i,j))
                    end do
                    rwork(k) = rwork(k) + s
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
                 call la_wsytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
                    ! multiply by diag(w)*inv(a**t).
                    call la_wsytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                 else if (kase == 2) then
                    ! multiply by inv(a)*diag(w).
                    do i = 1,n
                       work(i) = rwork(i)*work(i)
                    end do
                    call la_wsytrs(uplo,n,1,af,ldaf,ipiv,work,n,info)
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
     end subroutine la_wsyrfs
#endif

end module la_lapack_solve_ldl_comp
