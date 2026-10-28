!> SVD components: bidiagonal reduction and its orthogonal factors, Jacobi sweeps, generalized SVD
module la_lapack_svd_comp
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level3_gen
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_l1
     use la_lapack_blas_like_l2
     use la_lapack_givens_jacobi_rot
     use la_lapack_householder_reflectors
     use la_lapack_orthogonal_factors_ql
     use la_lapack_orthogonal_factors_qr
     use la_lapack_svd_comp2
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sgbbrd
     public :: la_sgebd2
     public :: la_sgsvj0
     public :: la_sgsvj1
     public :: la_stgsja
     public :: la_sgebrd
     public :: la_sorgbr
     public :: la_sormbr
     public :: la_dgbbrd
     public :: la_dgebd2
     public :: la_dgsvj0
     public :: la_dgsvj1
     public :: la_dtgsja
     public :: la_dgebrd
     public :: la_dorgbr
     public :: la_dormbr
#ifdef LA_WITH_XDP
     public :: la_xgbbrd
     public :: la_xgebd2
     public :: la_xgsvj0
     public :: la_xgsvj1
     public :: la_xtgsja
     public :: la_xgebrd
     public :: la_xorgbr
     public :: la_xormbr
#endif
#ifdef LA_WITH_QP
     public :: la_qgbbrd
     public :: la_qgebd2
     public :: la_qgsvj0
     public :: la_qgsvj1
     public :: la_qtgsja
     public :: la_qgebrd
     public :: la_qorgbr
     public :: la_qormbr
#endif
     public :: la_cgebd2
     public :: la_ctgsja
     public :: la_cgbbrd
     public :: la_cgebrd
     public :: la_cungbr
     public :: la_cunmbr
     public :: la_cgsvj0
     public :: la_cgsvj1
     public :: la_zgebd2
     public :: la_ztgsja
     public :: la_zgbbrd
     public :: la_zgebrd
     public :: la_zungbr
     public :: la_zunmbr
     public :: la_zgsvj0
     public :: la_zgsvj1
#ifdef LA_WITH_XDP
     public :: la_ygebd2
     public :: la_ytgsja
     public :: la_ygbbrd
     public :: la_ygebrd
     public :: la_yungbr
     public :: la_yunmbr
     public :: la_ygsvj0
     public :: la_ygsvj1
#endif
#ifdef LA_WITH_QP
     public :: la_wgebd2
     public :: la_wtgsja
     public :: la_wgbbrd
     public :: la_wgebrd
     public :: la_wungbr
     public :: la_wunmbr
     public :: la_wgsvj0
     public :: la_wgsvj1
#endif

     contains

     !> SGBBRD: reduces a real general m-by-n band matrix A to upper
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> The routine computes B, and optionally forms Q or P**T, or computes
     !> Q**T*C for a given matrix C.

     pure subroutine la_sgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(sp),intent(inout) :: ab(ldab,*),c(ldc,*)
           real(sp),intent(out) :: d(*),e(*),pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mn, &
                      mu,mu0,nr,nrt
           real(sp) :: ra,rb,rc,rs
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('SGBBRD',-info)
              return
           end if
           ! initialize q and p**t to the unit matrix, if needed
           if (wantq) call la_slaset('FULL',m,m,zero,one,q,ldq)
           if (wantpt) call la_slaset('FULL',n,n,zero,one,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the sines of the plane rotations are stored in work(1:max(m,n))
              ! and the cosines in work(max(m,n)+1:2*max(m,n)).
              mn = max(m,n)
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_slargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              work(mn + j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_slartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,work(mn + j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_slartg(ab(ku + ml - 1,i),ab(ku + ml,i),work(mn + i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_srot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,work(mn + i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_srot(m,q(1,j - 1),1,q(1,j),1,work(mn + j),work(j &
                                    ))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_srot(ncc,c(j - 1,1),ldc,c(j,1),ldc,work(mn + j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = work(mn + j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_slargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               work(mn + j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_slartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,work(mn + j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_slartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),work( &
                                    mn + i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_srot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,work(mn + i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**t
                       do j = j1,j2,kb1
                          call la_srot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,work( &
                                    mn + j + kun),work(j + kun))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = work(mn + j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, storing diagonal elements in d
              ! and off-diagonal elements in e
              do i = 1,min(m - 1,n)
                 call la_slartg(ab(1,i),ab(2,i),rc,rs,ra)
                 d(i) = ra
                 if (i < n) then
                    e(i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_srot(m,q(1,i),1,q(1,i + 1),1,rc,rs)
                 if (wantc) call la_srot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
              if (m <= n) d(m) = ab(1,m)
           else if (ku > 0) then
              ! a has been reduced to upper bidiagonal form
              if (m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right, storing diagonal elements in d and off-diagonal
                 ! elements in e
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_slartg(ab(ku + 1,i),rb,rc,rs,ra)
                    d(i) = ra
                    if (i > 1) then
                       rb = -rs*ab(ku,i)
                       e(i - 1) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_srot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc,rs)

                 end do
              else
                 ! copy off-diagonal elements to e and diagonal elements to d
                 do i = 1,minmn - 1
                    e(i) = ab(ku,i + 1)
                 end do
                 do i = 1,minmn
                    d(i) = ab(ku + 1,i)
                 end do
              end if
           else
              ! a is diagonal. set elements of e to zero and copy diagonal
              ! elements to d.
              do i = 1,minmn - 1
                 e(i) = zero
              end do
              do i = 1,minmn
                 d(i) = ab(1,i)
              end do
           end if
           return
     end subroutine la_sgbbrd
     !> DGBBRD: reduces a real general m-by-n band matrix A to upper
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> The routine computes B, and optionally forms Q or P**T, or computes
     !> Q**T*C for a given matrix C.

     pure subroutine la_dgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(dp),intent(inout) :: ab(ldab,*),c(ldc,*)
           real(dp),intent(out) :: d(*),e(*),pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mn, &
                      mu,mu0,nr,nrt
           real(dp) :: ra,rb,rc,rs
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('DGBBRD',-info)
              return
           end if
           ! initialize q and p**t to the unit matrix, if needed
           if (wantq) call la_dlaset('FULL',m,m,zero,one,q,ldq)
           if (wantpt) call la_dlaset('FULL',n,n,zero,one,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the sines of the plane rotations are stored in work(1:max(m,n))
              ! and the cosines in work(max(m,n)+1:2*max(m,n)).
              mn = max(m,n)
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_dlargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              work(mn + j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_dlartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,work(mn + j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_dlartg(ab(ku + ml - 1,i),ab(ku + ml,i),work(mn + i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_drot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,work(mn + i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_drot(m,q(1,j - 1),1,q(1,j),1,work(mn + j),work(j &
                                    ))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_drot(ncc,c(j - 1,1),ldc,c(j,1),ldc,work(mn + j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = work(mn + j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_dlargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               work(mn + j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_dlartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,work(mn + j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_dlartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),work( &
                                    mn + i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_drot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,work(mn + i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**t
                       do j = j1,j2,kb1
                          call la_drot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,work( &
                                    mn + j + kun),work(j + kun))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = work(mn + j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, storing diagonal elements in d
              ! and off-diagonal elements in e
              do i = 1,min(m - 1,n)
                 call la_dlartg(ab(1,i),ab(2,i),rc,rs,ra)
                 d(i) = ra
                 if (i < n) then
                    e(i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_drot(m,q(1,i),1,q(1,i + 1),1,rc,rs)
                 if (wantc) call la_drot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
              if (m <= n) d(m) = ab(1,m)
           else if (ku > 0) then
              ! a has been reduced to upper bidiagonal form
              if (m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right, storing diagonal elements in d and off-diagonal
                 ! elements in e
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_dlartg(ab(ku + 1,i),rb,rc,rs,ra)
                    d(i) = ra
                    if (i > 1) then
                       rb = -rs*ab(ku,i)
                       e(i - 1) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_drot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc,rs)

                 end do
              else
                 ! copy off-diagonal elements to e and diagonal elements to d
                 do i = 1,minmn - 1
                    e(i) = ab(ku,i + 1)
                 end do
                 do i = 1,minmn
                    d(i) = ab(ku + 1,i)
                 end do
              end if
           else
              ! a is diagonal. set elements of e to zero and copy diagonal
              ! elements to d.
              do i = 1,minmn - 1
                 e(i) = zero
              end do
              do i = 1,minmn
                 d(i) = ab(1,i)
              end do
           end if
           return
     end subroutine la_dgbbrd
#ifdef LA_WITH_XDP
     !> XGBBRD: reduces a real general m-by-n band matrix A to upper
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> The routine computes B, and optionally forms Q or P**T, or computes
     !> Q**T*C for a given matrix C.

     pure subroutine la_xgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(xdp),intent(inout) :: ab(ldab,*),c(ldc,*)
           real(xdp),intent(out) :: d(*),e(*),pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mn, &
                      mu,mu0,nr,nrt
           real(xdp) :: ra,rb,rc,rs
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('XGBBRD',-info)
              return
           end if
           ! initialize q and p**t to the unit matrix, if needed
           if (wantq) call la_xlaset('FULL',m,m,zero,one,q,ldq)
           if (wantpt) call la_xlaset('FULL',n,n,zero,one,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the sines of the plane rotations are stored in work(1:max(m,n))
              ! and the cosines in work(max(m,n)+1:2*max(m,n)).
              mn = max(m,n)
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_xlargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              work(mn + j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_xlartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,work(mn + j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_xlartg(ab(ku + ml - 1,i),ab(ku + ml,i),work(mn + i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_xrot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,work(mn + i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_xrot(m,q(1,j - 1),1,q(1,j),1,work(mn + j),work(j &
                                    ))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_xrot(ncc,c(j - 1,1),ldc,c(j,1),ldc,work(mn + j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = work(mn + j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_xlargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               work(mn + j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_xlartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,work(mn + j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_xlartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),work( &
                                    mn + i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_xrot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,work(mn + i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**t
                       do j = j1,j2,kb1
                          call la_xrot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,work( &
                                    mn + j + kun),work(j + kun))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = work(mn + j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, storing diagonal elements in d
              ! and off-diagonal elements in e
              do i = 1,min(m - 1,n)
                 call la_xlartg(ab(1,i),ab(2,i),rc,rs,ra)
                 d(i) = ra
                 if (i < n) then
                    e(i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_xrot(m,q(1,i),1,q(1,i + 1),1,rc,rs)
                 if (wantc) call la_xrot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
              if (m <= n) d(m) = ab(1,m)
           else if (ku > 0) then
              ! a has been reduced to upper bidiagonal form
              if (m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right, storing diagonal elements in d and off-diagonal
                 ! elements in e
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_xlartg(ab(ku + 1,i),rb,rc,rs,ra)
                    d(i) = ra
                    if (i > 1) then
                       rb = -rs*ab(ku,i)
                       e(i - 1) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_xrot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc,rs)

                 end do
              else
                 ! copy off-diagonal elements to e and diagonal elements to d
                 do i = 1,minmn - 1
                    e(i) = ab(ku,i + 1)
                 end do
                 do i = 1,minmn
                    d(i) = ab(ku + 1,i)
                 end do
              end if
           else
              ! a is diagonal. set elements of e to zero and copy diagonal
              ! elements to d.
              do i = 1,minmn - 1
                 e(i) = zero
              end do
              do i = 1,minmn
                 d(i) = ab(1,i)
              end do
           end if
           return
     end subroutine la_xgbbrd
#endif
#ifdef LA_WITH_QP
     !> QGBBRD: reduces a real general m-by-n band matrix A to upper
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> The routine computes B, and optionally forms Q or P**T, or computes
     !> Q**T*C for a given matrix C.

     pure subroutine la_qgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(qp),intent(inout) :: ab(ldab,*),c(ldc,*)
           real(qp),intent(out) :: d(*),e(*),pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mn, &
                      mu,mu0,nr,nrt
           real(qp) :: ra,rb,rc,rs
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('QGBBRD',-info)
              return
           end if
           ! initialize q and p**t to the unit matrix, if needed
           if (wantq) call la_qlaset('FULL',m,m,zero,one,q,ldq)
           if (wantpt) call la_qlaset('FULL',n,n,zero,one,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the sines of the plane rotations are stored in work(1:max(m,n))
              ! and the cosines in work(max(m,n)+1:2*max(m,n)).
              mn = max(m,n)
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_qlargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              work(mn + j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_qlartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,work(mn + j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_qlartg(ab(ku + ml - 1,i),ab(ku + ml,i),work(mn + i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_qrot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,work(mn + i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_qrot(m,q(1,j - 1),1,q(1,j),1,work(mn + j),work(j &
                                    ))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_qrot(ncc,c(j - 1,1),ldc,c(j,1),ldc,work(mn + j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = work(mn + j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_qlargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               work(mn + j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_qlartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,work(mn + j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_qlartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),work( &
                                    mn + i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_qrot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,work(mn + i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**t
                       do j = j1,j2,kb1
                          call la_qrot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,work( &
                                    mn + j + kun),work(j + kun))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = work(mn + j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, storing diagonal elements in d
              ! and off-diagonal elements in e
              do i = 1,min(m - 1,n)
                 call la_qlartg(ab(1,i),ab(2,i),rc,rs,ra)
                 d(i) = ra
                 if (i < n) then
                    e(i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_qrot(m,q(1,i),1,q(1,i + 1),1,rc,rs)
                 if (wantc) call la_qrot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
              if (m <= n) d(m) = ab(1,m)
           else if (ku > 0) then
              ! a has been reduced to upper bidiagonal form
              if (m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right, storing diagonal elements in d and off-diagonal
                 ! elements in e
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_qlartg(ab(ku + 1,i),rb,rc,rs,ra)
                    d(i) = ra
                    if (i > 1) then
                       rb = -rs*ab(ku,i)
                       e(i - 1) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_qrot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc,rs)

                 end do
              else
                 ! copy off-diagonal elements to e and diagonal elements to d
                 do i = 1,minmn - 1
                    e(i) = ab(ku,i + 1)
                 end do
                 do i = 1,minmn
                    d(i) = ab(ku + 1,i)
                 end do
              end if
           else
              ! a is diagonal. set elements of e to zero and copy diagonal
              ! elements to d.
              do i = 1,minmn - 1
                 e(i) = zero
              end do
              do i = 1,minmn
                 d(i) = ab(1,i)
              end do
           end if
           return
     end subroutine la_qgbbrd
#endif

     !> SGEBD2: reduces a real general m by n matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_sgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('SGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 call la_slarfg(m - i + 1,a(i,i),a(min(i + 1,m),i),1,tauq(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply h(i) to a(i:m,i+1:n) from the left
                 if (i < n) call la_slarf('LEFT',m - i + 1,n - i,a(i,i),1,tauq(i),a(i,i + &
                           1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_slarfg(n - i,a(i,i + 1),a(i,min(i + 2,n)),lda,taup(i))

                    e(i) = a(i,i + 1)
                    a(i,i + 1) = one
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_slarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = zero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_slarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,taup(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_slarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    call la_slarfg(m - i,a(i + 1,i),a(min(i + 2,m),i),1,tauq(i))

                    e(i) = a(i + 1,i)
                    a(i + 1,i) = one
                    ! apply h(i) to a(i+1:m,i+1:n) from the left
                    call la_slarf('LEFT',m - i,n - i,a(i + 1,i),1,tauq(i),a(i + 1,i + 1), &
                              lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = zero
                 end if
              end do
           end if
           return
     end subroutine la_sgebd2
     !> DGEBD2: reduces a real general m by n matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_dgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('DGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 call la_dlarfg(m - i + 1,a(i,i),a(min(i + 1,m),i),1,tauq(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply h(i) to a(i:m,i+1:n) from the left
                 if (i < n) call la_dlarf('LEFT',m - i + 1,n - i,a(i,i),1,tauq(i),a(i,i + &
                           1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_dlarfg(n - i,a(i,i + 1),a(i,min(i + 2,n)),lda,taup(i))

                    e(i) = a(i,i + 1)
                    a(i,i + 1) = one
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_dlarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = zero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_dlarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,taup(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_dlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    call la_dlarfg(m - i,a(i + 1,i),a(min(i + 2,m),i),1,tauq(i))

                    e(i) = a(i + 1,i)
                    a(i + 1,i) = one
                    ! apply h(i) to a(i+1:m,i+1:n) from the left
                    call la_dlarf('LEFT',m - i,n - i,a(i + 1,i),1,tauq(i),a(i + 1,i + 1), &
                              lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = zero
                 end if
              end do
           end if
           return
     end subroutine la_dgebd2
#ifdef LA_WITH_XDP
     !> XGEBD2: reduces a real general m by n matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_xgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('XGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 call la_xlarfg(m - i + 1,a(i,i),a(min(i + 1,m),i),1,tauq(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply h(i) to a(i:m,i+1:n) from the left
                 if (i < n) call la_xlarf('LEFT',m - i + 1,n - i,a(i,i),1,tauq(i),a(i,i + &
                           1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_xlarfg(n - i,a(i,i + 1),a(i,min(i + 2,n)),lda,taup(i))

                    e(i) = a(i,i + 1)
                    a(i,i + 1) = one
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_xlarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = zero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_xlarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,taup(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_xlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    call la_xlarfg(m - i,a(i + 1,i),a(min(i + 2,m),i),1,tauq(i))

                    e(i) = a(i + 1,i)
                    a(i + 1,i) = one
                    ! apply h(i) to a(i+1:m,i+1:n) from the left
                    call la_xlarf('LEFT',m - i,n - i,a(i + 1,i),1,tauq(i),a(i + 1,i + 1), &
                              lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = zero
                 end if
              end do
           end if
           return
     end subroutine la_xgebd2
#endif
#ifdef LA_WITH_QP
     !> QGEBD2: reduces a real general m by n matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_qgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('QGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 call la_qlarfg(m - i + 1,a(i,i),a(min(i + 1,m),i),1,tauq(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply h(i) to a(i:m,i+1:n) from the left
                 if (i < n) call la_qlarf('LEFT',m - i + 1,n - i,a(i,i),1,tauq(i),a(i,i + &
                           1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_qlarfg(n - i,a(i,i + 1),a(i,min(i + 2,n)),lda,taup(i))

                    e(i) = a(i,i + 1)
                    a(i,i + 1) = one
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_qlarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = zero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_qlarfg(n - i + 1,a(i,i),a(i,min(i + 1,n)),lda,taup(i))

                 d(i) = a(i,i)
                 a(i,i) = one
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_qlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    call la_qlarfg(m - i,a(i + 1,i),a(min(i + 2,m),i),1,tauq(i))

                    e(i) = a(i + 1,i)
                    a(i + 1,i) = one
                    ! apply h(i) to a(i+1:m,i+1:n) from the left
                    call la_qlarf('LEFT',m - i,n - i,a(i + 1,i),1,tauq(i),a(i + 1,i + 1), &
                              lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = zero
                 end if
              end do
           end if
           return
     end subroutine la_qgebd2
#endif

     !> SGSVJ0: is called from SGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as SGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_sgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_sp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(sp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),sva(n),d(n),v(ldv,*)
           real(sp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(sp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(sp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,float,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('SGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! .. row-cyclic jacobi svd algorithm with column pivoting ..
           emptsw = (n*(n - 1))/2
           notrot = 0
           fastr(1) = zero
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_sgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_sgesvj. for sweeps i=1:swband the procedure
           ! ......
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           swband = 0
           pskipped = 0
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_isamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_sswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_sswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          temp1 = d(p)
                          d(p) = d(q)
                          d(q) = temp1
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! some blas implementations compute la_snrm2(m,a(1,p),1)
              ! as sqrt(la_sdot(m,a(1,p),1,a(1,p),1)), which may result in
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_snrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_snrm2 is available, the if-then-else
              ! below should read "aapp = la_snrm2( m, a(1,p), 1 ) * d(p)".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_snrm2(m,a(1,p),1)*d(p)
                          else
                             temp1 = zero
                             aapp = one
                             call la_slassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)*d(p)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_sdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_scopy(m,a(1,p),1,work,1)
                                      call la_slascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_sdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_sdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_scopy(m,a(1,q),1,work,1)
                                      call la_slascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_sdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                 ! Rotate
                 ! rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_srotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_srotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_srotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_srotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_saxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_saxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                               if (rsvec) then
                                                  call la_saxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_saxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_saxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_saxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                               if (rsvec) then
                                                  call la_saxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_saxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_saxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_saxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_saxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_saxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_saxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_saxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_saxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_saxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_scopy(m,a(1,p),1,work,1)
                                      call la_slascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_slascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      temp1 = -aapq*d(p)/d(q)
                                      call la_saxpy(m,temp1,work,1,a(1,q),1)
                                      call la_slascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_snrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_slassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_snrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_slassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! Safe Gram Matrix Computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_sdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_scopy(m,a(1,p),1,work,1)
                                      call la_slascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_sdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_sdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_scopy(m,a(1,q),1,work,1)
                                      call la_slascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_sdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_srotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_srotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_srotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_srotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_saxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_saxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_saxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_saxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_saxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_saxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_saxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_saxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_saxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_saxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_saxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_saxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_saxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_saxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_saxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_saxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_scopy(m,a(1,p),1,work,1)
                                         call la_slascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_slascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_saxpy(m,temp1,work,1,a(1,q),1)

                                         call la_slascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_scopy(m,a(1,q),1,work,1)
                                         call la_slascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_slascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_saxpy(m,temp1,work,1,a(1,p),1)

                                         call la_slascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_snrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_slassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_snrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_slassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_snrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_slassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=sp)*tol) .and. (real(n,KIND=sp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of iterations.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given tolerance, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d.
           do p = 1,n - 1
              q = la_isamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_sswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_sswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_sgsvj0
     !> DGSVJ0: is called from DGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as DGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_dgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_dp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(dp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),sva(n),d(n),v(ldv,*)
           real(dp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(dp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(dp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('DGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! -#- row-cyclic jacobi svd algorithm with column pivoting -#-
           emptsw = (n*(n - 1))/2
           notrot = 0
           fastr(1) = zero
           ! -#- row-cyclic pivot strategy with de rijk's pivoting -#-
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_sgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_sgesvj. for sweeps i=1:swband the procedure
           ! ......
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           swband = 0
           pskipped = 0
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_idamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_dswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_dswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          temp1 = d(p)
                          d(p) = d(q)
                          d(q) = temp1
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! some blas implementations compute la_dnrm2(m,a(1,p),1)
              ! as sqrt(la_ddot(m,a(1,p),1,a(1,p),1)), which may result in
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_dnrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_dnrm2 is available, the if-then-else
              ! below should read "aapp = la_dnrm2( m, a(1,p), 1 ) * d(p)".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_dnrm2(m,a(1,p),1)*d(p)
                          else
                             temp1 = zero
                             aapp = one
                             call la_dlassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)*d(p)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_ddot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_dcopy(m,a(1,p),1,work,1)
                                      call la_dlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_ddot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_ddot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_dcopy(m,a(1,q),1,work,1)
                                      call la_dlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_ddot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                 ! Rotate
                 ! rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_drotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_drotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_drotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_drotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_daxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_daxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                               if (rsvec) then
                                                  call la_daxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_daxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_daxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_daxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                               if (rsvec) then
                                                  call la_daxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_daxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_daxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_daxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_daxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_daxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_daxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_daxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_daxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_daxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_dcopy(m,a(1,p),1,work,1)
                                      call la_dlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_dlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      temp1 = -aapq*d(p)/d(q)
                                      call la_daxpy(m,temp1,work,1,a(1,q),1)
                                      call la_dlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_dnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_dlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_dnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_dlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! -#- m x 2 jacobi svd -#-
              ! -#- safe gram matrix computation -#-
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_ddot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_dcopy(m,a(1,p),1,work,1)
                                      call la_dlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_ddot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_ddot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_dcopy(m,a(1,q),1,work,1)
                                      call la_dlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_ddot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_drotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_drotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_drotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_drotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_daxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_daxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_daxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_daxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_daxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_daxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_daxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_daxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_daxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_daxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_daxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_daxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_daxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_daxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_daxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_daxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_dcopy(m,a(1,p),1,work,1)
                                         call la_dlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_dlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_daxpy(m,temp1,work,1,a(1,q),1)

                                         call la_dlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_dcopy(m,a(1,q),1,work,1)
                                         call la_dlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_dlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_daxpy(m,temp1,work,1,a(1,p),1)

                                         call la_dlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_dnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_dlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_dnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_dlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_dnrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_dlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=dp)*tol) .and. (real(n,KIND=dp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of iterations.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given tolerance, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d.
           do p = 1,n - 1
              q = la_idamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_dswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_dswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_dgsvj0
#ifdef LA_WITH_XDP
     !> XGSVJ0: is called from XGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as XGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_xgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_xdp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(xdp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),sva(n),d(n),v(ldv,*)
           real(xdp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(xdp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('XGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! -#- row-cyclic jacobi svd algorithm with column pivoting -#-
           emptsw = (n*(n - 1))/2
           notrot = 0
           fastr(1) = zero
           ! -#- row-cyclic pivot strategy with de rijk's pivoting -#-
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_dgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_dgesvj. for sweeps i=1:swband the procedure
           ! ......
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           swband = 0
           pskipped = 0
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_ixamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_xswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_xswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          temp1 = d(p)
                          d(p) = d(q)
                          d(q) = temp1
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! some blas implementations compute la_xnrm2(m,a(1,p),1)
              ! as sqrt(la_xdot(m,a(1,p),1,a(1,p),1)), which may result in
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_xnrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_xnrm2 is available, the if-then-else
              ! below should read "aapp = la_xnrm2( m, a(1,p), 1 ) * d(p)".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_xnrm2(m,a(1,p),1)*d(p)
                          else
                             temp1 = zero
                             aapp = one
                             call la_xlassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)*d(p)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_xdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_xcopy(m,a(1,p),1,work,1)
                                      call la_xlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_xdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_xdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_xcopy(m,a(1,q),1,work,1)
                                      call la_xlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_xdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                 ! Rotate
                 ! rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_xrotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_xrotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_xrotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_xrotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_xaxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_xaxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                               if (rsvec) then
                                                  call la_xaxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_xaxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_xaxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_xaxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                               if (rsvec) then
                                                  call la_xaxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_xaxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_xaxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_xaxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_xaxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_xaxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_xaxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_xaxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_xaxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_xaxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_xcopy(m,a(1,p),1,work,1)
                                      call la_xlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_xlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      temp1 = -aapq*d(p)/d(q)
                                      call la_xaxpy(m,temp1,work,1,a(1,q),1)
                                      call la_xlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_xnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_xlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_xnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_xlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! -#- m x 2 jacobi svd -#-
              ! -#- safe gram matrix computation -#-
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_xdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_xcopy(m,a(1,p),1,work,1)
                                      call la_xlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_xdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_xdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_xcopy(m,a(1,q),1,work,1)
                                      call la_xlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_xdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_xrotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_xrotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_xrotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_xrotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_xaxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_xaxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_xaxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_xaxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_xaxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_xaxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_xaxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_xaxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_xaxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_xaxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_xaxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_xaxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_xaxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_xaxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_xaxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_xaxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_xcopy(m,a(1,p),1,work,1)
                                         call la_xlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_xlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_xaxpy(m,temp1,work,1,a(1,q),1)

                                         call la_xlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_xcopy(m,a(1,q),1,work,1)
                                         call la_xlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_xlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_xaxpy(m,temp1,work,1,a(1,p),1)

                                         call la_xlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_xnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_xlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_xnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_xlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_xnrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_xlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=xdp)*tol) .and. (real(n,KIND=xdp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of iterations.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given tolerance, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d.
           do p = 1,n - 1
              q = la_ixamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_xswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_xswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_xgsvj0
#endif
#ifdef LA_WITH_QP
     !> QGSVJ0: is called from QGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as QGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_qgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_qp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(qp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),sva(n),d(n),v(ldv,*)
           real(qp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(qp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(qp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('QGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! -#- row-cyclic jacobi svd algorithm with column pivoting -#-
           emptsw = (n*(n - 1))/2
           notrot = 0
           fastr(1) = zero
           ! -#- row-cyclic pivot strategy with de rijk's pivoting -#-
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_dgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_dgesvj. for sweeps i=1:swband the procedure
           ! ......
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           swband = 0
           pskipped = 0
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_iqamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_qswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_qswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          temp1 = d(p)
                          d(p) = d(q)
                          d(q) = temp1
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! some blas implementations compute la_qnrm2(m,a(1,p),1)
              ! as sqrt(la_qdot(m,a(1,p),1,a(1,p),1)), which may result in
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_qnrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_qnrm2 is available, the if-then-else
              ! below should read "aapp = la_qnrm2( m, a(1,p), 1 ) * d(p)".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_qnrm2(m,a(1,p),1)*d(p)
                          else
                             temp1 = zero
                             aapp = one
                             call la_qlassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)*d(p)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_qdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_qcopy(m,a(1,p),1,work,1)
                                      call la_qlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_qdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_qdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_qcopy(m,a(1,q),1,work,1)
                                      call la_qlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_qdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                 ! Rotate
                 ! rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_qrotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_qrotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_qrotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_qrotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_qaxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_qaxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                               if (rsvec) then
                                                  call la_qaxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_qaxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_qaxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_qaxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                               if (rsvec) then
                                                  call la_qaxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_qaxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_qaxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_qaxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_qaxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_qaxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_qaxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_qaxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_qaxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_qaxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_qcopy(m,a(1,p),1,work,1)
                                      call la_qlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_qlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      temp1 = -aapq*d(p)/d(q)
                                      call la_qaxpy(m,temp1,work,1,a(1,q),1)
                                      call la_qlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_qnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_qlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_qnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_qlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! -#- m x 2 jacobi svd -#-
              ! -#- safe gram matrix computation -#-
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_qdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_qcopy(m,a(1,p),1,work,1)
                                      call la_qlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_qdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_qdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_qcopy(m,a(1,q),1,work,1)
                                      call la_qlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_qdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_qrotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_qrotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_qrotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_qrotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_qaxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_qaxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_qaxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_qaxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_qaxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_qaxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_qaxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_qaxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_qaxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_qaxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_qaxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_qaxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_qaxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_qaxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_qaxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_qaxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_qcopy(m,a(1,p),1,work,1)
                                         call la_qlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_qlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_qaxpy(m,temp1,work,1,a(1,q),1)

                                         call la_qlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_qcopy(m,a(1,q),1,work,1)
                                         call la_qlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_qlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_qaxpy(m,temp1,work,1,a(1,p),1)

                                         call la_qlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_qnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_qlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_qnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_qlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_qnrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_qlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=qp)*tol) .and. (real(n,KIND=qp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of iterations.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given tolerance, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d.
           do p = 1,n - 1
              q = la_iqamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_qswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_qswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_qgsvj0
#endif

     !> SGSVJ1: is called from SGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as SGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> SGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_sgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_sp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),d(n),sva(n),v(ldv,*)
           real(sp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(sp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,large,mxaapq, &
           mxsinj,rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta, &
                     thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(sp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,float,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('SGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           large = big/sqrt(real(m*n,KIND=sp))
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           fastr(1) = zero
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_sgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_sgesvj.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = 1,nblc
                    jgl = n1 + (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! Safe Gram Matrix Computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_sdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_scopy(m,a(1,p),1,work,1)
                                      call la_slascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_sdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_sdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_scopy(m,a(1,q),1,work,1)
                                      call la_slascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_sdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_srotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_srotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_srotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_srotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_saxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_saxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_saxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_saxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_saxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_saxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_saxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_saxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_saxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_saxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_saxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_saxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_saxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_saxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_saxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_saxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_scopy(m,a(1,p),1,work,1)
                                         call la_slascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_slascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_saxpy(m,temp1,work,1,a(1,q),1)

                                         call la_slascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_scopy(m,a(1,q),1,work,1)
                                         call la_slascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_slascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_saxpy(m,temp1,work,1,a(1,p),1)

                                         call la_slascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_snrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_slassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_snrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_slassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                 ! skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
            ! if ( notrot >= emptsw )  go to 2011
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
      ! **      if ( notrot >= emptsw )  go to 2011
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **   if ( notrot >= emptsw ) go to 1994
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_snrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_slassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=sp)*tol) .and. (real(n,KIND=sp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of sweeps.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given threshold, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d
           do p = 1,n - 1
              q = la_isamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_sswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_sswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_sgsvj1
     !> DGSVJ1: is called from DGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as DGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> DGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_dgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_dp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),d(n),sva(n),v(ldv,*)
           real(dp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(dp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,large,mxaapq, &
           mxsinj,rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta, &
                     thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(dp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('DGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           large = big/sqrt(real(m*n,KIND=dp))
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           fastr(1) = zero
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_sgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_sgesvj.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = 1,nblc
                    jgl = n1 + (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! Safe Gram Matrix Computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_ddot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_dcopy(m,a(1,p),1,work,1)
                                      call la_dlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_ddot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_ddot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_dcopy(m,a(1,q),1,work,1)
                                      call la_dlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_ddot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_drotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_drotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_drotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_drotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_daxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_daxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_daxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_daxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_daxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_daxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_daxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_daxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_daxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_daxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_daxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_daxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_daxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_daxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_daxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_daxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_dcopy(m,a(1,p),1,work,1)
                                         call la_dlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_dlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_daxpy(m,temp1,work,1,a(1,q),1)

                                         call la_dlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_dcopy(m,a(1,q),1,work,1)
                                         call la_dlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_dlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_daxpy(m,temp1,work,1,a(1,p),1)

                                         call la_dlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_dnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_dlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_dnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_dlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                 ! skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
            ! if ( notrot >= emptsw )  go to 2011
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
      ! **      if ( notrot >= emptsw )  go to 2011
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **   if ( notrot >= emptsw ) go to 1994
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_dnrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_dlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=dp)*tol) .and. (real(n,KIND=dp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of sweeps.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given threshold, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d
           do p = 1,n - 1
              q = la_idamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_dswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_dswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_dgsvj1
#ifdef LA_WITH_XDP
     !> XGSVJ1: is called from XGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as XGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> XGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_xgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_xdp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),d(n),sva(n),v(ldv,*)
           real(xdp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(xdp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,large,mxaapq, &
           mxsinj,rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta, &
                     thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(xdp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('XGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           large = big/sqrt(real(m*n,KIND=xdp))
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           fastr(1) = zero
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_dgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_dgesvj.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = 1,nblc
                    jgl = n1 + (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! Safe Gram Matrix Computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_xdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_xcopy(m,a(1,p),1,work,1)
                                      call la_xlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_xdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_xdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_xcopy(m,a(1,q),1,work,1)
                                      call la_xlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_xdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_xrotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_xrotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_xrotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_xrotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_xaxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_xaxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_xaxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_xaxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_xaxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_xaxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_xaxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_xaxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_xaxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_xaxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_xaxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_xaxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_xaxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_xaxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_xaxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_xaxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_xcopy(m,a(1,p),1,work,1)
                                         call la_xlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_xlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_xaxpy(m,temp1,work,1,a(1,q),1)

                                         call la_xlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_xcopy(m,a(1,q),1,work,1)
                                         call la_xlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_xlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_xaxpy(m,temp1,work,1,a(1,p),1)

                                         call la_xlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_xnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_xlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_xnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_xlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                 ! skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
            ! if ( notrot >= emptsw )  go to 2011
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
      ! **      if ( notrot >= emptsw )  go to 2011
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **   if ( notrot >= emptsw ) go to 1994
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_xnrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_xlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=xdp)*tol) .and. (real(n,KIND=xdp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of sweeps.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given threshold, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d
           do p = 1,n - 1
              q = la_ixamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_xswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_xswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_xgsvj1
#endif
#ifdef LA_WITH_QP
     !> QGSVJ1: is called from QGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as QGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> QGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_qgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_qp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),d(n),sva(n),v(ldv,*)
           real(qp),intent(out) :: work(lwork)
        ! =====================================================================

           ! Local Scalars
           real(qp) :: aapp,aapp0,aapq,aaqq,apoaq,aqoap,big,bigtheta,cs,large,mxaapq, &
           mxsinj,rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta, &
                     thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Local Arrays
           real(qp) :: fastr(5)
           ! Intrinsic Functions
           intrinsic :: abs,max,real,min,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('QGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           large = big/sqrt(real(m*n,KIND=qp))
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           fastr(1) = zero
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_dgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_dgesvj.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
      ! ........................................................
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = 1,nblc
                    jgl = n1 + (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! Safe Gram Matrix Computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_qdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_qcopy(m,a(1,p),1,work,1)
                                      call la_qlascl('G',0,0,aapp,d(p),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_qdot(m,work,1,a(1,q),1)*d(q)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_qdot(m,a(1,p),1,a(1,q),1)*d(p) &
                                                *d(q)/aaqq)/aapp
                                   else
                                      call la_qcopy(m,a(1,q),1,work,1)
                                      call la_qlascl('G',0,0,aaqq,d(q),m,1,work,lda, &
                                                 ierr)
                                      aapq = la_qdot(m,work,1,a(1,p),1)*d(p)/ &
                                                aapp
                                   end if
                                end if
                                mxaapq = max(mxaapq,abs(aapq))
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq) > tol) then
                                   notrot = 0
                 ! rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         fastr(3) = t*d(p)/d(q)
                                         fastr(4) = -t*d(q)/d(p)
                                         call la_qrotm(m,a(1,p),1,a(1,q),1,fastr)

                                         if (rsvec) call la_qrotm(mvl,v(1,p),1,v(1,q), &
                                                    1,fastr)
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq))
                                         apoaq = d(p)/d(q)
                                         aqoap = d(q)/d(p)
                                         if (d(p) >= one) then
                                            if (d(q) >= one) then
                                               fastr(3) = t*apoaq
                                               fastr(4) = -t*aqoap
                                               d(p) = d(p)*cs
                                               d(q) = d(q)*cs
                                               call la_qrotm(m,a(1,p),1,a(1,q),1, &
                                                         fastr)
                                               if (rsvec) call la_qrotm(mvl,v(1,p),1,v( &
                                                         1,q),1,fastr)
                                            else
                                               call la_qaxpy(m,-t*aqoap,a(1,q),1,a(1, &
                                                         p),1)
                                               call la_qaxpy(m,cs*sn*apoaq,a(1,p),1,a( &
                                                         1,q),1)
                                               if (rsvec) then
                                                  call la_qaxpy(mvl,-t*aqoap,v(1,q),1,v( &
                                                             1,p),1)
                                                  call la_qaxpy(mvl,cs*sn*apoaq,v(1,p),1, &
                                                            v(1,q),1)
                                               end if
                                               d(p) = d(p)*cs
                                               d(q) = d(q)/cs
                                            end if
                                         else
                                            if (d(q) >= one) then
                                               call la_qaxpy(m,t*apoaq,a(1,p),1,a(1,q &
                                                         ),1)
                                               call la_qaxpy(m,-cs*sn*aqoap,a(1,q),1,a( &
                                                         1,p),1)
                                               if (rsvec) then
                                                  call la_qaxpy(mvl,t*apoaq,v(1,p),1,v( &
                                                            1,q),1)
                                                  call la_qaxpy(mvl,-cs*sn*aqoap,v(1,q), &
                                                            1,v(1,p),1)
                                               end if
                                               d(p) = d(p)/cs
                                               d(q) = d(q)*cs
                                            else
                                               if (d(p) >= d(q)) then
                                                  call la_qaxpy(m,-t*aqoap,a(1,q),1,a( &
                                                            1,p),1)
                                                  call la_qaxpy(m,cs*sn*apoaq,a(1,p),1, &
                                                            a(1,q),1)
                                                  d(p) = d(p)*cs
                                                  d(q) = d(q)/cs
                                                  if (rsvec) then
                                                     call la_qaxpy(mvl,-t*aqoap,v(1,q),1, &
                                                               v(1,p),1)
                                                     call la_qaxpy(mvl,cs*sn*apoaq,v(1,p), &
                                                                1,v(1,q),1)
                                                  end if
                                               else
                                                  call la_qaxpy(m,t*apoaq,a(1,p),1,a(1, &
                                                             q),1)
                                                  call la_qaxpy(m,-cs*sn*aqoap,a(1,q),1, &
                                                            a(1,p),1)
                                                  d(p) = d(p)/cs
                                                  d(q) = d(q)*cs
                                                  if (rsvec) then
                                                     call la_qaxpy(mvl,t*apoaq,v(1,p),1, &
                                                               v(1,q),1)
                                                     call la_qaxpy(mvl,-cs*sn*aqoap,v(1,q) &
                                                               ,1,v(1,p),1)
                                                  end if
                                               end if
                                            end if
                                         end if
                                      end if
                                   else
                                      if (aapp > aaqq) then
                                         call la_qcopy(m,a(1,p),1,work,1)
                                         call la_qlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                    ierr)
                                         call la_qlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(p)/d(q)
                                         call la_qaxpy(m,temp1,work,1,a(1,q),1)

                                         call la_qlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      else
                                         call la_qcopy(m,a(1,q),1,work,1)
                                         call la_qlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                    ierr)
                                         call la_qlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         temp1 = -aapq*d(q)/d(p)
                                         call la_qaxpy(m,temp1,work,1,a(1,p),1)

                                         call la_qlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq*aapq))
                                         mxsinj = max(mxsinj,sfmin)
                                      end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q)
                 ! .. recompute sva(q)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_qnrm2(m,a(1,q),1)*d(q)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_qlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)*d(q)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_qnrm2(m,a(1,p),1)*d(p)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_qlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)*d(p)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
                 ! skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
            ! if ( notrot >= emptsw )  go to 2011
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
      ! **      if ( notrot >= emptsw )  go to 2011
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **   if ( notrot >= emptsw ) go to 1994
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_qnrm2(m,a(1,n),1)*d(n)
              else
                 t = zero
                 aapp = one
                 call la_qlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)*d(n)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < real(n,KIND=qp)*tol) .and. (real(n,KIND=qp) &
                        *mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:) reaching this point means that the procedure has completed the given
           ! number of sweeps.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means that during the i-th sweep all pivots were
           ! below the given threshold, causing early exit.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector d
           do p = 1,n - 1
              q = la_iqamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 temp1 = d(p)
                 d(p) = d(q)
                 d(q) = temp1
                 call la_qswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_qswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_qgsvj1
#endif

     !> STGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two real upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine SGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ),
     !> where U, V and Q are orthogonal matrices.
     !> R is a nonsingular upper triangular matrix, and D1 and D2 are
     !> ``diagonal'' matrices, which are of the following structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the orthogonal transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_stgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(sp),intent(in) :: tola,tolb
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           real(sp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(sp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(sp) :: a1,a2,a3,b1,b2,b3,csq,csu,csv,error,gamma,rwk,snq,snu,snv, &
                     ssmin
           ! Intrinsic Functions
           intrinsic :: abs,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('STGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_slaset('FULL',m,m,zero,one,u,ldu)
           if (initv) call la_slaset('FULL',p,p,zero,one,v,ldv)
           if (initq) call la_slaset('FULL',n,n,zero,one,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = zero
                    a3 = zero
                    if (k + i <= m) a1 = a(k + i,n - l + i)
                    if (k + j <= m) a3 = a(k + j,n - l + j)
                    b1 = b(i,n - l + i)
                    b3 = b(j,n - l + j)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_slags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**t *a
                    if (k + j <= m) call la_srot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,snu)
                    ! update i-th and j-th rows of matrix b: v**t *b
                    call la_srot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,snv)
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_srot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_srot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = zero
                       b(i,n - l + j) = zero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = zero
                       b(j,n - l + i) = zero
                    end if
                    ! update orthogonal matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_srot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_srot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_srot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_scopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_scopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_slapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = a(k + i,n - l + i)
              b1 = b(i,n - l + i)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 ! change sign if necessary
                 if (gamma < zero) then
                    call la_sscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_sscal(p,-one,v(1,i),1)
                 end if
                 call la_slartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_sscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_sscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_scopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_scopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_stgsja
     !> DTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two real upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine DGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ),
     !> where U, V and Q are orthogonal matrices.
     !> R is a nonsingular upper triangular matrix, and D1 and D2 are
     !> ``diagonal'' matrices, which are of the following structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the orthogonal transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_dtgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(dp),intent(in) :: tola,tolb
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           real(dp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(dp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(dp) :: a1,a2,a3,b1,b2,b3,csq,csu,csv,error,gamma,rwk,snq,snu,snv, &
                     ssmin
           ! Intrinsic Functions
           intrinsic :: abs,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('DTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_dlaset('FULL',m,m,zero,one,u,ldu)
           if (initv) call la_dlaset('FULL',p,p,zero,one,v,ldv)
           if (initq) call la_dlaset('FULL',n,n,zero,one,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = zero
                    a3 = zero
                    if (k + i <= m) a1 = a(k + i,n - l + i)
                    if (k + j <= m) a3 = a(k + j,n - l + j)
                    b1 = b(i,n - l + i)
                    b3 = b(j,n - l + j)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_dlags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**t *a
                    if (k + j <= m) call la_drot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,snu)
                    ! update i-th and j-th rows of matrix b: v**t *b
                    call la_drot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,snv)
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_drot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_drot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = zero
                       b(i,n - l + j) = zero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = zero
                       b(j,n - l + i) = zero
                    end if
                    ! update orthogonal matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_drot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_drot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_drot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_dcopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_dcopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_dlapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = a(k + i,n - l + i)
              b1 = b(i,n - l + i)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 ! change sign if necessary
                 if (gamma < zero) then
                    call la_dscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_dscal(p,-one,v(1,i),1)
                 end if
                 call la_dlartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_dscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_dscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_dcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_dcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_dtgsja
#ifdef LA_WITH_XDP
     !> XTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two real upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine DGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ),
     !> where U, V and Q are orthogonal matrices.
     !> R is a nonsingular upper triangular matrix, and D1 and D2 are
     !> ``diagonal'' matrices, which are of the following structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the orthogonal transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_xtgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(xdp),intent(in) :: tola,tolb
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           real(xdp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(xdp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(xdp) :: a1,a2,a3,b1,b2,b3,csq,csu,csv,error,gamma,rwk,snq,snu,snv, &
                     ssmin
           ! Intrinsic Functions
           intrinsic :: abs,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('XTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_xlaset('FULL',m,m,zero,one,u,ldu)
           if (initv) call la_xlaset('FULL',p,p,zero,one,v,ldv)
           if (initq) call la_xlaset('FULL',n,n,zero,one,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = zero
                    a3 = zero
                    if (k + i <= m) a1 = a(k + i,n - l + i)
                    if (k + j <= m) a3 = a(k + j,n - l + j)
                    b1 = b(i,n - l + i)
                    b3 = b(j,n - l + j)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_xlags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**t *a
                    if (k + j <= m) call la_xrot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,snu)
                    ! update i-th and j-th rows of matrix b: v**t *b
                    call la_xrot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,snv)
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_xrot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_xrot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = zero
                       b(i,n - l + j) = zero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = zero
                       b(j,n - l + i) = zero
                    end if
                    ! update orthogonal matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_xrot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_xrot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_xrot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_xcopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_xcopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_xlapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = a(k + i,n - l + i)
              b1 = b(i,n - l + i)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 ! change sign if necessary
                 if (gamma < zero) then
                    call la_xscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_xscal(p,-one,v(1,i),1)
                 end if
                 call la_xlartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_xscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_xscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_xcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_xcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_xtgsja
#endif
#ifdef LA_WITH_QP
     !> QTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two real upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine DGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**T *A*Q = D1*( 0 R ),    V**T *B*Q = D2*( 0 R ),
     !> where U, V and Q are orthogonal matrices.
     !> R is a nonsingular upper triangular matrix, and D1 and D2 are
     !> ``diagonal'' matrices, which are of the following structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the orthogonal transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_qtgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(qp),intent(in) :: tola,tolb
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           real(qp),intent(out) :: alpha(*),beta(*),work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(qp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(qp) :: a1,a2,a3,b1,b2,b3,csq,csu,csv,error,gamma,rwk,snq,snu,snv, &
                     ssmin
           ! Intrinsic Functions
           intrinsic :: abs,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('QTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_qlaset('FULL',m,m,zero,one,u,ldu)
           if (initv) call la_qlaset('FULL',p,p,zero,one,v,ldv)
           if (initq) call la_qlaset('FULL',n,n,zero,one,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = zero
                    a3 = zero
                    if (k + i <= m) a1 = a(k + i,n - l + i)
                    if (k + j <= m) a3 = a(k + j,n - l + j)
                    b1 = b(i,n - l + i)
                    b3 = b(j,n - l + j)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_qlags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**t *a
                    if (k + j <= m) call la_qrot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,snu)
                    ! update i-th and j-th rows of matrix b: v**t *b
                    call la_qrot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,snv)
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_qrot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_qrot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = zero
                       b(i,n - l + j) = zero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = zero
                       b(j,n - l + i) = zero
                    end if
                    ! update orthogonal matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_qrot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_qrot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_qrot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_qcopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_qcopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_qlapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = a(k + i,n - l + i)
              b1 = b(i,n - l + i)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 ! change sign if necessary
                 if (gamma < zero) then
                    call la_qscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_qscal(p,-one,v(1,i),1)
                 end if
                 call la_qlartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_qscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_qscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_qcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_qcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_qtgsja
#endif

     !> SGEBRD: reduces a general real M-by-N matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_sgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: max,min,real
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'SGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=sp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('SGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'SGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'SGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+nb-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_slabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+nb:m,i+nb:n), using an update
              ! of the form  a := a - v*y**t - x*u**t
              call la_sgemm('NO TRANSPOSE','TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one,a(i + &
                        nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,one,a(i + nb,i + nb),lda)
              call la_sgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,one,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_sgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_sgebrd
     !> DGEBRD: reduces a general real M-by-N matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_dgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'DGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=dp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('DGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'DGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'DGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+nb-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_dlabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+nb:m,i+nb:n), using an update
              ! of the form  a := a - v*y**t - x*u**t
              call la_dgemm('NO TRANSPOSE','TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one,a(i + &
                        nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,one,a(i + nb,i + nb),lda)
              call la_dgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,one,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_dgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_dgebrd
#ifdef LA_WITH_XDP
     !> XGEBRD: reduces a general real M-by-N matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_xgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'XGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=xdp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('XGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'XGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'XGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+nb-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_xlabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+nb:m,i+nb:n), using an update
              ! of the form  a := a - v*y**t - x*u**t
              call la_xgemm('NO TRANSPOSE','TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one,a(i + &
                        nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,one,a(i + nb,i + nb),lda)
              call la_xgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,one,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_xgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_xgebrd
#endif
#ifdef LA_WITH_QP
     !> QGEBRD: reduces a general real M-by-N matrix A to upper or lower
     !> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_qgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(out) :: d(*),e(*),taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'QGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=qp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('QGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'QGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'QGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+nb-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_qlabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+nb:m,i+nb:n), using an update
              ! of the form  a := a - v*y**t - x*u**t
              call la_qgemm('NO TRANSPOSE','TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one,a(i + &
                        nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,one,a(i + nb,i + nb),lda)
              call la_qgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-one, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,one,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_qgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_qgebrd
#endif

     !> SORGBR: generates one of the real orthogonal matrices Q or P**T
     !> determined by SGEBRD when reducing a real matrix A to bidiagonal
     !> form: A = Q * B * P**T.  Q and P**T are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and SORGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and SORGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
     !> is of order N:
     !> if k < n, P**T = G(k) . . . G(2) G(1) and SORGBR returns the first m
     !> rows of P**T, where n >= m >= k;
     !> if k >= n, P**T = G(n-1) . . . G(2) G(1) and SORGBR returns P**T as
     !> an N-by-N matrix.

     pure subroutine la_sorgbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_sorgqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_sorgqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_sorglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_sorglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = work(1)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('SORGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_sgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_sorgqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = zero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = one
                 do i = 2,m
                    a(i,1) = zero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_sorgqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**t, determined by a call to la_sgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_sorglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! row downward, and set the first row and column of p**t to
                 ! those of the unit matrix
                 a(1,1) = one
                 do i = 2,n
                    a(i,1) = zero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = zero
                 end do
                 if (n > 1) then
                    ! form p**t(2:n,2:n)
                    call la_sorglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_sorgbr
     !> DORGBR: generates one of the real orthogonal matrices Q or P**T
     !> determined by DGEBRD when reducing a real matrix A to bidiagonal
     !> form: A = Q * B * P**T.  Q and P**T are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
     !> is of order N:
     !> if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
     !> rows of P**T, where n >= m >= k;
     !> if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
     !> an N-by-N matrix.

     pure subroutine la_dorgbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_dorgqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_dorgqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_dorglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_dorglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = work(1)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('DORGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_dgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_dorgqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = zero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = one
                 do i = 2,m
                    a(i,1) = zero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_dorgqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**t, determined by a call to la_dgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_dorglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! row downward, and set the first row and column of p**t to
                 ! those of the unit matrix
                 a(1,1) = one
                 do i = 2,n
                    a(i,1) = zero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = zero
                 end do
                 if (n > 1) then
                    ! form p**t(2:n,2:n)
                    call la_dorglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_dorgbr
#ifdef LA_WITH_XDP
     !> XORGBR: generates one of the real orthogonal matrices Q or P**T
     !> determined by XGEBRD when reducing a real matrix A to bidiagonal
     !> form: A = Q * B * P**T.  Q and P**T are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and XORGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and XORGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
     !> is of order N:
     !> if k < n, P**T = G(k) . . . G(2) G(1) and XORGBR returns the first m
     !> rows of P**T, where n >= m >= k;
     !> if k >= n, P**T = G(n-1) . . . G(2) G(1) and XORGBR returns P**T as
     !> an N-by-N matrix.

     pure subroutine la_xorgbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_xorgqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_xorgqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_xorglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_xorglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = work(1)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('XORGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_xgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_xorgqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = zero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = one
                 do i = 2,m
                    a(i,1) = zero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_xorgqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**t, determined by a call to la_xgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_xorglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! row downward, and set the first row and column of p**t to
                 ! those of the unit matrix
                 a(1,1) = one
                 do i = 2,n
                    a(i,1) = zero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = zero
                 end do
                 if (n > 1) then
                    ! form p**t(2:n,2:n)
                    call la_xorglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_xorgbr
#endif
#ifdef LA_WITH_QP
     !> QORGBR: generates one of the real orthogonal matrices Q or P**T
     !> determined by QGEBRD when reducing a real matrix A to bidiagonal
     !> form: A = Q * B * P**T.  Q and P**T are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and QORGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and QORGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
     !> is of order N:
     !> if k < n, P**T = G(k) . . . G(2) G(1) and QORGBR returns the first m
     !> rows of P**T, where n >= m >= k;
     !> if k >= n, P**T = G(n-1) . . . G(2) G(1) and QORGBR returns P**T as
     !> an N-by-N matrix.

     pure subroutine la_qorgbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_qorgqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_qorgqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_qorglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_qorglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = work(1)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('QORGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_qgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_qorgqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = zero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = one
                 do i = 2,m
                    a(i,1) = zero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_qorgqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**t, determined by a call to la_qgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_qorglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors one
                 ! row downward, and set the first row and column of p**t to
                 ! those of the unit matrix
                 a(1,1) = one
                 do i = 2,n
                    a(i,1) = zero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = zero
                 end do
                 if (n > 1) then
                    ! form p**t(2:n,2:n)
                    call la_qorglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_qorgbr
#endif

     !> If VECT = 'Q', SORMBR: overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'T':      P**T * C       C * P**T
     !> Here Q and P**T are the orthogonal matrices determined by SGEBRD when
     !> reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
     !> P**T are defined as products of elementary reflectors H(i) and G(i)
     !> respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the orthogonal matrix Q or P**T that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_sormbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(sp),intent(inout) :: a(lda,*),c(ldc,*)
           real(sp),intent(in) :: tau(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (applyq) then
                 if (left) then
                    nb = la_ilaenv(1,'SORMQR',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'SORMQR',side//trans,m,n - 1,n - 1,-1)
                 end if
              else
                 if (left) then
                    nb = la_ilaenv(1,'SORMLQ',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'SORMLQ',side//trans,m,n - 1,n - 1,-1)
                 end if
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('SORMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           work(1) = 1
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_sgebrd with nq >= k
                 call la_sormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_sgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_sormqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_sgebrd with nq > k
                 call la_sormlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_sgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_sormlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_sormbr
     !> If VECT = 'Q', DORMBR: overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'T':      P**T * C       C * P**T
     !> Here Q and P**T are the orthogonal matrices determined by DGEBRD when
     !> reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
     !> P**T are defined as products of elementary reflectors H(i) and G(i)
     !> respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the orthogonal matrix Q or P**T that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_dormbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(dp),intent(inout) :: a(lda,*),c(ldc,*)
           real(dp),intent(in) :: tau(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (applyq) then
                 if (left) then
                    nb = la_ilaenv(1,'DORMQR',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'DORMQR',side//trans,m,n - 1,n - 1,-1)
                 end if
              else
                 if (left) then
                    nb = la_ilaenv(1,'DORMLQ',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'DORMLQ',side//trans,m,n - 1,n - 1,-1)
                 end if
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('DORMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           work(1) = 1
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_dgebrd with nq >= k
                 call la_dormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_dgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_dormqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_dgebrd with nq > k
                 call la_dormlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_dgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_dormlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_dormbr
#ifdef LA_WITH_XDP
     !> If VECT = 'Q', XORMBR: overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> If VECT = 'P', XORMBR overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'T':      P**T * C       C * P**T
     !> Here Q and P**T are the orthogonal matrices determined by XGEBRD when
     !> reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
     !> P**T are defined as products of elementary reflectors H(i) and G(i)
     !> respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the orthogonal matrix Q or P**T that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_xormbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           real(xdp),intent(in) :: tau(*)
           real(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (applyq) then
                 if (left) then
                    nb = la_ilaenv(1,'XORMQR',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'XORMQR',side//trans,m,n - 1,n - 1,-1)
                 end if
              else
                 if (left) then
                    nb = la_ilaenv(1,'XORMLQ',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'XORMLQ',side//trans,m,n - 1,n - 1,-1)
                 end if
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('XORMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           work(1) = 1
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_xgebrd with nq >= k
                 call la_xormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_xgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_xormqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_xgebrd with nq > k
                 call la_xormlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_xgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_xormlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_xormbr
#endif
#ifdef LA_WITH_QP
     !> If VECT = 'Q', QORMBR: overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'T':      Q**T * C       C * Q**T
     !> If VECT = 'P', QORMBR overwrites the general real M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'T':      P**T * C       C * P**T
     !> Here Q and P**T are the orthogonal matrices determined by QGEBRD when
     !> reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
     !> P**T are defined as products of elementary reflectors H(i) and G(i)
     !> respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the orthogonal matrix Q or P**T that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_qormbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           real(qp),intent(inout) :: a(lda,*),c(ldc,*)
           real(qp),intent(in) :: tau(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'T')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (applyq) then
                 if (left) then
                    nb = la_ilaenv(1,'QORMQR',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'QORMQR',side//trans,m,n - 1,n - 1,-1)
                 end if
              else
                 if (left) then
                    nb = la_ilaenv(1,'QORMLQ',side//trans,m - 1,n,m - 1,-1)
                 else
                    nb = la_ilaenv(1,'QORMLQ',side//trans,m,n - 1,n - 1,-1)
                 end if
              end if
              lwkopt = nw*nb
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('QORMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           work(1) = 1
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_qgebrd with nq >= k
                 call la_qormqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_qgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_qormqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'T'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_qgebrd with nq > k
                 call la_qormlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_qgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_qormlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_qormbr
#endif

     !> CGEBD2: reduces a complex general m by n matrix A to upper or lower
     !> real bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_cgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(sp),intent(out) :: d(*),e(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(sp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('CGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 alpha = a(i,i)
                 call la_clarfg(m - i + 1,alpha,a(min(i + 1,m),i),1,tauq(i))
                 d(i) = real(alpha,KIND=sp)
                 a(i,i) = cone
                 ! apply h(i)**h to a(i:m,i+1:n) from the left
                 if (i < n) call la_clarf('LEFT',m - i + 1,n - i,a(i,i),1,conjg(tauq(i)), &
                           a(i,i + 1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_clacgv(n - i,a(i,i + 1),lda)
                    alpha = a(i,i + 1)
                    call la_clarfg(n - i,alpha,a(i,min(i + 2,n)),lda,taup(i))
                    e(i) = real(alpha,KIND=sp)
                    a(i,i + 1) = cone
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_clarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    call la_clacgv(n - i,a(i,i + 1),lda)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = czero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_clacgv(n - i + 1,a(i,i),lda)
                 alpha = a(i,i)
                 call la_clarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,taup(i))
                 d(i) = real(alpha,KIND=sp)
                 a(i,i) = cone
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_clarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 call la_clacgv(n - i + 1,a(i,i),lda)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    alpha = a(i + 1,i)
                    call la_clarfg(m - i,alpha,a(min(i + 2,m),i),1,tauq(i))
                    e(i) = real(alpha,KIND=sp)
                    a(i + 1,i) = cone
                    ! apply h(i)**h to a(i+1:m,i+1:n) from the left
                    call la_clarf('LEFT',m - i,n - i,a(i + 1,i),1,conjg(tauq(i)),a(i + &
                              1,i + 1),lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = czero
                 end if
              end do
           end if
           return
     end subroutine la_cgebd2
     !> ZGEBD2: reduces a complex general m by n matrix A to upper or lower
     !> real bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_zgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(dp),intent(out) :: d(*),e(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(dp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('ZGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 alpha = a(i,i)
                 call la_zlarfg(m - i + 1,alpha,a(min(i + 1,m),i),1,tauq(i))
                 d(i) = real(alpha,KIND=dp)
                 a(i,i) = cone
                 ! apply h(i)**h to a(i:m,i+1:n) from the left
                 if (i < n) call la_zlarf('LEFT',m - i + 1,n - i,a(i,i),1,conjg(tauq(i)), &
                           a(i,i + 1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_zlacgv(n - i,a(i,i + 1),lda)
                    alpha = a(i,i + 1)
                    call la_zlarfg(n - i,alpha,a(i,min(i + 2,n)),lda,taup(i))
                    e(i) = real(alpha,KIND=dp)
                    a(i,i + 1) = cone
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_zlarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    call la_zlacgv(n - i,a(i,i + 1),lda)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = czero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_zlacgv(n - i + 1,a(i,i),lda)
                 alpha = a(i,i)
                 call la_zlarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,taup(i))
                 d(i) = real(alpha,KIND=dp)
                 a(i,i) = cone
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_zlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 call la_zlacgv(n - i + 1,a(i,i),lda)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    alpha = a(i + 1,i)
                    call la_zlarfg(m - i,alpha,a(min(i + 2,m),i),1,tauq(i))
                    e(i) = real(alpha,KIND=dp)
                    a(i + 1,i) = cone
                    ! apply h(i)**h to a(i+1:m,i+1:n) from the left
                    call la_zlarf('LEFT',m - i,n - i,a(i + 1,i),1,conjg(tauq(i)),a(i + &
                              1,i + 1),lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = czero
                 end if
              end do
           end if
           return
     end subroutine la_zgebd2
#ifdef LA_WITH_XDP
     !> YGEBD2: reduces a complex general m by n matrix A to upper or lower
     !> real bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_ygebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(xdp),intent(out) :: d(*),e(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(xdp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('YGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 alpha = a(i,i)
                 call la_ylarfg(m - i + 1,alpha,a(min(i + 1,m),i),1,tauq(i))
                 d(i) = real(alpha,KIND=xdp)
                 a(i,i) = cone
                 ! apply h(i)**h to a(i:m,i+1:n) from the left
                 if (i < n) call la_ylarf('LEFT',m - i + 1,n - i,a(i,i),1,conjg(tauq(i)), &
                           a(i,i + 1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_ylacgv(n - i,a(i,i + 1),lda)
                    alpha = a(i,i + 1)
                    call la_ylarfg(n - i,alpha,a(i,min(i + 2,n)),lda,taup(i))
                    e(i) = real(alpha,KIND=xdp)
                    a(i,i + 1) = cone
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_ylarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    call la_ylacgv(n - i,a(i,i + 1),lda)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = czero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_ylacgv(n - i + 1,a(i,i),lda)
                 alpha = a(i,i)
                 call la_ylarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,taup(i))
                 d(i) = real(alpha,KIND=xdp)
                 a(i,i) = cone
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_ylarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 call la_ylacgv(n - i + 1,a(i,i),lda)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    alpha = a(i + 1,i)
                    call la_ylarfg(m - i,alpha,a(min(i + 2,m),i),1,tauq(i))
                    e(i) = real(alpha,KIND=xdp)
                    a(i + 1,i) = cone
                    ! apply h(i)**h to a(i+1:m,i+1:n) from the left
                    call la_ylarf('LEFT',m - i,n - i,a(i + 1,i),1,conjg(tauq(i)),a(i + &
                              1,i + 1),lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = czero
                 end if
              end do
           end if
           return
     end subroutine la_ygebd2
#endif
#ifdef LA_WITH_QP
     !> WGEBD2: reduces a complex general m by n matrix A to upper or lower
     !> real bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_wgebd2(m,n,a,lda,d,e,tauq,taup,work,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,m,n
           ! Array Arguments
           real(qp),intent(out) :: d(*),e(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           complex(qp) :: alpha
           ! Intrinsic Functions
           intrinsic :: conjg,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           end if
           if (info < 0) then
              call la_xerbla('WGEBD2',-info)
              return
           end if
           if (m >= n) then
              ! reduce to upper bidiagonal form
              do i = 1,n
                 ! generate elementary reflector h(i) to annihilate a(i+1:m,i)
                 alpha = a(i,i)
                 call la_wlarfg(m - i + 1,alpha,a(min(i + 1,m),i),1,tauq(i))
                 d(i) = real(alpha,KIND=qp)
                 a(i,i) = cone
                 ! apply h(i)**h to a(i:m,i+1:n) from the left
                 if (i < n) call la_wlarf('LEFT',m - i + 1,n - i,a(i,i),1,conjg(tauq(i)), &
                           a(i,i + 1),lda,work)
                 a(i,i) = d(i)
                 if (i < n) then
                    ! generate elementary reflector g(i) to annihilate
                    ! a(i,i+2:n)
                    call la_wlacgv(n - i,a(i,i + 1),lda)
                    alpha = a(i,i + 1)
                    call la_wlarfg(n - i,alpha,a(i,min(i + 2,n)),lda,taup(i))
                    e(i) = real(alpha,KIND=qp)
                    a(i,i + 1) = cone
                    ! apply g(i) to a(i+1:m,i+1:n) from the right
                    call la_wlarf('RIGHT',m - i,n - i,a(i,i + 1),lda,taup(i),a(i + 1,i + 1 &
                              ),lda,work)
                    call la_wlacgv(n - i,a(i,i + 1),lda)
                    a(i,i + 1) = e(i)
                 else
                    taup(i) = czero
                 end if
              end do
           else
              ! reduce to lower bidiagonal form
              do i = 1,m
                 ! generate elementary reflector g(i) to annihilate a(i,i+1:n)
                 call la_wlacgv(n - i + 1,a(i,i),lda)
                 alpha = a(i,i)
                 call la_wlarfg(n - i + 1,alpha,a(i,min(i + 1,n)),lda,taup(i))
                 d(i) = real(alpha,KIND=qp)
                 a(i,i) = cone
                 ! apply g(i) to a(i+1:m,i:n) from the right
                 if (i < m) call la_wlarf('RIGHT',m - i,n - i + 1,a(i,i),lda,taup(i),a(i + &
                           1,i),lda,work)
                 call la_wlacgv(n - i + 1,a(i,i),lda)
                 a(i,i) = d(i)
                 if (i < m) then
                    ! generate elementary reflector h(i) to annihilate
                    ! a(i+2:m,i)
                    alpha = a(i + 1,i)
                    call la_wlarfg(m - i,alpha,a(min(i + 2,m),i),1,tauq(i))
                    e(i) = real(alpha,KIND=qp)
                    a(i + 1,i) = cone
                    ! apply h(i)**h to a(i+1:m,i+1:n) from the left
                    call la_wlarf('LEFT',m - i,n - i,a(i + 1,i),1,conjg(tauq(i)),a(i + &
                              1,i + 1),lda,work)
                    a(i + 1,i) = e(i)
                 else
                    tauq(i) = czero
                 end if
              end do
           end if
           return
     end subroutine la_wgebd2
#endif

     !> CTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two complex upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine CGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ),
     !> where U, V and Q are unitary matrices.
     !> R is a nonsingular upper triangular matrix, and D1
     !> and D2 are ``diagonal'' matrices, which are of the following
     !> structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the unitary transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_ctgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(sp),intent(in) :: tola,tolb
           ! Array Arguments
           real(sp),intent(out) :: alpha(*),beta(*)
           complex(sp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(sp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(sp) :: a1,a3,b1,b3,csq,csu,csv,error,gamma,rwk,ssmin
           complex(sp) :: a2,b2,snq,snu,snv
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,min,real,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('CTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_claset('FULL',m,m,czero,cone,u,ldu)
           if (initv) call la_claset('FULL',p,p,czero,cone,v,ldv)
           if (initq) call la_claset('FULL',n,n,czero,cone,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = czero
                    a3 = zero
                    if (k + i <= m) a1 = real(a(k + i,n - l + i),KIND=sp)
                    if (k + j <= m) a3 = real(a(k + j,n - l + j),KIND=sp)
                    b1 = real(b(i,n - l + i),KIND=sp)
                    b3 = real(b(j,n - l + j),KIND=sp)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_clags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**h *a
                    if (k + j <= m) call la_crot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,conjg(snu))
                    ! update i-th and j-th rows of matrix b: v**h *b
                    call la_crot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,conjg(snv) &
                              )
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_crot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_crot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = czero
                       b(i,n - l + j) = czero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = czero
                       b(j,n - l + i) = czero
                    end if
                    ! ensure that the diagonal elements of a and b are real.
                    if (k + i <= m) a(k + i,n - l + i) = real(a(k + i,n - l + i),KIND=sp)
                    if (k + j <= m) a(k + j,n - l + j) = real(a(k + j,n - l + j),KIND=sp)
                    b(i,n - l + i) = real(b(i,n - l + i),KIND=sp)
                    b(j,n - l + j) = real(b(j,n - l + j),KIND=sp)
                    ! update unitary matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_crot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_crot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_crot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_ccopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_ccopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_clapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = real(a(k + i,n - l + i),KIND=sp)
              b1 = real(b(i,n - l + i),KIND=sp)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 if (gamma < zero) then
                    call la_csscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_csscal(p,-one,v(1,i),1)
                 end if
                 call la_slartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_csscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_csscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_ccopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_ccopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_ctgsja
     !> ZTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two complex upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine ZGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ),
     !> where U, V and Q are unitary matrices.
     !> R is a nonsingular upper triangular matrix, and D1
     !> and D2 are ``diagonal'' matrices, which are of the following
     !> structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the unitary transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_ztgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(dp),intent(in) :: tola,tolb
           ! Array Arguments
           real(dp),intent(out) :: alpha(*),beta(*)
           complex(dp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(dp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(dp) :: a1,a3,b1,b3,csq,csu,csv,error,gamma,rwk,ssmin
           complex(dp) :: a2,b2,snq,snu,snv
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('ZTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_zlaset('FULL',m,m,czero,cone,u,ldu)
           if (initv) call la_zlaset('FULL',p,p,czero,cone,v,ldv)
           if (initq) call la_zlaset('FULL',n,n,czero,cone,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = czero
                    a3 = zero
                    if (k + i <= m) a1 = real(a(k + i,n - l + i),KIND=dp)
                    if (k + j <= m) a3 = real(a(k + j,n - l + j),KIND=dp)
                    b1 = real(b(i,n - l + i),KIND=dp)
                    b3 = real(b(j,n - l + j),KIND=dp)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_zlags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**h *a
                    if (k + j <= m) call la_zrot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,conjg(snu))
                    ! update i-th and j-th rows of matrix b: v**h *b
                    call la_zrot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,conjg(snv) &
                              )
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_zrot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_zrot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = czero
                       b(i,n - l + j) = czero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = czero
                       b(j,n - l + i) = czero
                    end if
                    ! ensure that the diagonal elements of a and b are real.
                    if (k + i <= m) a(k + i,n - l + i) = real(a(k + i,n - l + i),KIND=dp)
                    if (k + j <= m) a(k + j,n - l + j) = real(a(k + j,n - l + j),KIND=dp)
                    b(i,n - l + i) = real(b(i,n - l + i),KIND=dp)
                    b(j,n - l + j) = real(b(j,n - l + j),KIND=dp)
                    ! update unitary matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_zrot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_zrot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_zrot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_zcopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_zcopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_zlapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = real(a(k + i,n - l + i),KIND=dp)
              b1 = real(b(i,n - l + i),KIND=dp)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 if (gamma < zero) then
                    call la_zdscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_zdscal(p,-one,v(1,i),1)
                 end if
                 call la_dlartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_zdscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_zdscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_zcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_zcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_ztgsja
#ifdef LA_WITH_XDP
     !> YTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two complex upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine ZGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ),
     !> where U, V and Q are unitary matrices.
     !> R is a nonsingular upper triangular matrix, and D1
     !> and D2 are ``diagonal'' matrices, which are of the following
     !> structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the unitary transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_ytgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(xdp),intent(in) :: tola,tolb
           ! Array Arguments
           real(xdp),intent(out) :: alpha(*),beta(*)
           complex(xdp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(xdp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(xdp) :: a1,a3,b1,b3,csq,csu,csv,error,gamma,rwk,ssmin
           complex(xdp) :: a2,b2,snq,snu,snv
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('YTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_ylaset('FULL',m,m,czero,cone,u,ldu)
           if (initv) call la_ylaset('FULL',p,p,czero,cone,v,ldv)
           if (initq) call la_ylaset('FULL',n,n,czero,cone,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = czero
                    a3 = zero
                    if (k + i <= m) a1 = real(a(k + i,n - l + i),KIND=xdp)
                    if (k + j <= m) a3 = real(a(k + j,n - l + j),KIND=xdp)
                    b1 = real(b(i,n - l + i),KIND=xdp)
                    b3 = real(b(j,n - l + j),KIND=xdp)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_ylags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**h *a
                    if (k + j <= m) call la_yrot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,conjg(snu))
                    ! update i-th and j-th rows of matrix b: v**h *b
                    call la_yrot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,conjg(snv) &
                              )
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_yrot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_yrot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = czero
                       b(i,n - l + j) = czero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = czero
                       b(j,n - l + i) = czero
                    end if
                    ! ensure that the diagonal elements of a and b are real.
                    if (k + i <= m) a(k + i,n - l + i) = real(a(k + i,n - l + i),KIND=xdp)
                    if (k + j <= m) a(k + j,n - l + j) = real(a(k + j,n - l + j),KIND=xdp)
                    b(i,n - l + i) = real(b(i,n - l + i),KIND=xdp)
                    b(j,n - l + j) = real(b(j,n - l + j),KIND=xdp)
                    ! update unitary matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_yrot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_yrot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_yrot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_ycopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_ycopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_ylapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = real(a(k + i,n - l + i),KIND=xdp)
              b1 = real(b(i,n - l + i),KIND=xdp)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 if (gamma < zero) then
                    call la_yxscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_yxscal(p,-one,v(1,i),1)
                 end if
                 call la_xlartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_yxscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_yxscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_ycopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_ycopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_ytgsja
#endif
#ifdef LA_WITH_QP
     !> WTGSJA: computes the generalized singular value decomposition (GSVD)
     !> of two complex upper triangular (or trapezoidal) matrices A and B.
     !> On entry, it is assumed that matrices A and B have the following
     !> forms, which may be obtained by the preprocessing subroutine ZGGSVP
     !> from a general M-by-N matrix A and P-by-N matrix B:
     !> N-K-L  K    L
     !> A =    K ( 0    A12  A13 ) if M-K-L >= 0;
     !> L ( 0     0   A23 )
     !> M-K-L ( 0     0    0  )
     !> N-K-L  K    L
     !> A =  K ( 0    A12  A13 ) if M-K-L < 0;
     !> M-K ( 0     0   A23 )
     !> N-K-L  K    L
     !> B =  L ( 0     0   B13 )
     !> P-L ( 0     0    0  )
     !> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
     !> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
     !> otherwise A23 is (M-K)-by-L upper trapezoidal.
     !> On exit,
     !> U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ),
     !> where U, V and Q are unitary matrices.
     !> R is a nonsingular upper triangular matrix, and D1
     !> and D2 are ``diagonal'' matrices, which are of the following
     !> structures:
     !> If M-K-L >= 0,
     !> K  L
     !> D1 =     K ( I  0 )
     !> L ( 0  C )
     !> M-K-L ( 0  0 )
     !> K  L
     !> D2 = L   ( 0  S )
     !> P-L ( 0  0 )
     !> N-K-L  K    L
     !> ( 0 R ) = K (  0   R11  R12 ) K
     !> L (  0    0   R22 ) L
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
     !> S = diag( BETA(K+1),  ... , BETA(K+L) ),
     !> C**2 + S**2 = I.
     !> R is stored in A(1:K+L,N-K-L+1:N) on exit.
     !> If M-K-L < 0,
     !> K M-K K+L-M
     !> D1 =   K ( I  0    0   )
     !> M-K ( 0  C    0   )
     !> K M-K K+L-M
     !> D2 =   M-K ( 0  S    0   )
     !> K+L-M ( 0  0    I   )
     !> P-L ( 0  0    0   )
     !> N-K-L  K   M-K  K+L-M
     !> ( 0 R ) =    K ( 0    R11  R12  R13  )
     !> M-K ( 0     0   R22  R23  )
     !> K+L-M ( 0     0    0   R33  )
     !> where
     !> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
     !> S = diag( BETA(K+1),  ... , BETA(M) ),
     !> C**2 + S**2 = I.
     !> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
     !> (  0  R22 R23 )
     !> in B(M-K+1:L,N+M-K-L+1:N) on exit.
     !> The computation of the unitary transformation matrices U, V or Q
     !> is optional.  These matrices may either be formed explicitly, or they
     !> may be postmultiplied into input matrices U1, V1, or Q1.

     pure subroutine la_wtgsja(jobu,jobv,jobq,m,p,n,k,l,a,lda,b,ldb,tola,tolb, &
               alpha,beta,u,ldu,v,ldv,q,ldq,work,ncycle,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobq,jobu,jobv
           integer(ilp),intent(out) :: info,ncycle
           integer(ilp),intent(in) :: k,l,lda,ldb,ldq,ldu,ldv,m,n,p
           real(qp),intent(in) :: tola,tolb
           ! Array Arguments
           real(qp),intent(out) :: alpha(*),beta(*)
           complex(qp),intent(inout) :: a(lda,*),b(ldb,*),q(ldq,*),u(ldu,*),v(ldv,*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 40
           real(qp),parameter :: hugenum = huge(zero)

           ! Local Scalars
           logical(lk) :: initq,initu,initv,upper,wantq,wantu,wantv
           integer(ilp) :: i,j,kcycle
           real(qp) :: a1,a3,b1,b3,csq,csu,csv,error,gamma,rwk,ssmin
           complex(qp) :: a2,b2,snq,snu,snv
           ! Intrinsic Functions
           intrinsic :: abs,real,conjg,max,min,huge
           ! Executable Statements
           ! decode and test the input parameters
           initu = la_lsame(jobu,'I')
           wantu = initu .or. la_lsame(jobu,'U')
           initv = la_lsame(jobv,'I')
           wantv = initv .or. la_lsame(jobv,'V')
           initq = la_lsame(jobq,'I')
           wantq = initq .or. la_lsame(jobq,'Q')
           info = 0
           if (.not. (initu .or. wantu .or. la_lsame(jobu,'N'))) then
              info = -1
           else if (.not. (initv .or. wantv .or. la_lsame(jobv,'N'))) then
              info = -2
           else if (.not. (initq .or. wantq .or. la_lsame(jobq,'N'))) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (p < 0) then
              info = -5
           else if (n < 0) then
              info = -6
           else if (lda < max(1,m)) then
              info = -10
           else if (ldb < max(1,p)) then
              info = -12
           else if (ldu < 1 .or. (wantu .and. ldu < m)) then
              info = -18
           else if (ldv < 1 .or. (wantv .and. ldv < p)) then
              info = -20
           else if (ldq < 1 .or. (wantq .and. ldq < n)) then
              info = -22
           end if
           if (info /= 0) then
              call la_xerbla('WTGSJA',-info)
              return
           end if
           ! initialize u, v and q, if necessary
           if (initu) call la_wlaset('FULL',m,m,czero,cone,u,ldu)
           if (initv) call la_wlaset('FULL',p,p,czero,cone,v,ldv)
           if (initq) call la_wlaset('FULL',n,n,czero,cone,q,ldq)
           ! loop until convergence
           upper = .false.
           loop_40: do kcycle = 1,maxit
              upper = .not. upper
              loop_20: do i = 1,l - 1
                 loop_10: do j = i + 1,l
                    a1 = zero
                    a2 = czero
                    a3 = zero
                    if (k + i <= m) a1 = real(a(k + i,n - l + i),KIND=qp)
                    if (k + j <= m) a3 = real(a(k + j,n - l + j),KIND=qp)
                    b1 = real(b(i,n - l + i),KIND=qp)
                    b3 = real(b(j,n - l + j),KIND=qp)
                    if (upper) then
                       if (k + i <= m) a2 = a(k + i,n - l + j)
                       b2 = b(i,n - l + j)
                    else
                       if (k + j <= m) a2 = a(k + j,n - l + i)
                       b2 = b(j,n - l + i)
                    end if
                    call la_wlags2(upper,a1,a2,a3,b1,b2,b3,csu,snu,csv,snv,csq, &
                              snq)
                    ! update (k+i)-th and (k+j)-th rows of matrix a: u**h *a
                    if (k + j <= m) call la_wrot(l,a(k + j,n - l + 1),lda,a(k + i,n - l + 1),lda, &
                              csu,conjg(snu))
                    ! update i-th and j-th rows of matrix b: v**h *b
                    call la_wrot(l,b(j,n - l + 1),ldb,b(i,n - l + 1),ldb,csv,conjg(snv) &
                              )
                    ! update (n-l+i)-th and (n-l+j)-th columns of matrices
                    ! a and b: a*q and b*q
                    call la_wrot(min(k + l,m),a(1,n - l + j),1,a(1,n - l + i),1,csq,snq)

                    call la_wrot(l,b(1,n - l + j),1,b(1,n - l + i),1,csq,snq)
                    if (upper) then
                       if (k + i <= m) a(k + i,n - l + j) = czero
                       b(i,n - l + j) = czero
                    else
                       if (k + j <= m) a(k + j,n - l + i) = czero
                       b(j,n - l + i) = czero
                    end if
                    ! ensure that the diagonal elements of a and b are real.
                    if (k + i <= m) a(k + i,n - l + i) = real(a(k + i,n - l + i),KIND=qp)
                    if (k + j <= m) a(k + j,n - l + j) = real(a(k + j,n - l + j),KIND=qp)
                    b(i,n - l + i) = real(b(i,n - l + i),KIND=qp)
                    b(j,n - l + j) = real(b(j,n - l + j),KIND=qp)
                    ! update unitary matrices u, v, q, if desired.
                    if (wantu .and. k + j <= m) call la_wrot(m,u(1,k + j),1,u(1,k + i),1, &
                              csu,snu)
                    if (wantv) call la_wrot(p,v(1,j),1,v(1,i),1,csv,snv)
                    if (wantq) call la_wrot(n,q(1,n - l + j),1,q(1,n - l + i),1,csq,snq)

                 end do loop_10
              end do loop_20
              if (.not. upper) then
                 ! the matrices a13 and b13 were lower triangular at the start
                 ! of the cycle, and are now upper triangular.
                 ! convergence test: test the parallelism of the corresponding
                 ! rows of a and b.
                 error = zero
                 do i = 1,min(l,m - k)
                    call la_wcopy(l - i + 1,a(k + i,n - l + i),lda,work,1)
                    call la_wcopy(l - i + 1,b(i,n - l + i),ldb,work(l + 1),1)
                    call la_wlapll(l - i + 1,work,1,work(l + 1),1,ssmin)
                    error = max(error,ssmin)
                 end do
                 if (abs(error) <= min(tola,tolb)) go to 50
              end if
              ! end of cycle loop
           end do loop_40
           ! the algorithm has not converged after maxit cycles.
           info = 1
           go to 100
           50 continue
           ! if error <= min(tola,tolb), then the algorithm has converged.
           ! compute the generalized singular value pairs (alpha, beta), and
           ! set the triangular matrix r to array a.
           do i = 1,k
              alpha(i) = one
              beta(i) = zero
           end do
           do i = 1,min(l,m - k)
              a1 = real(a(k + i,n - l + i),KIND=qp)
              b1 = real(b(i,n - l + i),KIND=qp)
              gamma = b1/a1
              if ((gamma <= hugenum) .and. (gamma >= -hugenum)) then
                 if (gamma < zero) then
                    call la_wqscal(l - i + 1,-one,b(i,n - l + i),ldb)
                    if (wantv) call la_wqscal(p,-one,v(1,i),1)
                 end if
                 call la_qlartg(abs(gamma),one,beta(k + i),alpha(k + i),rwk)
                 if (alpha(k + i) >= beta(k + i)) then
                    call la_wqscal(l - i + 1,one/alpha(k + i),a(k + i,n - l + i),lda)
                 else
                    call la_wqscal(l - i + 1,one/beta(k + i),b(i,n - l + i),ldb)
                    call la_wcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
                 end if
              else
                 alpha(k + i) = zero
                 beta(k + i) = one
                 call la_wcopy(l - i + 1,b(i,n - l + i),ldb,a(k + i,n - l + i),lda)
              end if
           end do
           ! post-assignment
           do i = m + 1,k + l
              alpha(i) = zero
              beta(i) = one
           end do
           if (k + l < n) then
              do i = k + l + 1,n
                 alpha(i) = zero
                 beta(i) = zero
              end do
           end if
           100 continue
           ncycle = kcycle
           return
     end subroutine la_wtgsja
#endif

     !> CGBBRD: reduces a complex general m-by-n band matrix A to real upper
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> The routine computes B, and optionally forms Q or P**H, or computes
     !> Q**H*C for a given matrix C.

     pure subroutine la_cgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,rwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(sp),intent(out) :: d(*),e(*),rwork(*)
           complex(sp),intent(inout) :: ab(ldab,*),c(ldc,*)
           complex(sp),intent(out) :: pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mu, &
                      mu0,nr,nrt
           real(sp) :: abst,rc
           complex(sp) :: ra,rb,rs,t
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('CGBBRD',-info)
              return
           end if
           ! initialize q and p**h to the unit matrix, if needed
           if (wantq) call la_claset('FULL',m,m,czero,cone,q,ldq)
           if (wantpt) call la_claset('FULL',n,n,czero,cone,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the complex sines of the plane rotations are stored in work,
              ! and the real cosines in rwork.
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_clargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              rwork(j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_clartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,rwork(j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_clartg(ab(ku + ml - 1,i),ab(ku + ml,i),rwork(i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_crot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,rwork(i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_crot(m,q(1,j - 1),1,q(1,j),1,rwork(j),conjg( &
                                    work(j)))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_crot(ncc,c(j - 1,1),ldc,c(j,1),ldc,rwork(j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = rwork(j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_clargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               rwork(j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_clartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,rwork(j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_clartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),rwork( &
                                    i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_crot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,rwork(i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**h
                       do j = j1,j2,kb1
                          call la_crot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,rwork( &
                                     j + kun),conjg(work(j + kun)))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = rwork(j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to complex lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, overwriting superdiagonal
              ! elements on subdiagonal elements
              do i = 1,min(m - 1,n)
                 call la_clartg(ab(1,i),ab(2,i),rc,rs,ra)
                 ab(1,i) = ra
                 if (i < n) then
                    ab(2,i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_crot(m,q(1,i),1,q(1,i + 1),1,rc,conjg(rs))

                 if (wantc) call la_crot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
           else
              ! a has been reduced to complex upper bidiagonal form or is
              ! diagonal
              if (ku > 0 .and. m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_clartg(ab(ku + 1,i),rb,rc,rs,ra)
                    ab(ku + 1,i) = ra
                    if (i > 1) then
                       rb = -conjg(rs)*ab(ku,i)
                       ab(ku,i) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_crot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc, &
                              conjg(rs))
                 end do
              end if
           end if
           ! make diagonal and superdiagonal elements real, storing them in d
           ! and e
           t = ab(ku + 1,1)
           loop_120: do i = 1,minmn
              abst = abs(t)
              d(i) = abst
              if (abst /= zero) then
                 t = t/abst
              else
                 t = cone
              end if
              if (wantq) call la_cscal(m,t,q(1,i),1)
              if (wantc) call la_cscal(ncc,conjg(t),c(i,1),ldc)
              if (i < minmn) then
                 if (ku == 0 .and. kl == 0) then
                    e(i) = zero
                    t = ab(1,i + 1)
                 else
                    if (ku == 0) then
                       t = ab(2,i)*conjg(t)
                    else
                       t = ab(ku,i + 1)*conjg(t)
                    end if
                    abst = abs(t)
                    e(i) = abst
                    if (abst /= zero) then
                       t = t/abst
                    else
                       t = cone
                    end if
                    if (wantpt) call la_cscal(n,t,pt(i + 1,1),ldpt)
                    t = ab(ku + 1,i + 1)*conjg(t)
                 end if
              end if
           end do loop_120
           return
     end subroutine la_cgbbrd
     !> ZGBBRD: reduces a complex general m-by-n band matrix A to real upper
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> The routine computes B, and optionally forms Q or P**H, or computes
     !> Q**H*C for a given matrix C.

     pure subroutine la_zgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,rwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(dp),intent(out) :: d(*),e(*),rwork(*)
           complex(dp),intent(inout) :: ab(ldab,*),c(ldc,*)
           complex(dp),intent(out) :: pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mu, &
                      mu0,nr,nrt
           real(dp) :: abst,rc
           complex(dp) :: ra,rb,rs,t
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('ZGBBRD',-info)
              return
           end if
           ! initialize q and p**h to the unit matrix, if needed
           if (wantq) call la_zlaset('FULL',m,m,czero,cone,q,ldq)
           if (wantpt) call la_zlaset('FULL',n,n,czero,cone,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the complex sines of the plane rotations are stored in work,
              ! and the real cosines in rwork.
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_zlargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              rwork(j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_zlartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,rwork(j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_zlartg(ab(ku + ml - 1,i),ab(ku + ml,i),rwork(i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_zrot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,rwork(i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_zrot(m,q(1,j - 1),1,q(1,j),1,rwork(j),conjg( &
                                    work(j)))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_zrot(ncc,c(j - 1,1),ldc,c(j,1),ldc,rwork(j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = rwork(j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_zlargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               rwork(j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_zlartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,rwork(j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_zlartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),rwork( &
                                    i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_zrot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,rwork(i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**h
                       do j = j1,j2,kb1
                          call la_zrot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,rwork( &
                                     j + kun),conjg(work(j + kun)))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = rwork(j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to complex lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, overwriting superdiagonal
              ! elements on subdiagonal elements
              do i = 1,min(m - 1,n)
                 call la_zlartg(ab(1,i),ab(2,i),rc,rs,ra)
                 ab(1,i) = ra
                 if (i < n) then
                    ab(2,i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_zrot(m,q(1,i),1,q(1,i + 1),1,rc,conjg(rs))

                 if (wantc) call la_zrot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
           else
              ! a has been reduced to complex upper bidiagonal form or is
              ! diagonal
              if (ku > 0 .and. m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_zlartg(ab(ku + 1,i),rb,rc,rs,ra)
                    ab(ku + 1,i) = ra
                    if (i > 1) then
                       rb = -conjg(rs)*ab(ku,i)
                       ab(ku,i) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_zrot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc, &
                              conjg(rs))
                 end do
              end if
           end if
           ! make diagonal and superdiagonal elements real, storing them in d
           ! and e
           t = ab(ku + 1,1)
           loop_120: do i = 1,minmn
              abst = abs(t)
              d(i) = abst
              if (abst /= zero) then
                 t = t/abst
              else
                 t = cone
              end if
              if (wantq) call la_zscal(m,t,q(1,i),1)
              if (wantc) call la_zscal(ncc,conjg(t),c(i,1),ldc)
              if (i < minmn) then
                 if (ku == 0 .and. kl == 0) then
                    e(i) = zero
                    t = ab(1,i + 1)
                 else
                    if (ku == 0) then
                       t = ab(2,i)*conjg(t)
                    else
                       t = ab(ku,i + 1)*conjg(t)
                    end if
                    abst = abs(t)
                    e(i) = abst
                    if (abst /= zero) then
                       t = t/abst
                    else
                       t = cone
                    end if
                    if (wantpt) call la_zscal(n,t,pt(i + 1,1),ldpt)
                    t = ab(ku + 1,i + 1)*conjg(t)
                 end if
              end if
           end do loop_120
           return
     end subroutine la_zgbbrd
#ifdef LA_WITH_XDP
     !> YGBBRD: reduces a complex general m-by-n band matrix A to real upper
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> The routine computes B, and optionally forms Q or P**H, or computes
     !> Q**H*C for a given matrix C.

     pure subroutine la_ygbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,rwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(xdp),intent(out) :: d(*),e(*),rwork(*)
           complex(xdp),intent(inout) :: ab(ldab,*),c(ldc,*)
           complex(xdp),intent(out) :: pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mu, &
                      mu0,nr,nrt
           real(xdp) :: abst,rc
           complex(xdp) :: ra,rb,rs,t
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('YGBBRD',-info)
              return
           end if
           ! initialize q and p**h to the unit matrix, if needed
           if (wantq) call la_ylaset('FULL',m,m,czero,cone,q,ldq)
           if (wantpt) call la_ylaset('FULL',n,n,czero,cone,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the complex sines of the plane rotations are stored in work,
              ! and the real cosines in rwork.
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_ylargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              rwork(j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_ylartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,rwork(j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_ylartg(ab(ku + ml - 1,i),ab(ku + ml,i),rwork(i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_yrot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,rwork(i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_yrot(m,q(1,j - 1),1,q(1,j),1,rwork(j),conjg( &
                                    work(j)))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_yrot(ncc,c(j - 1,1),ldc,c(j,1),ldc,rwork(j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = rwork(j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_ylargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               rwork(j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_ylartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,rwork(j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_ylartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),rwork( &
                                    i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_yrot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,rwork(i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**h
                       do j = j1,j2,kb1
                          call la_yrot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,rwork( &
                                     j + kun),conjg(work(j + kun)))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = rwork(j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to complex lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, overwriting superdiagonal
              ! elements on subdiagonal elements
              do i = 1,min(m - 1,n)
                 call la_ylartg(ab(1,i),ab(2,i),rc,rs,ra)
                 ab(1,i) = ra
                 if (i < n) then
                    ab(2,i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_yrot(m,q(1,i),1,q(1,i + 1),1,rc,conjg(rs))

                 if (wantc) call la_yrot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
           else
              ! a has been reduced to complex upper bidiagonal form or is
              ! diagonal
              if (ku > 0 .and. m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_ylartg(ab(ku + 1,i),rb,rc,rs,ra)
                    ab(ku + 1,i) = ra
                    if (i > 1) then
                       rb = -conjg(rs)*ab(ku,i)
                       ab(ku,i) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_yrot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc, &
                              conjg(rs))
                 end do
              end if
           end if
           ! make diagonal and superdiagonal elements real, storing them in d
           ! and e
           t = ab(ku + 1,1)
           loop_120: do i = 1,minmn
              abst = abs(t)
              d(i) = abst
              if (abst /= zero) then
                 t = t/abst
              else
                 t = cone
              end if
              if (wantq) call la_yscal(m,t,q(1,i),1)
              if (wantc) call la_yscal(ncc,conjg(t),c(i,1),ldc)
              if (i < minmn) then
                 if (ku == 0 .and. kl == 0) then
                    e(i) = zero
                    t = ab(1,i + 1)
                 else
                    if (ku == 0) then
                       t = ab(2,i)*conjg(t)
                    else
                       t = ab(ku,i + 1)*conjg(t)
                    end if
                    abst = abs(t)
                    e(i) = abst
                    if (abst /= zero) then
                       t = t/abst
                    else
                       t = cone
                    end if
                    if (wantpt) call la_yscal(n,t,pt(i + 1,1),ldpt)
                    t = ab(ku + 1,i + 1)*conjg(t)
                 end if
              end if
           end do loop_120
           return
     end subroutine la_ygbbrd
#endif
#ifdef LA_WITH_QP
     !> WGBBRD: reduces a complex general m-by-n band matrix A to real upper
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> The routine computes B, and optionally forms Q or P**H, or computes
     !> Q**H*C for a given matrix C.

     pure subroutine la_wgbbrd(vect,m,n,ncc,kl,ku,ab,ldab,d,e,q,ldq,pt,ldpt,c, &
               ldc,work,rwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: kl,ku,ldab,ldc,ldpt,ldq,m,n,ncc
           ! Array Arguments
           real(qp),intent(out) :: d(*),e(*),rwork(*)
           complex(qp),intent(inout) :: ab(ldab,*),c(ldc,*)
           complex(qp),intent(out) :: pt(ldpt,*),q(ldq,*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantb,wantc,wantpt,wantq
           integer(ilp) :: i,inca,j,j1,j2,kb,kb1,kk,klm,klu1,kun,l,minmn,ml,ml0,mu, &
                      mu0,nr,nrt
           real(qp) :: abst,rc
           complex(qp) :: ra,rb,rs,t
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,min
           ! Executable Statements
           ! test the input parameters
           wantb = la_lsame(vect,'B')
           wantq = la_lsame(vect,'Q') .or. wantb
           wantpt = la_lsame(vect,'P') .or. wantb
           wantc = ncc > 0
           klu1 = kl + ku + 1
           info = 0
           if (.not. wantq .and. .not. wantpt .and. .not. la_lsame(vect,'N')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (ncc < 0) then
              info = -4
           else if (kl < 0) then
              info = -5
           else if (ku < 0) then
              info = -6
           else if (ldab < klu1) then
              info = -8
           else if (ldq < 1 .or. wantq .and. ldq < max(1,m)) then
              info = -12
           else if (ldpt < 1 .or. wantpt .and. ldpt < max(1,n)) then
              info = -14
           else if (ldc < 1 .or. wantc .and. ldc < max(1,m)) then
              info = -16
           end if
           if (info /= 0) then
              call la_xerbla('WGBBRD',-info)
              return
           end if
           ! initialize q and p**h to the unit matrix, if needed
           if (wantq) call la_wlaset('FULL',m,m,czero,cone,q,ldq)
           if (wantpt) call la_wlaset('FULL',n,n,czero,cone,pt,ldpt)
           ! quick return if possible.
           if (m == 0 .or. n == 0) return
           minmn = min(m,n)
           if (kl + ku > 1) then
              ! reduce to upper bidiagonal form if ku > 0; if ku = 0, reduce
              ! first to lower bidiagonal form and then transform to upper
              ! bidiagonal
              if (ku > 0) then
                 ml0 = 1
                 mu0 = 2
              else
                 ml0 = 2
                 mu0 = 1
              end if
              ! wherever possible, plane rotations are generated and applied in
              ! vector operations of length nr over the index set j1:j2:klu1.
              ! the complex sines of the plane rotations are stored in work,
              ! and the real cosines in rwork.
              klm = min(m - 1,kl)
              kun = min(n - 1,ku)
              kb = klm + kun
              kb1 = kb + 1
              inca = kb1*ldab
              nr = 0
              j1 = klm + 2
              j2 = 1 - kun
              loop_90: do i = 1,minmn
                 ! reduce i-th column and i-th row of matrix to bidiagonal form
                 ml = klm + 1
                 mu = kun + 1
                 loop_80: do kk = 1,kb
                    j1 = j1 + kb
                    j2 = j2 + kb
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been created below the band
                    if (nr > 0) call la_wlargv(nr,ab(klu1,j1 - klm - 1),inca,work(j1),kb1, &
                              rwork(j1),kb1)
                    ! apply plane rotations from the left
                    do l = 1,kb
                       if (j2 - klm + l - 1 > n) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_wlartv(nrt,ab(klu1 - l,j1 - klm + l - 1),inca,ab( &
                                 klu1 - l + 1,j1 - klm + l - 1),inca,rwork(j1),work(j1),kb1)
                    end do
                    if (ml > ml0) then
                       if (ml <= m - i + 1) then
                          ! generate plane rotation to annihilate a(i+ml-1,i)
                          ! within the band, and apply rotation from the left
                          call la_wlartg(ab(ku + ml - 1,i),ab(ku + ml,i),rwork(i + ml - 1), &
                                    work(i + ml - 1),ra)
                          ab(ku + ml - 1,i) = ra
                          if (i < n) call la_wrot(min(ku + ml - 2,n - i),ab(ku + ml - 2,i + 1),ldab - &
                                    1,ab(ku + ml - 1,i + 1),ldab - 1,rwork(i + ml - 1),work(i + ml - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantq) then
                       ! accumulate product of plane rotations in q
                       do j = j1,j2,kb1
                          call la_wrot(m,q(1,j - 1),1,q(1,j),1,rwork(j),conjg( &
                                    work(j)))
                       end do
                    end if
                    if (wantc) then
                       ! apply plane rotations to c
                       do j = j1,j2,kb1
                          call la_wrot(ncc,c(j - 1,1),ldc,c(j,1),ldc,rwork(j), &
                                    work(j))
                       end do
                    end if
                    if (j2 + kun > n) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j-1,j+ku) above the band
                       ! and store it in work(n+1:2*n)
                       work(j + kun) = work(j)*ab(1,j + kun)
                       ab(1,j + kun) = rwork(j)*ab(1,j + kun)
                    end do
                    ! generate plane rotations to annihilate nonzero elements
                    ! which have been generated above the band
                    if (nr > 0) call la_wlargv(nr,ab(1,j1 + kun - 1),inca,work(j1 + kun),kb1, &
                               rwork(j1 + kun),kb1)
                    ! apply plane rotations from the right
                    do l = 1,kb
                       if (j2 + l - 1 > m) then
                          nrt = nr - 1
                       else
                          nrt = nr
                       end if
                       if (nrt > 0) call la_wlartv(nrt,ab(l + 1,j1 + kun - 1),inca,ab(l,j1 + &
                                 kun),inca,rwork(j1 + kun),work(j1 + kun),kb1)
                    end do
                    if (ml == ml0 .and. mu > mu0) then
                       if (mu <= n - i + 1) then
                          ! generate plane rotation to annihilate a(i,i+mu-1)
                          ! within the band, and apply rotation from the right
                          call la_wlartg(ab(ku - mu + 3,i + mu - 2),ab(ku - mu + 2,i + mu - 1),rwork( &
                                    i + mu - 1),work(i + mu - 1),ra)
                          ab(ku - mu + 3,i + mu - 2) = ra
                          call la_wrot(min(kl + mu - 2,m - i),ab(ku - mu + 4,i + mu - 2),1,ab(ku - &
                                    mu + 3,i + mu - 1),1,rwork(i + mu - 1),work(i + mu - 1))
                       end if
                       nr = nr + 1
                       j1 = j1 - kb1
                    end if
                    if (wantpt) then
                       ! accumulate product of plane rotations in p**h
                       do j = j1,j2,kb1
                          call la_wrot(n,pt(j + kun - 1,1),ldpt,pt(j + kun,1),ldpt,rwork( &
                                     j + kun),conjg(work(j + kun)))
                       end do
                    end if
                    if (j2 + kb > m) then
                       ! adjust j2 to keep within the bounds of the matrix
                       nr = nr - 1
                       j2 = j2 - kb1
                    end if
                    do j = j1,j2,kb1
                       ! create nonzero element a(j+kl+ku,j+ku-1) below the
                       ! band and store it in work(1:n)
                       work(j + kb) = work(j + kun)*ab(klu1,j + kun)
                       ab(klu1,j + kun) = rwork(j + kun)*ab(klu1,j + kun)
                    end do
                    if (ml > ml0) then
                       ml = ml - 1
                    else
                       mu = mu - 1
                    end if
                 end do loop_80
              end do loop_90
           end if
           if (ku == 0 .and. kl > 0) then
              ! a has been reduced to complex lower bidiagonal form
              ! transform lower bidiagonal form to upper bidiagonal by applying
              ! plane rotations from the left, overwriting superdiagonal
              ! elements on subdiagonal elements
              do i = 1,min(m - 1,n)
                 call la_wlartg(ab(1,i),ab(2,i),rc,rs,ra)
                 ab(1,i) = ra
                 if (i < n) then
                    ab(2,i) = rs*ab(1,i + 1)
                    ab(1,i + 1) = rc*ab(1,i + 1)
                 end if
                 if (wantq) call la_wrot(m,q(1,i),1,q(1,i + 1),1,rc,conjg(rs))

                 if (wantc) call la_wrot(ncc,c(i,1),ldc,c(i + 1,1),ldc,rc,rs)

              end do
           else
              ! a has been reduced to complex upper bidiagonal form or is
              ! diagonal
              if (ku > 0 .and. m < n) then
                 ! annihilate a(m,m+1) by applying plane rotations from the
                 ! right
                 rb = ab(ku,m + 1)
                 do i = m,1,-1
                    call la_wlartg(ab(ku + 1,i),rb,rc,rs,ra)
                    ab(ku + 1,i) = ra
                    if (i > 1) then
                       rb = -conjg(rs)*ab(ku,i)
                       ab(ku,i) = rc*ab(ku,i)
                    end if
                    if (wantpt) call la_wrot(n,pt(i,1),ldpt,pt(m + 1,1),ldpt,rc, &
                              conjg(rs))
                 end do
              end if
           end if
           ! make diagonal and superdiagonal elements real, storing them in d
           ! and e
           t = ab(ku + 1,1)
           loop_120: do i = 1,minmn
              abst = abs(t)
              d(i) = abst
              if (abst /= zero) then
                 t = t/abst
              else
                 t = cone
              end if
              if (wantq) call la_wscal(m,t,q(1,i),1)
              if (wantc) call la_wscal(ncc,conjg(t),c(i,1),ldc)
              if (i < minmn) then
                 if (ku == 0 .and. kl == 0) then
                    e(i) = zero
                    t = ab(1,i + 1)
                 else
                    if (ku == 0) then
                       t = ab(2,i)*conjg(t)
                    else
                       t = ab(ku,i + 1)*conjg(t)
                    end if
                    abst = abs(t)
                    e(i) = abst
                    if (abst /= zero) then
                       t = t/abst
                    else
                       t = cone
                    end if
                    if (wantpt) call la_wscal(n,t,pt(i + 1,1),ldpt)
                    t = ab(ku + 1,i + 1)*conjg(t)
                 end if
              end if
           end do loop_120
           return
     end subroutine la_wgbbrd
#endif

     !> CGEBRD: reduces a general complex M-by-N matrix A to upper or lower
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_cgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(sp),intent(out) :: d(*),e(*)
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: max,min,real
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'CGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=sp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('CGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'CGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'CGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+ib-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_clabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+ib:m,i+ib:n), using
              ! an update of the form  a := a - v*y**h - x*u**h
              call la_cgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,- &
              cone,a(i + nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,cone,a(i + nb,i + nb),lda)

              call la_cgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-cone, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,cone,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_cgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_cgebrd
     !> ZGEBRD: reduces a general complex M-by-N matrix A to upper or lower
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_zgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(dp),intent(out) :: d(*),e(*)
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'ZGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=dp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('ZGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'ZGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'ZGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+ib-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_zlabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+ib:m,i+ib:n), using
              ! an update of the form  a := a - v*y**h - x*u**h
              call la_zgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,- &
              cone,a(i + nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,cone,a(i + nb,i + nb),lda)

              call la_zgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-cone, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,cone,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_zgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_zgebrd
#ifdef LA_WITH_XDP
     !> YGEBRD: reduces a general complex M-by-N matrix A to upper or lower
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_ygebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(xdp),intent(out) :: d(*),e(*)
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'YGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=xdp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('YGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'YGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'YGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+ib-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_ylabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+ib:m,i+ib:n), using
              ! an update of the form  a := a - v*y**h - x*u**h
              call la_ygemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,- &
              cone,a(i + nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,cone,a(i + nb,i + nb),lda)

              call la_ygemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-cone, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,cone,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_ygebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_ygebrd
#endif
#ifdef LA_WITH_QP
     !> WGEBRD: reduces a general complex M-by-N matrix A to upper or lower
     !> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
     !> If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

     pure subroutine la_wgebrd(m,n,a,lda,d,e,tauq,taup,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,lwork,m,n
           ! Array Arguments
           real(qp),intent(out) :: d(*),e(*)
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(out) :: taup(*),tauq(*),work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: i,iinfo,j,ldwrkx,ldwrky,lwkopt,minmn,nb,nbmin,nx,ws
           ! Intrinsic Functions
           intrinsic :: real,max,min
           ! Executable Statements
           ! test the input parameters
           info = 0
           nb = max(1,la_ilaenv(1,'WGEBRD',' ',m,n,-1,-1))
           lwkopt = (m + n)*nb
           work(1) = real(lwkopt,KIND=qp)
           lquery = (lwork == -1)
           if (m < 0) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (lda < max(1,m)) then
              info = -4
           else if (lwork < max(1,m,n) .and. .not. lquery) then
              info = -10
           end if
           if (info < 0) then
              call la_xerbla('WGEBRD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           minmn = min(m,n)
           if (minmn == 0) then
              work(1) = 1
              return
           end if
           ws = max(m,n)
           ldwrkx = m
           ldwrky = n
           if (nb > 1 .and. nb < minmn) then
              ! set the crossover point nx.
              nx = max(nb,la_ilaenv(3,'WGEBRD',' ',m,n,-1,-1))
              ! determine when to switch from blocked to unblocked code.
              if (nx < minmn) then
                 ws = (m + n)*nb
                 if (lwork < ws) then
                    ! not enough work space for the optimal nb, consider using
                    ! a smaller block size.
                    nbmin = la_ilaenv(2,'WGEBRD',' ',m,n,-1,-1)
                    if (lwork >= (m + n)*nbmin) then
                       nb = lwork/(m + n)
                    else
                       nb = 1
                       nx = minmn
                    end if
                 end if
              end if
           else
              nx = minmn
           end if
           do i = 1,minmn - nx,nb
              ! reduce rows and columns i:i+ib-1 to bidiagonal form and return
              ! the matrices x and y which are needed to update the unreduced
              ! part of the matrix
              call la_wlabrd(m - i + 1,n - i + 1,nb,a(i,i),lda,d(i),e(i),tauq(i), &
                        taup(i),work,ldwrkx,work(ldwrkx*nb + 1),ldwrky)
              ! update the trailing submatrix a(i+ib:m,i+ib:n), using
              ! an update of the form  a := a - v*y**h - x*u**h
              call la_wgemm('NO TRANSPOSE','CONJUGATE TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,- &
              cone,a(i + nb,i),lda,work(ldwrkx*nb + nb + 1),ldwrky,cone,a(i + nb,i + nb),lda)

              call la_wgemm('NO TRANSPOSE','NO TRANSPOSE',m - i - nb + 1,n - i - nb + 1,nb,-cone, &
                        work(nb + 1),ldwrkx,a(i,i + nb),lda,cone,a(i + nb,i + nb),lda)
              ! copy diagonal and off-diagonal elements of b back into a
              if (m >= n) then
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j,j + 1) = e(j)
                 end do
              else
                 do j = i,i + nb - 1
                    a(j,j) = d(j)
                    a(j + 1,j) = e(j)
                 end do
              end if
           end do
           ! use unblocked code to reduce the remainder of the matrix
           call la_wgebd2(m - i + 1,n - i + 1,a(i,i),lda,d(i),e(i),tauq(i),taup(i), &
                     work,iinfo)
           work(1) = ws
           return
     end subroutine la_wgebrd
#endif

     !> CUNGBR: generates one of the complex unitary matrices Q or P**H
     !> determined by CGEBRD when reducing a complex matrix A to bidiagonal
     !> form: A = Q * B * P**H.  Q and P**H are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and CUNGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and CUNGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H
     !> is of order N:
     !> if k < n, P**H = G(k) . . . G(2) G(1) and CUNGBR returns the first m
     !> rows of P**H, where n >= m >= k;
     !> if k >= n, P**H = G(n-1) . . . G(2) G(1) and CUNGBR returns P**H as
     !> an N-by-N matrix.

     pure subroutine la_cungbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_cungqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_cungqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_cunglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_cunglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = real(work(1),KIND=sp)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('CUNGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_cgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_cungqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = czero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = cone
                 do i = 2,m
                    a(i,1) = czero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_cungqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**h, determined by a call to la_cgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_cunglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! row downward, and set the first row and column of p**h to
                 ! those of the unit matrix
                 a(1,1) = cone
                 do i = 2,n
                    a(i,1) = czero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = czero
                 end do
                 if (n > 1) then
                    ! form p**h(2:n,2:n)
                    call la_cunglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_cungbr
     !> ZUNGBR: generates one of the complex unitary matrices Q or P**H
     !> determined by ZGEBRD when reducing a complex matrix A to bidiagonal
     !> form: A = Q * B * P**H.  Q and P**H are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and ZUNGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and ZUNGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H
     !> is of order N:
     !> if k < n, P**H = G(k) . . . G(2) G(1) and ZUNGBR returns the first m
     !> rows of P**H, where n >= m >= k;
     !> if k >= n, P**H = G(n-1) . . . G(2) G(1) and ZUNGBR returns P**H as
     !> an N-by-N matrix.

     pure subroutine la_zungbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_zungqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_zungqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_zunglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_zunglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = real(work(1),KIND=dp)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('ZUNGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_zgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_zungqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = czero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = cone
                 do i = 2,m
                    a(i,1) = czero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_zungqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**h, determined by a call to la_zgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_zunglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! row downward, and set the first row and column of p**h to
                 ! those of the unit matrix
                 a(1,1) = cone
                 do i = 2,n
                    a(i,1) = czero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = czero
                 end do
                 if (n > 1) then
                    ! form p**h(2:n,2:n)
                    call la_zunglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_zungbr
#ifdef LA_WITH_XDP
     !> YUNGBR: generates one of the complex unitary matrices Q or P**H
     !> determined by YGEBRD when reducing a complex matrix A to bidiagonal
     !> form: A = Q * B * P**H.  Q and P**H are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and YUNGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and YUNGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H
     !> is of order N:
     !> if k < n, P**H = G(k) . . . G(2) G(1) and YUNGBR returns the first m
     !> rows of P**H, where n >= m >= k;
     !> if k >= n, P**H = G(n-1) . . . G(2) G(1) and YUNGBR returns P**H as
     !> an N-by-N matrix.

     pure subroutine la_yungbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_xdp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_yungqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_yungqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_yunglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_yunglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = real(work(1),KIND=xdp)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('YUNGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_ygebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_yungqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = czero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = cone
                 do i = 2,m
                    a(i,1) = czero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_yungqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**h, determined by a call to la_ygebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_yunglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! row downward, and set the first row and column of p**h to
                 ! those of the unit matrix
                 a(1,1) = cone
                 do i = 2,n
                    a(i,1) = czero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = czero
                 end do
                 if (n > 1) then
                    ! form p**h(2:n,2:n)
                    call la_yunglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_yungbr
#endif
#ifdef LA_WITH_QP
     !> WUNGBR: generates one of the complex unitary matrices Q or P**H
     !> determined by WGEBRD when reducing a complex matrix A to bidiagonal
     !> form: A = Q * B * P**H.  Q and P**H are defined as products of
     !> elementary reflectors H(i) or G(i) respectively.
     !> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
     !> is of order M:
     !> if m >= k, Q = H(1) H(2) . . . H(k) and WUNGBR returns the first n
     !> columns of Q, where m >= n >= k;
     !> if m < k, Q = H(1) H(2) . . . H(m-1) and WUNGBR returns Q as an
     !> M-by-M matrix.
     !> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H
     !> is of order N:
     !> if k < n, P**H = G(k) . . . G(2) G(1) and WUNGBR returns the first m
     !> rows of P**H, where n >= m >= k;
     !> if k >= n, P**H = G(n-1) . . . G(2) G(1) and WUNGBR returns P**H as
     !> an N-by-N matrix.

     pure subroutine la_wungbr(vect,m,n,k,a,lda,tau,work,lwork,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantq
           integer(ilp) :: i,iinfo,j,lwkopt,mn
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           wantq = la_lsame(vect,'Q')
           mn = min(m,n)
           lquery = (lwork == -1)
           if (.not. wantq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (m < 0) then
              info = -2
           else if (n < 0 .or. (wantq .and. (n > m .or. n < min(m,k))) .or. (.not. wantq .and. ( &
                     m > n .or. m < min(n,k)))) then
              info = -3
           else if (k < 0) then
              info = -4
           else if (lda < max(1,m)) then
              info = -6
           else if (lwork < max(1,mn) .and. .not. lquery) then
              info = -9
           end if
           if (info == 0) then
              work(1) = 1
              if (wantq) then
                 if (m >= k) then
                    call la_wungqr(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (m > 1) then
                       call la_wungqr(m - 1,m - 1,m - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              else
                 if (k < n) then
                    call la_wunglq(m,n,k,a,lda,tau,work,-1,iinfo)
                 else
                    if (n > 1) then
                       call la_wunglq(n - 1,n - 1,n - 1,a,lda,tau,work,-1,iinfo)
                    end if
                 end if
              end if
              lwkopt = real(work(1),KIND=qp)
              lwkopt = max(lwkopt,mn)
           end if
           if (info /= 0) then
              call la_xerbla('WUNGBR',-info)
              return
           else if (lquery) then
              work(1) = lwkopt
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) then
              work(1) = 1
              return
           end if
           if (wantq) then
              ! form q, determined by a call to la_wgebrd to reduce an m-by-k
              ! matrix
              if (m >= k) then
                 ! if m >= k, assume m >= n >= k
                 call la_wungqr(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if m < k, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! column to the right, and set the first row and column of q
                 ! to those of the unit matrix
                 do j = m,2,-1
                    a(1,j) = czero
                    do i = j + 1,m
                       a(i,j) = a(i,j - 1)
                    end do
                 end do
                 a(1,1) = cone
                 do i = 2,m
                    a(i,1) = czero
                 end do
                 if (m > 1) then
                    ! form q(2:m,2:m)
                    call la_wungqr(m - 1,m - 1,m - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           else
              ! form p**h, determined by a call to la_wgebrd to reduce a k-by-n
              ! matrix
              if (k < n) then
                 ! if k < n, assume k <= m <= n
                 call la_wunglq(m,n,k,a,lda,tau,work,lwork,iinfo)
              else
                 ! if k >= n, assume m = n
                 ! shift the vectors which define the elementary reflectors cone
                 ! row downward, and set the first row and column of p**h to
                 ! those of the unit matrix
                 a(1,1) = cone
                 do i = 2,n
                    a(i,1) = czero
                 end do
                 do j = 2,n
                    do i = j - 1,2,-1
                       a(i,j) = a(i - 1,j)
                    end do
                    a(1,j) = czero
                 end do
                 if (n > 1) then
                    ! form p**h(2:n,2:n)
                    call la_wunglq(n - 1,n - 1,n - 1,a(2,2),lda,tau,work,lwork,iinfo)

                 end if
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_wungbr
#endif

     !> If VECT = 'Q', CUNMBR: overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> If VECT = 'P', CUNMBR overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'C':      P**H * C       C * P**H
     !> Here Q and P**H are the unitary matrices determined by CGEBRD when
     !> reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q
     !> and P**H are defined as products of elementary reflectors H(i) and
     !> G(i) respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the unitary matrix Q or P**H that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_cunmbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(sp),intent(in) :: tau(*)
           complex(sp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (m > 0 .and. n > 0) then
                 if (applyq) then
                    if (left) then
                       nb = la_ilaenv(1,'CUNMQR',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'CUNMQR',side//trans,m,n - 1,n - 1,-1)
                    end if
                 else
                    if (left) then
                       nb = la_ilaenv(1,'CUNMLQ',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'CUNMLQ',side//trans,m,n - 1,n - 1,-1)
                    end if
                 end if
                 lwkopt = nw*nb
              else
                 lwkopt = 1
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('CUNMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_cgebrd with nq >= k
                 call la_cunmqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_cgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_cunmqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_cgebrd with nq > k
                 call la_cunmlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_cgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_cunmlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_cunmbr
     !> If VECT = 'Q', ZUNMBR: overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> If VECT = 'P', ZUNMBR overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'C':      P**H * C       C * P**H
     !> Here Q and P**H are the unitary matrices determined by ZGEBRD when
     !> reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q
     !> and P**H are defined as products of elementary reflectors H(i) and
     !> G(i) respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the unitary matrix Q or P**H that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_zunmbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(dp),intent(in) :: tau(*)
           complex(dp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (m > 0 .and. n > 0) then
                 if (applyq) then
                    if (left) then
                       nb = la_ilaenv(1,'ZUNMQR',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'ZUNMQR',side//trans,m,n - 1,n - 1,-1)
                    end if
                 else
                    if (left) then
                       nb = la_ilaenv(1,'ZUNMLQ',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'ZUNMLQ',side//trans,m,n - 1,n - 1,-1)
                    end if
                 end if
                 lwkopt = nw*nb
              else
                 lwkopt = 1
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('ZUNMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_zgebrd with nq >= k
                 call la_zunmqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_zgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_zunmqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_zgebrd with nq > k
                 call la_zunmlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_zgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_zunmlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_zunmbr
#ifdef LA_WITH_XDP
     !> If VECT = 'Q', YUNMBR: overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> If VECT = 'P', YUNMBR overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'C':      P**H * C       C * P**H
     !> Here Q and P**H are the unitary matrices determined by YGEBRD when
     !> reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q
     !> and P**H are defined as products of elementary reflectors H(i) and
     !> G(i) respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the unitary matrix Q or P**H that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_yunmbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(xdp),intent(in) :: tau(*)
           complex(xdp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (m > 0 .and. n > 0) then
                 if (applyq) then
                    if (left) then
                       nb = la_ilaenv(1,'YUNMQR',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'YUNMQR',side//trans,m,n - 1,n - 1,-1)
                    end if
                 else
                    if (left) then
                       nb = la_ilaenv(1,'YUNMLQ',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'YUNMLQ',side//trans,m,n - 1,n - 1,-1)
                    end if
                 end if
                 lwkopt = nw*nb
              else
                 lwkopt = 1
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('YUNMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_ygebrd with nq >= k
                 call la_yunmqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_ygebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_yunmqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_ygebrd with nq > k
                 call la_yunmlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_ygebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_yunmlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_yunmbr
#endif
#ifdef LA_WITH_QP
     !> If VECT = 'Q', WUNMBR: overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      Q * C          C * Q
     !> TRANS = 'C':      Q**H * C       C * Q**H
     !> If VECT = 'P', WUNMBR overwrites the general complex M-by-N matrix C
     !> with
     !> SIDE = 'L'     SIDE = 'R'
     !> TRANS = 'N':      P * C          C * P
     !> TRANS = 'C':      P**H * C       C * P**H
     !> Here Q and P**H are the unitary matrices determined by WGEBRD when
     !> reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q
     !> and P**H are defined as products of elementary reflectors H(i) and
     !> G(i) respectively.
     !> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
     !> order of the unitary matrix Q or P**H that is applied.
     !> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
     !> if nq >= k, Q = H(1) H(2) . . . H(k);
     !> if nq < k, Q = H(1) H(2) . . . H(nq-1).
     !> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
     !> if k < nq, P = G(1) G(2) . . . G(k);
     !> if k >= nq, P = G(1) G(2) . . . G(nq-1).

     pure subroutine la_wunmbr(vect,side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork, &
               info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: side,trans,vect
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: k,lda,ldc,lwork,m,n
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),c(ldc,*)
           complex(qp),intent(in) :: tau(*)
           complex(qp),intent(out) :: work(*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: applyq,left,lquery,notran
           character :: transt
           integer(ilp) :: i1,i2,iinfo,lwkopt,mi,nb,ni,nq,nw
           ! Intrinsic Functions
           intrinsic :: max,min
           ! Executable Statements
           ! test the input arguments
           info = 0
           applyq = la_lsame(vect,'Q')
           left = la_lsame(side,'L')
           notran = la_lsame(trans,'N')
           lquery = (lwork == -1)
           ! nq is the order of q or p and nw is the minimum dimension of work
           if (left) then
              nq = m
              nw = max(1,n)
           else
              nq = n
              nw = max(1,m)
           end if
           if (.not. applyq .and. .not. la_lsame(vect,'P')) then
              info = -1
           else if (.not. left .and. .not. la_lsame(side,'R')) then
              info = -2
           else if (.not. notran .and. .not. la_lsame(trans,'C')) then
              info = -3
           else if (m < 0) then
              info = -4
           else if (n < 0) then
              info = -5
           else if (k < 0) then
              info = -6
           else if ((applyq .and. lda < max(1,nq)) .or. (.not. applyq .and. lda < max(1,min(nq, &
                      k)))) then
              info = -8
           else if (ldc < max(1,m)) then
              info = -11
           else if (lwork < nw .and. .not. lquery) then
              info = -13
           end if
           if (info == 0) then
              if (m > 0 .and. n > 0) then
                 if (applyq) then
                    if (left) then
                       nb = la_ilaenv(1,'WUNMQR',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'WUNMQR',side//trans,m,n - 1,n - 1,-1)
                    end if
                 else
                    if (left) then
                       nb = la_ilaenv(1,'WUNMLQ',side//trans,m - 1,n,m - 1,-1)
                    else
                       nb = la_ilaenv(1,'WUNMLQ',side//trans,m,n - 1,n - 1,-1)
                    end if
                 end if
                 lwkopt = nw*nb
              else
                 lwkopt = 1
              end if
              work(1) = lwkopt
           end if
           if (info /= 0) then
              call la_xerbla('WUNMBR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (m == 0 .or. n == 0) return
           if (applyq) then
              ! apply q
              if (nq >= k) then
                 ! q was determined by a call to la_wgebrd with nq >= k
                 call la_wunmqr(side,trans,m,n,k,a,lda,tau,c,ldc,work,lwork,iinfo &
                           )
              else if (nq > 1) then
                 ! q was determined by a call to la_wgebrd with nq < k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_wunmqr(side,trans,mi,ni,nq - 1,a(2,1),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           else
              ! apply p
              if (notran) then
                 transt = 'C'
              else
                 transt = 'N'
              end if
              if (nq > k) then
                 ! p was determined by a call to la_wgebrd with nq > k
                 call la_wunmlq(side,transt,m,n,k,a,lda,tau,c,ldc,work,lwork, &
                           iinfo)
              else if (nq > 1) then
                 ! p was determined by a call to la_wgebrd with nq <= k
                 if (left) then
                    mi = m - 1
                    ni = n
                    i1 = 2
                    i2 = 1
                 else
                    mi = m
                    ni = n - 1
                    i1 = 1
                    i2 = 2
                 end if
                 call la_wunmlq(side,transt,mi,ni,nq - 1,a(1,2),lda,tau,c(i1,i2), &
                           ldc,work,lwork,iinfo)
              end if
           end if
           work(1) = lwkopt
           return
     end subroutine la_wunmbr
#endif

     !> CGSVJ0: is called from CGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as CGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_cgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_sp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(sp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(sp),intent(out) :: work(lwork)
           real(sp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(sp) :: aapq,ompq
           real(sp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,max,conjg,real,min,sign,sqrt
           ! from lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('CGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! .. row-cyclic jacobi svd algorithm with column pivoting ..
           emptsw = (n*(n - 1))/2
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           swband = 0
      ! [tp] swband is a tuning parameter [tp]. it is meaningful and effective
           ! if la_cgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_cgejsv. for sweeps i=1:swband the procedure
           ! works on pivots inside a band-like region around the diagonal.
           ! the boundaries are determined dynamically, based on the number of
           ! pivots above a threshold.
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = kbl**2
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           ! quasi block transformations, using the lower (upper) triangular
           ! structure of the input matrix. the quasi-block-cycling usually
           ! invokes cubic convergence. big part of this cycle is done inside
           ! canonical subspaces of dimensions less than m.
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_isamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_cswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_cswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          aapq = d(p)
                          d(p) = d(q)
                          d(q) = aapq
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! unfortunately, some blas implementations compute sncrm2(m,a(1,p),1)
              ! as sqrt(s=la_cdotc(m,a(1,p),1,a(1,p),1)), which may cause the result to
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and to
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_scnrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_scnrm2 is available, the if-then-else-end if
              ! below should be replaced with "aapp = la_scnrm2( m, a(1,p), 1 )".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_scnrm2(m,a(1,p),1)
                          else
                             temp1 = zero
                             aapp = one
                             call la_classq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_cdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_ccopy(m,a(1,p),1,work,1)
                                      call la_clascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_cdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_cdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aapp)/aaqq
                                   else
                                      call la_ccopy(m,a(1,q),1,work,1)
                                      call la_clascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_cdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg( cwork(p) ) * cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                 ! Rotate
      ! [rtd]      rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_crot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_crot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_crot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_crot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                      else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_ccopy(m,a(1,p),1,work,1)
                                      call la_clascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_clascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      call la_caxpy(m,-aapq,work,1,a(1,q),1)
                                      call la_clascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_scnrm2(m,a(1,q),1)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_classq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_scnrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_classq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_cdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_ccopy(m,a(1,p),1,work,1)
                                      call la_clascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_cdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_cdotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_ccopy(m,a(1,q),1,work,1)
                                      call la_clascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_cdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_crot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_crot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_crot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_crot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_ccopy(m,a(1,p),1,work,1)
                                         call la_clascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_clascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_caxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_clascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_ccopy(m,a(1,q),1,work,1)
                                         call la_clascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_clascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_caxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_clascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_scnrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_classq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_scnrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_classq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_scnrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_classq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=sp))*tol) .and. (real(n, &
                        KIND=sp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_isamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_cswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_cswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_cgsvj0
     !> ZGSVJ0: is called from ZGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as ZGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_zgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_dp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(dp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(dp),intent(out) :: work(lwork)
           real(dp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(dp) :: aapq,ompq
           real(dp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,max,conjg,real,min,sign,sqrt
           ! from lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('ZGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! .. row-cyclic jacobi svd algorithm with column pivoting ..
           emptsw = (n*(n - 1))/2
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           swband = 0
      ! [tp] swband is a tuning parameter [tp]. it is meaningful and effective
           ! if la_zgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_zgejsv. for sweeps i=1:swband the procedure
           ! works on pivots inside a band-like region around the diagonal.
           ! the boundaries are determined dynamically, based on the number of
           ! pivots above a threshold.
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = kbl**2
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           ! quasi block transformations, using the lower (upper) triangular
           ! structure of the input matrix. the quasi-block-cycling usually
           ! invokes cubic convergence. big part of this cycle is done inside
           ! canonical subspaces of dimensions less than m.
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_idamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_zswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_zswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          aapq = d(p)
                          d(p) = d(q)
                          d(q) = aapq
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! unfortunately, some blas implementations compute sncrm2(m,a(1,p),1)
              ! as sqrt(s=la_zdotc(m,a(1,p),1,a(1,p),1)), which may cause the result to
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and to
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_dznrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_dznrm2 is available, the if-then-else-end if
              ! below should be replaced with "aapp = la_dznrm2( m, a(1,p), 1 )".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_dznrm2(m,a(1,p),1)
                          else
                             temp1 = zero
                             aapp = one
                             call la_zlassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_zdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_zcopy(m,a(1,p),1,work,1)
                                      call la_zlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_zdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_zdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aapp)/aaqq
                                   else
                                      call la_zcopy(m,a(1,q),1,work,1)
                                      call la_zlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_zdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg( cwork(p) ) * cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                 ! Rotate
      ! [rtd]      rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_zrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_zrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_zrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_zrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                      else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_zcopy(m,a(1,p),1,work,1)
                                      call la_zlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_zlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      call la_zaxpy(m,-aapq,work,1,a(1,q),1)
                                      call la_zlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_dznrm2(m,a(1,q),1)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_zlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_dznrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_zlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_zdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_zcopy(m,a(1,p),1,work,1)
                                      call la_zlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_zdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_zdotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_zcopy(m,a(1,q),1,work,1)
                                      call la_zlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_zdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_zrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_zrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_zrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_zrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_zcopy(m,a(1,p),1,work,1)
                                         call la_zlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_zlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_zaxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_zlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_zcopy(m,a(1,q),1,work,1)
                                         call la_zlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_zlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_zaxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_zlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_dznrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_zlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_dznrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_zlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_dznrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_zlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=dp))*tol) .and. (real(n, &
                        KIND=dp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_idamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_zswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_zswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_zgsvj0
#ifdef LA_WITH_XDP
     !> YGSVJ0: is called from YGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as YGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_ygsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_xdp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(xdp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(xdp),intent(out) :: work(lwork)
           real(xdp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(xdp) :: aapq,ompq
           real(xdp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,max,conjg,real,min,sign,sqrt
           ! from lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('YGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! .. row-cyclic jacobi svd algorithm with column pivoting ..
           emptsw = (n*(n - 1))/2
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           swband = 0
      ! [tp] swband is a tuning parameter [tp]. it is meaningful and effective
           ! if la_ygesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_ygejsv. for sweeps i=1:swband the procedure
           ! works on pivots inside a band-like region around the diagonal.
           ! the boundaries are determined dynamically, based on the number of
           ! pivots above a threshold.
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = kbl**2
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           ! quasi block transformations, using the lower (upper) triangular
           ! structure of the input matrix. the quasi-block-cycling usually
           ! invokes cubic convergence. big part of this cycle is done inside
           ! canonical subspaces of dimensions less than m.
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_ixamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_yswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_yswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          aapq = d(p)
                          d(p) = d(q)
                          d(q) = aapq
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! unfortunately, some blas implementations compute sncrm2(m,a(1,p),1)
              ! as sqrt(s=la_ydotc(m,a(1,p),1,a(1,p),1)), which may cause the result to
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and to
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_xynrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_xynrm2 is available, the if-then-else-end if
              ! below should be replaced with "aapp = la_xynrm2( m, a(1,p), 1 )".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_xynrm2(m,a(1,p),1)
                          else
                             temp1 = zero
                             aapp = one
                             call la_ylassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_ydotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_ycopy(m,a(1,p),1,work,1)
                                      call la_ylascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_ydotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_ydotc(m,a(1,p),1,a(1,q),1)/ &
                                                aapp)/aaqq
                                   else
                                      call la_ycopy(m,a(1,q),1,work,1)
                                      call la_ylascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_ydotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg( cwork(p) ) * cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                 ! Rotate
      ! [rtd]      rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_yrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_yrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_yrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_yrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                      else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_ycopy(m,a(1,p),1,work,1)
                                      call la_ylascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_ylascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      call la_yaxpy(m,-aapq,work,1,a(1,q),1)
                                      call la_ylascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_xynrm2(m,a(1,q),1)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_ylassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_xynrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_ylassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_ydotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_ycopy(m,a(1,p),1,work,1)
                                      call la_ylascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_ydotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_ydotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_ycopy(m,a(1,q),1,work,1)
                                      call la_ylascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_ydotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_yrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_yrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_yrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_yrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_ycopy(m,a(1,p),1,work,1)
                                         call la_ylascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_ylascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_yaxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_ylascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_ycopy(m,a(1,q),1,work,1)
                                         call la_ylascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_ylascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_yaxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_ylascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_xynrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_ylassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_xynrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_ylassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_xynrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_ylassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=xdp))*tol) .and. (real(n, &
                        KIND=xdp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_ixamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_yswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_yswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_ygsvj0
#endif
#ifdef LA_WITH_QP
     !> WGSVJ0: is called from WGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as WGESVJ does, but
     !> it does not check convergence (stopping criterion). Few tuning
     !> parameters (marked by [TP]) are available for the implementer.

     pure subroutine la_wgsvj0(jobv,m,n,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_qp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,nsweep
           real(qp),intent(in) :: eps,sfmin,tol
           character,intent(in) :: jobv
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(qp),intent(out) :: work(lwork)
           real(qp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(qp) :: aapq,ompq
           real(qp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,ierr,igl,ijblsk,ir1,iswrot,jbc,jgl,kbl, &
                     lkahead,mvl,nbl,notrot,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,max,conjg,real,min,sign,sqrt
           ! from lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (lda < m) then
              info = -5
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -8
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -10
           else if (tol <= eps) then
              info = -13
           else if (nsweep < 0) then
              info = -14
           else if (lwork < m) then
              info = -16
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('WGSVJ0',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! .. row-cyclic jacobi svd algorithm with column pivoting ..
           emptsw = (n*(n - 1))/2
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           swband = 0
      ! [tp] swband is a tuning parameter [tp]. it is meaningful and effective
           ! if la_wgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_wgejsv. for sweeps i=1:swband the procedure
           ! works on pivots inside a band-like region around the diagonal.
           ! the boundaries are determined dynamically, based on the number of
           ! pivots above a threshold.
           kbl = min(8,n)
      ! [tp] kbl is a tuning parameter that defines the tile size in the
           ! tiling of the p-q loops of pivot pairs. in general, an optimal
           ! value of kbl depends on the matrix dimensions and on the
           ! parameters of the computer's memory.
           nbl = n/kbl
           if ((nbl*kbl) /= n) nbl = nbl + 1
           blskip = kbl**2
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           lkahead = 1
      ! [tp] lkahead is a tuning parameter.
           ! quasi block transformations, using the lower (upper) triangular
           ! structure of the input matrix. the quasi-block-cycling usually
           ! invokes cubic convergence. big part of this cycle is done inside
           ! canonical subspaces of dimensions less than m.
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nbl
                 igl = (ibr - 1)*kbl + 1
                 loop_1002: do ir1 = 0,min(lkahead,nbl - ibr)
                    igl = igl + ir1*kbl
                    loop_2001: do p = igl,min(igl + kbl - 1,n - 1)
           ! .. de rijk's pivoting
                       q = la_iqamax(n - p + 1,sva(p),1) + p - 1
                       if (p /= q) then
                          call la_wswap(m,a(1,p),1,a(1,q),1)
                          if (rsvec) call la_wswap(mvl,v(1,p),1,v(1,q),1)
                          temp1 = sva(p)
                          sva(p) = sva(q)
                          sva(q) = temp1
                          aapq = d(p)
                          d(p) = d(q)
                          d(q) = aapq
                       end if
                       if (ir1 == 0) then
              ! column norms are periodically updated by explicit
              ! norm computation.
              ! caveat:
              ! unfortunately, some blas implementations compute sncrm2(m,a(1,p),1)
              ! as sqrt(s=la_wdotc(m,a(1,p),1,a(1,p),1)), which may cause the result to
              ! overflow for ||a(:,p)||_2 > sqrt(overflow_threshold), and to
              ! underflow for ||a(:,p)||_2 < sqrt(underflow_threshold).
              ! hence, la_qwnrm2 cannot be trusted, not even in the case when
              ! the true norm is far from the under(over)flow boundaries.
              ! if properly implemented la_qwnrm2 is available, the if-then-else-end if
              ! below should be replaced with "aapp = la_qwnrm2( m, a(1,p), 1 )".
                          if ((sva(p) < rootbig) .and. (sva(p) > rootsfmin)) then
                             sva(p) = la_qwnrm2(m,a(1,p),1)
                          else
                             temp1 = zero
                             aapp = one
                             call la_wlassq(m,a(1,p),1,temp1,aapp)
                             sva(p) = temp1*sqrt(aapp)
                          end if
                          aapp = sva(p)
                       else
                          aapp = sva(p)
                       end if
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2002: do q = p + 1,min(igl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
                                if (aaqq >= one) then
                                   rotok = (small*aapp) <= aaqq
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_wdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_wcopy(m,a(1,p),1,work,1)
                                      call la_wlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_wdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   rotok = aapp <= (aaqq/small)
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_wdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aapp)/aaqq
                                   else
                                      call la_wcopy(m,a(1,q),1,work,1)
                                      call la_wlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_wdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg( cwork(p) ) * cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                 ! Rotate
      ! [rtd]      rotated = rotated + one
                                   if (ir1 == 0) then
                                      notrot = 0
                                      pskipped = 0
                                      iswrot = iswrot + 1
                                   end if
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_wrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_wrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_wrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_wrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                      else
                    ! .. have to use modified gram-schmidt like transformation
                                      call la_wcopy(m,a(1,p),1,work,1)
                                      call la_wlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      call la_wlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                lda,ierr)
                                      call la_waxpy(m,-aapq,work,1,a(1,q),1)
                                      call la_wlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                lda,ierr)
                                      sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))
                                      mxsinj = max(mxsinj,sfmin)
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! recompute sva(q), sva(p).
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_qwnrm2(m,a(1,q),1)
                                      else
                                         t = zero
                                         aaqq = one
                                         call la_wlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0) <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_qwnrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_wlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                                else
              ! a(:,p) and a(:,q) already numerically orthogonal
                                   if (ir1 == 0) notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                end if
                             else
              ! a(:,q) is zero column
                                if (ir1 == 0) notrot = notrot + 1
                                pskipped = pskipped + 1
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                if (ir1 == 0) aapp = -aapp
                                notrot = 0
                                go to 2103
                             end if
                          end do loop_2002
           ! end q-loop
           2103 continue
           ! bailed out of q-loop
                          sva(p) = aapp
                       else
                          sva(p) = aapp
                          if ((ir1 == 0) .and. (aapp == zero)) notrot = notrot + min(igl + kbl - 1, &
                                    n) - p
                       end if
                    end do loop_2001
           ! end of the p-loop
           ! end of doing the block ( ibr, ibr )
                 end do loop_1002
           ! end of ir1-loop
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                 loop_2010: do jbc = ibr + 1,nbl
                    jgl = (jbc - 1)*kbl + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_wdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_wcopy(m,a(1,p),1,work,1)
                                      call la_wlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_wdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_wdotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_wcopy(m,a(1,q),1,work,1)
                                      call la_wlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_wdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_wrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_wrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_wrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_wrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_wcopy(m,a(1,p),1,work,1)
                                         call la_wlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_wlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_waxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_wlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_wcopy(m,a(1,q),1,work,1)
                                         call la_wlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_wlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_waxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_wlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_qwnrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_wlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_qwnrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_wlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_qwnrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_wlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=qp))*tol) .and. (real(n, &
                        KIND=qp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_iqamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_wswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_wswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_wgsvj0
#endif

     !> CGSVJ1: is called from CGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as CGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> CGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_cgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_sp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(sp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           complex(sp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(sp),intent(out) :: work(lwork)
           real(sp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(sp) :: aapq,ompq
           real(sp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,max,conjg,real,min,sign,sqrt
           ! From Lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('CGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           ! large = big / sqrt( real( m*n,KIND=sp) )
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_cgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_cgejsv.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                  ! do 2010 jbc = ibr + 1, nbl
                 loop_2010: do jbc = 1,nblc
                    jgl = (jbc - 1)*kbl + n1 + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_cdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_ccopy(m,a(1,p),1,work,1)
                                      call la_clascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_cdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_cdotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_ccopy(m,a(1,q),1,work,1)
                                      call la_clascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_cdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_crot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_crot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_crot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_crot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_ccopy(m,a(1,p),1,work,1)
                                         call la_clascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_clascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_caxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_clascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_ccopy(m,a(1,q),1,work,1)
                                         call la_clascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_clascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_caxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_clascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_scnrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_classq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_scnrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_classq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_scnrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_classq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=sp))*tol) .and. (real(n, &
                        KIND=sp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_isamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_cswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_cswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_cgsvj1
     !> ZGSVJ1: is called from ZGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as ZGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> ZGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_zgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_dp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(dp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           complex(dp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(dp),intent(out) :: work(lwork)
           real(dp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(dp) :: aapq,ompq
           real(dp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real,min,sign,sqrt
           ! From Lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('ZGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           ! large = big / sqrt( real( m*n,KIND=dp) )
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_zgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_zgejsv.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                  ! do 2010 jbc = ibr + 1, nbl
                 loop_2010: do jbc = 1,nblc
                    jgl = (jbc - 1)*kbl + n1 + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_zdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_zcopy(m,a(1,p),1,work,1)
                                      call la_zlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_zdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_zdotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_zcopy(m,a(1,q),1,work,1)
                                      call la_zlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_zdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_zrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_zrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_zrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_zrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_zcopy(m,a(1,p),1,work,1)
                                         call la_zlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_zlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_zaxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_zlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_zcopy(m,a(1,q),1,work,1)
                                         call la_zlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_zlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_zaxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_zlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_dznrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_zlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_dznrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_zlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_dznrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_zlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=dp))*tol) .and. (real(n, &
                        KIND=dp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_idamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_zswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_zswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_zgsvj1
#ifdef LA_WITH_XDP
     !> YGSVJ1: is called from YGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as YGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> YGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_ygsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_xdp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(xdp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           complex(xdp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(xdp),intent(out) :: work(lwork)
           real(xdp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(xdp) :: aapq,ompq
           real(xdp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real,min,sign,sqrt
           ! From Lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('YGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           ! large = big / sqrt( real( m*n,KIND=xdp) )
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_ygesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_ygejsv.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                  ! do 2010 jbc = ibr + 1, nbl
                 loop_2010: do jbc = 1,nblc
                    jgl = (jbc - 1)*kbl + n1 + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_ydotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_ycopy(m,a(1,p),1,work,1)
                                      call la_ylascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_ydotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_ydotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_ycopy(m,a(1,q),1,work,1)
                                      call la_ylascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_ydotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_yrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_yrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_yrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_yrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_ycopy(m,a(1,p),1,work,1)
                                         call la_ylascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_ylascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_yaxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_ylascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_ycopy(m,a(1,q),1,work,1)
                                         call la_ylascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_ylascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_yaxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_ylascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_xynrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_ylassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_xynrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_ylassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_xynrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_ylassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=xdp))*tol) .and. (real(n, &
                        KIND=xdp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_ixamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_yswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_yswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_ygsvj1
#endif
#ifdef LA_WITH_QP
     !> WGSVJ1: is called from WGESVJ as a pre-processor and that is its main
     !> purpose. It applies Jacobi rotations in the same way as WGESVJ does, but
     !> it targets only particular pivots and it does not check convergence
     !> (stopping criterion). Few tuning parameters (marked by [TP]) are
     !> available for the implementer.
     !> Further Details
     !>
     !> WGSVJ1 applies few sweeps of Jacobi rotations in the column space of
     !> the input M-by-N matrix A. The pivot pairs are taken from the (1,2)
     !> off-diagonal block in the corresponding N-by-N Gram matrix A^T * A. The
     !> block-entries (tiles) of the (1,2) off-diagonal block are marked by the
     !> [x]'s in the following scheme:
     !> | *  *  * [x] [x] [x]|
     !> | *  *  * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
     !> | *  *  * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> |[x] [x] [x] *  *  * |
     !> In terms of the columns of A, the first N1 columns are rotated 'against'
     !> the remaining N-N1 columns, trying to increase the angle between the
     !> corresponding subspaces. The off-diagonal block is N1-by(N-N1) and it is
     !> tiled using quadratic tiles of side KBL. Here, KBL is a tuning parameter.
     !> The number of sweeps is given in NSWEEP and the orthogonality threshold
     !> is given in TOL.

     pure subroutine la_wgsvj1(jobv,m,n,n1,a,lda,d,sva,mv,v,ldv,eps,sfmin,tol, &
               nsweep,work,lwork,info)
        use la_constants_qp,only:zero,half,one
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           real(qp),intent(in) :: eps,sfmin,tol
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: lda,ldv,lwork,m,mv,n,n1,nsweep
           character,intent(in) :: jobv
           ! Array Arguments
           complex(qp),intent(inout) :: a(lda,*),d(n),v(ldv,*)
           complex(qp),intent(out) :: work(lwork)
           real(qp),intent(inout) :: sva(n)
        ! =====================================================================

           ! Local Scalars
           complex(qp) :: aapq,ompq
           real(qp) :: aapp,aapp0,aapq1,aaqq,apoaq,aqoap,big,bigtheta,cs,mxaapq,mxsinj, &
                     rootbig,rooteps,rootsfmin,roottol,small,sn,t,temp1,theta,thsign
           integer(ilp) :: blskip,emptsw,i,ibr,igl,ierr,ijblsk,iswrot,jbc,jgl,kbl,mvl, &
                     notrot,nblc,nblr,p,pskipped,q,rowskip,swband
           logical(lk) :: applv,rotok,rsvec
           ! Intrinsic Functions
           intrinsic :: abs,conjg,max,real,min,sign,sqrt
           ! From Lapack
           ! Executable Statements
           ! test the input parameters.
           applv = la_lsame(jobv,'A')
           rsvec = la_lsame(jobv,'V')
           if (.not. (rsvec .or. applv .or. la_lsame(jobv,'N'))) then
              info = -1
           else if (m < 0) then
              info = -2
           else if ((n < 0) .or. (n > m)) then
              info = -3
           else if (n1 < 0) then
              info = -4
           else if (lda < m) then
              info = -6
           else if ((rsvec .or. applv) .and. (mv < 0)) then
              info = -9
           else if ((rsvec .and. (ldv < n)) .or. (applv .and. (ldv < mv))) then
              info = -11
           else if (tol <= eps) then
              info = -14
           else if (nsweep < 0) then
              info = -15
           else if (lwork < m) then
              info = -17
           else
              info = 0
           end if
           ! #:(
           if (info /= 0) then
              call la_xerbla('WGSVJ1',-info)
              return
           end if
           if (rsvec) then
              mvl = n
           else if (applv) then
              mvl = mv
           end if
           rsvec = rsvec .or. applv
           rooteps = sqrt(eps)
           rootsfmin = sqrt(sfmin)
           small = sfmin/eps
           big = one/sfmin
           rootbig = one/rootsfmin
           ! large = big / sqrt( real( m*n,KIND=qp) )
           bigtheta = one/rooteps
           roottol = sqrt(tol)
           ! Initialize The Right Singular Vector Matrix
           ! rsvec = la_lsame( jobv, 'y' )
           emptsw = n1*(n - n1)
           notrot = 0
           ! .. row-cyclic pivot strategy with de rijk's pivoting ..
           kbl = min(8,n)
           nblr = n1/kbl
           if ((nblr*kbl) /= n1) nblr = nblr + 1
           ! .. the tiling is nblr-by-nblc [tiles]
           nblc = (n - n1)/kbl
           if ((nblc*kbl) /= (n - n1)) nblc = nblc + 1
           blskip = (kbl**2) + 1
      ! [tp] blkskip is a tuning parameter that depends on swband and kbl.
           rowskip = min(5,kbl)
      ! [tp] rowskip is a tuning parameter.
           swband = 0
      ! [tp] swband is a tuning parameter. it is meaningful and effective
           ! if la_wgesvj is used as a computational routine in the preconditioned
           ! jacobi svd algorithm la_wgejsv.
           ! | *   *   * [x] [x] [x]|
           ! | *   *   * [x] [x] [x]|    row-cycling in the nblr-by-nblc [x] blocks.
           ! | *   *   * [x] [x] [x]|    row-cyclic pivoting inside each [x] block.
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           ! |[x] [x] [x] *   *   * |
           loop_1993: do i = 1,nsweep
           ! .. go go go ...
              mxaapq = zero
              mxsinj = zero
              iswrot = 0
              notrot = 0
              pskipped = 0
           ! each sweep is unrolled using kbl-by-kbl tiles over the pivot pairs
           ! 1 <= p < q <= n. this is the first step toward a blocked implementation
           ! of the rotations. new implementation, based on block transformations,
           ! is under development.
              loop_2000: do ibr = 1,nblr
                 igl = (ibr - 1)*kbl + 1
       ! ... go to the off diagonal blocks
                 igl = (ibr - 1)*kbl + 1
                  ! do 2010 jbc = ibr + 1, nbl
                 loop_2010: do jbc = 1,nblc
                    jgl = (jbc - 1)*kbl + n1 + 1
              ! doing the block at ( ibr, jbc )
                    ijblsk = 0
                    loop_2100: do p = igl,min(igl + kbl - 1,n1)
                       aapp = sva(p)
                       if (aapp > zero) then
                          pskipped = 0
                          loop_2200: do q = jgl,min(jgl + kbl - 1,n)
                             aaqq = sva(q)
                             if (aaqq > zero) then
                                aapp0 = aapp
           ! M X 2 Jacobi Svd
              ! safe gram matrix computation
                                if (aaqq >= one) then
                                   if (aapp >= aaqq) then
                                      rotok = (small*aapp) <= aaqq
                                   else
                                      rotok = (small*aaqq) <= aapp
                                   end if
                                   if (aapp < (big/aaqq)) then
                                      aapq = (la_wdotc(m,a(1,p),1,a(1,q),1)/ &
                                                aaqq)/aapp
                                   else
                                      call la_wcopy(m,a(1,p),1,work,1)
                                      call la_wlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_wdotc(m,work,1,a(1,q),1)/ &
                                                aaqq
                                   end if
                                else
                                   if (aapp >= aaqq) then
                                      rotok = aapp <= (aaqq/small)
                                   else
                                      rotok = aaqq <= (aapp/small)
                                   end if
                                   if (aapp > (small/aaqq)) then
                                      aapq = (la_wdotc(m,a(1,p),1,a(1,q),1)/max( &
                                                aaqq,aapp))/min(aaqq,aapp)
                                   else
                                      call la_wcopy(m,a(1,q),1,work,1)
                                      call la_wlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                ierr)
                                      aapq = la_wdotc(m,a(1,p),1,work,1)/ &
                                                aapp
                                   end if
                                end if
                                 ! aapq = aapq * conjg(cwork(p))*cwork(q)
                                aapq1 = -abs(aapq)
                                mxaapq = max(mxaapq,-aapq1)
              ! to rotate or not to rotate, that is the question ...
                                if (abs(aapq1) > tol) then
                                   ompq = aapq/abs(aapq)
                                   notrot = 0
      ! [rtd]      rotated  = rotated + 1
                                   pskipped = 0
                                   iswrot = iswrot + 1
                                   if (rotok) then
                                      aqoap = aaqq/aapp
                                      apoaq = aapp/aaqq
                                      theta = -half*abs(aqoap - apoaq)/aapq1
                                      if (aaqq > aapp0) theta = -theta
                                      if (abs(theta) > bigtheta) then
                                         t = half/theta
                                         cs = one
                                         call la_wrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *t)
                                         if (rsvec) then
                                             call la_wrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*t)
                                         end if
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         mxsinj = max(mxsinj,abs(t))
                                      else
                       ! Choose Correct Signum For Theta And Rotate
                                         thsign = -sign(one,aapq1)
                                         if (aaqq > aapp0) thsign = -thsign
                                         t = one/(theta + thsign*sqrt(one + theta*theta))

                                         cs = sqrt(one/(one + t*t))
                                         sn = t*cs
                                         mxsinj = max(mxsinj,abs(sn))
                                         sva(q) = aaqq*sqrt(max(zero,one + t*apoaq*aapq1))

                                         aapp = aapp*sqrt(max(zero,one - t*aqoap*aapq1))
                                         call la_wrot(m,a(1,p),1,a(1,q),1,cs,conjg(ompq) &
                                                   *sn)
                                         if (rsvec) then
                                             call la_wrot(mvl,v(1,p),1,v(1,q),1,cs, &
                                                       conjg(ompq)*sn)
                                         end if
                                      end if
                                      d(p) = -d(q)*ompq
                                   else
                    ! .. have to use modified gram-schmidt like transformation
                                    if (aapp > aaqq) then
                                         call la_wcopy(m,a(1,p),1,work,1)
                                         call la_wlascl('G',0,0,aapp,one,m,1,work,lda, &
                                                   ierr)
                                         call la_wlascl('G',0,0,aaqq,one,m,1,a(1,q), &
                                                    lda,ierr)
                                         call la_waxpy(m,-aapq,work,1,a(1,q),1)

                                         call la_wlascl('G',0,0,one,aaqq,m,1,a(1,q), &
                                                    lda,ierr)
                                         sva(q) = aaqq*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    else
                                        call la_wcopy(m,a(1,q),1,work,1)
                                         call la_wlascl('G',0,0,aaqq,one,m,1,work,lda, &
                                                   ierr)
                                         call la_wlascl('G',0,0,aapp,one,m,1,a(1,p), &
                                                    lda,ierr)
                                         call la_waxpy(m,-conjg(aapq),work,1,a(1,p),1 &
                                                   )
                                         call la_wlascl('G',0,0,one,aapp,m,1,a(1,p), &
                                                    lda,ierr)
                                         sva(p) = aapp*sqrt(max(zero,one - aapq1*aapq1))

                                         mxsinj = max(mxsinj,sfmin)
                                    end if
                                   end if
                 ! end if rotok then ... else
                 ! in the case of cancellation in updating sva(q), sva(p)
                 ! .. recompute sva(q), sva(p)
                                   if ((sva(q)/aaqq)**2 <= rooteps) then
                                      if ((aaqq < rootbig) .and. (aaqq > rootsfmin)) then
                                         sva(q) = la_qwnrm2(m,a(1,q),1)
                                       else
                                         t = zero
                                         aaqq = one
                                         call la_wlassq(m,a(1,q),1,t,aaqq)
                                         sva(q) = t*sqrt(aaqq)
                                      end if
                                   end if
                                   if ((aapp/aapp0)**2 <= rooteps) then
                                      if ((aapp < rootbig) .and. (aapp > rootsfmin)) then
                                         aapp = la_qwnrm2(m,a(1,p),1)
                                      else
                                         t = zero
                                         aapp = one
                                         call la_wlassq(m,a(1,p),1,t,aapp)
                                         aapp = t*sqrt(aapp)
                                      end if
                                      sva(p) = aapp
                                   end if
                    ! end of ok rotation
                                else
                                   notrot = notrot + 1
      ! [rtd]      skipped  = skipped  + 1
                                   pskipped = pskipped + 1
                                   ijblsk = ijblsk + 1
                                end if
                             else
                                notrot = notrot + 1
                                pskipped = pskipped + 1
                                ijblsk = ijblsk + 1
                             end if
                             if ((i <= swband) .and. (ijblsk >= blskip)) then
                                sva(p) = aapp
                                notrot = 0
                                go to 2011
                             end if
                             if ((i <= swband) .and. (pskipped > rowskip)) then
                                aapp = -aapp
                                notrot = 0
                                go to 2203
                             end if
                          end do loop_2200
              ! end of the q-loop
              2203 continue
                          sva(p) = aapp
                       else
                          if (aapp == zero) notrot = notrot + min(jgl + kbl - 1,n) - jgl + 1
                          if (aapp < zero) notrot = 0
                       end if
                    end do loop_2100
           ! end of the p-loop
                 end do loop_2010
           ! end of the jbc-loop
           2011 continue
      ! 2011 bailed out of the jbc-loop
                 do p = igl,min(igl + kbl - 1,n)
                    sva(p) = abs(sva(p))
                 end do
      ! **
              end do loop_2000
      ! 2000 :: end of the ibr-loop
           ! .. update sva(n)
              if ((sva(n) < rootbig) .and. (sva(n) > rootsfmin)) then
                 sva(n) = la_qwnrm2(m,a(1,n),1)
              else
                 t = zero
                 aapp = one
                 call la_wlassq(m,a(1,n),1,t,aapp)
                 sva(n) = t*sqrt(aapp)
              end if
           ! additional steering devices
              if ((i < swband) .and. ((mxaapq <= roottol) .or. (iswrot <= n))) swband = i
              if ((i > swband + 1) .and. (mxaapq < sqrt(real(n,KIND=qp))*tol) .and. (real(n, &
                        KIND=qp)*mxaapq*mxsinj < tol)) then
                 go to 1994
              end if
              if (notrot >= emptsw) go to 1994
           end do loop_1993
           ! end i=1:nsweep loop
       ! #:( reaching this point means that the procedure has not converged.
           info = nsweep - 1
           go to 1995
           1994 continue
       ! #:) reaching this point means numerical convergence after the i-th
           ! sweep.
           info = 0
       ! #:) info = 0 confirms successful iterations.
       1995 continue
           ! sort the vector sva() of column norms.
           do p = 1,n - 1
              q = la_iqamax(n - p + 1,sva(p),1) + p - 1
              if (p /= q) then
                 temp1 = sva(p)
                 sva(p) = sva(q)
                 sva(q) = temp1
                 aapq = d(p)
                 d(p) = d(q)
                 d(q) = aapq
                 call la_wswap(m,a(1,p),1,a(1,q),1)
                 if (rsvec) call la_wswap(mvl,v(1,p),1,v(1,q),1)
              end if
           end do
           return
     end subroutine la_wgsvj1
#endif

end module la_lapack_svd_comp
