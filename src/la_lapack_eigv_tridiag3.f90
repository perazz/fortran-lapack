!> Symmetric tridiagonal eigenvalue drivers: divide and conquer, MRRR, bisection and inverse iteration
module la_lapack_eigv_tridiag3
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
     use la_lapack_eigv_tridiag
     use la_lapack_eigv_tridiag2
     use la_lapack_solve_chol_comp
     use la_lapack_svd_bidiag_qr
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_sstebz
     public :: la_sstein
     public :: la_ssterf
     public :: la_sstev
     public :: la_sstevx
     public :: la_sstedc
     public :: la_sstevd
     public :: la_spteqr
     public :: la_sstegr
     public :: la_sstemr
     public :: la_sstevr
     public :: la_dstebz
     public :: la_dstein
     public :: la_dsterf
     public :: la_dstev
     public :: la_dstevx
     public :: la_dstedc
     public :: la_dstevd
     public :: la_dpteqr
     public :: la_dstegr
     public :: la_dstemr
     public :: la_dstevr
     public :: la_qstebz
     public :: la_qstein
     public :: la_qsterf
     public :: la_qstev
     public :: la_qstevx
     public :: la_qstedc
     public :: la_qstevd
     public :: la_qpteqr
     public :: la_qstegr
     public :: la_qstemr
     public :: la_qstevr
     public :: la_cstein
     public :: la_cpteqr
     public :: la_cstemr
     public :: la_cstedc
     public :: la_cstegr
     public :: la_zstein
     public :: la_zpteqr
     public :: la_zstemr
     public :: la_zstedc
     public :: la_zstegr
     public :: la_wstein
     public :: la_wpteqr
     public :: la_wstemr
     public :: la_wstedc
     public :: la_wstegr

     contains

     !> SSTEBZ: computes the eigenvalues of a symmetric tridiagonal
     !> matrix T.  The user may ask for all eigenvalues, all eigenvalues
     !> in the half-open interval (VL, VU], or the IL-th through IU-th
     !> eigenvalues.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_sstebz(range,order,n,vl,vu,il,iu,abstol,d,e,m,nsplit,w, &
               iblock,isplit,work,iwork,info)
        use la_constants_sp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: order,range
           integer(ilp),intent(in) :: il,iu,n
           integer(ilp),intent(out) :: info,m,nsplit
           real(sp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),isplit(*),iwork(*)
           real(sp),intent(in) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: fudge = 2.1_sp
           real(sp),parameter :: relfac = 2.0_sp

           ! Local Scalars
           logical(lk) :: ncnvrg,toofew
           integer(ilp) :: ib,ibegin,idiscl,idiscu,ie,iend,iinfo,im,in,ioff,iorder, &
                     iout,irange,itmax,itmp1,iw,iwoff,j,jb,jdisc,je,nb,nwl,nwu
           real(sp) :: atoli,bnorm,gl,gu,pivmin,rtoli,safemn,tmp1,tmp2,tnorm,ulp,wkill, &
                      wl,wlu,wu,wul
           ! Local Arrays
           integer(ilp) :: idumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           ! decode range
           if (la_lsame(range,'A')) then
              irange = 1
           else if (la_lsame(range,'V')) then
              irange = 2
           else if (la_lsame(range,'I')) then
              irange = 3
           else
              irange = 0
           end if
           ! decode order
           if (la_lsame(order,'B')) then
              iorder = 2
           else if (la_lsame(order,'E')) then
              iorder = 1
           else
              iorder = 0
           end if
           ! check for errors
           if (irange <= 0) then
              info = -1
           else if (iorder <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (irange == 2) then
              if (vl >= vu) info = -5
           else if (irange == 3 .and. (il < 1 .or. il > max(1,n))) then
              info = -6
           else if (irange == 3 .and. (iu < min(n,il) .or. iu > n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('SSTEBZ',-info)
              return
           end if
           ! initialize error flags
           info = 0
           ncnvrg = .false.
           toofew = .false.
           ! quick return if possible
           m = 0
           if (n == 0) return
           ! simplifications:
           if (irange == 3 .and. il == 1 .and. iu == n) irange = 1
           ! get machine constants
           ! nb is the minimum vector length for vector bisection, or 0
           ! if only scalar is to be done.
           safemn = la_slamch('S')
           ulp = la_slamch('P')
           rtoli = ulp*relfac
           nb = la_ilaenv(1,'SSTEBZ',' ',n,-1,-1,-1)
           if (nb <= 1) nb = 0
           ! special case when n=1
           if (n == 1) then
              nsplit = 1
              isplit(1) = 1
              if (irange == 2 .and. (vl >= d(1) .or. vu < d(1))) then
                 m = 0
              else
                 w(1) = d(1)
                 iblock(1) = 1
                 m = 1
              end if
              return
           end if
           ! compute splitting points
           nsplit = 1
           work(n) = zero
           pivmin = one
           do j = 2,n
              tmp1 = e(j - 1)**2
              if (abs(d(j)*d(j - 1))*ulp**2 + safemn > tmp1) then
                 isplit(nsplit) = j - 1
                 nsplit = nsplit + 1
                 work(j - 1) = zero
              else
                 work(j - 1) = tmp1
                 pivmin = max(pivmin,tmp1)
              end if
           end do
           isplit(nsplit) = n
           pivmin = pivmin*safemn
           ! compute interval and atoli
           if (irange == 3) then
              ! range='i': compute the interval containing eigenvalues
                         ! il through iu.
              ! compute gershgorin interval for entire (split) matrix
              ! and use it as the initial interval
              gu = d(1)
              gl = d(1)
              tmp1 = zero
              do j = 1,n - 1
                 tmp2 = sqrt(work(j))
                 gu = max(gu,d(j) + tmp1 + tmp2)
                 gl = min(gl,d(j) - tmp1 - tmp2)
                 tmp1 = tmp2
              end do
              gu = max(gu,d(n) + tmp1)
              gl = min(gl,d(n) - tmp1)
              tnorm = max(abs(gl),abs(gu))
              gl = gl - fudge*tnorm*ulp*n - fudge*two*pivmin
              gu = gu + fudge*tnorm*ulp*n + fudge*pivmin
              ! compute iteration parameters
              itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
              if (abstol <= zero) then
                 atoli = ulp*tnorm
              else
                 atoli = abstol
              end if
              work(n + 1) = gl
              work(n + 2) = gl
              work(n + 3) = gu
              work(n + 4) = gu
              work(n + 5) = gl
              work(n + 6) = gu
              iwork(1) = -1
              iwork(2) = -1
              iwork(3) = n + 1
              iwork(4) = n + 1
              iwork(5) = il - 1
              iwork(6) = iu
              call la_slaebz(3,itmax,n,2,2,nb,atoli,rtoli,pivmin,d,e,work,iwork( &
                        5),work(n + 1),work(n + 5),iout,iwork,w,iblock,iinfo)
              if (iwork(6) == iu) then
                 wl = work(n + 1)
                 wlu = work(n + 3)
                 nwl = iwork(1)
                 wu = work(n + 4)
                 wul = work(n + 2)
                 nwu = iwork(4)
              else
                 wl = work(n + 2)
                 wlu = work(n + 4)
                 nwl = iwork(2)
                 wu = work(n + 3)
                 wul = work(n + 1)
                 nwu = iwork(3)
              end if
              if (nwl < 0 .or. nwl >= n .or. nwu < 1 .or. nwu > n) then
                 info = 4
                 return
              end if
           else
              ! range='a' or 'v' -- set atoli
              tnorm = max(abs(d(1)) + abs(e(1)),abs(d(n)) + abs(e(n - 1)))
              do j = 2,n - 1
                 tnorm = max(tnorm,abs(d(j)) + abs(e(j - 1)) + abs(e(j)))
              end do
              if (abstol <= zero) then
                 atoli = ulp*tnorm
              else
                 atoli = abstol
              end if
              if (irange == 2) then
                 wl = vl
                 wu = vu
              else
                 wl = zero
                 wu = zero
              end if
           end if
           ! find eigenvalues -- loop over blocks and recompute nwl and nwu.
           ! nwl accumulates the number of eigenvalues .le. wl,
           ! nwu accumulates the number of eigenvalues .le. wu
           m = 0
           iend = 0
           info = 0
           nwl = 0
           nwu = 0
           loop_70: do jb = 1,nsplit
              ioff = iend
              ibegin = ioff + 1
              iend = isplit(jb)
              in = iend - ioff
              if (in == 1) then
                 ! special case -- in=1
                 if (irange == 1 .or. wl >= d(ibegin) - pivmin) nwl = nwl + 1
                 if (irange == 1 .or. wu >= d(ibegin) - pivmin) nwu = nwu + 1
                 if (irange == 1 .or. (wl < d(ibegin) - pivmin .and. wu >= d(ibegin) - pivmin)) &
                           then
                    m = m + 1
                    w(m) = d(ibegin)
                    iblock(m) = jb
                 end if
              else
                 ! general case -- in > 1
                 ! compute gershgorin interval
                 ! and use it as the initial interval
                 gu = d(ibegin)
                 gl = d(ibegin)
                 tmp1 = zero
                 do j = ibegin,iend - 1
                    tmp2 = abs(e(j))
                    gu = max(gu,d(j) + tmp1 + tmp2)
                    gl = min(gl,d(j) - tmp1 - tmp2)
                    tmp1 = tmp2
                 end do
                 gu = max(gu,d(iend) + tmp1)
                 gl = min(gl,d(iend) - tmp1)
                 bnorm = max(abs(gl),abs(gu))
                 gl = gl - fudge*bnorm*ulp*in - fudge*pivmin
                 gu = gu + fudge*bnorm*ulp*in + fudge*pivmin
                 ! compute atoli for the current submatrix
                 if (abstol <= zero) then
                    atoli = ulp*max(abs(gl),abs(gu))
                 else
                    atoli = abstol
                 end if
                 if (irange > 1) then
                    if (gu < wl) then
                       nwl = nwl + in
                       nwu = nwu + in
                       cycle loop_70
                    end if
                    gl = max(gl,wl)
                    gu = min(gu,wu)
                    if (gl >= gu) cycle loop_70
                 end if
                 ! set up initial interval
                 work(n + 1) = gl
                 work(n + in + 1) = gu
                 call la_slaebz(1,0,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                 ibegin),work(ibegin),idumma,work(n + 1),work(n + 2*in + 1),im,iwork,w(m + 1 &
                           ),iblock(m + 1),iinfo)
                 nwl = nwl + iwork(1)
                 nwu = nwu + iwork(in + 1)
                 iwoff = m - iwork(1)
                 ! compute eigenvalues
                 itmax = int((log(gu - gl + pivmin) - log(pivmin))/log(two),KIND=ilp) + &
                           2
                 call la_slaebz(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                  ibegin),work(ibegin),idumma,work(n + 1),work(n + 2*in + 1),iout,iwork,w( &
                            m + 1),iblock(m + 1),iinfo)
                 ! copy eigenvalues into w and iblock
                 ! use -jb for block number for unconverged eigenvalues.
                 do j = 1,iout
                    tmp1 = half*(work(j + n) + work(j + in + n))
                    ! flag non-convergence.
                    if (j > iout - iinfo) then
                       ncnvrg = .true.
                       ib = -jb
                    else
                       ib = jb
                    end if
                    do je = iwork(j) + 1 + iwoff,iwork(j + in) + iwoff
                       w(je) = tmp1
                       iblock(je) = ib
                    end do
                 end do
                 m = m + im
              end if
           end do loop_70
           ! if range='i', then (wl,wu) contains eigenvalues nwl+1,...,nwu
           ! if nwl+1 < il or nwu > iu, discard extra eigenvalues.
           if (irange == 3) then
              im = 0
              idiscl = il - 1 - nwl
              idiscu = nwu - iu
              if (idiscl > 0 .or. idiscu > 0) then
                 do je = 1,m
                    if (w(je) <= wlu .and. idiscl > 0) then
                       idiscl = idiscl - 1
                    else if (w(je) >= wul .and. idiscu > 0) then
                       idiscu = idiscu - 1
                    else
                       im = im + 1
                       w(im) = w(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl > 0 .or. idiscu > 0) then
                 ! code to deal with effects of bad arithmetic:
                 ! some low eigenvalues to be discarded are not in (wl,wlu],
                 ! or high eigenvalues to be discarded are not in (wul,wu]
                 ! so just kill off the smallest idiscl/largest idiscu
                 ! eigenvalues, by simply finding the smallest/largest
                 ! eigenvalue(s).
                 ! (if n(w) is monotone non-decreasing, this should never
                     ! happen.)
                 if (idiscl > 0) then
                    wkill = wu
                    do jdisc = 1,idiscl
                       iw = 0
                       do je = 1,m
                          if (iblock(je) /= 0 .and. (w(je) < wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 if (idiscu > 0) then
                    wkill = wl
                    do jdisc = 1,idiscu
                       iw = 0
                       do je = 1,m
                          if (iblock(je) /= 0 .and. (w(je) > wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 im = 0
                 do je = 1,m
                    if (iblock(je) /= 0) then
                       im = im + 1
                       w(im) = w(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl < 0 .or. idiscu < 0) then
                 toofew = .true.
              end if
           end if
           ! if order='b', do nothing -- the eigenvalues are already sorted
              ! by block.
           ! if order='e', sort the eigenvalues from smallest to largest
           if (iorder == 1 .and. nsplit > 1) then
              do je = 1,m - 1
                 ie = 0
                 tmp1 = w(je)
                 do j = je + 1,m
                    if (w(j) < tmp1) then
                       ie = j
                       tmp1 = w(j)
                    end if
                 end do
                 if (ie /= 0) then
                    itmp1 = iblock(ie)
                    w(ie) = w(je)
                    iblock(ie) = iblock(je)
                    w(je) = tmp1
                    iblock(je) = itmp1
                 end if
              end do
           end if
           info = 0
           if (ncnvrg) info = info + 1
           if (toofew) info = info + 2
           return
     end subroutine la_sstebz
     !> DSTEBZ: computes the eigenvalues of a symmetric tridiagonal
     !> matrix T.  The user may ask for all eigenvalues, all eigenvalues
     !> in the half-open interval (VL, VU], or the IL-th through IU-th
     !> eigenvalues.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_dstebz(range,order,n,vl,vu,il,iu,abstol,d,e,m,nsplit,w, &
               iblock,isplit,work,iwork,info)
        use la_constants_dp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: order,range
           integer(ilp),intent(in) :: il,iu,n
           integer(ilp),intent(out) :: info,m,nsplit
           real(dp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),isplit(*),iwork(*)
           real(dp),intent(in) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: fudge = 2.1_dp
           real(dp),parameter :: relfac = 2.0_dp

           ! Local Scalars
           logical(lk) :: ncnvrg,toofew
           integer(ilp) :: ib,ibegin,idiscl,idiscu,ie,iend,iinfo,im,in,ioff,iorder, &
                     iout,irange,itmax,itmp1,iw,iwoff,j,jb,jdisc,je,nb,nwl,nwu
           real(dp) :: atoli,bnorm,gl,gu,pivmin,rtoli,safemn,tmp1,tmp2,tnorm,ulp,wkill, &
                      wl,wlu,wu,wul
           ! Local Arrays
           integer(ilp) :: idumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           ! decode range
           if (la_lsame(range,'A')) then
              irange = 1
           else if (la_lsame(range,'V')) then
              irange = 2
           else if (la_lsame(range,'I')) then
              irange = 3
           else
              irange = 0
           end if
           ! decode order
           if (la_lsame(order,'B')) then
              iorder = 2
           else if (la_lsame(order,'E')) then
              iorder = 1
           else
              iorder = 0
           end if
           ! check for errors
           if (irange <= 0) then
              info = -1
           else if (iorder <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (irange == 2) then
              if (vl >= vu) info = -5
           else if (irange == 3 .and. (il < 1 .or. il > max(1,n))) then
              info = -6
           else if (irange == 3 .and. (iu < min(n,il) .or. iu > n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('DSTEBZ',-info)
              return
           end if
           ! initialize error flags
           info = 0
           ncnvrg = .false.
           toofew = .false.
           ! quick return if possible
           m = 0
           if (n == 0) return
           ! simplifications:
           if (irange == 3 .and. il == 1 .and. iu == n) irange = 1
           ! get machine constants
           ! nb is the minimum vector length for vector bisection, or 0
           ! if only scalar is to be done.
           safemn = la_dlamch('S')
           ulp = la_dlamch('P')
           rtoli = ulp*relfac
           nb = la_ilaenv(1,'DSTEBZ',' ',n,-1,-1,-1)
           if (nb <= 1) nb = 0
           ! special case when n=1
           if (n == 1) then
              nsplit = 1
              isplit(1) = 1
              if (irange == 2 .and. (vl >= d(1) .or. vu < d(1))) then
                 m = 0
              else
                 w(1) = d(1)
                 iblock(1) = 1
                 m = 1
              end if
              return
           end if
           ! compute splitting points
           nsplit = 1
           work(n) = zero
           pivmin = one
           do j = 2,n
              tmp1 = e(j - 1)**2
              if (abs(d(j)*d(j - 1))*ulp**2 + safemn > tmp1) then
                 isplit(nsplit) = j - 1
                 nsplit = nsplit + 1
                 work(j - 1) = zero
              else
                 work(j - 1) = tmp1
                 pivmin = max(pivmin,tmp1)
              end if
           end do
           isplit(nsplit) = n
           pivmin = pivmin*safemn
           ! compute interval and atoli
           if (irange == 3) then
              ! range='i': compute the interval containing eigenvalues
                         ! il through iu.
              ! compute gershgorin interval for entire (split) matrix
              ! and use it as the initial interval
              gu = d(1)
              gl = d(1)
              tmp1 = zero
              do j = 1,n - 1
                 tmp2 = sqrt(work(j))
                 gu = max(gu,d(j) + tmp1 + tmp2)
                 gl = min(gl,d(j) - tmp1 - tmp2)
                 tmp1 = tmp2
              end do
              gu = max(gu,d(n) + tmp1)
              gl = min(gl,d(n) - tmp1)
              tnorm = max(abs(gl),abs(gu))
              gl = gl - fudge*tnorm*ulp*n - fudge*two*pivmin
              gu = gu + fudge*tnorm*ulp*n + fudge*pivmin
              ! compute iteration parameters
              itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
              if (abstol <= zero) then
                 atoli = ulp*tnorm
              else
                 atoli = abstol
              end if
              work(n + 1) = gl
              work(n + 2) = gl
              work(n + 3) = gu
              work(n + 4) = gu
              work(n + 5) = gl
              work(n + 6) = gu
              iwork(1) = -1
              iwork(2) = -1
              iwork(3) = n + 1
              iwork(4) = n + 1
              iwork(5) = il - 1
              iwork(6) = iu
              call la_dlaebz(3,itmax,n,2,2,nb,atoli,rtoli,pivmin,d,e,work,iwork( &
                        5),work(n + 1),work(n + 5),iout,iwork,w,iblock,iinfo)
              if (iwork(6) == iu) then
                 wl = work(n + 1)
                 wlu = work(n + 3)
                 nwl = iwork(1)
                 wu = work(n + 4)
                 wul = work(n + 2)
                 nwu = iwork(4)
              else
                 wl = work(n + 2)
                 wlu = work(n + 4)
                 nwl = iwork(2)
                 wu = work(n + 3)
                 wul = work(n + 1)
                 nwu = iwork(3)
              end if
              if (nwl < 0 .or. nwl >= n .or. nwu < 1 .or. nwu > n) then
                 info = 4
                 return
              end if
           else
              ! range='a' or 'v' -- set atoli
              tnorm = max(abs(d(1)) + abs(e(1)),abs(d(n)) + abs(e(n - 1)))
              do j = 2,n - 1
                 tnorm = max(tnorm,abs(d(j)) + abs(e(j - 1)) + abs(e(j)))
              end do
              if (abstol <= zero) then
                 atoli = ulp*tnorm
              else
                 atoli = abstol
              end if
              if (irange == 2) then
                 wl = vl
                 wu = vu
              else
                 wl = zero
                 wu = zero
              end if
           end if
           ! find eigenvalues -- loop over blocks and recompute nwl and nwu.
           ! nwl accumulates the number of eigenvalues .le. wl,
           ! nwu accumulates the number of eigenvalues .le. wu
           m = 0
           iend = 0
           info = 0
           nwl = 0
           nwu = 0
           loop_70: do jb = 1,nsplit
              ioff = iend
              ibegin = ioff + 1
              iend = isplit(jb)
              in = iend - ioff
              if (in == 1) then
                 ! special case -- in=1
                 if (irange == 1 .or. wl >= d(ibegin) - pivmin) nwl = nwl + 1
                 if (irange == 1 .or. wu >= d(ibegin) - pivmin) nwu = nwu + 1
                 if (irange == 1 .or. (wl < d(ibegin) - pivmin .and. wu >= d(ibegin) - pivmin)) &
                           then
                    m = m + 1
                    w(m) = d(ibegin)
                    iblock(m) = jb
                 end if
              else
                 ! general case -- in > 1
                 ! compute gershgorin interval
                 ! and use it as the initial interval
                 gu = d(ibegin)
                 gl = d(ibegin)
                 tmp1 = zero
                 do j = ibegin,iend - 1
                    tmp2 = abs(e(j))
                    gu = max(gu,d(j) + tmp1 + tmp2)
                    gl = min(gl,d(j) - tmp1 - tmp2)
                    tmp1 = tmp2
                 end do
                 gu = max(gu,d(iend) + tmp1)
                 gl = min(gl,d(iend) - tmp1)
                 bnorm = max(abs(gl),abs(gu))
                 gl = gl - fudge*bnorm*ulp*in - fudge*pivmin
                 gu = gu + fudge*bnorm*ulp*in + fudge*pivmin
                 ! compute atoli for the current submatrix
                 if (abstol <= zero) then
                    atoli = ulp*max(abs(gl),abs(gu))
                 else
                    atoli = abstol
                 end if
                 if (irange > 1) then
                    if (gu < wl) then
                       nwl = nwl + in
                       nwu = nwu + in
                       cycle loop_70
                    end if
                    gl = max(gl,wl)
                    gu = min(gu,wu)
                    if (gl >= gu) cycle loop_70
                 end if
                 ! set up initial interval
                 work(n + 1) = gl
                 work(n + in + 1) = gu
                 call la_dlaebz(1,0,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                 ibegin),work(ibegin),idumma,work(n + 1),work(n + 2*in + 1),im,iwork,w(m + 1 &
                           ),iblock(m + 1),iinfo)
                 nwl = nwl + iwork(1)
                 nwu = nwu + iwork(in + 1)
                 iwoff = m - iwork(1)
                 ! compute eigenvalues
                 itmax = int((log(gu - gl + pivmin) - log(pivmin))/log(two),KIND=ilp) + &
                           2
                 call la_dlaebz(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                  ibegin),work(ibegin),idumma,work(n + 1),work(n + 2*in + 1),iout,iwork,w( &
                            m + 1),iblock(m + 1),iinfo)
                 ! copy eigenvalues into w and iblock
                 ! use -jb for block number for unconverged eigenvalues.
                 do j = 1,iout
                    tmp1 = half*(work(j + n) + work(j + in + n))
                    ! flag non-convergence.
                    if (j > iout - iinfo) then
                       ncnvrg = .true.
                       ib = -jb
                    else
                       ib = jb
                    end if
                    do je = iwork(j) + 1 + iwoff,iwork(j + in) + iwoff
                       w(je) = tmp1
                       iblock(je) = ib
                    end do
                 end do
                 m = m + im
              end if
           end do loop_70
           ! if range='i', then (wl,wu) contains eigenvalues nwl+1,...,nwu
           ! if nwl+1 < il or nwu > iu, discard extra eigenvalues.
           if (irange == 3) then
              im = 0
              idiscl = il - 1 - nwl
              idiscu = nwu - iu
              if (idiscl > 0 .or. idiscu > 0) then
                 do je = 1,m
                    if (w(je) <= wlu .and. idiscl > 0) then
                       idiscl = idiscl - 1
                    else if (w(je) >= wul .and. idiscu > 0) then
                       idiscu = idiscu - 1
                    else
                       im = im + 1
                       w(im) = w(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl > 0 .or. idiscu > 0) then
                 ! code to deal with effects of bad arithmetic:
                 ! some low eigenvalues to be discarded are not in (wl,wlu],
                 ! or high eigenvalues to be discarded are not in (wul,wu]
                 ! so just kill off the smallest idiscl/largest idiscu
                 ! eigenvalues, by simply finding the smallest/largest
                 ! eigenvalue(s).
                 ! (if n(w) is monotone non-decreasing, this should never
                     ! happen.)
                 if (idiscl > 0) then
                    wkill = wu
                    do jdisc = 1,idiscl
                       iw = 0
                       do je = 1,m
                          if (iblock(je) /= 0 .and. (w(je) < wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 if (idiscu > 0) then
                    wkill = wl
                    do jdisc = 1,idiscu
                       iw = 0
                       do je = 1,m
                          if (iblock(je) /= 0 .and. (w(je) > wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 im = 0
                 do je = 1,m
                    if (iblock(je) /= 0) then
                       im = im + 1
                       w(im) = w(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl < 0 .or. idiscu < 0) then
                 toofew = .true.
              end if
           end if
           ! if order='b', do nothing -- the eigenvalues are already sorted
              ! by block.
           ! if order='e', sort the eigenvalues from smallest to largest
           if (iorder == 1 .and. nsplit > 1) then
              do je = 1,m - 1
                 ie = 0
                 tmp1 = w(je)
                 do j = je + 1,m
                    if (w(j) < tmp1) then
                       ie = j
                       tmp1 = w(j)
                    end if
                 end do
                 if (ie /= 0) then
                    itmp1 = iblock(ie)
                    w(ie) = w(je)
                    iblock(ie) = iblock(je)
                    w(je) = tmp1
                    iblock(je) = itmp1
                 end if
              end do
           end if
           info = 0
           if (ncnvrg) info = info + 1
           if (toofew) info = info + 2
           return
     end subroutine la_dstebz
     !> QSTEBZ: computes the eigenvalues of a symmetric tridiagonal
     !> matrix T.  The user may ask for all eigenvalues, all eigenvalues
     !> in the half-open interval (VL, VU], or the IL-th through IU-th
     !> eigenvalues.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_qstebz(range,order,n,vl,vu,il,iu,abstol,d,e,m,nsplit,w, &
               iblock,isplit,work,iwork,info)
        use la_constants_qp,only:zero,half,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: order,range
           integer(ilp),intent(in) :: il,iu,n
           integer(ilp),intent(out) :: info,m,nsplit
           real(qp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),isplit(*),iwork(*)
           real(qp),intent(in) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: fudge = 2.1_qp
           real(qp),parameter :: relfac = 2.0_qp

           ! Local Scalars
           logical(lk) :: ncnvrg,toofew
           integer(ilp) :: ib,ibegin,idiscl,idiscu,ie,iend,iinfo,im,in,ioff,iorder, &
                     iout,irange,itmax,itmp1,iw,iwoff,j,jb,jdisc,je,nb,nwl,nwu
           real(qp) :: atoli,bnorm,gl,gu,pivmin,rtoli,safemn,tmp1,tmp2,tnorm,ulp,wkill, &
                      wl,wlu,wu,wul
           ! Local Arrays
           integer(ilp) :: idumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min,sqrt
           ! Executable Statements
           info = 0
           ! decode range
           if (la_lsame(range,'A')) then
              irange = 1
           else if (la_lsame(range,'V')) then
              irange = 2
           else if (la_lsame(range,'I')) then
              irange = 3
           else
              irange = 0
           end if
           ! decode order
           if (la_lsame(order,'B')) then
              iorder = 2
           else if (la_lsame(order,'E')) then
              iorder = 1
           else
              iorder = 0
           end if
           ! check for errors
           if (irange <= 0) then
              info = -1
           else if (iorder <= 0) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (irange == 2) then
              if (vl >= vu) info = -5
           else if (irange == 3 .and. (il < 1 .or. il > max(1,n))) then
              info = -6
           else if (irange == 3 .and. (iu < min(n,il) .or. iu > n)) then
              info = -7
           end if
           if (info /= 0) then
              call la_xerbla('QSTEBZ',-info)
              return
           end if
           ! initialize error flags
           info = 0
           ncnvrg = .false.
           toofew = .false.
           ! quick return if possible
           m = 0
           if (n == 0) return
           ! simplifications:
           if (irange == 3 .and. il == 1 .and. iu == n) irange = 1
           ! get machine constants
           ! nb is the minimum vector length for vector bisection, or 0
           ! if only scalar is to be done.
           safemn = la_qlamch('S')
           ulp = la_qlamch('P')
           rtoli = ulp*relfac
           nb = la_ilaenv(1,'QSTEBZ',' ',n,-1,-1,-1)
           if (nb <= 1) nb = 0
           ! special case when n=1
           if (n == 1) then
              nsplit = 1
              isplit(1) = 1
              if (irange == 2 .and. (vl >= d(1) .or. vu < d(1))) then
                 m = 0
              else
                 w(1) = d(1)
                 iblock(1) = 1
                 m = 1
              end if
              return
           end if
           ! compute splitting points
           nsplit = 1
           work(n) = zero
           pivmin = one
           do j = 2,n
              tmp1 = e(j - 1)**2
              if (abs(d(j)*d(j - 1))*ulp**2 + safemn > tmp1) then
                 isplit(nsplit) = j - 1
                 nsplit = nsplit + 1
                 work(j - 1) = zero
              else
                 work(j - 1) = tmp1
                 pivmin = max(pivmin,tmp1)
              end if
           end do
           isplit(nsplit) = n
           pivmin = pivmin*safemn
           ! compute interval and atoli
           if (irange == 3) then
              ! range='i': compute the interval containing eigenvalues
                         ! il through iu.
              ! compute gershgorin interval for entire (split) matrix
              ! and use it as the initial interval
              gu = d(1)
              gl = d(1)
              tmp1 = zero
              do j = 1,n - 1
                 tmp2 = sqrt(work(j))
                 gu = max(gu,d(j) + tmp1 + tmp2)
                 gl = min(gl,d(j) - tmp1 - tmp2)
                 tmp1 = tmp2
              end do
              gu = max(gu,d(n) + tmp1)
              gl = min(gl,d(n) - tmp1)
              tnorm = max(abs(gl),abs(gu))
              gl = gl - fudge*tnorm*ulp*n - fudge*two*pivmin
              gu = gu + fudge*tnorm*ulp*n + fudge*pivmin
              ! compute iteration parameters
              itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
              if (abstol <= zero) then
                 atoli = ulp*tnorm
              else
                 atoli = abstol
              end if
              work(n + 1) = gl
              work(n + 2) = gl
              work(n + 3) = gu
              work(n + 4) = gu
              work(n + 5) = gl
              work(n + 6) = gu
              iwork(1) = -1
              iwork(2) = -1
              iwork(3) = n + 1
              iwork(4) = n + 1
              iwork(5) = il - 1
              iwork(6) = iu
              call la_qlaebz(3,itmax,n,2,2,nb,atoli,rtoli,pivmin,d,e,work,iwork( &
                        5),work(n + 1),work(n + 5),iout,iwork,w,iblock,iinfo)
              if (iwork(6) == iu) then
                 wl = work(n + 1)
                 wlu = work(n + 3)
                 nwl = iwork(1)
                 wu = work(n + 4)
                 wul = work(n + 2)
                 nwu = iwork(4)
              else
                 wl = work(n + 2)
                 wlu = work(n + 4)
                 nwl = iwork(2)
                 wu = work(n + 3)
                 wul = work(n + 1)
                 nwu = iwork(3)
              end if
              if (nwl < 0 .or. nwl >= n .or. nwu < 1 .or. nwu > n) then
                 info = 4
                 return
              end if
           else
              ! range='a' or 'v' -- set atoli
              tnorm = max(abs(d(1)) + abs(e(1)),abs(d(n)) + abs(e(n - 1)))
              do j = 2,n - 1
                 tnorm = max(tnorm,abs(d(j)) + abs(e(j - 1)) + abs(e(j)))
              end do
              if (abstol <= zero) then
                 atoli = ulp*tnorm
              else
                 atoli = abstol
              end if
              if (irange == 2) then
                 wl = vl
                 wu = vu
              else
                 wl = zero
                 wu = zero
              end if
           end if
           ! find eigenvalues -- loop over blocks and recompute nwl and nwu.
           ! nwl accumulates the number of eigenvalues .le. wl,
           ! nwu accumulates the number of eigenvalues .le. wu
           m = 0
           iend = 0
           info = 0
           nwl = 0
           nwu = 0
           loop_70: do jb = 1,nsplit
              ioff = iend
              ibegin = ioff + 1
              iend = isplit(jb)
              in = iend - ioff
              if (in == 1) then
                 ! special case -- in=1
                 if (irange == 1 .or. wl >= d(ibegin) - pivmin) nwl = nwl + 1
                 if (irange == 1 .or. wu >= d(ibegin) - pivmin) nwu = nwu + 1
                 if (irange == 1 .or. (wl < d(ibegin) - pivmin .and. wu >= d(ibegin) - pivmin)) &
                           then
                    m = m + 1
                    w(m) = d(ibegin)
                    iblock(m) = jb
                 end if
              else
                 ! general case -- in > 1
                 ! compute gershgorin interval
                 ! and use it as the initial interval
                 gu = d(ibegin)
                 gl = d(ibegin)
                 tmp1 = zero
                 do j = ibegin,iend - 1
                    tmp2 = abs(e(j))
                    gu = max(gu,d(j) + tmp1 + tmp2)
                    gl = min(gl,d(j) - tmp1 - tmp2)
                    tmp1 = tmp2
                 end do
                 gu = max(gu,d(iend) + tmp1)
                 gl = min(gl,d(iend) - tmp1)
                 bnorm = max(abs(gl),abs(gu))
                 gl = gl - fudge*bnorm*ulp*in - fudge*pivmin
                 gu = gu + fudge*bnorm*ulp*in + fudge*pivmin
                 ! compute atoli for the current submatrix
                 if (abstol <= zero) then
                    atoli = ulp*max(abs(gl),abs(gu))
                 else
                    atoli = abstol
                 end if
                 if (irange > 1) then
                    if (gu < wl) then
                       nwl = nwl + in
                       nwu = nwu + in
                       cycle loop_70
                    end if
                    gl = max(gl,wl)
                    gu = min(gu,wu)
                    if (gl >= gu) cycle loop_70
                 end if
                 ! set up initial interval
                 work(n + 1) = gl
                 work(n + in + 1) = gu
                 call la_qlaebz(1,0,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                 ibegin),work(ibegin),idumma,work(n + 1),work(n + 2*in + 1),im,iwork,w(m + 1 &
                           ),iblock(m + 1),iinfo)
                 nwl = nwl + iwork(1)
                 nwu = nwu + iwork(in + 1)
                 iwoff = m - iwork(1)
                 ! compute eigenvalues
                 itmax = int((log(gu - gl + pivmin) - log(pivmin))/log(two),KIND=ilp) + &
                           2
                 call la_qlaebz(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                  ibegin),work(ibegin),idumma,work(n + 1),work(n + 2*in + 1),iout,iwork,w( &
                            m + 1),iblock(m + 1),iinfo)
                 ! copy eigenvalues into w and iblock
                 ! use -jb for block number for unconverged eigenvalues.
                 do j = 1,iout
                    tmp1 = half*(work(j + n) + work(j + in + n))
                    ! flag non-convergence.
                    if (j > iout - iinfo) then
                       ncnvrg = .true.
                       ib = -jb
                    else
                       ib = jb
                    end if
                    do je = iwork(j) + 1 + iwoff,iwork(j + in) + iwoff
                       w(je) = tmp1
                       iblock(je) = ib
                    end do
                 end do
                 m = m + im
              end if
           end do loop_70
           ! if range='i', then (wl,wu) contains eigenvalues nwl+1,...,nwu
           ! if nwl+1 < il or nwu > iu, discard extra eigenvalues.
           if (irange == 3) then
              im = 0
              idiscl = il - 1 - nwl
              idiscu = nwu - iu
              if (idiscl > 0 .or. idiscu > 0) then
                 do je = 1,m
                    if (w(je) <= wlu .and. idiscl > 0) then
                       idiscl = idiscl - 1
                    else if (w(je) >= wul .and. idiscu > 0) then
                       idiscu = idiscu - 1
                    else
                       im = im + 1
                       w(im) = w(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl > 0 .or. idiscu > 0) then
                 ! code to deal with effects of bad arithmetic:
                 ! some low eigenvalues to be discarded are not in (wl,wlu],
                 ! or high eigenvalues to be discarded are not in (wul,wu]
                 ! so just kill off the smallest idiscl/largest idiscu
                 ! eigenvalues, by simply finding the smallest/largest
                 ! eigenvalue(s).
                 ! (if n(w) is monotone non-decreasing, this should never
                     ! happen.)
                 if (idiscl > 0) then
                    wkill = wu
                    do jdisc = 1,idiscl
                       iw = 0
                       do je = 1,m
                          if (iblock(je) /= 0 .and. (w(je) < wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 if (idiscu > 0) then
                    wkill = wl
                    do jdisc = 1,idiscu
                       iw = 0
                       do je = 1,m
                          if (iblock(je) /= 0 .and. (w(je) > wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 im = 0
                 do je = 1,m
                    if (iblock(je) /= 0) then
                       im = im + 1
                       w(im) = w(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl < 0 .or. idiscu < 0) then
                 toofew = .true.
              end if
           end if
           ! if order='b', do nothing -- the eigenvalues are already sorted
              ! by block.
           ! if order='e', sort the eigenvalues from smallest to largest
           if (iorder == 1 .and. nsplit > 1) then
              do je = 1,m - 1
                 ie = 0
                 tmp1 = w(je)
                 do j = je + 1,m
                    if (w(j) < tmp1) then
                       ie = j
                       tmp1 = w(j)
                    end if
                 end do
                 if (ie /= 0) then
                    itmp1 = iblock(ie)
                    w(ie) = w(je)
                    iblock(ie) = iblock(je)
                    w(je) = tmp1
                    iblock(je) = itmp1
                 end if
              end do
           end if
           info = 0
           if (ncnvrg) info = info + 1
           if (toofew) info = info + 2
           return
     end subroutine la_qstebz

     !> SSTEIN: computes the eigenvectors of a real symmetric tridiagonal
     !> matrix T corresponding to specified eigenvalues, using inverse
     !> iteration.
     !> The maximum number of iterations allowed for each eigenvector is
     !> specified by an internal parameter MAXITS (currently set to 5).

     pure subroutine la_sstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail, &
               info)
        use la_constants_sp,only:zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,m,n
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),isplit(*)
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(sp),intent(in) :: d(*),e(*),w(*)
           real(sp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: odm3 = 1.0e-3_sp
           real(sp),parameter :: odm1 = 1.0e-1_sp
           integer(ilp),parameter :: maxits = 5
           integer(ilp),parameter :: extra = 2

           ! Local Scalars
           integer(ilp) :: b1,blksiz,bn,gpind,i,iinfo,indrv1,indrv2,indrv3,indrv4, &
                     indrv5,its,j,j1,jblk,jmax,nblk,nrmchk
           real(sp) :: ctr,eps,eps1,nrm,onenrm,ortol,pertol,scl,sep,stpcrt,tol,xj, &
                     xjm
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           do i = 1,m
              ifail(i) = 0
           end do
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -4
           else if (ldz < max(1,n)) then
              info = -9
           else
              do j = 2,m
                 if (iblock(j) < iblock(j - 1)) then
                    info = -6
                    go to 30
                 end if
                 if (iblock(j) == iblock(j - 1) .and. w(j) < w(j - 1)) then
                    info = -5
                    go to 30
                 end if
              end do
              30 continue
           end if
           if (info /= 0) then
              call la_xerbla('SSTEIN',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) then
              return
           else if (n == 1) then
              z(1,1) = one
              return
           end if
           ! get machine constants.
           eps = la_slamch('PRECISION')
           ! initialize seed for random number generator la_slarnv.
           do i = 1,4
              iseed(i) = 1
           end do
           ! initialize pointers.
           indrv1 = 0
           indrv2 = indrv1 + n
           indrv3 = indrv2 + n
           indrv4 = indrv3 + n
           indrv5 = indrv4 + n
           ! compute eigenvectors of matrix blocks.
           j1 = 1
           loop_160: do nblk = 1,iblock(m)
              ! find starting and ending indices of block nblk.
              if (nblk == 1) then
                 b1 = 1
              else
                 b1 = isplit(nblk - 1) + 1
              end if
              bn = isplit(nblk)
              blksiz = bn - b1 + 1
              if (blksiz == 1) go to 60
              gpind = j1
              ! compute reorthogonalization criterion and stopping criterion.
              onenrm = abs(d(b1)) + abs(e(b1))
              onenrm = max(onenrm,abs(d(bn)) + abs(e(bn - 1)))
              do i = b1 + 1,bn - 1
                 onenrm = max(onenrm,abs(d(i)) + abs(e(i - 1)) + abs(e(i)))
              end do
              ortol = odm3*onenrm
              stpcrt = sqrt(odm1/blksiz)
              ! loop through eigenvalues of block nblk.
              60 continue
              jblk = 0
              loop_150: do j = j1,m
                 if (iblock(j) /= nblk) then
                    j1 = j
                    cycle loop_160
                 end if
                 jblk = jblk + 1
                 xj = w(j)
                 ! skip all the work if the block size is one.
                 if (blksiz == 1) then
                    work(indrv1 + 1) = one
                    go to 120
                 end if
                 ! if eigenvalues j and j-1 are too close, add a relatively
                 ! small perturbation.
                 if (jblk > 1) then
                    eps1 = abs(eps*xj)
                    pertol = ten*eps1
                    sep = xj - xjm
                    if (sep < pertol) xj = xjm + pertol
                 end if
                 its = 0
                 nrmchk = 0
                 ! get random starting vector.
                 call la_slarnv(2,iseed,blksiz,work(indrv1 + 1))
                 ! copy the matrix t so it won't be destroyed in factorization.
                 call la_scopy(blksiz,d(b1),1,work(indrv4 + 1),1)
                 call la_scopy(blksiz - 1,e(b1),1,work(indrv2 + 2),1)
                 call la_scopy(blksiz - 1,e(b1),1,work(indrv3 + 1),1)
                 ! compute lu factors with partial pivoting  ( pt = lu )
                 tol = zero
                 call la_slagtf(blksiz,work(indrv4 + 1),xj,work(indrv2 + 2),work(indrv3 + &
                           1),tol,work(indrv5 + 1),iwork,iinfo)
                 ! update iteration count.
                 70 continue
                 its = its + 1
                 if (its > maxits) go to 100
                 ! normalize and scale the righthand side vector pb.
                 jmax = la_isamax(blksiz,work(indrv1 + 1),1)
                 scl = blksiz*onenrm*max(eps,abs(work(indrv4 + blksiz)))/abs(work(indrv1 + &
                           jmax))
                 call la_sscal(blksiz,scl,work(indrv1 + 1),1)
                 ! solve the system lu = pb.
                 call la_slagts(-1,blksiz,work(indrv4 + 1),work(indrv2 + 2),work(indrv3 + &
                           1),work(indrv5 + 1),iwork,work(indrv1 + 1),tol,iinfo)
                 ! reorthogonalize by modified gram-schmidt if eigenvalues are
                 ! close enough.
                 if (jblk == 1) go to 90
                 if (abs(xj - xjm) > ortol) gpind = j
                 if (gpind /= j) then
                    do i = gpind,j - 1
                       ctr = -la_sdot(blksiz,work(indrv1 + 1),1,z(b1,i),1)
                       call la_saxpy(blksiz,ctr,z(b1,i),1,work(indrv1 + 1),1)
                    end do
                 end if
                 ! check the infinity norm of the iterate.
                 90 continue
                 jmax = la_isamax(blksiz,work(indrv1 + 1),1)
                 nrm = abs(work(indrv1 + jmax))
                 ! continue for additional iterations after norm reaches
                 ! stopping criterion.
                 if (nrm < stpcrt) go to 70
                 nrmchk = nrmchk + 1
                 if (nrmchk < extra + 1) go to 70
                 go to 110
                 ! if stopping criterion was not satisfied, update info and
                 ! store eigenvector number in array ifail.
                 100 continue
                 info = info + 1
                 ifail(info) = j
                 ! accept iterate as jth eigenvector.
                 110 continue
                 scl = one/la_snrm2(blksiz,work(indrv1 + 1),1)
                 jmax = la_isamax(blksiz,work(indrv1 + 1),1)
                 if (work(indrv1 + jmax) < zero) scl = -scl
                 call la_sscal(blksiz,scl,work(indrv1 + 1),1)
                 120 continue
                 do i = 1,n
                    z(i,j) = zero
                 end do
                 do i = 1,blksiz
                    z(b1 + i - 1,j) = work(indrv1 + i)
                 end do
                 ! save the shift to check eigenvalue spacing at next
                 ! iteration.
                 xjm = xj
              end do loop_150
           end do loop_160
           return
     end subroutine la_sstein
     !> DSTEIN: computes the eigenvectors of a real symmetric tridiagonal
     !> matrix T corresponding to specified eigenvalues, using inverse
     !> iteration.
     !> The maximum number of iterations allowed for each eigenvector is
     !> specified by an internal parameter MAXITS (currently set to 5).

     pure subroutine la_dstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail, &
               info)
        use la_constants_dp,only:zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,m,n
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),isplit(*)
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(dp),intent(in) :: d(*),e(*),w(*)
           real(dp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: odm3 = 1.0e-3_dp
           real(dp),parameter :: odm1 = 1.0e-1_dp
           integer(ilp),parameter :: maxits = 5
           integer(ilp),parameter :: extra = 2

           ! Local Scalars
           integer(ilp) :: b1,blksiz,bn,gpind,i,iinfo,indrv1,indrv2,indrv3,indrv4, &
                     indrv5,its,j,j1,jblk,jmax,nblk,nrmchk
           real(dp) :: dtpcrt,eps,eps1,nrm,onenrm,ortol,pertol,scl,sep,tol,xj,xjm, &
                     ztr
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           do i = 1,m
              ifail(i) = 0
           end do
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -4
           else if (ldz < max(1,n)) then
              info = -9
           else
              do j = 2,m
                 if (iblock(j) < iblock(j - 1)) then
                    info = -6
                    go to 30
                 end if
                 if (iblock(j) == iblock(j - 1) .and. w(j) < w(j - 1)) then
                    info = -5
                    go to 30
                 end if
              end do
              30 continue
           end if
           if (info /= 0) then
              call la_xerbla('DSTEIN',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) then
              return
           else if (n == 1) then
              z(1,1) = one
              return
           end if
           ! get machine constants.
           eps = la_dlamch('PRECISION')
           ! initialize seed for random number generator la_dlarnv.
           do i = 1,4
              iseed(i) = 1
           end do
           ! initialize pointers.
           indrv1 = 0
           indrv2 = indrv1 + n
           indrv3 = indrv2 + n
           indrv4 = indrv3 + n
           indrv5 = indrv4 + n
           ! compute eigenvectors of matrix blocks.
           j1 = 1
           loop_160: do nblk = 1,iblock(m)
              ! find starting and ending indices of block nblk.
              if (nblk == 1) then
                 b1 = 1
              else
                 b1 = isplit(nblk - 1) + 1
              end if
              bn = isplit(nblk)
              blksiz = bn - b1 + 1
              if (blksiz == 1) go to 60
              gpind = j1
              ! compute reorthogonalization criterion and stopping criterion.
              onenrm = abs(d(b1)) + abs(e(b1))
              onenrm = max(onenrm,abs(d(bn)) + abs(e(bn - 1)))
              do i = b1 + 1,bn - 1
                 onenrm = max(onenrm,abs(d(i)) + abs(e(i - 1)) + abs(e(i)))
              end do
              ortol = odm3*onenrm
              dtpcrt = sqrt(odm1/blksiz)
              ! loop through eigenvalues of block nblk.
              60 continue
              jblk = 0
              loop_150: do j = j1,m
                 if (iblock(j) /= nblk) then
                    j1 = j
                    cycle loop_160
                 end if
                 jblk = jblk + 1
                 xj = w(j)
                 ! skip all the work if the block size is one.
                 if (blksiz == 1) then
                    work(indrv1 + 1) = one
                    go to 120
                 end if
                 ! if eigenvalues j and j-1 are too close, add a relatively
                 ! small perturbation.
                 if (jblk > 1) then
                    eps1 = abs(eps*xj)
                    pertol = ten*eps1
                    sep = xj - xjm
                    if (sep < pertol) xj = xjm + pertol
                 end if
                 its = 0
                 nrmchk = 0
                 ! get random starting vector.
                 call la_dlarnv(2,iseed,blksiz,work(indrv1 + 1))
                 ! copy the matrix t so it won't be destroyed in factorization.
                 call la_dcopy(blksiz,d(b1),1,work(indrv4 + 1),1)
                 call la_dcopy(blksiz - 1,e(b1),1,work(indrv2 + 2),1)
                 call la_dcopy(blksiz - 1,e(b1),1,work(indrv3 + 1),1)
                 ! compute lu factors with partial pivoting  ( pt = lu )
                 tol = zero
                 call la_dlagtf(blksiz,work(indrv4 + 1),xj,work(indrv2 + 2),work(indrv3 + &
                           1),tol,work(indrv5 + 1),iwork,iinfo)
                 ! update iteration count.
                 70 continue
                 its = its + 1
                 if (its > maxits) go to 100
                 ! normalize and scale the righthand side vector pb.
                 jmax = la_idamax(blksiz,work(indrv1 + 1),1)
                 scl = blksiz*onenrm*max(eps,abs(work(indrv4 + blksiz)))/abs(work(indrv1 + &
                           jmax))
                 call la_dscal(blksiz,scl,work(indrv1 + 1),1)
                 ! solve the system lu = pb.
                 call la_dlagts(-1,blksiz,work(indrv4 + 1),work(indrv2 + 2),work(indrv3 + &
                           1),work(indrv5 + 1),iwork,work(indrv1 + 1),tol,iinfo)
                 ! reorthogonalize by modified gram-schmidt if eigenvalues are
                 ! close enough.
                 if (jblk == 1) go to 90
                 if (abs(xj - xjm) > ortol) gpind = j
                 if (gpind /= j) then
                    do i = gpind,j - 1
                       ztr = -la_ddot(blksiz,work(indrv1 + 1),1,z(b1,i),1)
                       call la_daxpy(blksiz,ztr,z(b1,i),1,work(indrv1 + 1),1)
                    end do
                 end if
                 ! check the infinity norm of the iterate.
                 90 continue
                 jmax = la_idamax(blksiz,work(indrv1 + 1),1)
                 nrm = abs(work(indrv1 + jmax))
                 ! continue for additional iterations after norm reaches
                 ! stopping criterion.
                 if (nrm < dtpcrt) go to 70
                 nrmchk = nrmchk + 1
                 if (nrmchk < extra + 1) go to 70
                 go to 110
                 ! if stopping criterion was not satisfied, update info and
                 ! store eigenvector number in array ifail.
                 100 continue
                 info = info + 1
                 ifail(info) = j
                 ! accept iterate as jth eigenvector.
                 110 continue
                 scl = one/la_dnrm2(blksiz,work(indrv1 + 1),1)
                 jmax = la_idamax(blksiz,work(indrv1 + 1),1)
                 if (work(indrv1 + jmax) < zero) scl = -scl
                 call la_dscal(blksiz,scl,work(indrv1 + 1),1)
                 120 continue
                 do i = 1,n
                    z(i,j) = zero
                 end do
                 do i = 1,blksiz
                    z(b1 + i - 1,j) = work(indrv1 + i)
                 end do
                 ! save the shift to check eigenvalue spacing at next
                 ! iteration.
                 xjm = xj
              end do loop_150
           end do loop_160
           return
     end subroutine la_dstein
     !> QSTEIN: computes the eigenvectors of a real symmetric tridiagonal
     !> matrix T corresponding to specified eigenvalues, using inverse
     !> iteration.
     !> The maximum number of iterations allowed for each eigenvector is
     !> specified by an internal parameter MAXITS (currently set to 5).

     pure subroutine la_qstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail, &
               info)
        use la_constants_qp,only:zero,one,ten
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,m,n
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),isplit(*)
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(qp),intent(in) :: d(*),e(*),w(*)
           real(qp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: odm3 = 1.0e-3_qp
           real(qp),parameter :: odm1 = 1.0e-1_qp
           integer(ilp),parameter :: maxits = 5
           integer(ilp),parameter :: extra = 2

           ! Local Scalars
           integer(ilp) :: b1,blksiz,bn,gpind,i,iinfo,indrv1,indrv2,indrv3,indrv4, &
                     indrv5,its,j,j1,jblk,jmax,nblk,nrmchk
           real(qp) :: qtpcrt,eps,eps1,nrm,onenrm,ortol,pertol,scl,sep,tol,xj,xjm, &
                     wtr
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           do i = 1,m
              ifail(i) = 0
           end do
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -4
           else if (ldz < max(1,n)) then
              info = -9
           else
              do j = 2,m
                 if (iblock(j) < iblock(j - 1)) then
                    info = -6
                    go to 30
                 end if
                 if (iblock(j) == iblock(j - 1) .and. w(j) < w(j - 1)) then
                    info = -5
                    go to 30
                 end if
              end do
              30 continue
           end if
           if (info /= 0) then
              call la_xerbla('QSTEIN',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) then
              return
           else if (n == 1) then
              z(1,1) = one
              return
           end if
           ! get machine constants.
           eps = la_qlamch('PRECISION')
           ! initialize seed for random number generator la_qlarnv.
           do i = 1,4
              iseed(i) = 1
           end do
           ! initialize pointers.
           indrv1 = 0
           indrv2 = indrv1 + n
           indrv3 = indrv2 + n
           indrv4 = indrv3 + n
           indrv5 = indrv4 + n
           ! compute eigenvectors of matrix blocks.
           j1 = 1
           loop_160: do nblk = 1,iblock(m)
              ! find starting and ending indices of block nblk.
              if (nblk == 1) then
                 b1 = 1
              else
                 b1 = isplit(nblk - 1) + 1
              end if
              bn = isplit(nblk)
              blksiz = bn - b1 + 1
              if (blksiz == 1) go to 60
              gpind = j1
              ! compute reorthogonalization criterion and stopping criterion.
              onenrm = abs(d(b1)) + abs(e(b1))
              onenrm = max(onenrm,abs(d(bn)) + abs(e(bn - 1)))
              do i = b1 + 1,bn - 1
                 onenrm = max(onenrm,abs(d(i)) + abs(e(i - 1)) + abs(e(i)))
              end do
              ortol = odm3*onenrm
              qtpcrt = sqrt(odm1/blksiz)
              ! loop through eigenvalues of block nblk.
              60 continue
              jblk = 0
              loop_150: do j = j1,m
                 if (iblock(j) /= nblk) then
                    j1 = j
                    cycle loop_160
                 end if
                 jblk = jblk + 1
                 xj = w(j)
                 ! skip all the work if the block size is one.
                 if (blksiz == 1) then
                    work(indrv1 + 1) = one
                    go to 120
                 end if
                 ! if eigenvalues j and j-1 are too close, add a relatively
                 ! small perturbation.
                 if (jblk > 1) then
                    eps1 = abs(eps*xj)
                    pertol = ten*eps1
                    sep = xj - xjm
                    if (sep < pertol) xj = xjm + pertol
                 end if
                 its = 0
                 nrmchk = 0
                 ! get random starting vector.
                 call la_qlarnv(2,iseed,blksiz,work(indrv1 + 1))
                 ! copy the matrix t so it won't be destroyed in factorization.
                 call la_qcopy(blksiz,d(b1),1,work(indrv4 + 1),1)
                 call la_qcopy(blksiz - 1,e(b1),1,work(indrv2 + 2),1)
                 call la_qcopy(blksiz - 1,e(b1),1,work(indrv3 + 1),1)
                 ! compute lu factors with partial pivoting  ( pt = lu )
                 tol = zero
                 call la_qlagtf(blksiz,work(indrv4 + 1),xj,work(indrv2 + 2),work(indrv3 + &
                           1),tol,work(indrv5 + 1),iwork,iinfo)
                 ! update iteration count.
                 70 continue
                 its = its + 1
                 if (its > maxits) go to 100
                 ! normalize and scale the righthand side vector pb.
                 jmax = la_iqamax(blksiz,work(indrv1 + 1),1)
                 scl = blksiz*onenrm*max(eps,abs(work(indrv4 + blksiz)))/abs(work(indrv1 + &
                           jmax))
                 call la_qscal(blksiz,scl,work(indrv1 + 1),1)
                 ! solve the system lu = pb.
                 call la_qlagts(-1,blksiz,work(indrv4 + 1),work(indrv2 + 2),work(indrv3 + &
                           1),work(indrv5 + 1),iwork,work(indrv1 + 1),tol,iinfo)
                 ! reorthogonalize by modified gram-schmidt if eigenvalues are
                 ! close enough.
                 if (jblk == 1) go to 90
                 if (abs(xj - xjm) > ortol) gpind = j
                 if (gpind /= j) then
                    do i = gpind,j - 1
                       wtr = -la_qdot(blksiz,work(indrv1 + 1),1,z(b1,i),1)
                       call la_qaxpy(blksiz,wtr,z(b1,i),1,work(indrv1 + 1),1)
                    end do
                 end if
                 ! check the infinity norm of the iterate.
                 90 continue
                 jmax = la_iqamax(blksiz,work(indrv1 + 1),1)
                 nrm = abs(work(indrv1 + jmax))
                 ! continue for additional iterations after norm reaches
                 ! stopping criterion.
                 if (nrm < qtpcrt) go to 70
                 nrmchk = nrmchk + 1
                 if (nrmchk < extra + 1) go to 70
                 go to 110
                 ! if stopping criterion was not satisfied, update info and
                 ! store eigenvector number in array ifail.
                 100 continue
                 info = info + 1
                 ifail(info) = j
                 ! accept iterate as jth eigenvector.
                 110 continue
                 scl = one/la_qnrm2(blksiz,work(indrv1 + 1),1)
                 jmax = la_iqamax(blksiz,work(indrv1 + 1),1)
                 if (work(indrv1 + jmax) < zero) scl = -scl
                 call la_qscal(blksiz,scl,work(indrv1 + 1),1)
                 120 continue
                 do i = 1,n
                    z(i,j) = zero
                 end do
                 do i = 1,blksiz
                    z(b1 + i - 1,j) = work(indrv1 + i)
                 end do
                 ! save the shift to check eigenvalue spacing at next
                 ! iteration.
                 xjm = xj
              end do loop_150
           end do loop_160
           return
     end subroutine la_qstein

     !> SSTERF: computes all eigenvalues of a symmetric tridiagonal matrix
     !> using the Pal-Walker-Kahan variant of the QL or QR algorithm.

     pure subroutine la_ssterf(n,d,e,info)
        use la_constants_sp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: d(*),e(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,iscale,jtot,l,l1,lend,lendsv,lsv,m,nmaxit
           real(sp) :: alpha,anorm,bb,c,eps,eps2,gamma,oldc,oldgam,p,r,rt1,rt2,rte, &
                     s,safmax,safmin,sigma,ssfmax,ssfmin
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! quick return if possible
           if (n < 0) then
              info = -1
              call la_xerbla('SSTERF',-info)
              return
           end if
           if (n <= 1) return
           ! determine the unit roundoff for this environment.
           eps = la_slamch('E')
           eps2 = eps**2
           safmin = la_slamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           ! compute the eigenvalues of the tridiagonal matrix.
           nmaxit = n*maxit
           sigma = zero
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           10 continue
           if (l1 > n) go to 170
           if (l1 > 1) e(l1 - 1) = zero
           do m = l1,n - 1
              if (abs(e(m)) <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) &
                        then
                 e(m) = zero
                 go to 30
              end if
           end do
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
           do i = l,lend - 1
              e(i) = e(i)**2
           end do
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend >= l) then
              ! ql iteration
              ! look for small subdiagonal element.
              50 continue
              if (l /= lend) then
                 do m = l,lend - 1
                    if (abs(e(m)) <= eps2*abs(d(m)*d(m + 1))) go to 70
                 end do
              end if
              m = lend
              70 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 90
              ! if remaining matrix is 2 by 2, use la_slae2 to compute its
              ! eigenvalues.
              if (m == l + 1) then
                 rte = sqrt(e(l))
                 call la_slae2(d(l),rte,d(l + 1),rt1,rt2)
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 50
                 go to 150
              end if
              if (jtot == nmaxit) go to 150
              jtot = jtot + 1
              ! form shift.
              rte = sqrt(e(l))
              sigma = (d(l + 1) - p)/(two*rte)
              r = la_slapy2(sigma,one)
              sigma = p - (rte/(sigma + sign(r,sigma)))
              c = one
              s = zero
              gamma = d(m) - sigma
              p = gamma*gamma
              ! inner loop
              do i = m - 1,l,-1
                 bb = e(i)
                 r = p + bb
                 if (i /= m - 1) e(i + 1) = s*r
                 oldc = c
                 c = p/r
                 s = bb/r
                 oldgam = gamma
                 alpha = d(i)
                 gamma = c*(alpha - sigma) - s*oldgam
                 d(i + 1) = oldgam + (alpha - gamma)
                 if (c /= zero) then
                    p = (gamma*gamma)/c
                 else
                    p = oldc*bb
                 end if
              end do
              e(l) = s*p
              d(l) = sigma + gamma
              go to 50
              ! eigenvalue found.
              90 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 50
              go to 150
           else
              ! qr iteration
              ! look for small superdiagonal element.
              100 continue
              do m = l,lend + 1,-1
                 if (abs(e(m - 1)) <= eps2*abs(d(m)*d(m - 1))) go to 120
              end do
              m = lend
              120 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 140
              ! if remaining matrix is 2 by 2, use la_slae2 to compute its
              ! eigenvalues.
              if (m == l - 1) then
                 rte = sqrt(e(l - 1))
                 call la_slae2(d(l),rte,d(l - 1),rt1,rt2)
                 d(l) = rt1
                 d(l - 1) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 100
                 go to 150
              end if
              if (jtot == nmaxit) go to 150
              jtot = jtot + 1
              ! form shift.
              rte = sqrt(e(l - 1))
              sigma = (d(l - 1) - p)/(two*rte)
              r = la_slapy2(sigma,one)
              sigma = p - (rte/(sigma + sign(r,sigma)))
              c = one
              s = zero
              gamma = d(m) - sigma
              p = gamma*gamma
              ! inner loop
              do i = m,l - 1
                 bb = e(i)
                 r = p + bb
                 if (i /= m) e(i - 1) = s*r
                 oldc = c
                 c = p/r
                 s = bb/r
                 oldgam = gamma
                 alpha = d(i + 1)
                 gamma = c*(alpha - sigma) - s*oldgam
                 d(i) = oldgam + (alpha - gamma)
                 if (c /= zero) then
                    p = (gamma*gamma)/c
                 else
                    p = oldc*bb
                 end if
              end do
              e(l - 1) = s*p
              d(l) = sigma + gamma
              go to 100
              ! eigenvalue found.
              140 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 100
              go to 150
           end if
           ! undo scaling if necessary
           150 continue
           if (iscale == 1) call la_slascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv), &
                     n,info)
           if (iscale == 2) call la_slascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv), &
                     n,info)
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 180
           ! sort eigenvalues in increasing order.
           170 continue
           call la_slasrt('I',n,d,info)
           180 continue
           return
     end subroutine la_ssterf
     !> DSTERF: computes all eigenvalues of a symmetric tridiagonal matrix
     !> using the Pal-Walker-Kahan variant of the QL or QR algorithm.

     pure subroutine la_dsterf(n,d,e,info)
        use la_constants_dp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: d(*),e(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,iscale,jtot,l,l1,lend,lendsv,lsv,m,nmaxit
           real(dp) :: alpha,anorm,bb,c,eps,eps2,gamma,oldc,oldgam,p,r,rt1,rt2,rte, &
                     s,safmax,safmin,sigma,ssfmax,ssfmin,rmax
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! quick return if possible
           if (n < 0) then
              info = -1
              call la_xerbla('DSTERF',-info)
              return
           end if
           if (n <= 1) return
           ! determine the unit roundoff for this environment.
           eps = la_dlamch('E')
           eps2 = eps**2
           safmin = la_dlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           rmax = la_dlamch('O')
           ! compute the eigenvalues of the tridiagonal matrix.
           nmaxit = n*maxit
           sigma = zero
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           10 continue
           if (l1 > n) go to 170
           if (l1 > 1) e(l1 - 1) = zero
           do m = l1,n - 1
              if (abs(e(m)) <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) &
                        then
                 e(m) = zero
                 go to 30
              end if
           end do
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
           if ((anorm > ssfmax)) then
              iscale = 1
              call la_dlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_dlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_dlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_dlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           do i = l,lend - 1
              e(i) = e(i)**2
           end do
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend >= l) then
              ! ql iteration
              ! look for small subdiagonal element.
              50 continue
              if (l /= lend) then
                 do m = l,lend - 1
                    if (abs(e(m)) <= eps2*abs(d(m)*d(m + 1))) go to 70
                 end do
              end if
              m = lend
              70 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 90
              ! if remaining matrix is 2 by 2, use la_dlae2 to compute its
              ! eigenvalues.
              if (m == l + 1) then
                 rte = sqrt(e(l))
                 call la_dlae2(d(l),rte,d(l + 1),rt1,rt2)
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 50
                 go to 150
              end if
              if (jtot == nmaxit) go to 150
              jtot = jtot + 1
              ! form shift.
              rte = sqrt(e(l))
              sigma = (d(l + 1) - p)/(two*rte)
              r = la_dlapy2(sigma,one)
              sigma = p - (rte/(sigma + sign(r,sigma)))
              c = one
              s = zero
              gamma = d(m) - sigma
              p = gamma*gamma
              ! inner loop
              do i = m - 1,l,-1
                 bb = e(i)
                 r = p + bb
                 if (i /= m - 1) e(i + 1) = s*r
                 oldc = c
                 c = p/r
                 s = bb/r
                 oldgam = gamma
                 alpha = d(i)
                 gamma = c*(alpha - sigma) - s*oldgam
                 d(i + 1) = oldgam + (alpha - gamma)
                 if (c /= zero) then
                    p = (gamma*gamma)/c
                 else
                    p = oldc*bb
                 end if
              end do
              e(l) = s*p
              d(l) = sigma + gamma
              go to 50
              ! eigenvalue found.
              90 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 50
              go to 150
           else
              ! qr iteration
              ! look for small superdiagonal element.
              100 continue
              do m = l,lend + 1,-1
                 if (abs(e(m - 1)) <= eps2*abs(d(m)*d(m - 1))) go to 120
              end do
              m = lend
              120 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 140
              ! if remaining matrix is 2 by 2, use la_dlae2 to compute its
              ! eigenvalues.
              if (m == l - 1) then
                 rte = sqrt(e(l - 1))
                 call la_dlae2(d(l),rte,d(l - 1),rt1,rt2)
                 d(l) = rt1
                 d(l - 1) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 100
                 go to 150
              end if
              if (jtot == nmaxit) go to 150
              jtot = jtot + 1
              ! form shift.
              rte = sqrt(e(l - 1))
              sigma = (d(l - 1) - p)/(two*rte)
              r = la_dlapy2(sigma,one)
              sigma = p - (rte/(sigma + sign(r,sigma)))
              c = one
              s = zero
              gamma = d(m) - sigma
              p = gamma*gamma
              ! inner loop
              do i = m,l - 1
                 bb = e(i)
                 r = p + bb
                 if (i /= m) e(i - 1) = s*r
                 oldc = c
                 c = p/r
                 s = bb/r
                 oldgam = gamma
                 alpha = d(i + 1)
                 gamma = c*(alpha - sigma) - s*oldgam
                 d(i) = oldgam + (alpha - gamma)
                 if (c /= zero) then
                    p = (gamma*gamma)/c
                 else
                    p = oldc*bb
                 end if
              end do
              e(l - 1) = s*p
              d(l) = sigma + gamma
              go to 100
              ! eigenvalue found.
              140 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 100
              go to 150
           end if
           ! undo scaling if necessary
           150 continue
           if (iscale == 1) call la_dlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv), &
                     n,info)
           if (iscale == 2) call la_dlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv), &
                     n,info)
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 180
           ! sort eigenvalues in increasing order.
           170 continue
           call la_dlasrt('I',n,d,info)
           180 continue
           return
     end subroutine la_dsterf
     !> QSTERF: computes all eigenvalues of a symmetric tridiagonal matrix
     !> using the Pal-Walker-Kahan variant of the QL or QR algorithm.

     pure subroutine la_qsterf(n,d,e,info)
        use la_constants_qp,only:zero,one,two,three
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: d(*),e(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxit = 30

           ! Local Scalars
           integer(ilp) :: i,iscale,jtot,l,l1,lend,lendsv,lsv,m,nmaxit
           real(qp) :: alpha,anorm,bb,c,eps,eps2,gamma,oldc,oldgam,p,r,rt1,rt2,rte, &
                     s,safmax,safmin,sigma,ssfmax,ssfmin,rmax
           ! Intrinsic Functions
           intrinsic :: abs,sign,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           ! quick return if possible
           if (n < 0) then
              info = -1
              call la_xerbla('QSTERF',-info)
              return
           end if
           if (n <= 1) return
           ! determine the unit roundoff for this environment.
           eps = la_qlamch('E')
           eps2 = eps**2
           safmin = la_qlamch('S')
           safmax = one/safmin
           ssfmax = sqrt(safmax)/three
           ssfmin = sqrt(safmin)/eps2
           rmax = la_qlamch('O')
           ! compute the eigenvalues of the tridiagonal matrix.
           nmaxit = n*maxit
           sigma = zero
           jtot = 0
           ! determine where the matrix splits and choose ql or qr iteration
           ! for each block, according to whether top or bottom diagonal
           ! element is smaller.
           l1 = 1
           10 continue
           if (l1 > n) go to 170
           if (l1 > 1) e(l1 - 1) = zero
           do m = l1,n - 1
              if (abs(e(m)) <= (sqrt(abs(d(m)))*sqrt(abs(d(m + 1))))*eps) &
                        then
                 e(m) = zero
                 go to 30
              end if
           end do
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
           if ((anorm > ssfmax)) then
              iscale = 1
              call la_qlascl('G',0,0,anorm,ssfmax,lend - l + 1,1,d(l),n,info)
              call la_qlascl('G',0,0,anorm,ssfmax,lend - l,1,e(l),n,info)
           else if (anorm < ssfmin) then
              iscale = 2
              call la_qlascl('G',0,0,anorm,ssfmin,lend - l + 1,1,d(l),n,info)
              call la_qlascl('G',0,0,anorm,ssfmin,lend - l,1,e(l),n,info)
           end if
           do i = l,lend - 1
              e(i) = e(i)**2
           end do
           ! choose between ql and qr iteration
           if (abs(d(lend)) < abs(d(l))) then
              lend = lsv
              l = lendsv
           end if
           if (lend >= l) then
              ! ql iteration
              ! look for small subdiagonal element.
              50 continue
              if (l /= lend) then
                 do m = l,lend - 1
                    if (abs(e(m)) <= eps2*abs(d(m)*d(m + 1))) go to 70
                 end do
              end if
              m = lend
              70 continue
              if (m < lend) e(m) = zero
              p = d(l)
              if (m == l) go to 90
              ! if remaining matrix is 2 by 2, use la_qlae2 to compute its
              ! eigenvalues.
              if (m == l + 1) then
                 rte = sqrt(e(l))
                 call la_qlae2(d(l),rte,d(l + 1),rt1,rt2)
                 d(l) = rt1
                 d(l + 1) = rt2
                 e(l) = zero
                 l = l + 2
                 if (l <= lend) go to 50
                 go to 150
              end if
              if (jtot == nmaxit) go to 150
              jtot = jtot + 1
              ! form shift.
              rte = sqrt(e(l))
              sigma = (d(l + 1) - p)/(two*rte)
              r = la_qlapy2(sigma,one)
              sigma = p - (rte/(sigma + sign(r,sigma)))
              c = one
              s = zero
              gamma = d(m) - sigma
              p = gamma*gamma
              ! inner loop
              do i = m - 1,l,-1
                 bb = e(i)
                 r = p + bb
                 if (i /= m - 1) e(i + 1) = s*r
                 oldc = c
                 c = p/r
                 s = bb/r
                 oldgam = gamma
                 alpha = d(i)
                 gamma = c*(alpha - sigma) - s*oldgam
                 d(i + 1) = oldgam + (alpha - gamma)
                 if (c /= zero) then
                    p = (gamma*gamma)/c
                 else
                    p = oldc*bb
                 end if
              end do
              e(l) = s*p
              d(l) = sigma + gamma
              go to 50
              ! eigenvalue found.
              90 continue
              d(l) = p
              l = l + 1
              if (l <= lend) go to 50
              go to 150
           else
              ! qr iteration
              ! look for small superdiagonal element.
              100 continue
              do m = l,lend + 1,-1
                 if (abs(e(m - 1)) <= eps2*abs(d(m)*d(m - 1))) go to 120
              end do
              m = lend
              120 continue
              if (m > lend) e(m - 1) = zero
              p = d(l)
              if (m == l) go to 140
              ! if remaining matrix is 2 by 2, use la_qlae2 to compute its
              ! eigenvalues.
              if (m == l - 1) then
                 rte = sqrt(e(l - 1))
                 call la_qlae2(d(l),rte,d(l - 1),rt1,rt2)
                 d(l) = rt1
                 d(l - 1) = rt2
                 e(l - 1) = zero
                 l = l - 2
                 if (l >= lend) go to 100
                 go to 150
              end if
              if (jtot == nmaxit) go to 150
              jtot = jtot + 1
              ! form shift.
              rte = sqrt(e(l - 1))
              sigma = (d(l - 1) - p)/(two*rte)
              r = la_qlapy2(sigma,one)
              sigma = p - (rte/(sigma + sign(r,sigma)))
              c = one
              s = zero
              gamma = d(m) - sigma
              p = gamma*gamma
              ! inner loop
              do i = m,l - 1
                 bb = e(i)
                 r = p + bb
                 if (i /= m) e(i - 1) = s*r
                 oldc = c
                 c = p/r
                 s = bb/r
                 oldgam = gamma
                 alpha = d(i + 1)
                 gamma = c*(alpha - sigma) - s*oldgam
                 d(i) = oldgam + (alpha - gamma)
                 if (c /= zero) then
                    p = (gamma*gamma)/c
                 else
                    p = oldc*bb
                 end if
              end do
              e(l - 1) = s*p
              d(l) = sigma + gamma
              go to 100
              ! eigenvalue found.
              140 continue
              d(l) = p
              l = l - 1
              if (l >= lend) go to 100
              go to 150
           end if
           ! undo scaling if necessary
           150 continue
           if (iscale == 1) call la_qlascl('G',0,0,ssfmax,anorm,lendsv - lsv + 1,1,d(lsv), &
                     n,info)
           if (iscale == 2) call la_qlascl('G',0,0,ssfmin,anorm,lendsv - lsv + 1,1,d(lsv), &
                     n,info)
           ! check for no convergence to an eigenvalue after a total
           ! of n*maxit iterations.
           if (jtot < nmaxit) go to 10
           do i = 1,n - 1
              if (e(i) /= zero) info = info + 1
           end do
           go to 180
           ! sort eigenvalues in increasing order.
           170 continue
           call la_qlasrt('I',n,d,info)
           180 continue
           return
     end subroutine la_qsterf

     !> SSTEV: computes all eigenvalues and, optionally, eigenvectors of a
     !> real symmetric tridiagonal matrix A.

     pure subroutine la_sstev(jobz,n,d,e,z,ldz,work,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantz
           integer(ilp) :: imax,iscale
           real(sp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tnrm
           ! Intrinsic Functions
           intrinsic :: sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('SSTEV ',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = sqrt(bignum)
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           tnrm = la_slanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_sscal(n,sigma,d,1)
              call la_sscal(n - 1,sigma,e(1),1)
           end if
           ! for eigenvalues only, call la_ssterf.  for eigenvalues and
           ! eigenvectors, call la_ssteqr.
           if (.not. wantz) then
              call la_ssterf(n,d,e,info)
           else
              call la_ssteqr('I',n,d,e,z,ldz,work,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           if (iscale == 1) then
              if (info == 0) then
                 imax = n
              else
                 imax = info - 1
              end if
              call la_sscal(imax,one/sigma,d,1)
           end if
           return
     end subroutine la_sstev
     !> DSTEV: computes all eigenvalues and, optionally, eigenvectors of a
     !> real symmetric tridiagonal matrix A.

     pure subroutine la_dstev(jobz,n,d,e,z,ldz,work,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantz
           integer(ilp) :: imax,iscale
           real(dp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tnrm
           ! Intrinsic Functions
           intrinsic :: sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('DSTEV ',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = sqrt(bignum)
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           tnrm = la_dlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_dscal(n,sigma,d,1)
              call la_dscal(n - 1,sigma,e(1),1)
           end if
           ! for eigenvalues only, call la_dsterf.  for eigenvalues and
           ! eigenvectors, call la_dsteqr.
           if (.not. wantz) then
              call la_dsterf(n,d,e,info)
           else
              call la_dsteqr('I',n,d,e,z,ldz,work,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           if (iscale == 1) then
              if (info == 0) then
                 imax = n
              else
                 imax = info - 1
              end if
              call la_dscal(imax,one/sigma,d,1)
           end if
           return
     end subroutine la_dstev
     !> QSTEV: computes all eigenvalues and, optionally, eigenvectors of a
     !> real symmetric tridiagonal matrix A.

     pure subroutine la_qstev(jobz,n,d,e,z,ldz,work,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,n
           ! Array Arguments
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: wantz
           integer(ilp) :: imax,iscale
           real(qp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tnrm
           ! Intrinsic Functions
           intrinsic :: sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -6
           end if
           if (info /= 0) then
              call la_xerbla('QSTEV ',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = sqrt(bignum)
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           tnrm = la_qlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_qscal(n,sigma,d,1)
              call la_qscal(n - 1,sigma,e(1),1)
           end if
           ! for eigenvalues only, call la_qsterf.  for eigenvalues and
           ! eigenvectors, call la_qsteqr.
           if (.not. wantz) then
              call la_qsterf(n,d,e,info)
           else
              call la_qsteqr('I',n,d,e,z,ldz,work,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           if (iscale == 1) then
              if (info == 0) then
                 imax = n
              else
                 imax = info - 1
              end if
              call la_qscal(imax,one/sigma,d,1)
           end if
           return
     end subroutine la_qstev

     !> SSTEVX: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix A.  Eigenvalues and
     !> eigenvectors can be selected by specifying either a range of values
     !> or a range of indices for the desired eigenvalues.

     pure subroutine la_sstevx(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               work,iwork,ifail,info)
        use la_constants_sp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: alleig,indeig,test,valeig,wantz
           character :: order
           integer(ilp) :: i,imax,indibl,indisp,indiwo,indwrk,iscale,itmp1,j,jj, &
                     nsplit
           real(sp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tmp1,tnrm,vll, &
                     vuu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else
              if (valeig) then
                 if (n > 0 .and. vu <= vl) info = -7
              else if (indeig) then
                 if (il < 1 .or. il > max(1,n)) then
                    info = -8
                 else if (iu < min(n,il) .or. iu > n) then
                    info = -9
                 end if
              end if
           end if
           if (info == 0) then
              if (ldz < 1 .or. (wantz .and. ldz < n)) info = -14
           end if
           if (info /= 0) then
              call la_xerbla('SSTEVX',-info)
              return
           end if
           ! quick return if possible
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (vl < d(1) .and. vu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           if (valeig) then
              vll = vl
              vuu = vu
           else
              vll = zero
              vuu = zero
           end if
           tnrm = la_slanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_sscal(n,sigma,d,1)
              call la_sscal(n - 1,sigma,e(1),1)
              if (valeig) then
                 vll = vl*sigma
                 vuu = vu*sigma
              end if
           end if
           ! if all eigenvalues are desired and abstol is less than zero, then
           ! call la_ssterf or la_ssteqr.  if this fails for some eigenvalue, then
           ! try la_sstebz.
           test = .false.
           if (indeig) then
              if (il == 1 .and. iu == n) then
                 test = .true.
              end if
           end if
           if ((alleig .or. test) .and. (abstol <= zero)) then
              call la_scopy(n,d,1,w,1)
              call la_scopy(n - 1,e(1),1,work(1),1)
              indwrk = n + 1
              if (.not. wantz) then
                 call la_ssterf(n,w,work,info)
              else
                 call la_ssteqr('I',n,w,work,z,ldz,work(indwrk),info)
                 if (info == 0) then
                    do i = 1,n
                       ifail(i) = 0
                    end do
                 end if
              end if
              if (info == 0) then
                 m = n
                 go to 20
              end if
              info = 0
           end if
           ! otherwise, call la_sstebz and, if eigenvectors are desired, la_sstein.
           if (wantz) then
              order = 'B'
           else
              order = 'E'
           end if
           indwrk = 1
           indibl = 1
           indisp = indibl + n
           indiwo = indisp + n
           call la_sstebz(range,order,n,vll,vuu,il,iu,abstol,d,e,m,nsplit,w, &
                     iwork(indibl),iwork(indisp),work(indwrk),iwork(indiwo),info)
           if (wantz) then
              call la_sstein(n,d,e,m,w,iwork(indibl),iwork(indisp),z,ldz,work( &
                        indwrk),iwork(indiwo),ifail,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           20 continue
           if (iscale == 1) then
              if (info == 0) then
                 imax = m
              else
                 imax = info - 1
              end if
              call la_sscal(imax,one/sigma,w,1)
           end if
           ! if eigenvalues are not in order, then sort them, along with
           ! eigenvectors.
           if (wantz) then
              do j = 1,m - 1
                 i = 0
                 tmp1 = w(j)
                 do jj = j + 1,m
                    if (w(jj) < tmp1) then
                       i = jj
                       tmp1 = w(jj)
                    end if
                 end do
                 if (i /= 0) then
                    itmp1 = iwork(indibl + i - 1)
                    w(i) = w(j)
                    iwork(indibl + i - 1) = iwork(indibl + j - 1)
                    w(j) = tmp1
                    iwork(indibl + j - 1) = itmp1
                    call la_sswap(n,z(1,i),1,z(1,j),1)
                    if (info /= 0) then
                       itmp1 = ifail(i)
                       ifail(i) = ifail(j)
                       ifail(j) = itmp1
                    end if
                 end if
              end do
           end if
           return
     end subroutine la_sstevx
     !> DSTEVX: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix A.  Eigenvalues and
     !> eigenvectors can be selected by specifying either a range of values
     !> or a range of indices for the desired eigenvalues.

     pure subroutine la_dstevx(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               work,iwork,ifail,info)
        use la_constants_dp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: alleig,indeig,test,valeig,wantz
           character :: order
           integer(ilp) :: i,imax,indibl,indisp,indiwo,indwrk,iscale,itmp1,j,jj, &
                     nsplit
           real(dp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tmp1,tnrm,vll, &
                     vuu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else
              if (valeig) then
                 if (n > 0 .and. vu <= vl) info = -7
              else if (indeig) then
                 if (il < 1 .or. il > max(1,n)) then
                    info = -8
                 else if (iu < min(n,il) .or. iu > n) then
                    info = -9
                 end if
              end if
           end if
           if (info == 0) then
              if (ldz < 1 .or. (wantz .and. ldz < n)) info = -14
           end if
           if (info /= 0) then
              call la_xerbla('DSTEVX',-info)
              return
           end if
           ! quick return if possible
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (vl < d(1) .and. vu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           if (valeig) then
              vll = vl
              vuu = vu
           else
              vll = zero
              vuu = zero
           end if
           tnrm = la_dlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_dscal(n,sigma,d,1)
              call la_dscal(n - 1,sigma,e(1),1)
              if (valeig) then
                 vll = vl*sigma
                 vuu = vu*sigma
              end if
           end if
           ! if all eigenvalues are desired and abstol is less than zero, then
           ! call la_dsterf or la_ssteqr.  if this fails for some eigenvalue, then
           ! try la_dstebz.
           test = .false.
           if (indeig) then
              if (il == 1 .and. iu == n) then
                 test = .true.
              end if
           end if
           if ((alleig .or. test) .and. (abstol <= zero)) then
              call la_dcopy(n,d,1,w,1)
              call la_dcopy(n - 1,e(1),1,work(1),1)
              indwrk = n + 1
              if (.not. wantz) then
                 call la_dsterf(n,w,work,info)
              else
                 call la_dsteqr('I',n,w,work,z,ldz,work(indwrk),info)
                 if (info == 0) then
                    do i = 1,n
                       ifail(i) = 0
                    end do
                 end if
              end if
              if (info == 0) then
                 m = n
                 go to 20
              end if
              info = 0
           end if
           ! otherwise, call la_dstebz and, if eigenvectors are desired, la_sstein.
           if (wantz) then
              order = 'B'
           else
              order = 'E'
           end if
           indwrk = 1
           indibl = 1
           indisp = indibl + n
           indiwo = indisp + n
           call la_dstebz(range,order,n,vll,vuu,il,iu,abstol,d,e,m,nsplit,w, &
                     iwork(indibl),iwork(indisp),work(indwrk),iwork(indiwo),info)
           if (wantz) then
              call la_dstein(n,d,e,m,w,iwork(indibl),iwork(indisp),z,ldz,work( &
                        indwrk),iwork(indiwo),ifail,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           20 continue
           if (iscale == 1) then
              if (info == 0) then
                 imax = m
              else
                 imax = info - 1
              end if
              call la_dscal(imax,one/sigma,w,1)
           end if
           ! if eigenvalues are not in order, then sort them, along with
           ! eigenvectors.
           if (wantz) then
              do j = 1,m - 1
                 i = 0
                 tmp1 = w(j)
                 do jj = j + 1,m
                    if (w(jj) < tmp1) then
                       i = jj
                       tmp1 = w(jj)
                    end if
                 end do
                 if (i /= 0) then
                    itmp1 = iwork(indibl + i - 1)
                    w(i) = w(j)
                    iwork(indibl + i - 1) = iwork(indibl + j - 1)
                    w(j) = tmp1
                    iwork(indibl + j - 1) = itmp1
                    call la_dswap(n,z(1,i),1,z(1,j),1)
                    if (info /= 0) then
                       itmp1 = ifail(i)
                       ifail(i) = ifail(j)
                       ifail(j) = itmp1
                    end if
                 end if
              end do
           end if
           return
     end subroutine la_dstevx
     !> QSTEVX: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix A.  Eigenvalues and
     !> eigenvectors can be selected by specifying either a range of values
     !> or a range of indices for the desired eigenvalues.

     pure subroutine la_qstevx(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               work,iwork,ifail,info)
        use la_constants_qp,only:zero,one
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: alleig,indeig,test,valeig,wantz
           character :: order
           integer(ilp) :: i,imax,indibl,indisp,indiwo,indwrk,iscale,itmp1,j,jj, &
                     nsplit
           real(qp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tmp1,tnrm,vll, &
                     vuu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else
              if (valeig) then
                 if (n > 0 .and. vu <= vl) info = -7
              else if (indeig) then
                 if (il < 1 .or. il > max(1,n)) then
                    info = -8
                 else if (iu < min(n,il) .or. iu > n) then
                    info = -9
                 end if
              end if
           end if
           if (info == 0) then
              if (ldz < 1 .or. (wantz .and. ldz < n)) info = -14
           end if
           if (info /= 0) then
              call la_xerbla('QSTEVX',-info)
              return
           end if
           ! quick return if possible
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (vl < d(1) .and. vu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           if (valeig) then
              vll = vl
              vuu = vu
           else
              vll = zero
              vuu = zero
           end if
           tnrm = la_qlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_qscal(n,sigma,d,1)
              call la_qscal(n - 1,sigma,e(1),1)
              if (valeig) then
                 vll = vl*sigma
                 vuu = vu*sigma
              end if
           end if
           ! if all eigenvalues are desired and abstol is less than zero, then
           ! call la_qsterf or la_dsteqr.  if this fails for some eigenvalue, then
           ! try la_qstebz.
           test = .false.
           if (indeig) then
              if (il == 1 .and. iu == n) then
                 test = .true.
              end if
           end if
           if ((alleig .or. test) .and. (abstol <= zero)) then
              call la_qcopy(n,d,1,w,1)
              call la_qcopy(n - 1,e(1),1,work(1),1)
              indwrk = n + 1
              if (.not. wantz) then
                 call la_qsterf(n,w,work,info)
              else
                 call la_qsteqr('I',n,w,work,z,ldz,work(indwrk),info)
                 if (info == 0) then
                    do i = 1,n
                       ifail(i) = 0
                    end do
                 end if
              end if
              if (info == 0) then
                 m = n
                 go to 20
              end if
              info = 0
           end if
           ! otherwise, call la_qstebz and, if eigenvectors are desired, la_dstein.
           if (wantz) then
              order = 'B'
           else
              order = 'E'
           end if
           indwrk = 1
           indibl = 1
           indisp = indibl + n
           indiwo = indisp + n
           call la_qstebz(range,order,n,vll,vuu,il,iu,abstol,d,e,m,nsplit,w, &
                     iwork(indibl),iwork(indisp),work(indwrk),iwork(indiwo),info)
           if (wantz) then
              call la_qstein(n,d,e,m,w,iwork(indibl),iwork(indisp),z,ldz,work( &
                        indwrk),iwork(indiwo),ifail,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           20 continue
           if (iscale == 1) then
              if (info == 0) then
                 imax = m
              else
                 imax = info - 1
              end if
              call la_qscal(imax,one/sigma,w,1)
           end if
           ! if eigenvalues are not in order, then sort them, along with
           ! eigenvectors.
           if (wantz) then
              do j = 1,m - 1
                 i = 0
                 tmp1 = w(j)
                 do jj = j + 1,m
                    if (w(jj) < tmp1) then
                       i = jj
                       tmp1 = w(jj)
                    end if
                 end do
                 if (i /= 0) then
                    itmp1 = iwork(indibl + i - 1)
                    w(i) = w(j)
                    iwork(indibl + i - 1) = iwork(indibl + j - 1)
                    w(j) = tmp1
                    iwork(indibl + j - 1) = itmp1
                    call la_qswap(n,z(1,i),1,z(1,j),1)
                    if (info /= 0) then
                       itmp1 = ifail(i)
                       ifail(i) = ifail(j)
                       ifail(j) = itmp1
                    end if
                 end if
              end do
           end if
           return
     end subroutine la_qstevx

     !> SSTEDC: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.
     !> The eigenvectors of a full or band real symmetric matrix can also be
     !> found if SSYTRD or SSPTRD or SSBTRD has been used to reduce this
     !> matrix to tridiagonal form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See SLAED3 for details.

     pure subroutine la_sstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one,two

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: finish,i,icompz,ii,j,k,lgn,liwmin,lwmin,m,smlsiz,start, &
                     storez,strtrw
           real(sp) :: eps,orgnrm,p,tiny
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,mod,real,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
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
           if (info == 0) then
              ! compute the workspace requirements
              smlsiz = la_ilaenv(9,'SSTEDC',' ',0,0,0,0)
              if (n <= 1 .or. icompz == 0) then
                 liwmin = 1
                 lwmin = 1
              else if (n <= smlsiz) then
                 liwmin = 1
                 lwmin = 2*(n - 1)
              else
                 lgn = int(log(real(n,KIND=sp))/log(two),KIND=ilp)
                 if (2**lgn < n) lgn = lgn + 1
                 if (2**lgn < n) lgn = lgn + 1
                 if (icompz == 1) then
                    lwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                    liwmin = 6 + 6*n + 5*n*lgn
                 else if (icompz == 2) then
                    lwmin = 1 + 4*n + n**2
                    liwmin = 3 + 5*n
                 end if
              end if
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -10
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SSTEDC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz /= 0) z(1,1) = one
              return
           end if
           ! if the following conditional clause is removed, then the routine
           ! will use the divide and conquer routine to compute only the
           ! eigenvalues, which requires (3n + 3n**2) real workspace and
           ! (2 + 5n + 2n lg(n)) integer workspace.
           ! since on many architectures la_ssterf is much faster than any other
           ! algorithm for finding eigenvalues only, it is used here
           ! as the default. if the conditional clause is removed, then
           ! information on the size of workspace needs to be changed.
           ! if compz = 'n', use la_ssterf to compute the eigenvalues.
           if (icompz == 0) then
              call la_ssterf(n,d,e,info)
              go to 50
           end if
           ! if n is smaller than the minimum divide size (smlsiz+1), then
           ! solve the problem with another solver.
           if (n <= smlsiz) then
              call la_ssteqr(compz,n,d,e,z,ldz,work,info)
           else
              ! if compz = 'v', the z matrix must be stored elsewhere for later
              ! use.
              if (icompz == 1) then
                 storez = 1 + n*n
              else
                 storez = 1
              end if
              if (icompz == 2) then
                 call la_slaset('FULL',n,n,zero,one,z,ldz)
              end if
              ! scale.
              orgnrm = la_slanst('M',n,d,e)
              if (orgnrm == zero) go to 50
              eps = la_slamch('EPSILON')
              start = 1
              ! while ( start <= n )
              10 continue
              if (start <= n) then
                 ! let finish be the position of the next subdiagonal entry
                 ! such that e( finish ) <= tiny or finish = n if no such
                 ! subdiagonal exists.  the matrix identified by the elements
                 ! between start and finish constitutes an independent
                 ! sub-problem.
                 finish = start
                 20 continue
                 if (finish < n) then
                    tiny = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                    if (abs(e(finish)) > tiny) then
                       finish = finish + 1
                       go to 20
                    end if
                 end if
                 ! (sub) problem determined.  compute its size and solve it.
                 m = finish - start + 1
                 if (m == 1) then
                    start = finish + 1
                    go to 10
                 end if
                 if (m > smlsiz) then
                    ! scale.
                    orgnrm = la_slanst('M',m,d(start),e(start))
                    call la_slascl('G',0,0,orgnrm,one,m,1,d(start),m,info)
                    call la_slascl('G',0,0,orgnrm,one,m - 1,1,e(start),m - 1,info)

                    if (icompz == 1) then
                       strtrw = 1
                    else
                       strtrw = start
                    end if
                    call la_slaed0(icompz,n,m,d(start),e(start),z(strtrw,start), &
                              ldz,work(1),n,work(storez),iwork,info)
                    if (info /= 0) then
                       info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, (m + 1)) + start - &
                                 1
                       go to 50
                    end if
                    ! scale back.
                    call la_slascl('G',0,0,one,orgnrm,m,1,d(start),m,info)
                 else
                    if (icompz == 1) then
                       ! since qr won't update a z matrix which is larger than
                       ! the length of d, we must solve the sub-problem in a
                       ! workspace and then multiply back into z.
                       call la_ssteqr('I',m,d(start),e(start),work,m,work(m*m + 1), &
                                 info)
                       call la_slacpy('A',n,m,z(1,start),ldz,work(storez),n)

                       call la_sgemm('N','N',n,m,m,one,work(storez),n,work,m,zero, &
                                 z(1,start),ldz)
                    else if (icompz == 2) then
                       call la_ssteqr('I',m,d(start),e(start),z(start,start),ldz, &
                                 work,info)
                    else
                       call la_ssterf(m,d(start),e(start),info)
                    end if
                    if (info /= 0) then
                       info = start*(n + 1) + finish
                       go to 50
                    end if
                 end if
                 start = finish + 1
                 go to 10
              end if
              ! endwhile
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
           end if
           50 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_sstedc
     !> DSTEDC: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.
     !> The eigenvectors of a full or band real symmetric matrix can also be
     !> found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
     !> matrix to tridiagonal form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See DLAED3 for details.

     pure subroutine la_dstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one,two

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: finish,i,icompz,ii,j,k,lgn,liwmin,lwmin,m,smlsiz,start, &
                     storez,strtrw
           real(dp) :: eps,orgnrm,p,tiny
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max,mod,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
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
           if (info == 0) then
              ! compute the workspace requirements
              smlsiz = la_ilaenv(9,'DSTEDC',' ',0,0,0,0)
              if (n <= 1 .or. icompz == 0) then
                 liwmin = 1
                 lwmin = 1
              else if (n <= smlsiz) then
                 liwmin = 1
                 lwmin = 2*(n - 1)
              else
                 lgn = int(log(real(n,KIND=dp))/log(two),KIND=ilp)
                 if (2**lgn < n) lgn = lgn + 1
                 if (2**lgn < n) lgn = lgn + 1
                 if (icompz == 1) then
                    lwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                    liwmin = 6 + 6*n + 5*n*lgn
                 else if (icompz == 2) then
                    lwmin = 1 + 4*n + n**2
                    liwmin = 3 + 5*n
                 end if
              end if
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -10
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DSTEDC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz /= 0) z(1,1) = one
              return
           end if
           ! if the following conditional clause is removed, then the routine
           ! will use the divide and conquer routine to compute only the
           ! eigenvalues, which requires (3n + 3n**2) real workspace and
           ! (2 + 5n + 2n lg(n)) integer workspace.
           ! since on many architectures la_dsterf is much faster than any other
           ! algorithm for finding eigenvalues only, it is used here
           ! as the default. if the conditional clause is removed, then
           ! information on the size of workspace needs to be changed.
           ! if compz = 'n', use la_dsterf to compute the eigenvalues.
           if (icompz == 0) then
              call la_dsterf(n,d,e,info)
              go to 50
           end if
           ! if n is smaller than the minimum divide size (smlsiz+1), then
           ! solve the problem with another solver.
           if (n <= smlsiz) then
              call la_dsteqr(compz,n,d,e,z,ldz,work,info)
           else
              ! if compz = 'v', the z matrix must be stored elsewhere for later
              ! use.
              if (icompz == 1) then
                 storez = 1 + n*n
              else
                 storez = 1
              end if
              if (icompz == 2) then
                 call la_dlaset('FULL',n,n,zero,one,z,ldz)
              end if
              ! scale.
              orgnrm = la_dlanst('M',n,d,e)
              if (orgnrm == zero) go to 50
              eps = la_dlamch('EPSILON')
              start = 1
              ! while ( start <= n )
              10 continue
              if (start <= n) then
                 ! let finish be the position of the next subdiagonal entry
                 ! such that e( finish ) <= tiny or finish = n if no such
                 ! subdiagonal exists.  the matrix identified by the elements
                 ! between start and finish constitutes an independent
                 ! sub-problem.
                 finish = start
                 20 continue
                 if (finish < n) then
                    tiny = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                    if (abs(e(finish)) > tiny) then
                       finish = finish + 1
                       go to 20
                    end if
                 end if
                 ! (sub) problem determined.  compute its size and solve it.
                 m = finish - start + 1
                 if (m == 1) then
                    start = finish + 1
                    go to 10
                 end if
                 if (m > smlsiz) then
                    ! scale.
                    orgnrm = la_dlanst('M',m,d(start),e(start))
                    call la_dlascl('G',0,0,orgnrm,one,m,1,d(start),m,info)
                    call la_dlascl('G',0,0,orgnrm,one,m - 1,1,e(start),m - 1,info)

                    if (icompz == 1) then
                       strtrw = 1
                    else
                       strtrw = start
                    end if
                    call la_dlaed0(icompz,n,m,d(start),e(start),z(strtrw,start), &
                              ldz,work(1),n,work(storez),iwork,info)
                    if (info /= 0) then
                       info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, (m + 1)) + start - &
                                 1
                       go to 50
                    end if
                    ! scale back.
                    call la_dlascl('G',0,0,one,orgnrm,m,1,d(start),m,info)
                 else
                    if (icompz == 1) then
                       ! since qr won't update a z matrix which is larger than
                       ! the length of d, we must solve the sub-problem in a
                       ! workspace and then multiply back into z.
                       call la_dsteqr('I',m,d(start),e(start),work,m,work(m*m + 1), &
                                 info)
                       call la_dlacpy('A',n,m,z(1,start),ldz,work(storez),n)

                       call la_dgemm('N','N',n,m,m,one,work(storez),n,work,m,zero, &
                                 z(1,start),ldz)
                    else if (icompz == 2) then
                       call la_dsteqr('I',m,d(start),e(start),z(start,start),ldz, &
                                 work,info)
                    else
                       call la_dsterf(m,d(start),e(start),info)
                    end if
                    if (info /= 0) then
                       info = start*(n + 1) + finish
                       go to 50
                    end if
                 end if
                 start = finish + 1
                 go to 10
              end if
              ! endwhile
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
           end if
           50 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_dstedc
     !> QSTEDC: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.
     !> The eigenvectors of a full or band real symmetric matrix can also be
     !> found if QSYTRD or QSPTRD or QSBTRD has been used to reduce this
     !> matrix to tridiagonal form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See QLAED3 for details.

     pure subroutine la_qstedc(compz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one,two

        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),e(*),z(ldz,*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: finish,i,icompz,ii,j,k,lgn,liwmin,lwmin,m,smlsiz,start, &
                     storez,strtrw
           real(qp) :: eps,orgnrm,p,tiny
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max,mod,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1 .or. liwork == -1)
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
           if (info == 0) then
              ! compute the workspace requirements
              smlsiz = la_ilaenv(9,'QSTEDC',' ',0,0,0,0)
              if (n <= 1 .or. icompz == 0) then
                 liwmin = 1
                 lwmin = 1
              else if (n <= smlsiz) then
                 liwmin = 1
                 lwmin = 2*(n - 1)
              else
                 lgn = int(log(real(n,KIND=qp))/log(two),KIND=ilp)
                 if (2**lgn < n) lgn = lgn + 1
                 if (2**lgn < n) lgn = lgn + 1
                 if (icompz == 1) then
                    lwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                    liwmin = 6 + 6*n + 5*n*lgn
                 else if (icompz == 2) then
                    lwmin = 1 + 4*n + n**2
                    liwmin = 3 + 5*n
                 end if
              end if
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -10
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QSTEDC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz /= 0) z(1,1) = one
              return
           end if
           ! if the following conditional clause is removed, then the routine
           ! will use the divide and conquer routine to compute only the
           ! eigenvalues, which requires (3n + 3n**2) real workspace and
           ! (2 + 5n + 2n lg(n)) integer workspace.
           ! since on many architectures la_qsterf is much faster than any other
           ! algorithm for finding eigenvalues only, it is used here
           ! as the default. if the conditional clause is removed, then
           ! information on the size of workspace needs to be changed.
           ! if compz = 'n', use la_qsterf to compute the eigenvalues.
           if (icompz == 0) then
              call la_qsterf(n,d,e,info)
              go to 50
           end if
           ! if n is smaller than the minimum divide size (smlsiz+1), then
           ! solve the problem with another solver.
           if (n <= smlsiz) then
              call la_qsteqr(compz,n,d,e,z,ldz,work,info)
           else
              ! if compz = 'v', the z matrix must be stored elsewhere for later
              ! use.
              if (icompz == 1) then
                 storez = 1 + n*n
              else
                 storez = 1
              end if
              if (icompz == 2) then
                 call la_qlaset('FULL',n,n,zero,one,z,ldz)
              end if
              ! scale.
              orgnrm = la_qlanst('M',n,d,e)
              if (orgnrm == zero) go to 50
              eps = la_qlamch('EPSILON')
              start = 1
              ! while ( start <= n )
              10 continue
              if (start <= n) then
                 ! let finish be the position of the next subdiagonal entry
                 ! such that e( finish ) <= tiny or finish = n if no such
                 ! subdiagonal exists.  the matrix identified by the elements
                 ! between start and finish constitutes an independent
                 ! sub-problem.
                 finish = start
                 20 continue
                 if (finish < n) then
                    tiny = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                    if (abs(e(finish)) > tiny) then
                       finish = finish + 1
                       go to 20
                    end if
                 end if
                 ! (sub) problem determined.  compute its size and solve it.
                 m = finish - start + 1
                 if (m == 1) then
                    start = finish + 1
                    go to 10
                 end if
                 if (m > smlsiz) then
                    ! scale.
                    orgnrm = la_qlanst('M',m,d(start),e(start))
                    call la_qlascl('G',0,0,orgnrm,one,m,1,d(start),m,info)
                    call la_qlascl('G',0,0,orgnrm,one,m - 1,1,e(start),m - 1,info)

                    if (icompz == 1) then
                       strtrw = 1
                    else
                       strtrw = start
                    end if
                    call la_qlaed0(icompz,n,m,d(start),e(start),z(strtrw,start), &
                              ldz,work(1),n,work(storez),iwork,info)
                    if (info /= 0) then
                       info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, (m + 1)) + start - &
                                 1
                       go to 50
                    end if
                    ! scale back.
                    call la_qlascl('G',0,0,one,orgnrm,m,1,d(start),m,info)
                 else
                    if (icompz == 1) then
                       ! since qr won't update a z matrix which is larger than
                       ! the length of d, we must solve the sub-problem in a
                       ! workspace and then multiply back into z.
                       call la_qsteqr('I',m,d(start),e(start),work,m,work(m*m + 1), &
                                 info)
                       call la_qlacpy('A',n,m,z(1,start),ldz,work(storez),n)

                       call la_qgemm('N','N',n,m,m,one,work(storez),n,work,m,zero, &
                                 z(1,start),ldz)
                    else if (icompz == 2) then
                       call la_qsteqr('I',m,d(start),e(start),z(start,start),ldz, &
                                 work,info)
                    else
                       call la_qsterf(m,d(start),e(start),info)
                    end if
                    if (info /= 0) then
                       info = start*(n + 1) + finish
                       go to 50
                    end if
                 end if
                 start = finish + 1
                 go to 10
              end if
              ! endwhile
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
           end if
           50 continue
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_qstedc

     !> SSTEVD: computes all eigenvalues and, optionally, eigenvectors of a
     !> real symmetric tridiagonal matrix. If eigenvectors are desired, it
     !> uses a divide and conquer algorithm.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_sstevd(jobz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantz
           integer(ilp) :: iscale,liwmin,lwmin
           real(sp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tnrm
           ! Intrinsic Functions
           intrinsic :: sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           lquery = (lwork == -1 .or. liwork == -1)
           info = 0
           liwmin = 1
           lwmin = 1
           if (n > 1 .and. wantz) then
              lwmin = 1 + 4*n + n**2
              liwmin = 3 + 5*n
           end if
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -6
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -10
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SSTEVD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = sqrt(bignum)
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           tnrm = la_slanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_sscal(n,sigma,d,1)
              call la_sscal(n - 1,sigma,e(1),1)
           end if
           ! for eigenvalues only, call la_ssterf.  for eigenvalues and
           ! eigenvectors, call la_sstedc.
           if (.not. wantz) then
              call la_ssterf(n,d,e,info)
           else
              call la_sstedc('I',n,d,e,z,ldz,work,lwork,iwork,liwork,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           if (iscale == 1) call la_sscal(n,one/sigma,d,1)
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_sstevd
     !> DSTEVD: computes all eigenvalues and, optionally, eigenvectors of a
     !> real symmetric tridiagonal matrix. If eigenvectors are desired, it
     !> uses a divide and conquer algorithm.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_dstevd(jobz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantz
           integer(ilp) :: iscale,liwmin,lwmin
           real(dp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tnrm
           ! Intrinsic Functions
           intrinsic :: sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           lquery = (lwork == -1 .or. liwork == -1)
           info = 0
           liwmin = 1
           lwmin = 1
           if (n > 1 .and. wantz) then
              lwmin = 1 + 4*n + n**2
              liwmin = 3 + 5*n
           end if
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -6
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -10
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DSTEVD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = sqrt(bignum)
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           tnrm = la_dlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_dscal(n,sigma,d,1)
              call la_dscal(n - 1,sigma,e(1),1)
           end if
           ! for eigenvalues only, call la_dsterf.  for eigenvalues and
           ! eigenvectors, call la_dstedc.
           if (.not. wantz) then
              call la_dsterf(n,d,e,info)
           else
              call la_dstedc('I',n,d,e,z,ldz,work,lwork,iwork,liwork,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           if (iscale == 1) call la_dscal(n,one/sigma,d,1)
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_dstevd
     !> QSTEVD: computes all eigenvalues and, optionally, eigenvectors of a
     !> real symmetric tridiagonal matrix. If eigenvectors are desired, it
     !> uses a divide and conquer algorithm.
     !> The divide and conquer algorithm makes very mild assumptions about
     !> floating point arithmetic. It will work on machines with a guard
     !> digit in add/subtract, or on those binary machines without guard
     !> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
     !> Cray-2. It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.

     pure subroutine la_qstevd(jobz,n,d,e,z,ldz,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one

        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery,wantz
           integer(ilp) :: iscale,liwmin,lwmin
           real(qp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tnrm
           ! Intrinsic Functions
           intrinsic :: sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           lquery = (lwork == -1 .or. liwork == -1)
           info = 0
           liwmin = 1
           lwmin = 1
           if (n > 1 .and. wantz) then
              lwmin = 1 + 4*n + n**2
              liwmin = 3 + 5*n
           end if
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (n < 0) then
              info = -2
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -6
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -10
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QSTEVD',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = sqrt(bignum)
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           tnrm = la_qlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_qscal(n,sigma,d,1)
              call la_qscal(n - 1,sigma,e(1),1)
           end if
           ! for eigenvalues only, call la_qsterf.  for eigenvalues and
           ! eigenvectors, call la_qstedc.
           if (.not. wantz) then
              call la_qsterf(n,d,e,info)
           else
              call la_qstedc('I',n,d,e,z,ldz,work,lwork,iwork,liwork,info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           if (iscale == 1) call la_qscal(n,one/sigma,d,1)
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_qstevd

     !> SPTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric positive definite tridiagonal matrix by first factoring the
     !> matrix using SPTTRF, and then calling SBDSQR to compute the singular
     !> values of the bidiagonal factor.
     !> This routine computes the eigenvalues of the positive definite
     !> tridiagonal matrix to high relative accuracy.  This means that if the
     !> eigenvalues range over many orders of magnitude in size, then the
     !> small eigenvalues and corresponding eigenvectors will be computed
     !> more accurately than, for example, with the standard QR method.
     !> The eigenvectors of a full or band symmetric positive definite matrix
     !> can also be found if SSYTRD, SSPTRD, or SSBTRD has been used to
     !> reduce this matrix to tridiagonal form. (The reduction to tridiagonal
     !> form, however, may preclude the possibility of obtaining high
     !> relative accuracy in the small eigenvalues of the original matrix, if
     !> these eigenvalues range over many orders of magnitude.)

     pure subroutine la_spteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_sp
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

           ! Local Arrays
           real(sp) :: c(1,1),vt(1,1)
           ! Local Scalars
           integer(ilp) :: i,icompz,nru
           ! Intrinsic Functions
           intrinsic :: max,sqrt
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
              call la_xerbla('SPTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz > 0) z(1,1) = one
              return
           end if
           if (icompz == 2) call la_slaset('FULL',n,n,zero,one,z,ldz)
           ! call la_spttrf to factor the matrix.
           call la_spttrf(n,d,e,info)
           if (info /= 0) return
           do i = 1,n
              d(i) = sqrt(d(i))
           end do
           do i = 1,n - 1
              e(i) = e(i)*d(i)
           end do
           ! call la_sbdsqr to compute the singular values/vectors of the
           ! bidiagonal factor.
           if (icompz > 0) then
              nru = n
           else
              nru = 0
           end if
           call la_sbdsqr('LOWER',n,0,nru,0,d,e,vt,1,z,ldz,c,1,work,info)

           ! square the singular values.
           if (info == 0) then
              do i = 1,n
                 d(i) = d(i)*d(i)
              end do
           else
              info = n + info
           end if
           return
     end subroutine la_spteqr
     !> DPTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric positive definite tridiagonal matrix by first factoring the
     !> matrix using DPTTRF, and then calling DBDSQR to compute the singular
     !> values of the bidiagonal factor.
     !> This routine computes the eigenvalues of the positive definite
     !> tridiagonal matrix to high relative accuracy.  This means that if the
     !> eigenvalues range over many orders of magnitude in size, then the
     !> small eigenvalues and corresponding eigenvectors will be computed
     !> more accurately than, for example, with the standard QR method.
     !> The eigenvectors of a full or band symmetric positive definite matrix
     !> can also be found if DSYTRD, DSPTRD, or DSBTRD has been used to
     !> reduce this matrix to tridiagonal form. (The reduction to tridiagonal
     !> form, however, may preclude the possibility of obtaining high
     !> relative accuracy in the small eigenvalues of the original matrix, if
     !> these eigenvalues range over many orders of magnitude.)

     pure subroutine la_dpteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_dp
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

           ! Local Arrays
           real(dp) :: c(1,1),vt(1,1)
           ! Local Scalars
           integer(ilp) :: i,icompz,nru
           ! Intrinsic Functions
           intrinsic :: max,sqrt
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
              call la_xerbla('DPTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz > 0) z(1,1) = one
              return
           end if
           if (icompz == 2) call la_dlaset('FULL',n,n,zero,one,z,ldz)
           ! call la_dpttrf to factor the matrix.
           call la_dpttrf(n,d,e,info)
           if (info /= 0) return
           do i = 1,n
              d(i) = sqrt(d(i))
           end do
           do i = 1,n - 1
              e(i) = e(i)*d(i)
           end do
           ! call la_dbdsqr to compute the singular values/vectors of the
           ! bidiagonal factor.
           if (icompz > 0) then
              nru = n
           else
              nru = 0
           end if
           call la_dbdsqr('LOWER',n,0,nru,0,d,e,vt,1,z,ldz,c,1,work,info)

           ! square the singular values.
           if (info == 0) then
              do i = 1,n
                 d(i) = d(i)*d(i)
              end do
           else
              info = n + info
           end if
           return
     end subroutine la_dpteqr
     !> QPTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric positive definite tridiagonal matrix by first factoring the
     !> matrix using QPTTRF, and then calling QBDSQR to compute the singular
     !> values of the bidiagonal factor.
     !> This routine computes the eigenvalues of the positive definite
     !> tridiagonal matrix to high relative accuracy.  This means that if the
     !> eigenvalues range over many orders of magnitude in size, then the
     !> small eigenvalues and corresponding eigenvectors will be computed
     !> more accurately than, for example, with the standard QR method.
     !> The eigenvectors of a full or band symmetric positive definite matrix
     !> can also be found if QSYTRD, QSPTRD, or QSBTRD has been used to
     !> reduce this matrix to tridiagonal form. (The reduction to tridiagonal
     !> form, however, may preclude the possibility of obtaining high
     !> relative accuracy in the small eigenvalues of the original matrix, if
     !> these eigenvalues range over many orders of magnitude.)

     pure subroutine la_qpteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_qp
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

           ! Local Arrays
           real(qp) :: c(1,1),vt(1,1)
           ! Local Scalars
           integer(ilp) :: i,icompz,nru
           ! Intrinsic Functions
           intrinsic :: max,sqrt
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
              call la_xerbla('QPTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz > 0) z(1,1) = one
              return
           end if
           if (icompz == 2) call la_qlaset('FULL',n,n,zero,one,z,ldz)
           ! call la_qpttrf to factor the matrix.
           call la_qpttrf(n,d,e,info)
           if (info /= 0) return
           do i = 1,n
              d(i) = sqrt(d(i))
           end do
           do i = 1,n - 1
              e(i) = e(i)*d(i)
           end do
           ! call la_qbdsqr to compute the singular values/vectors of the
           ! bidiagonal factor.
           if (icompz > 0) then
              nru = n
           else
              nru = 0
           end if
           call la_qbdsqr('LOWER',n,0,nru,0,d,e,vt,1,z,ldz,c,1,work,info)

           ! square the singular values.
           if (info == 0) then
              do i = 1,n
                 d(i) = d(i)*d(i)
              end do
           else
              info = n + info
           end if
           return
     end subroutine la_qpteqr

     !> SSTEGR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> SSTEGR is a compatibility wrapper around the improved SSTEMR routine.
     !> See SSTEMR for further details.
     !> One important change is that the ABSTOL parameter no longer provides any
     !> benefit and hence is no longer used.
     !> Note : SSTEGR and SSTEMR work only on machines which follow
     !> IEEE-754 floating-point standard in their handling of infinities and
     !> NaNs.  Normal execution may create these exceptiona values and hence
     !> may abort due to a floating point exception in environments which
     !> do not conform to the IEEE-754 standard.

     pure subroutine la_sstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*)
           real(sp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: tryrac
           ! Executable Statements
           info = 0
           tryrac = .false.
           call la_sstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,n,isuppz, &
                     tryrac,work,lwork,iwork,liwork,info)
     end subroutine la_sstegr
     !> DSTEGR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> DSTEGR is a compatibility wrapper around the improved DSTEMR routine.
     !> See DSTEMR for further details.
     !> One important change is that the ABSTOL parameter no longer provides any
     !> benefit and hence is no longer used.
     !> Note : DSTEGR and DSTEMR work only on machines which follow
     !> IEEE-754 floating-point standard in their handling of infinities and
     !> NaNs.  Normal execution may create these exceptiona values and hence
     !> may abort due to a floating point exception in environments which
     !> do not conform to the IEEE-754 standard.

     pure subroutine la_dstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*)
           real(dp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: tryrac
           ! Executable Statements
           info = 0
           tryrac = .false.
           call la_dstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,n,isuppz, &
                     tryrac,work,lwork,iwork,liwork,info)
     end subroutine la_dstegr
     !> QSTEGR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> QSTEGR is a compatibility wrapper around the improved QSTEMR routine.
     !> See QSTEMR for further details.
     !> One important change is that the ABSTOL parameter no longer provides any
     !> benefit and hence is no longer used.
     !> Note : QSTEGR and QSTEMR work only on machines which follow
     !> IEEE-754 floating-point standard in their handling of infinities and
     !> NaNs.  Normal execution may create these exceptiona values and hence
     !> may abort due to a floating point exception in environments which
     !> do not conform to the IEEE-754 standard.

     pure subroutine la_qstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*)
           real(qp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: tryrac
           ! Executable Statements
           info = 0
           tryrac = .false.
           call la_qstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,n,isuppz, &
                     tryrac,work,lwork,iwork,liwork,info)
     end subroutine la_qstegr

     !> SSTEMR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> Depending on the number of desired eigenvalues, these are computed either
     !> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
     !> computed by the use of various suitable L D L^T factorizations near clusters
     !> of close eigenvalues (referred to as RRRs, Relatively Robust
     !> Representations). An informal sketch of the algorithm follows.
     !> For each unreduced block (submatrix) of T,
     !> (a) Compute T - sigma I  = L D L^T, so that L and D
     !> define all the wanted eigenvalues to high relative accuracy.
     !> This means that small relative changes in the entries of D and L
     !> cause only small relative changes in the eigenvalues and
     !> eigenvectors. The standard (unfactored) representation of the
     !> tridiagonal matrix T does not have this property in general.
     !> (b) Compute the eigenvalues to suitable accuracy.
     !> If the eigenvectors are desired, the algorithm attains full
     !> accuracy of the computed eigenvalues only right before
     !> the corresponding vectors have to be computed, see steps c) and d).
     !> (c) For each cluster of close eigenvalues, select a new
     !> shift close to the cluster, find a new factorization, and refine
     !> the shifted eigenvalues to suitable accuracy.
     !> (d) For each eigenvalue with a large enough relative separation compute
     !> the corresponding eigenvector by forming a rank revealing twisted
     !> factorization. Go back to (c) for any clusters that remain.
     !> For more details, see:
     !> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
     !> to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
     !> Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     !> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
     !> Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
     !> 2004.  Also LAPACK Working Note 154.
     !> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem",
     !> Computer Science Division Technical Report No. UCB/CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Further Details
     !> 1.SSTEMR works only on machines which follow IEEE-754
     !> floating-point standard in their handling of infinities and NaNs.
     !> This permits the use of efficient inner loops avoiding a check for
     !> zero divisors.

     pure subroutine la_sstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc, &
               isuppz,tryrac,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           logical(lk),intent(inout) :: tryrac
           integer(ilp),intent(in) :: il,iu,ldz,nzc,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*)
           real(sp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: minrgp = 3.0e-3_sp

           ! Local Scalars
           logical(lk) :: alleig,indeig,lquery,valeig,wantz,zquery
           integer(ilp) :: i,ibegin,iend,ifirst,iil,iindbl,iindw,iindwk,iinfo,iinspl, &
           iiu,ilast,in,indd,inde2,inderr,indgp,indgrs,indwrk,itmp,itmp2,j,jblk,jj, &
                     liwmin,lwmin,nsplit,nzcmin,offset,wbegin,wend
           real(sp) :: bignum,cs,eps,pivmin,r1,r2,rmax,rmin,rtol1,rtol2,safmin,scale, &
                     smlnum,sn,thresh,tmp,tnrm,wl,wu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           zquery = (nzc == -1)
           ! la_sstemr needs work of size 6*n, iwork of size 3*n.
           ! in addition, la_slarre needs work of size 6*n, iwork of size 5*n.
           ! furthermore, la_slarrv needs work of size 12*n, iwork of size 7*n.
           if (wantz) then
              lwmin = 18*n
              liwmin = 10*n
           else
              ! need less workspace if only the eigenvalues are wanted
              lwmin = 12*n
              liwmin = 8*n
           end if
           wl = zero
           wu = zero
           iil = 0
           iiu = 0
           nsplit = 0
           if (valeig) then
              ! we do not reference vl, vu in the cases range = 'i','a'
              ! the interval (wl, wu] contains all the wanted eigenvalues.
              ! it is either given by the user or computed in la_slarre.
              wl = vl
              wu = vu
           elseif (indeig) then
              ! we do not reference il, iu in the cases range = 'v','a'
              iil = il
              iiu = iu
           end if
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (valeig .and. n > 0 .and. wu <= wl) then
              info = -7
           else if (indeig .and. (iil < 1 .or. iil > n)) then
              info = -8
           else if (indeig .and. (iiu < iil .or. iiu > n)) then
              info = -9
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -13
           else if (lwork < lwmin .and. .not. lquery) then
              info = -17
           else if (liwork < liwmin .and. .not. lquery) then
              info = -19
           end if
           ! get machine constants.
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (wantz .and. alleig) then
                 nzcmin = n
              else if (wantz .and. valeig) then
                 call la_slarrc('T',n,vl,vu,d,e,safmin,nzcmin,itmp,itmp2,info)

              else if (wantz .and. indeig) then
                 nzcmin = iiu - iil + 1
              else
                 ! wantz == false.
                 nzcmin = 0
              end if
              if (zquery .and. info == 0) then
                 z(1,1) = nzcmin
              else if (nzc < nzcmin .and. .not. zquery) then
                 info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SSTEMR',-info)
              return
           else if (lquery .or. zquery) then
              return
           end if
           ! handle n = 0, 1, and 2 cases immediately
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (wl < d(1) .and. wu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz .and. (.not. zquery)) then
                 z(1,1) = one
                 isuppz(1) = 1
                 isuppz(2) = 1
              end if
              return
           end if
           if (n == 2) then
              if (.not. wantz) then
                 call la_slae2(d(1),e(1),d(2),r1,r2)
              else if (wantz .and. (.not. zquery)) then
                 call la_slaev2(d(1),e(1),d(2),r1,r2,cs,sn)
              end if
              if (alleig .or. (valeig .and. (r2 > wl) .and. (r2 <= wu)) .or. (indeig .and. (iil == 1))) &
                        then
                 m = m + 1
                 w(m) = r2
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = -sn
                    z(2,m) = cs
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
              if (alleig .or. (valeig .and. (r1 > wl) .and. (r1 <= wu)) .or. (indeig .and. (iiu == 2))) &
                        then
                 m = m + 1
                 w(m) = r1
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = cs
                    z(2,m) = sn
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
           else
           ! continue with general n
              indgrs = 1
              inderr = 2*n + 1
              indgp = 3*n + 1
              indd = 4*n + 1
              inde2 = 5*n + 1
              indwrk = 6*n + 1
              iinspl = 1
              iindbl = n + 1
              iindw = 2*n + 1
              iindwk = 3*n + 1
              ! scale matrix to allowable range, if necessary.
              ! the allowable range is related to the pivmin parameter; see the
              ! comments in la_slarrd.  the preference for scaling small values
              ! up is heuristic; we expect users' matrices not to be close to the
              ! rmax threshold.
              scale = one
              tnrm = la_slanst('M',n,d,e)
              if (tnrm > zero .and. tnrm < rmin) then
                 scale = rmin/tnrm
              else if (tnrm > rmax) then
                 scale = rmax/tnrm
              end if
              if (scale /= one) then
                 call la_sscal(n,scale,d,1)
                 call la_sscal(n - 1,scale,e,1)
                 tnrm = tnrm*scale
                 if (valeig) then
                    ! if eigenvalues in interval have to be found,
                    ! scale (wl, wu] accordingly
                    wl = wl*scale
                    wu = wu*scale
                 end if
              end if
              ! compute the desired eigenvalues of the tridiagonal after splitting
              ! into smaller subblocks if the corresponding off-diagonal elements
              ! are small
              ! thresh is the splitting parameter for la_slarre
              ! a negative thresh forces the old splitting criterion based on the
              ! size of the off-diagonal. a positive thresh switches to splitting
              ! which preserves relative accuracy.
              if (tryrac) then
                 ! test whether the matrix warrants the more expensive relative approach.
                 call la_slarrr(n,d,e,iinfo)
              else
                 ! the user does not care about relative accurately eigenvalues
                 iinfo = -1
              end if
              ! set the splitting criterion
              if (iinfo == 0) then
                 thresh = eps
              else
                 thresh = -eps
                 ! relative accuracy is desired but t does not guarantee it
                 tryrac = .false.
              end if
              if (tryrac) then
                 ! copy original diagonal, needed to guarantee relative accuracy
                 call la_scopy(n,d,1,work(indd),1)
              end if
              ! store the squares of the offdiagonal values of t
              do j = 1,n - 1
                 work(inde2 + j - 1) = e(j)**2
              end do
              ! set the tolerance parameters for bisection
              if (.not. wantz) then
                 ! la_slarre computes the eigenvalues to full precision.
                 rtol1 = four*eps
                 rtol2 = four*eps
              else
                 ! la_slarre computes the eigenvalues to less than full precision.
                 ! la_slarrv will refine the eigenvalue approximations, and we can
                 ! need less accurate initial bisection in la_slarre.
                 ! note: these settings do only affect the subset case and la_slarre
                 rtol1 = max(sqrt(eps)*5.0e-2_sp,four*eps)
                 rtol2 = max(sqrt(eps)*5.0e-3_sp,four*eps)
              end if
              call la_slarre(range,n,wl,wu,iil,iiu,d,e,work(inde2),rtol1,rtol2, &
              thresh,nsplit,iwork(iinspl),m,w,work(inderr),work(indgp),iwork(iindbl), &
              iwork(iindw),work(indgrs),pivmin,work(indwrk),iwork(iindwk),iinfo)

              if (iinfo /= 0) then
                 info = 10 + abs(iinfo)
                 return
              end if
              ! note that if range /= 'v', la_slarre computes bounds on the desired
              ! part of the spectrum. all desired eigenvalues are contained in
              ! (wl,wu]
              if (wantz) then
                 ! compute the desired eigenvectors corresponding to the computed
                 ! eigenvalues
                 call la_slarrv(n,wl,wu,d,e,pivmin,iwork(iinspl),m,1,m,minrgp, &
                 rtol1,rtol2,w,work(inderr),work(indgp),iwork(iindbl),iwork(iindw), &
                           work(indgrs),z,ldz,isuppz,work(indwrk),iwork(iindwk),iinfo)
                 if (iinfo /= 0) then
                    info = 20 + abs(iinfo)
                    return
                 end if
              else
                 ! la_slarre computes eigenvalues of the (shifted) root representation
                 ! la_slarrv returns the eigenvalues of the unshifted matrix.
                 ! however, if the eigenvectors are not desired by the user, we need
                 ! to apply the corresponding shifts from la_slarre to obtain the
                 ! eigenvalues of the original matrix.
                 do j = 1,m
                    itmp = iwork(iindbl + j - 1)
                    w(j) = w(j) + e(iwork(iinspl + itmp - 1))
                 end do
              end if
              if (tryrac) then
                 ! refine computed eigenvalues so that they are relatively accurate
                 ! with respect to the original matrix t.
                 ibegin = 1
                 wbegin = 1
                 loop_39: do jblk = 1,iwork(iindbl + m - 1)
                    iend = iwork(iinspl + jblk - 1)
                    in = iend - ibegin + 1
                    wend = wbegin - 1
                    ! check if any eigenvalues have to be refined in this block
                    36 continue
                    if (wend < m) then
                       if (iwork(iindbl + wend) == jblk) then
                          wend = wend + 1
                          go to 36
                       end if
                    end if
                    if (wend < wbegin) then
                       ibegin = iend + 1
                       cycle loop_39
                    end if
                    offset = iwork(iindw + wbegin - 1) - 1
                    ifirst = iwork(iindw + wbegin - 1)
                    ilast = iwork(iindw + wend - 1)
                    rtol2 = four*eps
                    call la_slarrj(in,work(indd + ibegin - 1),work(inde2 + ibegin - 1),ifirst, &
                    ilast,rtol2,offset,w(wbegin),work(inderr + wbegin - 1),work(indwrk),iwork( &
                               iindwk),pivmin,tnrm,iinfo)
                    ibegin = iend + 1
                    wbegin = wend + 1
                 end do loop_39
              end if
              ! if matrix was scaled, then rescale eigenvalues appropriately.
              if (scale /= one) then
                 call la_sscal(m,one/scale,w,1)
              end if
           end if
           ! if eigenvalues are not in increasing order, then sort them,
           ! possibly along with eigenvectors.
           if (nsplit > 1 .or. n == 2) then
              if (.not. wantz) then
                 call la_slasrt('I',m,w,iinfo)
                 if (iinfo /= 0) then
                    info = 3
                    return
                 end if
              else
                 do j = 1,m - 1
                    i = 0
                    tmp = w(j)
                    do jj = j + 1,m
                       if (w(jj) < tmp) then
                          i = jj
                          tmp = w(jj)
                       end if
                    end do
                    if (i /= 0) then
                       w(i) = w(j)
                       w(j) = tmp
                       if (wantz) then
                          call la_sswap(n,z(1,i),1,z(1,j),1)
                          itmp = isuppz(2*i - 1)
                          isuppz(2*i - 1) = isuppz(2*j - 1)
                          isuppz(2*j - 1) = itmp
                          itmp = isuppz(2*i)
                          isuppz(2*i) = isuppz(2*j)
                          isuppz(2*j) = itmp
                       end if
                    end if
                 end do
              end if
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_sstemr
     !> DSTEMR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> Depending on the number of desired eigenvalues, these are computed either
     !> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
     !> computed by the use of various suitable L D L^T factorizations near clusters
     !> of close eigenvalues (referred to as RRRs, Relatively Robust
     !> Representations). An informal sketch of the algorithm follows.
     !> For each unreduced block (submatrix) of T,
     !> (a) Compute T - sigma I  = L D L^T, so that L and D
     !> define all the wanted eigenvalues to high relative accuracy.
     !> This means that small relative changes in the entries of D and L
     !> cause only small relative changes in the eigenvalues and
     !> eigenvectors. The standard (unfactored) representation of the
     !> tridiagonal matrix T does not have this property in general.
     !> (b) Compute the eigenvalues to suitable accuracy.
     !> If the eigenvectors are desired, the algorithm attains full
     !> accuracy of the computed eigenvalues only right before
     !> the corresponding vectors have to be computed, see steps c) and d).
     !> (c) For each cluster of close eigenvalues, select a new
     !> shift close to the cluster, find a new factorization, and refine
     !> the shifted eigenvalues to suitable accuracy.
     !> (d) For each eigenvalue with a large enough relative separation compute
     !> the corresponding eigenvector by forming a rank revealing twisted
     !> factorization. Go back to (c) for any clusters that remain.
     !> For more details, see:
     !> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
     !> to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
     !> Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     !> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
     !> Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
     !> 2004.  Also LAPACK Working Note 154.
     !> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem",
     !> Computer Science Division Technical Report No. UCB/CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Further Details
     !> 1.DSTEMR works only on machines which follow IEEE-754
     !> floating-point standard in their handling of infinities and NaNs.
     !> This permits the use of efficient inner loops avoiding a check for
     !> zero divisors.

     pure subroutine la_dstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc, &
               isuppz,tryrac,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           logical(lk),intent(inout) :: tryrac
           integer(ilp),intent(in) :: il,iu,ldz,nzc,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*)
           real(dp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: minrgp = 1.0e-3_dp

           ! Local Scalars
           logical(lk) :: alleig,indeig,lquery,valeig,wantz,zquery
           integer(ilp) :: i,ibegin,iend,ifirst,iil,iindbl,iindw,iindwk,iinfo,iinspl, &
           iiu,ilast,in,indd,inde2,inderr,indgp,indgrs,indwrk,itmp,itmp2,j,jblk,jj, &
                     liwmin,lwmin,nsplit,nzcmin,offset,wbegin,wend
           real(dp) :: bignum,cs,eps,pivmin,r1,r2,rmax,rmin,rtol1,rtol2,safmin,scale, &
                     smlnum,sn,thresh,tmp,tnrm,wl,wu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           zquery = (nzc == -1)
           ! la_dstemr needs work of size 6*n, iwork of size 3*n.
           ! in addition, la_dlarre needs work of size 6*n, iwork of size 5*n.
           ! furthermore, la_dlarrv needs work of size 12*n, iwork of size 7*n.
           if (wantz) then
              lwmin = 18*n
              liwmin = 10*n
           else
              ! need less workspace if only the eigenvalues are wanted
              lwmin = 12*n
              liwmin = 8*n
           end if
           wl = zero
           wu = zero
           iil = 0
           iiu = 0
           nsplit = 0
           if (valeig) then
              ! we do not reference vl, vu in the cases range = 'i','a'
              ! the interval (wl, wu] contains all the wanted eigenvalues.
              ! it is either given by the user or computed in la_dlarre.
              wl = vl
              wu = vu
           elseif (indeig) then
              ! we do not reference il, iu in the cases range = 'v','a'
              iil = il
              iiu = iu
           end if
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (valeig .and. n > 0 .and. wu <= wl) then
              info = -7
           else if (indeig .and. (iil < 1 .or. iil > n)) then
              info = -8
           else if (indeig .and. (iiu < iil .or. iiu > n)) then
              info = -9
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -13
           else if (lwork < lwmin .and. .not. lquery) then
              info = -17
           else if (liwork < liwmin .and. .not. lquery) then
              info = -19
           end if
           ! get machine constants.
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (wantz .and. alleig) then
                 nzcmin = n
              else if (wantz .and. valeig) then
                 call la_dlarrc('T',n,vl,vu,d,e,safmin,nzcmin,itmp,itmp2,info)

              else if (wantz .and. indeig) then
                 nzcmin = iiu - iil + 1
              else
                 ! wantz == false.
                 nzcmin = 0
              end if
              if (zquery .and. info == 0) then
                 z(1,1) = nzcmin
              else if (nzc < nzcmin .and. .not. zquery) then
                 info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DSTEMR',-info)
              return
           else if (lquery .or. zquery) then
              return
           end if
           ! handle n = 0, 1, and 2 cases immediately
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (wl < d(1) .and. wu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz .and. (.not. zquery)) then
                 z(1,1) = one
                 isuppz(1) = 1
                 isuppz(2) = 1
              end if
              return
           end if
           if (n == 2) then
              if (.not. wantz) then
                 call la_dlae2(d(1),e(1),d(2),r1,r2)
              else if (wantz .and. (.not. zquery)) then
                 call la_dlaev2(d(1),e(1),d(2),r1,r2,cs,sn)
              end if
              if (alleig .or. (valeig .and. (r2 > wl) .and. (r2 <= wu)) .or. (indeig .and. (iil == 1))) &
                        then
                 m = m + 1
                 w(m) = r2
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = -sn
                    z(2,m) = cs
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
              if (alleig .or. (valeig .and. (r1 > wl) .and. (r1 <= wu)) .or. (indeig .and. (iiu == 2))) &
                        then
                 m = m + 1
                 w(m) = r1
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = cs
                    z(2,m) = sn
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
           else
           ! continue with general n
              indgrs = 1
              inderr = 2*n + 1
              indgp = 3*n + 1
              indd = 4*n + 1
              inde2 = 5*n + 1
              indwrk = 6*n + 1
              iinspl = 1
              iindbl = n + 1
              iindw = 2*n + 1
              iindwk = 3*n + 1
              ! scale matrix to allowable range, if necessary.
              ! the allowable range is related to the pivmin parameter; see the
              ! comments in la_dlarrd.  the preference for scaling small values
              ! up is heuristic; we expect users' matrices not to be close to the
              ! rmax threshold.
              scale = one
              tnrm = la_dlanst('M',n,d,e)
              if (tnrm > zero .and. tnrm < rmin) then
                 scale = rmin/tnrm
              else if (tnrm > rmax) then
                 scale = rmax/tnrm
              end if
              if (scale /= one) then
                 call la_dscal(n,scale,d,1)
                 call la_dscal(n - 1,scale,e,1)
                 tnrm = tnrm*scale
                 if (valeig) then
                    ! if eigenvalues in interval have to be found,
                    ! scale (wl, wu] accordingly
                    wl = wl*scale
                    wu = wu*scale
                 end if
              end if
              ! compute the desired eigenvalues of the tridiagonal after splitting
              ! into smaller subblocks if the corresponding off-diagonal elements
              ! are small
              ! thresh is the splitting parameter for la_dlarre
              ! a negative thresh forces the old splitting criterion based on the
              ! size of the off-diagonal. a positive thresh switches to splitting
              ! which preserves relative accuracy.
              if (tryrac) then
                 ! test whether the matrix warrants the more expensive relative approach.
                 call la_dlarrr(n,d,e,iinfo)
              else
                 ! the user does not care about relative accurately eigenvalues
                 iinfo = -1
              end if
              ! set the splitting criterion
              if (iinfo == 0) then
                 thresh = eps
              else
                 thresh = -eps
                 ! relative accuracy is desired but t does not guarantee it
                 tryrac = .false.
              end if
              if (tryrac) then
                 ! copy original diagonal, needed to guarantee relative accuracy
                 call la_dcopy(n,d,1,work(indd),1)
              end if
              ! store the squares of the offdiagonal values of t
              do j = 1,n - 1
                 work(inde2 + j - 1) = e(j)**2
              end do
              ! set the tolerance parameters for bisection
              if (.not. wantz) then
                 ! la_dlarre computes the eigenvalues to full precision.
                 rtol1 = four*eps
                 rtol2 = four*eps
              else
                 ! la_dlarre computes the eigenvalues to less than full precision.
                 ! la_dlarrv will refine the eigenvalue approximations, and we can
                 ! need less accurate initial bisection in la_dlarre.
                 ! note: these settings do only affect the subset case and la_dlarre
                 rtol1 = sqrt(eps)
                 rtol2 = max(sqrt(eps)*5.0e-3_dp,four*eps)
              end if
              call la_dlarre(range,n,wl,wu,iil,iiu,d,e,work(inde2),rtol1,rtol2, &
              thresh,nsplit,iwork(iinspl),m,w,work(inderr),work(indgp),iwork(iindbl), &
              iwork(iindw),work(indgrs),pivmin,work(indwrk),iwork(iindwk),iinfo)

              if (iinfo /= 0) then
                 info = 10 + abs(iinfo)
                 return
              end if
              ! note that if range /= 'v', la_dlarre computes bounds on the desired
              ! part of the spectrum. all desired eigenvalues are contained in
              ! (wl,wu]
              if (wantz) then
                 ! compute the desired eigenvectors corresponding to the computed
                 ! eigenvalues
                 call la_dlarrv(n,wl,wu,d,e,pivmin,iwork(iinspl),m,1,m,minrgp, &
                 rtol1,rtol2,w,work(inderr),work(indgp),iwork(iindbl),iwork(iindw), &
                           work(indgrs),z,ldz,isuppz,work(indwrk),iwork(iindwk),iinfo)
                 if (iinfo /= 0) then
                    info = 20 + abs(iinfo)
                    return
                 end if
              else
                 ! la_dlarre computes eigenvalues of the (shifted) root representation
                 ! la_dlarrv returns the eigenvalues of the unshifted matrix.
                 ! however, if the eigenvectors are not desired by the user, we need
                 ! to apply the corresponding shifts from la_dlarre to obtain the
                 ! eigenvalues of the original matrix.
                 do j = 1,m
                    itmp = iwork(iindbl + j - 1)
                    w(j) = w(j) + e(iwork(iinspl + itmp - 1))
                 end do
              end if
              if (tryrac) then
                 ! refine computed eigenvalues so that they are relatively accurate
                 ! with respect to the original matrix t.
                 ibegin = 1
                 wbegin = 1
                 loop_39: do jblk = 1,iwork(iindbl + m - 1)
                    iend = iwork(iinspl + jblk - 1)
                    in = iend - ibegin + 1
                    wend = wbegin - 1
                    ! check if any eigenvalues have to be refined in this block
                    36 continue
                    if (wend < m) then
                       if (iwork(iindbl + wend) == jblk) then
                          wend = wend + 1
                          go to 36
                       end if
                    end if
                    if (wend < wbegin) then
                       ibegin = iend + 1
                       cycle loop_39
                    end if
                    offset = iwork(iindw + wbegin - 1) - 1
                    ifirst = iwork(iindw + wbegin - 1)
                    ilast = iwork(iindw + wend - 1)
                    rtol2 = four*eps
                    call la_dlarrj(in,work(indd + ibegin - 1),work(inde2 + ibegin - 1),ifirst, &
                    ilast,rtol2,offset,w(wbegin),work(inderr + wbegin - 1),work(indwrk),iwork( &
                               iindwk),pivmin,tnrm,iinfo)
                    ibegin = iend + 1
                    wbegin = wend + 1
                 end do loop_39
              end if
              ! if matrix was scaled, then rescale eigenvalues appropriately.
              if (scale /= one) then
                 call la_dscal(m,one/scale,w,1)
              end if
           end if
           ! if eigenvalues are not in increasing order, then sort them,
           ! possibly along with eigenvectors.
           if (nsplit > 1 .or. n == 2) then
              if (.not. wantz) then
                 call la_dlasrt('I',m,w,iinfo)
                 if (iinfo /= 0) then
                    info = 3
                    return
                 end if
              else
                 do j = 1,m - 1
                    i = 0
                    tmp = w(j)
                    do jj = j + 1,m
                       if (w(jj) < tmp) then
                          i = jj
                          tmp = w(jj)
                       end if
                    end do
                    if (i /= 0) then
                       w(i) = w(j)
                       w(j) = tmp
                       if (wantz) then
                          call la_dswap(n,z(1,i),1,z(1,j),1)
                          itmp = isuppz(2*i - 1)
                          isuppz(2*i - 1) = isuppz(2*j - 1)
                          isuppz(2*j - 1) = itmp
                          itmp = isuppz(2*i)
                          isuppz(2*i) = isuppz(2*j)
                          isuppz(2*j) = itmp
                       end if
                    end if
                 end do
              end if
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_dstemr
     !> QSTEMR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> Depending on the number of desired eigenvalues, these are computed either
     !> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
     !> computed by the use of various suitable L D L^T factorizations near clusters
     !> of close eigenvalues (referred to as RRRs, Relatively Robust
     !> Representations). An informal sketch of the algorithm follows.
     !> For each unreduced block (submatrix) of T,
     !> (a) Compute T - sigma I  = L D L^T, so that L and D
     !> define all the wanted eigenvalues to high relative accuracy.
     !> This means that small relative changes in the entries of D and L
     !> cause only small relative changes in the eigenvalues and
     !> eigenvectors. The standard (unfactored) representation of the
     !> tridiagonal matrix T does not have this property in general.
     !> (b) Compute the eigenvalues to suitable accuracy.
     !> If the eigenvectors are desired, the algorithm attains full
     !> accuracy of the computed eigenvalues only right before
     !> the corresponding vectors have to be computed, see steps c) and d).
     !> (c) For each cluster of close eigenvalues, select a new
     !> shift close to the cluster, find a new factorization, and refine
     !> the shifted eigenvalues to suitable accuracy.
     !> (d) For each eigenvalue with a large enough relative separation compute
     !> the corresponding eigenvector by forming a rank revealing twisted
     !> factorization. Go back to (c) for any clusters that remain.
     !> For more details, see:
     !> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
     !> to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
     !> Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     !> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
     !> Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
     !> 2004.  Also LAPACK Working Note 154.
     !> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem",
     !> Computer Science Division Technical Report No. UCB/CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Further Details
     !> 1.QSTEMR works only on machines which follow IEEE-754
     !> floating-point standard in their handling of infinities and NaNs.
     !> This permits the use of efficient inner loops avoiding a check for
     !> zero divisors.

     pure subroutine la_qstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc, &
               isuppz,tryrac,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           logical(lk),intent(inout) :: tryrac
           integer(ilp),intent(in) :: il,iu,ldz,nzc,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*)
           real(qp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: minrgp = 1.0e-3_qp

           ! Local Scalars
           logical(lk) :: alleig,indeig,lquery,valeig,wantz,zquery
           integer(ilp) :: i,ibegin,iend,ifirst,iil,iindbl,iindw,iindwk,iinfo,iinspl, &
           iiu,ilast,in,indd,inde2,inderr,indgp,indgrs,indwrk,itmp,itmp2,j,jblk,jj, &
                     liwmin,lwmin,nsplit,nzcmin,offset,wbegin,wend
           real(qp) :: bignum,cs,eps,pivmin,r1,r2,rmax,rmin,rtol1,rtol2,safmin,scale, &
                     smlnum,sn,thresh,tmp,tnrm,wl,wu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           zquery = (nzc == -1)
           ! la_qstemr needs work of size 6*n, iwork of size 3*n.
           ! in addition, la_qlarre needs work of size 6*n, iwork of size 5*n.
           ! furthermore, la_qlarrv needs work of size 12*n, iwork of size 7*n.
           if (wantz) then
              lwmin = 18*n
              liwmin = 10*n
           else
              ! need less workspace if only the eigenvalues are wanted
              lwmin = 12*n
              liwmin = 8*n
           end if
           wl = zero
           wu = zero
           iil = 0
           iiu = 0
           nsplit = 0
           if (valeig) then
              ! we do not reference vl, vu in the cases range = 'i','a'
              ! the interval (wl, wu] contains all the wanted eigenvalues.
              ! it is either given by the user or computed in la_qlarre.
              wl = vl
              wu = vu
           elseif (indeig) then
              ! we do not reference il, iu in the cases range = 'v','a'
              iil = il
              iiu = iu
           end if
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (valeig .and. n > 0 .and. wu <= wl) then
              info = -7
           else if (indeig .and. (iil < 1 .or. iil > n)) then
              info = -8
           else if (indeig .and. (iiu < iil .or. iiu > n)) then
              info = -9
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -13
           else if (lwork < lwmin .and. .not. lquery) then
              info = -17
           else if (liwork < liwmin .and. .not. lquery) then
              info = -19
           end if
           ! get machine constants.
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (wantz .and. alleig) then
                 nzcmin = n
              else if (wantz .and. valeig) then
                 call la_qlarrc('T',n,vl,vu,d,e,safmin,nzcmin,itmp,itmp2,info)

              else if (wantz .and. indeig) then
                 nzcmin = iiu - iil + 1
              else
                 ! wantz == false.
                 nzcmin = 0
              end if
              if (zquery .and. info == 0) then
                 z(1,1) = nzcmin
              else if (nzc < nzcmin .and. .not. zquery) then
                 info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QSTEMR',-info)
              return
           else if (lquery .or. zquery) then
              return
           end if
           ! handle n = 0, 1, and 2 cases immediately
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (wl < d(1) .and. wu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz .and. (.not. zquery)) then
                 z(1,1) = one
                 isuppz(1) = 1
                 isuppz(2) = 1
              end if
              return
           end if
           if (n == 2) then
              if (.not. wantz) then
                 call la_qlae2(d(1),e(1),d(2),r1,r2)
              else if (wantz .and. (.not. zquery)) then
                 call la_qlaev2(d(1),e(1),d(2),r1,r2,cs,sn)
              end if
              if (alleig .or. (valeig .and. (r2 > wl) .and. (r2 <= wu)) .or. (indeig .and. (iil == 1))) &
                        then
                 m = m + 1
                 w(m) = r2
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = -sn
                    z(2,m) = cs
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
              if (alleig .or. (valeig .and. (r1 > wl) .and. (r1 <= wu)) .or. (indeig .and. (iiu == 2))) &
                        then
                 m = m + 1
                 w(m) = r1
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = cs
                    z(2,m) = sn
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
           else
           ! continue with general n
              indgrs = 1
              inderr = 2*n + 1
              indgp = 3*n + 1
              indd = 4*n + 1
              inde2 = 5*n + 1
              indwrk = 6*n + 1
              iinspl = 1
              iindbl = n + 1
              iindw = 2*n + 1
              iindwk = 3*n + 1
              ! scale matrix to allowable range, if necessary.
              ! the allowable range is related to the pivmin parameter; see the
              ! comments in la_qlarrd.  the preference for scaling small values
              ! up is heuristic; we expect users' matrices not to be close to the
              ! rmax threshold.
              scale = one
              tnrm = la_qlanst('M',n,d,e)
              if (tnrm > zero .and. tnrm < rmin) then
                 scale = rmin/tnrm
              else if (tnrm > rmax) then
                 scale = rmax/tnrm
              end if
              if (scale /= one) then
                 call la_qscal(n,scale,d,1)
                 call la_qscal(n - 1,scale,e,1)
                 tnrm = tnrm*scale
                 if (valeig) then
                    ! if eigenvalues in interval have to be found,
                    ! scale (wl, wu] accordingly
                    wl = wl*scale
                    wu = wu*scale
                 end if
              end if
              ! compute the desired eigenvalues of the tridiagonal after splitting
              ! into smaller subblocks if the corresponding off-diagonal elements
              ! are small
              ! thresh is the splitting parameter for la_qlarre
              ! a negative thresh forces the old splitting criterion based on the
              ! size of the off-diagonal. a positive thresh switches to splitting
              ! which preserves relative accuracy.
              if (tryrac) then
                 ! test whether the matrix warrants the more expensive relative approach.
                 call la_qlarrr(n,d,e,iinfo)
              else
                 ! the user does not care about relative accurately eigenvalues
                 iinfo = -1
              end if
              ! set the splitting criterion
              if (iinfo == 0) then
                 thresh = eps
              else
                 thresh = -eps
                 ! relative accuracy is desired but t does not guarantee it
                 tryrac = .false.
              end if
              if (tryrac) then
                 ! copy original diagonal, needed to guarantee relative accuracy
                 call la_qcopy(n,d,1,work(indd),1)
              end if
              ! store the squares of the offdiagonal values of t
              do j = 1,n - 1
                 work(inde2 + j - 1) = e(j)**2
              end do
              ! set the tolerance parameters for bisection
              if (.not. wantz) then
                 ! la_qlarre computes the eigenvalues to full precision.
                 rtol1 = four*eps
                 rtol2 = four*eps
              else
                 ! la_qlarre computes the eigenvalues to less than full precision.
                 ! la_qlarrv will refine the eigenvalue approximations, and we can
                 ! need less accurate initial bisection in la_qlarre.
                 ! note: these settings do only affect the subset case and la_qlarre
                 rtol1 = sqrt(eps)
                 rtol2 = max(sqrt(eps)*5.0e-3_qp,four*eps)
              end if
              call la_qlarre(range,n,wl,wu,iil,iiu,d,e,work(inde2),rtol1,rtol2, &
              thresh,nsplit,iwork(iinspl),m,w,work(inderr),work(indgp),iwork(iindbl), &
              iwork(iindw),work(indgrs),pivmin,work(indwrk),iwork(iindwk),iinfo)

              if (iinfo /= 0) then
                 info = 10 + abs(iinfo)
                 return
              end if
              ! note that if range /= 'v', la_qlarre computes bounds on the desired
              ! part of the spectrum. all desired eigenvalues are contained in
              ! (wl,wu]
              if (wantz) then
                 ! compute the desired eigenvectors corresponding to the computed
                 ! eigenvalues
                 call la_qlarrv(n,wl,wu,d,e,pivmin,iwork(iinspl),m,1,m,minrgp, &
                 rtol1,rtol2,w,work(inderr),work(indgp),iwork(iindbl),iwork(iindw), &
                           work(indgrs),z,ldz,isuppz,work(indwrk),iwork(iindwk),iinfo)
                 if (iinfo /= 0) then
                    info = 20 + abs(iinfo)
                    return
                 end if
              else
                 ! la_qlarre computes eigenvalues of the (shifted) root representation
                 ! la_qlarrv returns the eigenvalues of the unshifted matrix.
                 ! however, if the eigenvectors are not desired by the user, we need
                 ! to apply the corresponding shifts from la_qlarre to obtain the
                 ! eigenvalues of the original matrix.
                 do j = 1,m
                    itmp = iwork(iindbl + j - 1)
                    w(j) = w(j) + e(iwork(iinspl + itmp - 1))
                 end do
              end if
              if (tryrac) then
                 ! refine computed eigenvalues so that they are relatively accurate
                 ! with respect to the original matrix t.
                 ibegin = 1
                 wbegin = 1
                 loop_39: do jblk = 1,iwork(iindbl + m - 1)
                    iend = iwork(iinspl + jblk - 1)
                    in = iend - ibegin + 1
                    wend = wbegin - 1
                    ! check if any eigenvalues have to be refined in this block
                    36 continue
                    if (wend < m) then
                       if (iwork(iindbl + wend) == jblk) then
                          wend = wend + 1
                          go to 36
                       end if
                    end if
                    if (wend < wbegin) then
                       ibegin = iend + 1
                       cycle loop_39
                    end if
                    offset = iwork(iindw + wbegin - 1) - 1
                    ifirst = iwork(iindw + wbegin - 1)
                    ilast = iwork(iindw + wend - 1)
                    rtol2 = four*eps
                    call la_qlarrj(in,work(indd + ibegin - 1),work(inde2 + ibegin - 1),ifirst, &
                    ilast,rtol2,offset,w(wbegin),work(inderr + wbegin - 1),work(indwrk),iwork( &
                               iindwk),pivmin,tnrm,iinfo)
                    ibegin = iend + 1
                    wbegin = wend + 1
                 end do loop_39
              end if
              ! if matrix was scaled, then rescale eigenvalues appropriately.
              if (scale /= one) then
                 call la_qscal(m,one/scale,w,1)
              end if
           end if
           ! if eigenvalues are not in increasing order, then sort them,
           ! possibly along with eigenvectors.
           if (nsplit > 1 .or. n == 2) then
              if (.not. wantz) then
                 call la_qlasrt('I',m,w,iinfo)
                 if (iinfo /= 0) then
                    info = 3
                    return
                 end if
              else
                 do j = 1,m - 1
                    i = 0
                    tmp = w(j)
                    do jj = j + 1,m
                       if (w(jj) < tmp) then
                          i = jj
                          tmp = w(jj)
                       end if
                    end do
                    if (i /= 0) then
                       w(i) = w(j)
                       w(j) = tmp
                       if (wantz) then
                          call la_qswap(n,z(1,i),1,z(1,j),1)
                          itmp = isuppz(2*i - 1)
                          isuppz(2*i - 1) = isuppz(2*j - 1)
                          isuppz(2*j - 1) = itmp
                          itmp = isuppz(2*i)
                          isuppz(2*i) = isuppz(2*j)
                          isuppz(2*j) = itmp
                       end if
                    end if
                 end do
              end if
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_qstemr

     !> SSTEVR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T.  Eigenvalues and
     !> eigenvectors can be selected by specifying either a range of values
     !> or a range of indices for the desired eigenvalues.
     !> Whenever possible, SSTEVR calls SSTEMR to compute the
     !> eigenspectrum using Relatively Robust Representations.  SSTEMR
     !> computes eigenvalues by the dqds algorithm, while orthogonal
     !> eigenvectors are computed from various "good" L D L^T representations
     !> (also known as Relatively Robust Representations). Gram-Schmidt
     !> orthogonalization is avoided as far as possible. More specifically,
     !> the various steps of the algorithm are as follows. For the i-th
     !> unreduced block of T,
     !> (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
     !> is a relatively robust representation,
     !> (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
     !> relative accuracy by the dqds algorithm,
     !> (c) If there is a cluster of close eigenvalues, "choose" sigma_i
     !> close to the cluster, and go to step (a),
     !> (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
     !> compute the corresponding eigenvector by forming a
     !> rank-revealing twisted factorization.
     !> The desired accuracy of the output can be specified by the input
     !> parameter ABSTOL.
     !> For more details, see "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
     !> Computer Science Division Technical Report No. UCB//CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Note 1 : SSTEVR calls SSTEMR when the full spectrum is requested
     !> on machines which conform to the ieee-754 floating point standard.
     !> SSTEVR calls SSTEBZ and SSTEIN on non-ieee machines and
     !> when partial spectrum requests are made.
     !> Normal execution of SSTEMR may create NaNs and infinities and
     !> hence may abort due to a floating point exception in environments
     !> which do not handle NaNs and infinities in the ieee standard default
     !> manner.

     pure subroutine la_sstevr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: alleig,indeig,test,lquery,valeig,wantz,tryrac
           character :: order
           integer(ilp) :: i,ieeeok,imax,indibl,indifl,indisp,indiwo,iscale,j,jj,liwmin, &
                      lwmin,nsplit
           real(sp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tmp1,tnrm,vll, &
                     vuu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           ieeeok = la_ilaenv(10,'SSTEVR','N',1,2,3,4)
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           lwmin = max(1,20*n)
           liwmin = max(1,10*n)
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else
              if (valeig) then
                 if (n > 0 .and. vu <= vl) info = -7
              else if (indeig) then
                 if (il < 1 .or. il > max(1,n)) then
                    info = -8
                 else if (iu < min(n,il) .or. iu > n) then
                    info = -9
                 end if
              end if
           end if
           if (info == 0) then
              if (ldz < 1 .or. (wantz .and. ldz < n)) then
                 info = -14
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -17
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -19
              end if
           end if
           if (info /= 0) then
              call la_xerbla('SSTEVR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (vl < d(1) .and. vu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           if (valeig) then
              vll = vl
              vuu = vu
           end if
           tnrm = la_slanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_sscal(n,sigma,d,1)
              call la_sscal(n - 1,sigma,e(1),1)
              if (valeig) then
                 vll = vl*sigma
                 vuu = vu*sigma
              end if
           end if
           ! initialize indices into workspaces.  note: these indices are used only
           ! if la_ssterf or la_sstemr fail.
           ! iwork(indibl:indibl+m-1) corresponds to iblock in la_sstebz and
           ! stores the block indices of each of the m<=n eigenvalues.
           indibl = 1
           ! iwork(indisp:indisp+nsplit-1) corresponds to isplit in la_sstebz and
           ! stores the starting and finishing indices of each block.
           indisp = indibl + n
           ! iwork(indifl:indifl+n-1) stores the indices of eigenvectors
           ! that corresponding to eigenvectors that fail to converge in
           ! la_sstein.  this information is discarded; if any fail, the driver
           ! returns info > 0.
           indifl = indisp + n
           ! indiwo is the offset of the remaining integer workspace.
           indiwo = indisp + n
           ! if all eigenvalues are desired, then
           ! call la_ssterf or la_sstemr.  if this fails for some eigenvalue, then
           ! try la_sstebz.
           test = .false.
           if (indeig) then
              if (il == 1 .and. iu == n) then
                 test = .true.
              end if
           end if
           if ((alleig .or. test) .and. ieeeok == 1) then
              call la_scopy(n - 1,e(1),1,work(1),1)
              if (.not. wantz) then
                 call la_scopy(n,d,1,w,1)
                 call la_ssterf(n,w,work,info)
              else
                 call la_scopy(n,d,1,work(n + 1),1)
                 if (abstol <= two*n*eps) then
                    tryrac = .true.
                 else
                    tryrac = .false.
                 end if
                 call la_sstemr(jobz,'A',n,work(n + 1),work,vl,vu,il,iu,m,w,z,ldz, &
                            n,isuppz,tryrac,work(2*n + 1),lwork - 2*n,iwork,liwork,info)
              end if
              if (info == 0) then
                 m = n
                 go to 10
              end if
              info = 0
           end if
           ! otherwise, call la_sstebz and, if eigenvectors are desired, la_sstein.
           if (wantz) then
              order = 'B'
           else
              order = 'E'
           end if
           call la_sstebz(range,order,n,vll,vuu,il,iu,abstol,d,e,m,nsplit,w, &
                     iwork(indibl),iwork(indisp),work,iwork(indiwo),info)
           if (wantz) then
              call la_sstein(n,d,e,m,w,iwork(indibl),iwork(indisp),z,ldz,work, &
                        iwork(indiwo),iwork(indifl),info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           10 continue
           if (iscale == 1) then
              if (info == 0) then
                 imax = m
              else
                 imax = info - 1
              end if
              call la_sscal(imax,one/sigma,w,1)
           end if
           ! if eigenvalues are not in order, then sort them, along with
           ! eigenvectors.
           if (wantz) then
              do j = 1,m - 1
                 i = 0
                 tmp1 = w(j)
                 do jj = j + 1,m
                    if (w(jj) < tmp1) then
                       i = jj
                       tmp1 = w(jj)
                    end if
                 end do
                 if (i /= 0) then
                    w(i) = w(j)
                    w(j) = tmp1
                    call la_sswap(n,z(1,i),1,z(1,j),1)
                 end if
              end do
           end if
            ! causes problems with tests 19
            ! if (wantz .and. indeig ) z( 1,1) = z(1,1) / 1.002_sp + .002
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_sstevr
     !> DSTEVR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T.  Eigenvalues and
     !> eigenvectors can be selected by specifying either a range of values
     !> or a range of indices for the desired eigenvalues.
     !> Whenever possible, DSTEVR calls DSTEMR to compute the
     !> eigenspectrum using Relatively Robust Representations.  DSTEMR
     !> computes eigenvalues by the dqds algorithm, while orthogonal
     !> eigenvectors are computed from various "good" L D L^T representations
     !> (also known as Relatively Robust Representations). Gram-Schmidt
     !> orthogonalization is avoided as far as possible. More specifically,
     !> the various steps of the algorithm are as follows. For the i-th
     !> unreduced block of T,
     !> (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
     !> is a relatively robust representation,
     !> (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
     !> relative accuracy by the dqds algorithm,
     !> (c) If there is a cluster of close eigenvalues, "choose" sigma_i
     !> close to the cluster, and go to step (a),
     !> (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
     !> compute the corresponding eigenvector by forming a
     !> rank-revealing twisted factorization.
     !> The desired accuracy of the output can be specified by the input
     !> parameter ABSTOL.
     !> For more details, see "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
     !> Computer Science Division Technical Report No. UCB//CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested
     !> on machines which conform to the ieee-754 floating point standard.
     !> DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and
     !> when partial spectrum requests are made.
     !> Normal execution of DSTEMR may create NaNs and infinities and
     !> hence may abort due to a floating point exception in environments
     !> which do not handle NaNs and infinities in the ieee standard default
     !> manner.

     pure subroutine la_dstevr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: alleig,indeig,test,lquery,valeig,wantz,tryrac
           character :: order
           integer(ilp) :: i,ieeeok,imax,indibl,indifl,indisp,indiwo,iscale,itmp1,j,jj, &
                     liwmin,lwmin,nsplit
           real(dp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tmp1,tnrm,vll, &
                     vuu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           ieeeok = la_ilaenv(10,'DSTEVR','N',1,2,3,4)
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           lwmin = max(1,20*n)
           liwmin = max(1,10*n)
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else
              if (valeig) then
                 if (n > 0 .and. vu <= vl) info = -7
              else if (indeig) then
                 if (il < 1 .or. il > max(1,n)) then
                    info = -8
                 else if (iu < min(n,il) .or. iu > n) then
                    info = -9
                 end if
              end if
           end if
           if (info == 0) then
              if (ldz < 1 .or. (wantz .and. ldz < n)) then
                 info = -14
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -17
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -19
              end if
           end if
           if (info /= 0) then
              call la_xerbla('DSTEVR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (vl < d(1) .and. vu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           if (valeig) then
              vll = vl
              vuu = vu
           end if
           tnrm = la_dlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_dscal(n,sigma,d,1)
              call la_dscal(n - 1,sigma,e(1),1)
              if (valeig) then
                 vll = vl*sigma
                 vuu = vu*sigma
              end if
           end if
           ! initialize indices into workspaces.  note: these indices are used only
           ! if la_dsterf or la_dstemr fail.
           ! iwork(indibl:indibl+m-1) corresponds to iblock in la_dstebz and
           ! stores the block indices of each of the m<=n eigenvalues.
           indibl = 1
           ! iwork(indisp:indisp+nsplit-1) corresponds to isplit in la_dstebz and
           ! stores the starting and finishing indices of each block.
           indisp = indibl + n
           ! iwork(indifl:indifl+n-1) stores the indices of eigenvectors
           ! that corresponding to eigenvectors that fail to converge in
           ! la_dstein.  this information is discarded; if any fail, the driver
           ! returns info > 0.
           indifl = indisp + n
           ! indiwo is the offset of the remaining integer workspace.
           indiwo = indisp + n
           ! if all eigenvalues are desired, then
           ! call la_dsterf or la_dstemr.  if this fails for some eigenvalue, then
           ! try la_dstebz.
           test = .false.
           if (indeig) then
              if (il == 1 .and. iu == n) then
                 test = .true.
              end if
           end if
           if ((alleig .or. test) .and. ieeeok == 1) then
              call la_dcopy(n - 1,e(1),1,work(1),1)
              if (.not. wantz) then
                 call la_dcopy(n,d,1,w,1)
                 call la_dsterf(n,w,work,info)
              else
                 call la_dcopy(n,d,1,work(n + 1),1)
                 if (abstol <= two*n*eps) then
                    tryrac = .true.
                 else
                    tryrac = .false.
                 end if
                 call la_dstemr(jobz,'A',n,work(n + 1),work,vl,vu,il,iu,m,w,z,ldz, &
                            n,isuppz,tryrac,work(2*n + 1),lwork - 2*n,iwork,liwork,info)
              end if
              if (info == 0) then
                 m = n
                 go to 10
              end if
              info = 0
           end if
           ! otherwise, call la_dstebz and, if eigenvectors are desired, la_dstein.
           if (wantz) then
              order = 'B'
           else
              order = 'E'
           end if
           call la_dstebz(range,order,n,vll,vuu,il,iu,abstol,d,e,m,nsplit,w, &
                     iwork(indibl),iwork(indisp),work,iwork(indiwo),info)
           if (wantz) then
              call la_dstein(n,d,e,m,w,iwork(indibl),iwork(indisp),z,ldz,work, &
                        iwork(indiwo),iwork(indifl),info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           10 continue
           if (iscale == 1) then
              if (info == 0) then
                 imax = m
              else
                 imax = info - 1
              end if
              call la_dscal(imax,one/sigma,w,1)
           end if
           ! if eigenvalues are not in order, then sort them, along with
           ! eigenvectors.
           if (wantz) then
              do j = 1,m - 1
                 i = 0
                 tmp1 = w(j)
                 do jj = j + 1,m
                    if (w(jj) < tmp1) then
                       i = jj
                       tmp1 = w(jj)
                    end if
                 end do
                 if (i /= 0) then
                    itmp1 = iwork(i)
                    w(i) = w(j)
                    iwork(i) = iwork(j)
                    w(j) = tmp1
                    iwork(j) = itmp1
                    call la_dswap(n,z(1,i),1,z(1,j),1)
                 end if
              end do
           end if
            ! causes problems with tests 19
            ! if (wantz .and. indeig ) z( 1,1) = z(1,1) / 1.002_dp + .002
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_dstevr
     !> QSTEVR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T.  Eigenvalues and
     !> eigenvectors can be selected by specifying either a range of values
     !> or a range of indices for the desired eigenvalues.
     !> Whenever possible, QSTEVR calls QSTEMR to compute the
     !> eigenspectrum using Relatively Robust Representations.  QSTEMR
     !> computes eigenvalues by the dqds algorithm, while orthogonal
     !> eigenvectors are computed from various "good" L D L^T representations
     !> (also known as Relatively Robust Representations). Gram-Schmidt
     !> orthogonalization is avoided as far as possible. More specifically,
     !> the various steps of the algorithm are as follows. For the i-th
     !> unreduced block of T,
     !> (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
     !> is a relatively robust representation,
     !> (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
     !> relative accuracy by the dqds algorithm,
     !> (c) If there is a cluster of close eigenvalues, "choose" sigma_i
     !> close to the cluster, and go to step (a),
     !> (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
     !> compute the corresponding eigenvector by forming a
     !> rank-revealing twisted factorization.
     !> The desired accuracy of the output can be specified by the input
     !> parameter ABSTOL.
     !> For more details, see "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
     !> Computer Science Division Technical Report No. UCB//CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Note 1 : QSTEVR calls QSTEMR when the full spectrum is requested
     !> on machines which conform to the ieee-754 floating point standard.
     !> QSTEVR calls QSTEBZ and QSTEIN on non-ieee machines and
     !> when partial spectrum requests are made.
     !> Normal execution of QSTEMR may create NaNs and infinities and
     !> hence may abort due to a floating point exception in environments
     !> which do not handle NaNs and infinities in the ieee standard default
     !> manner.

     pure subroutine la_qstevr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one,two
        ! -- lapack driver routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*),z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: alleig,indeig,test,lquery,valeig,wantz,tryrac
           character :: order
           integer(ilp) :: i,ieeeok,imax,indibl,indifl,indisp,indiwo,iscale,itmp1,j,jj, &
                     liwmin,lwmin,nsplit
           real(qp) :: bignum,eps,rmax,rmin,safmin,sigma,smlnum,tmp1,tnrm,vll, &
                     vuu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           ieeeok = la_ilaenv(10,'QSTEVR','N',1,2,3,4)
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           lwmin = max(1,20*n)
           liwmin = max(1,10*n)
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else
              if (valeig) then
                 if (n > 0 .and. vu <= vl) info = -7
              else if (indeig) then
                 if (il < 1 .or. il > max(1,n)) then
                    info = -8
                 else if (iu < min(n,il) .or. iu > n) then
                    info = -9
                 end if
              end if
           end if
           if (info == 0) then
              if (ldz < 1 .or. (wantz .and. ldz < n)) then
                 info = -14
              end if
           end if
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -17
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -19
              end if
           end if
           if (info /= 0) then
              call la_xerbla('QSTEVR',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (vl < d(1) .and. vu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz) z(1,1) = one
              return
           end if
           ! get machine constants.
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           ! scale matrix to allowable range, if necessary.
           iscale = 0
           if (valeig) then
              vll = vl
              vuu = vu
           end if
           tnrm = la_qlanst('M',n,d,e)
           if (tnrm > zero .and. tnrm < rmin) then
              iscale = 1
              sigma = rmin/tnrm
           else if (tnrm > rmax) then
              iscale = 1
              sigma = rmax/tnrm
           end if
           if (iscale == 1) then
              call la_qscal(n,sigma,d,1)
              call la_qscal(n - 1,sigma,e(1),1)
              if (valeig) then
                 vll = vl*sigma
                 vuu = vu*sigma
              end if
           end if
           ! initialize indices into workspaces.  note: these indices are used only
           ! if la_qsterf or la_qstemr fail.
           ! iwork(indibl:indibl+m-1) corresponds to iblock in la_qstebz and
           ! stores the block indices of each of the m<=n eigenvalues.
           indibl = 1
           ! iwork(indisp:indisp+nsplit-1) corresponds to isplit in la_qstebz and
           ! stores the starting and finishing indices of each block.
           indisp = indibl + n
           ! iwork(indifl:indifl+n-1) stores the indices of eigenvectors
           ! that corresponding to eigenvectors that fail to converge in
           ! la_qstein.  this information is discarded; if any fail, the driver
           ! returns info > 0.
           indifl = indisp + n
           ! indiwo is the offset of the remaining integer workspace.
           indiwo = indisp + n
           ! if all eigenvalues are desired, then
           ! call la_qsterf or la_qstemr.  if this fails for some eigenvalue, then
           ! try la_qstebz.
           test = .false.
           if (indeig) then
              if (il == 1 .and. iu == n) then
                 test = .true.
              end if
           end if
           if ((alleig .or. test) .and. ieeeok == 1) then
              call la_qcopy(n - 1,e(1),1,work(1),1)
              if (.not. wantz) then
                 call la_qcopy(n,d,1,w,1)
                 call la_qsterf(n,w,work,info)
              else
                 call la_qcopy(n,d,1,work(n + 1),1)
                 if (abstol <= two*n*eps) then
                    tryrac = .true.
                 else
                    tryrac = .false.
                 end if
                 call la_qstemr(jobz,'A',n,work(n + 1),work,vl,vu,il,iu,m,w,z,ldz, &
                            n,isuppz,tryrac,work(2*n + 1),lwork - 2*n,iwork,liwork,info)
              end if
              if (info == 0) then
                 m = n
                 go to 10
              end if
              info = 0
           end if
           ! otherwise, call la_qstebz and, if eigenvectors are desired, la_qstein.
           if (wantz) then
              order = 'B'
           else
              order = 'E'
           end if
           call la_qstebz(range,order,n,vll,vuu,il,iu,abstol,d,e,m,nsplit,w, &
                     iwork(indibl),iwork(indisp),work,iwork(indiwo),info)
           if (wantz) then
              call la_qstein(n,d,e,m,w,iwork(indibl),iwork(indisp),z,ldz,work, &
                        iwork(indiwo),iwork(indifl),info)
           end if
           ! if matrix was scaled, then rescale eigenvalues appropriately.
           10 continue
           if (iscale == 1) then
              if (info == 0) then
                 imax = m
              else
                 imax = info - 1
              end if
              call la_qscal(imax,one/sigma,w,1)
           end if
           ! if eigenvalues are not in order, then sort them, along with
           ! eigenvectors.
           if (wantz) then
              do j = 1,m - 1
                 i = 0
                 tmp1 = w(j)
                 do jj = j + 1,m
                    if (w(jj) < tmp1) then
                       i = jj
                       tmp1 = w(jj)
                    end if
                 end do
                 if (i /= 0) then
                    itmp1 = iwork(i)
                    w(i) = w(j)
                    iwork(i) = iwork(j)
                    w(j) = tmp1
                    iwork(j) = itmp1
                    call la_qswap(n,z(1,i),1,z(1,j),1)
                 end if
              end do
           end if
            ! causes problems with tests 19
            ! if (wantz .and. indeig ) z( 1,1) = z(1,1) / 1.002_qp + .002
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_qstevr

     !> CSTEIN: computes the eigenvectors of a real symmetric tridiagonal
     !> matrix T corresponding to specified eigenvalues, using inverse
     !> iteration.
     !> The maximum number of iterations allowed for each eigenvector is
     !> specified by an internal parameter MAXITS (currently set to 5).
     !> Although the eigenvectors are real, they are stored in a complex
     !> array, which may be passed to CUNMTR or CUPMTR for back
     !> transformation to the eigenvectors of a complex Hermitian matrix
     !> which was reduced to tridiagonal form.

     pure subroutine la_cstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail, &
               info)
        use la_constants_sp,only:zero,one,ten,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,m,n
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),isplit(*)
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(sp),intent(in) :: d(*),e(*),w(*)
           real(sp),intent(out) :: work(*)
           complex(sp),intent(out) :: z(ldz,*)
       ! =====================================================================
           ! Parameters
           real(sp),parameter :: odm3 = 1.0e-3_sp
           real(sp),parameter :: odm1 = 1.0e-1_sp
           integer(ilp),parameter :: maxits = 5
           integer(ilp),parameter :: extra = 2

           ! Local Scalars
           integer(ilp) :: b1,blksiz,bn,gpind,i,iinfo,indrv1,indrv2,indrv3,indrv4, &
                     indrv5,its,j,j1,jblk,jmax,jr,nblk,nrmchk
           real(sp) :: ctr,eps,eps1,nrm,onenrm,ortol,pertol,scl,sep,stpcrt,tol,xj, &
                     xjm
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,cmplx,max,real,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           do i = 1,m
              ifail(i) = 0
           end do
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -4
           else if (ldz < max(1,n)) then
              info = -9
           else
              do j = 2,m
                 if (iblock(j) < iblock(j - 1)) then
                    info = -6
                    go to 30
                 end if
                 if (iblock(j) == iblock(j - 1) .and. w(j) < w(j - 1)) then
                    info = -5
                    go to 30
                 end if
              end do
              30 continue
           end if
           if (info /= 0) then
              call la_xerbla('CSTEIN',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) then
              return
           else if (n == 1) then
              z(1,1) = cone
              return
           end if
           ! get machine constants.
           eps = la_slamch('PRECISION')
           ! initialize seed for random number generator la_slarnv.
           do i = 1,4
              iseed(i) = 1
           end do
           ! initialize pointers.
           indrv1 = 0
           indrv2 = indrv1 + n
           indrv3 = indrv2 + n
           indrv4 = indrv3 + n
           indrv5 = indrv4 + n
           ! compute eigenvectors of matrix blocks.
           j1 = 1
           loop_180: do nblk = 1,iblock(m)
              ! find starting and ending indices of block nblk.
              if (nblk == 1) then
                 b1 = 1
              else
                 b1 = isplit(nblk - 1) + 1
              end if
              bn = isplit(nblk)
              blksiz = bn - b1 + 1
              if (blksiz == 1) go to 60
              gpind = j1
              ! compute reorthogonalization criterion and stopping criterion.
              onenrm = abs(d(b1)) + abs(e(b1))
              onenrm = max(onenrm,abs(d(bn)) + abs(e(bn - 1)))
              do i = b1 + 1,bn - 1
                 onenrm = max(onenrm,abs(d(i)) + abs(e(i - 1)) + abs(e(i)))
              end do
              ortol = odm3*onenrm
              stpcrt = sqrt(odm1/blksiz)
              ! loop through eigenvalues of block nblk.
              60 continue
              jblk = 0
              loop_170: do j = j1,m
                 if (iblock(j) /= nblk) then
                    j1 = j
                    cycle loop_180
                 end if
                 jblk = jblk + 1
                 xj = w(j)
                 ! skip all the work if the block size is one.
                 if (blksiz == 1) then
                    work(indrv1 + 1) = one
                    go to 140
                 end if
                 ! if eigenvalues j and j-1 are too close, add a relatively
                 ! small perturbation.
                 if (jblk > 1) then
                    eps1 = abs(eps*xj)
                    pertol = ten*eps1
                    sep = xj - xjm
                    if (sep < pertol) xj = xjm + pertol
                 end if
                 its = 0
                 nrmchk = 0
                 ! get random starting vector.
                 call la_slarnv(2,iseed,blksiz,work(indrv1 + 1))
                 ! copy the matrix t so it won't be destroyed in factorization.
                 call la_scopy(blksiz,d(b1),1,work(indrv4 + 1),1)
                 call la_scopy(blksiz - 1,e(b1),1,work(indrv2 + 2),1)
                 call la_scopy(blksiz - 1,e(b1),1,work(indrv3 + 1),1)
                 ! compute lu factors with partial pivoting  ( pt = lu )
                 tol = zero
                 call la_slagtf(blksiz,work(indrv4 + 1),xj,work(indrv2 + 2),work(indrv3 + &
                           1),tol,work(indrv5 + 1),iwork,iinfo)
                 ! update iteration count.
                 70 continue
                 its = its + 1
                 if (its > maxits) go to 120
                 ! normalize and scale the righthand side vector pb.
                 jmax = la_isamax(blksiz,work(indrv1 + 1),1)
                 scl = blksiz*onenrm*max(eps,abs(work(indrv4 + blksiz)))/abs(work(indrv1 + &
                           jmax))
                 call la_sscal(blksiz,scl,work(indrv1 + 1),1)
                 ! solve the system lu = pb.
                 call la_slagts(-1,blksiz,work(indrv4 + 1),work(indrv2 + 2),work(indrv3 + &
                           1),work(indrv5 + 1),iwork,work(indrv1 + 1),tol,iinfo)
                 ! reorthogonalize by modified gram-schmidt if eigenvalues are
                 ! close enough.
                 if (jblk == 1) go to 110
                 if (abs(xj - xjm) > ortol) gpind = j
                 if (gpind /= j) then
                    do i = gpind,j - 1
                       ctr = zero
                       do jr = 1,blksiz
                          ctr = ctr + work(indrv1 + jr)*real(z(b1 - 1 + jr,i),KIND=sp)
                       end do
                       do jr = 1,blksiz
                          work(indrv1 + jr) = work(indrv1 + jr) - ctr*real(z(b1 - 1 + jr,i), &
                                    KIND=sp)
                       end do
                    end do
                 end if
                 ! check the infinity norm of the iterate.
                 110 continue
                 jmax = la_isamax(blksiz,work(indrv1 + 1),1)
                 nrm = abs(work(indrv1 + jmax))
                 ! continue for additional iterations after norm reaches
                 ! stopping criterion.
                 if (nrm < stpcrt) go to 70
                 nrmchk = nrmchk + 1
                 if (nrmchk < extra + 1) go to 70
                 go to 130
                 ! if stopping criterion was not satisfied, update info and
                 ! store eigenvector number in array ifail.
                 120 continue
                 info = info + 1
                 ifail(info) = j
                 ! accept iterate as jth eigenvector.
                 130 continue
                 scl = one/la_snrm2(blksiz,work(indrv1 + 1),1)
                 jmax = la_isamax(blksiz,work(indrv1 + 1),1)
                 if (work(indrv1 + jmax) < zero) scl = -scl
                 call la_sscal(blksiz,scl,work(indrv1 + 1),1)
                 140 continue
                 do i = 1,n
                    z(i,j) = czero
                 end do
                 do i = 1,blksiz
                    z(b1 + i - 1,j) = cmplx(work(indrv1 + i),zero,KIND=sp)
                 end do
                 ! save the shift to check eigenvalue spacing at next
                 ! iteration.
                 xjm = xj
              end do loop_170
           end do loop_180
           return
     end subroutine la_cstein
     !> ZSTEIN: computes the eigenvectors of a real symmetric tridiagonal
     !> matrix T corresponding to specified eigenvalues, using inverse
     !> iteration.
     !> The maximum number of iterations allowed for each eigenvector is
     !> specified by an internal parameter MAXITS (currently set to 5).
     !> Although the eigenvectors are real, they are stored in a complex
     !> array, which may be passed to ZUNMTR or ZUPMTR for back
     !> transformation to the eigenvectors of a complex Hermitian matrix
     !> which was reduced to tridiagonal form.

     pure subroutine la_zstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail, &
               info)
        use la_constants_dp,only:zero,one,ten,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,m,n
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),isplit(*)
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(dp),intent(in) :: d(*),e(*),w(*)
           real(dp),intent(out) :: work(*)
           complex(dp),intent(out) :: z(ldz,*)
       ! =====================================================================
           ! Parameters
           real(dp),parameter :: odm3 = 1.0e-3_dp
           real(dp),parameter :: odm1 = 1.0e-1_dp
           integer(ilp),parameter :: maxits = 5
           integer(ilp),parameter :: extra = 2

           ! Local Scalars
           integer(ilp) :: b1,blksiz,bn,gpind,i,iinfo,indrv1,indrv2,indrv3,indrv4, &
                     indrv5,its,j,j1,jblk,jmax,jr,nblk,nrmchk
           real(dp) :: dtpcrt,eps,eps1,nrm,onenrm,ortol,pertol,scl,sep,tol,xj,xjm, &
                     ztr
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           do i = 1,m
              ifail(i) = 0
           end do
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -4
           else if (ldz < max(1,n)) then
              info = -9
           else
              do j = 2,m
                 if (iblock(j) < iblock(j - 1)) then
                    info = -6
                    go to 30
                 end if
                 if (iblock(j) == iblock(j - 1) .and. w(j) < w(j - 1)) then
                    info = -5
                    go to 30
                 end if
              end do
              30 continue
           end if
           if (info /= 0) then
              call la_xerbla('ZSTEIN',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) then
              return
           else if (n == 1) then
              z(1,1) = cone
              return
           end if
           ! get machine constants.
           eps = la_dlamch('PRECISION')
           ! initialize seed for random number generator la_dlarnv.
           do i = 1,4
              iseed(i) = 1
           end do
           ! initialize pointers.
           indrv1 = 0
           indrv2 = indrv1 + n
           indrv3 = indrv2 + n
           indrv4 = indrv3 + n
           indrv5 = indrv4 + n
           ! compute eigenvectors of matrix blocks.
           j1 = 1
           loop_180: do nblk = 1,iblock(m)
              ! find starting and ending indices of block nblk.
              if (nblk == 1) then
                 b1 = 1
              else
                 b1 = isplit(nblk - 1) + 1
              end if
              bn = isplit(nblk)
              blksiz = bn - b1 + 1
              if (blksiz == 1) go to 60
              gpind = j1
              ! compute reorthogonalization criterion and stopping criterion.
              onenrm = abs(d(b1)) + abs(e(b1))
              onenrm = max(onenrm,abs(d(bn)) + abs(e(bn - 1)))
              do i = b1 + 1,bn - 1
                 onenrm = max(onenrm,abs(d(i)) + abs(e(i - 1)) + abs(e(i)))
              end do
              ortol = odm3*onenrm
              dtpcrt = sqrt(odm1/blksiz)
              ! loop through eigenvalues of block nblk.
              60 continue
              jblk = 0
              loop_170: do j = j1,m
                 if (iblock(j) /= nblk) then
                    j1 = j
                    cycle loop_180
                 end if
                 jblk = jblk + 1
                 xj = w(j)
                 ! skip all the work if the block size is one.
                 if (blksiz == 1) then
                    work(indrv1 + 1) = one
                    go to 140
                 end if
                 ! if eigenvalues j and j-1 are too close, add a relatively
                 ! small perturbation.
                 if (jblk > 1) then
                    eps1 = abs(eps*xj)
                    pertol = ten*eps1
                    sep = xj - xjm
                    if (sep < pertol) xj = xjm + pertol
                 end if
                 its = 0
                 nrmchk = 0
                 ! get random starting vector.
                 call la_dlarnv(2,iseed,blksiz,work(indrv1 + 1))
                 ! copy the matrix t so it won't be destroyed in factorization.
                 call la_dcopy(blksiz,d(b1),1,work(indrv4 + 1),1)
                 call la_dcopy(blksiz - 1,e(b1),1,work(indrv2 + 2),1)
                 call la_dcopy(blksiz - 1,e(b1),1,work(indrv3 + 1),1)
                 ! compute lu factors with partial pivoting  ( pt = lu )
                 tol = zero
                 call la_dlagtf(blksiz,work(indrv4 + 1),xj,work(indrv2 + 2),work(indrv3 + &
                           1),tol,work(indrv5 + 1),iwork,iinfo)
                 ! update iteration count.
                 70 continue
                 its = its + 1
                 if (its > maxits) go to 120
                 ! normalize and scale the righthand side vector pb.
                 jmax = la_idamax(blksiz,work(indrv1 + 1),1)
                 scl = blksiz*onenrm*max(eps,abs(work(indrv4 + blksiz)))/abs(work(indrv1 + &
                           jmax))
                 call la_dscal(blksiz,scl,work(indrv1 + 1),1)
                 ! solve the system lu = pb.
                 call la_dlagts(-1,blksiz,work(indrv4 + 1),work(indrv2 + 2),work(indrv3 + &
                           1),work(indrv5 + 1),iwork,work(indrv1 + 1),tol,iinfo)
                 ! reorthogonalize by modified gram-schmidt if eigenvalues are
                 ! close enough.
                 if (jblk == 1) go to 110
                 if (abs(xj - xjm) > ortol) gpind = j
                 if (gpind /= j) then
                    do i = gpind,j - 1
                       ztr = zero
                       do jr = 1,blksiz
                          ztr = ztr + work(indrv1 + jr)*real(z(b1 - 1 + jr,i),KIND=dp)
                       end do
                       do jr = 1,blksiz
                          work(indrv1 + jr) = work(indrv1 + jr) - ztr*real(z(b1 - 1 + jr,i), &
                                    KIND=dp)
                       end do
                    end do
                 end if
                 ! check the infinity norm of the iterate.
                 110 continue
                 jmax = la_idamax(blksiz,work(indrv1 + 1),1)
                 nrm = abs(work(indrv1 + jmax))
                 ! continue for additional iterations after norm reaches
                 ! stopping criterion.
                 if (nrm < dtpcrt) go to 70
                 nrmchk = nrmchk + 1
                 if (nrmchk < extra + 1) go to 70
                 go to 130
                 ! if stopping criterion was not satisfied, update info and
                 ! store eigenvector number in array ifail.
                 120 continue
                 info = info + 1
                 ifail(info) = j
                 ! accept iterate as jth eigenvector.
                 130 continue
                 scl = one/la_dnrm2(blksiz,work(indrv1 + 1),1)
                 jmax = la_idamax(blksiz,work(indrv1 + 1),1)
                 if (work(indrv1 + jmax) < zero) scl = -scl
                 call la_dscal(blksiz,scl,work(indrv1 + 1),1)
                 140 continue
                 do i = 1,n
                    z(i,j) = czero
                 end do
                 do i = 1,blksiz
                    z(b1 + i - 1,j) = cmplx(work(indrv1 + i),zero,KIND=dp)
                 end do
                 ! save the shift to check eigenvalue spacing at next
                 ! iteration.
                 xjm = xj
              end do loop_170
           end do loop_180
           return
     end subroutine la_zstein
     !> WSTEIN: computes the eigenvectors of a real symmetric tridiagonal
     !> matrix T corresponding to specified eigenvalues, using inverse
     !> iteration.
     !> The maximum number of iterations allowed for each eigenvector is
     !> specified by an internal parameter MAXITS (currently set to 5).
     !> Although the eigenvectors are real, they are stored in a complex
     !> array, which may be passed to WUNMTR or WUPMTR for back
     !> transformation to the eigenvectors of a complex Hermitian matrix
     !> which was reduced to tridiagonal form.

     pure subroutine la_wstein(n,d,e,m,w,iblock,isplit,z,ldz,work,iwork,ifail, &
               info)
        use la_constants_qp,only:zero,one,ten,czero,cone
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,m,n
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),isplit(*)
           integer(ilp),intent(out) :: ifail(*),iwork(*)
           real(qp),intent(in) :: d(*),e(*),w(*)
           real(qp),intent(out) :: work(*)
           complex(qp),intent(out) :: z(ldz,*)
       ! =====================================================================
           ! Parameters
           real(qp),parameter :: odm3 = 1.0e-3_qp
           real(qp),parameter :: odm1 = 1.0e-1_qp
           integer(ilp),parameter :: maxits = 5
           integer(ilp),parameter :: extra = 2

           ! Local Scalars
           integer(ilp) :: b1,blksiz,bn,gpind,i,iinfo,indrv1,indrv2,indrv3,indrv4, &
                     indrv5,its,j,j1,jblk,jmax,jr,nblk,nrmchk
           real(qp) :: qtpcrt,eps,eps1,nrm,onenrm,ortol,pertol,scl,sep,tol,xj,xjm, &
                     wtr
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,real,cmplx,max,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           do i = 1,m
              ifail(i) = 0
           end do
           if (n < 0) then
              info = -1
           else if (m < 0 .or. m > n) then
              info = -4
           else if (ldz < max(1,n)) then
              info = -9
           else
              do j = 2,m
                 if (iblock(j) < iblock(j - 1)) then
                    info = -6
                    go to 30
                 end if
                 if (iblock(j) == iblock(j - 1) .and. w(j) < w(j - 1)) then
                    info = -5
                    go to 30
                 end if
              end do
              30 continue
           end if
           if (info /= 0) then
              call la_xerbla('WSTEIN',-info)
              return
           end if
           ! quick return if possible
           if (n == 0 .or. m == 0) then
              return
           else if (n == 1) then
              z(1,1) = cone
              return
           end if
           ! get machine constants.
           eps = la_qlamch('PRECISION')
           ! initialize seed for random number generator la_qlarnv.
           do i = 1,4
              iseed(i) = 1
           end do
           ! initialize pointers.
           indrv1 = 0
           indrv2 = indrv1 + n
           indrv3 = indrv2 + n
           indrv4 = indrv3 + n
           indrv5 = indrv4 + n
           ! compute eigenvectors of matrix blocks.
           j1 = 1
           loop_180: do nblk = 1,iblock(m)
              ! find starting and ending indices of block nblk.
              if (nblk == 1) then
                 b1 = 1
              else
                 b1 = isplit(nblk - 1) + 1
              end if
              bn = isplit(nblk)
              blksiz = bn - b1 + 1
              if (blksiz == 1) go to 60
              gpind = j1
              ! compute reorthogonalization criterion and stopping criterion.
              onenrm = abs(d(b1)) + abs(e(b1))
              onenrm = max(onenrm,abs(d(bn)) + abs(e(bn - 1)))
              do i = b1 + 1,bn - 1
                 onenrm = max(onenrm,abs(d(i)) + abs(e(i - 1)) + abs(e(i)))
              end do
              ortol = odm3*onenrm
              qtpcrt = sqrt(odm1/blksiz)
              ! loop through eigenvalues of block nblk.
              60 continue
              jblk = 0
              loop_170: do j = j1,m
                 if (iblock(j) /= nblk) then
                    j1 = j
                    cycle loop_180
                 end if
                 jblk = jblk + 1
                 xj = w(j)
                 ! skip all the work if the block size is one.
                 if (blksiz == 1) then
                    work(indrv1 + 1) = one
                    go to 140
                 end if
                 ! if eigenvalues j and j-1 are too close, add a relatively
                 ! small perturbation.
                 if (jblk > 1) then
                    eps1 = abs(eps*xj)
                    pertol = ten*eps1
                    sep = xj - xjm
                    if (sep < pertol) xj = xjm + pertol
                 end if
                 its = 0
                 nrmchk = 0
                 ! get random starting vector.
                 call la_qlarnv(2,iseed,blksiz,work(indrv1 + 1))
                 ! copy the matrix t so it won't be destroyed in factorization.
                 call la_qcopy(blksiz,d(b1),1,work(indrv4 + 1),1)
                 call la_qcopy(blksiz - 1,e(b1),1,work(indrv2 + 2),1)
                 call la_qcopy(blksiz - 1,e(b1),1,work(indrv3 + 1),1)
                 ! compute lu factors with partial pivoting  ( pt = lu )
                 tol = zero
                 call la_qlagtf(blksiz,work(indrv4 + 1),xj,work(indrv2 + 2),work(indrv3 + &
                           1),tol,work(indrv5 + 1),iwork,iinfo)
                 ! update iteration count.
                 70 continue
                 its = its + 1
                 if (its > maxits) go to 120
                 ! normalize and scale the righthand side vector pb.
                 jmax = la_iqamax(blksiz,work(indrv1 + 1),1)
                 scl = blksiz*onenrm*max(eps,abs(work(indrv4 + blksiz)))/abs(work(indrv1 + &
                           jmax))
                 call la_qscal(blksiz,scl,work(indrv1 + 1),1)
                 ! solve the system lu = pb.
                 call la_qlagts(-1,blksiz,work(indrv4 + 1),work(indrv2 + 2),work(indrv3 + &
                           1),work(indrv5 + 1),iwork,work(indrv1 + 1),tol,iinfo)
                 ! reorthogonalize by modified gram-schmidt if eigenvalues are
                 ! close enough.
                 if (jblk == 1) go to 110
                 if (abs(xj - xjm) > ortol) gpind = j
                 if (gpind /= j) then
                    do i = gpind,j - 1
                       wtr = zero
                       do jr = 1,blksiz
                          wtr = wtr + work(indrv1 + jr)*real(z(b1 - 1 + jr,i),KIND=qp)
                       end do
                       do jr = 1,blksiz
                          work(indrv1 + jr) = work(indrv1 + jr) - wtr*real(z(b1 - 1 + jr,i), &
                                    KIND=qp)
                       end do
                    end do
                 end if
                 ! check the infinity norm of the iterate.
                 110 continue
                 jmax = la_iqamax(blksiz,work(indrv1 + 1),1)
                 nrm = abs(work(indrv1 + jmax))
                 ! continue for additional iterations after norm reaches
                 ! stopping criterion.
                 if (nrm < qtpcrt) go to 70
                 nrmchk = nrmchk + 1
                 if (nrmchk < extra + 1) go to 70
                 go to 130
                 ! if stopping criterion was not satisfied, update info and
                 ! store eigenvector number in array ifail.
                 120 continue
                 info = info + 1
                 ifail(info) = j
                 ! accept iterate as jth eigenvector.
                 130 continue
                 scl = one/la_qnrm2(blksiz,work(indrv1 + 1),1)
                 jmax = la_iqamax(blksiz,work(indrv1 + 1),1)
                 if (work(indrv1 + jmax) < zero) scl = -scl
                 call la_qscal(blksiz,scl,work(indrv1 + 1),1)
                 140 continue
                 do i = 1,n
                    z(i,j) = czero
                 end do
                 do i = 1,blksiz
                    z(b1 + i - 1,j) = cmplx(work(indrv1 + i),zero,KIND=qp)
                 end do
                 ! save the shift to check eigenvalue spacing at next
                 ! iteration.
                 xjm = xj
              end do loop_170
           end do loop_180
           return
     end subroutine la_wstein

     !> CPTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric positive definite tridiagonal matrix by first factoring the
     !> matrix using SPTTRF and then calling CBDSQR to compute the singular
     !> values of the bidiagonal factor.
     !> This routine computes the eigenvalues of the positive definite
     !> tridiagonal matrix to high relative accuracy.  This means that if the
     !> eigenvalues range over many orders of magnitude in size, then the
     !> small eigenvalues and corresponding eigenvectors will be computed
     !> more accurately than, for example, with the standard QR method.
     !> The eigenvectors of a full or band positive definite Hermitian matrix
     !> can also be found if CHETRD, CHPTRD, or CHBTRD has been used to
     !> reduce this matrix to tridiagonal form.  (The reduction to
     !> tridiagonal form, however, may preclude the possibility of obtaining
     !> high relative accuracy in the small eigenvalues of the original
     !> matrix, if these eigenvalues range over many orders of magnitude.)

     pure subroutine la_cpteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_sp
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
        ! ====================================================================

           ! Local Arrays
           complex(sp) :: c(1,1),vt(1,1)
           ! Local Scalars
           integer(ilp) :: i,icompz,nru
           ! Intrinsic Functions
           intrinsic :: max,sqrt
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
              call la_xerbla('CPTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz > 0) z(1,1) = cone
              return
           end if
           if (icompz == 2) call la_claset('FULL',n,n,czero,cone,z,ldz)
           ! call la_spttrf to factor the matrix.
           call la_spttrf(n,d,e,info)
           if (info /= 0) return
           do i = 1,n
              d(i) = sqrt(d(i))
           end do
           do i = 1,n - 1
              e(i) = e(i)*d(i)
           end do
           ! call la_cbdsqr to compute the singular values/vectors of the
           ! bidiagonal factor.
           if (icompz > 0) then
              nru = n
           else
              nru = 0
           end if
           call la_cbdsqr('LOWER',n,0,nru,0,d,e,vt,1,z,ldz,c,1,work,info)

           ! square the singular values.
           if (info == 0) then
              do i = 1,n
                 d(i) = d(i)*d(i)
              end do
           else
              info = n + info
           end if
           return
     end subroutine la_cpteqr
     !> ZPTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric positive definite tridiagonal matrix by first factoring the
     !> matrix using DPTTRF and then calling ZBDSQR to compute the singular
     !> values of the bidiagonal factor.
     !> This routine computes the eigenvalues of the positive definite
     !> tridiagonal matrix to high relative accuracy.  This means that if the
     !> eigenvalues range over many orders of magnitude in size, then the
     !> small eigenvalues and corresponding eigenvectors will be computed
     !> more accurately than, for example, with the standard QR method.
     !> The eigenvectors of a full or band positive definite Hermitian matrix
     !> can also be found if ZHETRD, ZHPTRD, or ZHBTRD has been used to
     !> reduce this matrix to tridiagonal form.  (The reduction to
     !> tridiagonal form, however, may preclude the possibility of obtaining
     !> high relative accuracy in the small eigenvalues of the original
     !> matrix, if these eigenvalues range over many orders of magnitude.)

     pure subroutine la_zpteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_dp
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
        ! ====================================================================

           ! Local Arrays
           complex(dp) :: c(1,1),vt(1,1)
           ! Local Scalars
           integer(ilp) :: i,icompz,nru
           ! Intrinsic Functions
           intrinsic :: max,sqrt
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
              call la_xerbla('ZPTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz > 0) z(1,1) = cone
              return
           end if
           if (icompz == 2) call la_zlaset('FULL',n,n,czero,cone,z,ldz)
           ! call la_dpttrf to factor the matrix.
           call la_dpttrf(n,d,e,info)
           if (info /= 0) return
           do i = 1,n
              d(i) = sqrt(d(i))
           end do
           do i = 1,n - 1
              e(i) = e(i)*d(i)
           end do
           ! call la_zbdsqr to compute the singular values/vectors of the
           ! bidiagonal factor.
           if (icompz > 0) then
              nru = n
           else
              nru = 0
           end if
           call la_zbdsqr('LOWER',n,0,nru,0,d,e,vt,1,z,ldz,c,1,work,info)

           ! square the singular values.
           if (info == 0) then
              do i = 1,n
                 d(i) = d(i)*d(i)
              end do
           else
              info = n + info
           end if
           return
     end subroutine la_zpteqr
     !> WPTEQR: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric positive definite tridiagonal matrix by first factoring the
     !> matrix using QPTTRF and then calling WBDSQR to compute the singular
     !> values of the bidiagonal factor.
     !> This routine computes the eigenvalues of the positive definite
     !> tridiagonal matrix to high relative accuracy.  This means that if the
     !> eigenvalues range over many orders of magnitude in size, then the
     !> small eigenvalues and corresponding eigenvectors will be computed
     !> more accurately than, for example, with the standard QR method.
     !> The eigenvectors of a full or band positive definite Hermitian matrix
     !> can also be found if WHETRD, WHPTRD, or WHBTRD has been used to
     !> reduce this matrix to tridiagonal form.  (The reduction to
     !> tridiagonal form, however, may preclude the possibility of obtaining
     !> high relative accuracy in the small eigenvalues of the original
     !> matrix, if these eigenvalues range over many orders of magnitude.)

     pure subroutine la_wpteqr(compz,n,d,e,z,ldz,work,info)
        use la_constants_qp
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
        ! ====================================================================

           ! Local Arrays
           complex(qp) :: c(1,1),vt(1,1)
           ! Local Scalars
           integer(ilp) :: i,icompz,nru
           ! Intrinsic Functions
           intrinsic :: max,sqrt
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
              call la_xerbla('WPTEQR',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz > 0) z(1,1) = cone
              return
           end if
           if (icompz == 2) call la_wlaset('FULL',n,n,czero,cone,z,ldz)
           ! call la_qpttrf to factor the matrix.
           call la_qpttrf(n,d,e,info)
           if (info /= 0) return
           do i = 1,n
              d(i) = sqrt(d(i))
           end do
           do i = 1,n - 1
              e(i) = e(i)*d(i)
           end do
           ! call la_wbdsqr to compute the singular values/vectors of the
           ! bidiagonal factor.
           if (icompz > 0) then
              nru = n
           else
              nru = 0
           end if
           call la_wbdsqr('LOWER',n,0,nru,0,d,e,vt,1,z,ldz,c,1,work,info)

           ! square the singular values.
           if (info == 0) then
              do i = 1,n
                 d(i) = d(i)*d(i)
              end do
           else
              info = n + info
           end if
           return
     end subroutine la_wpteqr

     !> CSTEMR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> Depending on the number of desired eigenvalues, these are computed either
     !> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
     !> computed by the use of various suitable L D L^T factorizations near clusters
     !> of close eigenvalues (referred to as RRRs, Relatively Robust
     !> Representations). An informal sketch of the algorithm follows.
     !> For each unreduced block (submatrix) of T,
     !> (a) Compute T - sigma I  = L D L^T, so that L and D
     !> define all the wanted eigenvalues to high relative accuracy.
     !> This means that small relative changes in the entries of D and L
     !> cause only small relative changes in the eigenvalues and
     !> eigenvectors. The standard (unfactored) representation of the
     !> tridiagonal matrix T does not have this property in general.
     !> (b) Compute the eigenvalues to suitable accuracy.
     !> If the eigenvectors are desired, the algorithm attains full
     !> accuracy of the computed eigenvalues only right before
     !> the corresponding vectors have to be computed, see steps c) and d).
     !> (c) For each cluster of close eigenvalues, select a new
     !> shift close to the cluster, find a new factorization, and refine
     !> the shifted eigenvalues to suitable accuracy.
     !> (d) For each eigenvalue with a large enough relative separation compute
     !> the corresponding eigenvector by forming a rank revealing twisted
     !> factorization. Go back to (c) for any clusters that remain.
     !> For more details, see:
     !> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
     !> to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
     !> Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     !> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
     !> Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
     !> 2004.  Also LAPACK Working Note 154.
     !> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem",
     !> Computer Science Division Technical Report No. UCB/CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Further Details
     !> 1.CSTEMR works only on machines which follow IEEE-754
     !> floating-point standard in their handling of infinities and NaNs.
     !> This permits the use of efficient inner loops avoiding a check for
     !> zero divisors.
     !> 2. LAPACK routines can be used to reduce a complex Hermitean matrix to
     !> real symmetric tridiagonal form.
     !> (Any complex Hermitean tridiagonal matrix has real values on its diagonal
     !> and potentially complex numbers on its off-diagonals. By applying a
     !> similarity transform with an appropriate diagonal matrix
     !> diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex Hermitean
     !> matrix can be transformed into a real symmetric matrix and complex
     !> arithmetic can be entirely avoided.)
     !> While the eigenvectors of the real symmetric tridiagonal matrix are real,
     !> the eigenvectors of original complex Hermitean matrix have complex entries
     !> in general.
     !> Since LAPACK drivers overwrite the matrix data with the eigenvectors,
     !> CSTEMR accepts complex workspace to facilitate interoperability
     !> with CUNMTR or CUPMTR.

     pure subroutine la_cstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc, &
               isuppz,tryrac,work,lwork,iwork,liwork,info)
        use la_constants_sp,only:zero,one,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           logical(lk),intent(inout) :: tryrac
           integer(ilp),intent(in) :: il,iu,ldz,nzc,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*)
           complex(sp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: minrgp = 3.0e-3_sp

           ! Local Scalars
           logical(lk) :: alleig,indeig,lquery,valeig,wantz,zquery
           integer(ilp) :: i,ibegin,iend,ifirst,iil,iindbl,iindw,iindwk,iinfo,iinspl, &
           iiu,ilast,in,indd,inde2,inderr,indgp,indgrs,indwrk,itmp,itmp2,j,jblk,jj, &
                     liwmin,lwmin,nsplit,nzcmin,offset,wbegin,wend
           real(sp) :: bignum,cs,eps,pivmin,r1,r2,rmax,rmin,rtol1,rtol2,safmin,scale, &
                     smlnum,sn,thresh,tmp,tnrm,wl,wu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           zquery = (nzc == -1)
           ! la_sstemr needs work of size 6*n, iwork of size 3*n.
           ! in addition, la_slarre needs work of size 6*n, iwork of size 5*n.
           ! furthermore, la_clarrv needs work of size 12*n, iwork of size 7*n.
           if (wantz) then
              lwmin = 18*n
              liwmin = 10*n
           else
              ! need less workspace if only the eigenvalues are wanted
              lwmin = 12*n
              liwmin = 8*n
           end if
           wl = zero
           wu = zero
           iil = 0
           iiu = 0
           nsplit = 0
           if (valeig) then
              ! we do not reference vl, vu in the cases range = 'i','a'
              ! the interval (wl, wu] contains all the wanted eigenvalues.
              ! it is either given by the user or computed in la_slarre.
              wl = vl
              wu = vu
           elseif (indeig) then
              ! we do not reference il, iu in the cases range = 'v','a'
              iil = il
              iiu = iu
           end if
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (valeig .and. n > 0 .and. wu <= wl) then
              info = -7
           else if (indeig .and. (iil < 1 .or. iil > n)) then
              info = -8
           else if (indeig .and. (iiu < iil .or. iiu > n)) then
              info = -9
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -13
           else if (lwork < lwmin .and. .not. lquery) then
              info = -17
           else if (liwork < liwmin .and. .not. lquery) then
              info = -19
           end if
           ! get machine constants.
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (wantz .and. alleig) then
                 nzcmin = n
              else if (wantz .and. valeig) then
                 call la_slarrc('T',n,vl,vu,d,e,safmin,nzcmin,itmp,itmp2,info)

              else if (wantz .and. indeig) then
                 nzcmin = iiu - iil + 1
              else
                 ! wantz == false.
                 nzcmin = 0
              end if
              if (zquery .and. info == 0) then
                 z(1,1) = nzcmin
              else if (nzc < nzcmin .and. .not. zquery) then
                 info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CSTEMR',-info)
              return
           else if (lquery .or. zquery) then
              return
           end if
           ! handle n = 0, 1, and 2 cases immediately
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (wl < d(1) .and. wu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz .and. (.not. zquery)) then
                 z(1,1) = one
                 isuppz(1) = 1
                 isuppz(2) = 1
              end if
              return
           end if
           if (n == 2) then
              if (.not. wantz) then
                 call la_slae2(d(1),e(1),d(2),r1,r2)
              else if (wantz .and. (.not. zquery)) then
                 call la_slaev2(d(1),e(1),d(2),r1,r2,cs,sn)
              end if
              if (alleig .or. (valeig .and. (r2 > wl) .and. (r2 <= wu)) .or. (indeig .and. (iil == 1))) &
                        then
                 m = m + 1
                 w(m) = r2
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = -sn
                    z(2,m) = cs
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
              if (alleig .or. (valeig .and. (r1 > wl) .and. (r1 <= wu)) .or. (indeig .and. (iiu == 2))) &
                        then
                 m = m + 1
                 w(m) = r1
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = cs
                    z(2,m) = sn
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
           else
              ! continue with general n
              indgrs = 1
              inderr = 2*n + 1
              indgp = 3*n + 1
              indd = 4*n + 1
              inde2 = 5*n + 1
              indwrk = 6*n + 1
              iinspl = 1
              iindbl = n + 1
              iindw = 2*n + 1
              iindwk = 3*n + 1
              ! scale matrix to allowable range, if necessary.
              ! the allowable range is related to the pivmin parameter; see the
              ! comments in la_slarrd.  the preference for scaling small values
              ! up is heuristic; we expect users' matrices not to be close to the
              ! rmax threshold.
              scale = one
              tnrm = la_slanst('M',n,d,e)
              if (tnrm > zero .and. tnrm < rmin) then
                 scale = rmin/tnrm
              else if (tnrm > rmax) then
                 scale = rmax/tnrm
              end if
              if (scale /= one) then
                 call la_sscal(n,scale,d,1)
                 call la_sscal(n - 1,scale,e,1)
                 tnrm = tnrm*scale
                 if (valeig) then
                    ! if eigenvalues in interval have to be found,
                    ! scale (wl, wu] accordingly
                    wl = wl*scale
                    wu = wu*scale
                 end if
              end if
              ! compute the desired eigenvalues of the tridiagonal after splitting
              ! into smaller subblocks if the corresponding off-diagonal elements
              ! are small
              ! thresh is the splitting parameter for la_slarre
              ! a negative thresh forces the old splitting criterion based on the
              ! size of the off-diagonal. a positive thresh switches to splitting
              ! which preserves relative accuracy.
              if (tryrac) then
                 ! test whether the matrix warrants the more expensive relative approach.
                 call la_slarrr(n,d,e,iinfo)
              else
                 ! the user does not care about relative accurately eigenvalues
                 iinfo = -1
              end if
              ! set the splitting criterion
              if (iinfo == 0) then
                 thresh = eps
              else
                 thresh = -eps
                 ! relative accuracy is desired but t does not guarantee it
                 tryrac = .false.
              end if
              if (tryrac) then
                 ! copy original diagonal, needed to guarantee relative accuracy
                 call la_scopy(n,d,1,work(indd),1)
              end if
              ! store the squares of the offdiagonal values of t
              do j = 1,n - 1
                 work(inde2 + j - 1) = e(j)**2
              end do
              ! set the tolerance parameters for bisection
              if (.not. wantz) then
                 ! la_slarre computes the eigenvalues to full precision.
                 rtol1 = four*eps
                 rtol2 = four*eps
              else
                 ! la_slarre computes the eigenvalues to less than full precision.
                 ! la_clarrv will refine the eigenvalue approximations, and we only
                 ! need less accurate initial bisection in la_slarre.
                 ! note: these settings do only affect the subset case and la_slarre
                 rtol1 = max(sqrt(eps)*5.0e-2_sp,four*eps)
                 rtol2 = max(sqrt(eps)*5.0e-3_sp,four*eps)
              end if
              call la_slarre(range,n,wl,wu,iil,iiu,d,e,work(inde2),rtol1,rtol2, &
              thresh,nsplit,iwork(iinspl),m,w,work(inderr),work(indgp),iwork(iindbl), &
              iwork(iindw),work(indgrs),pivmin,work(indwrk),iwork(iindwk),iinfo)

              if (iinfo /= 0) then
                 info = 10 + abs(iinfo)
                 return
              end if
              ! note that if range /= 'v', la_slarre computes bounds on the desired
              ! part of the spectrum. all desired eigenvalues are contained in
              ! (wl,wu]
              if (wantz) then
                 ! compute the desired eigenvectors corresponding to the computed
                 ! eigenvalues
                 call la_clarrv(n,wl,wu,d,e,pivmin,iwork(iinspl),m,1,m,minrgp, &
                 rtol1,rtol2,w,work(inderr),work(indgp),iwork(iindbl),iwork(iindw), &
                           work(indgrs),z,ldz,isuppz,work(indwrk),iwork(iindwk),iinfo)
                 if (iinfo /= 0) then
                    info = 20 + abs(iinfo)
                    return
                 end if
              else
                 ! la_slarre computes eigenvalues of the (shifted) root representation
                 ! la_clarrv returns the eigenvalues of the unshifted matrix.
                 ! however, if the eigenvectors are not desired by the user, we need
                 ! to apply the corresponding shifts from la_slarre to obtain the
                 ! eigenvalues of the original matrix.
                 do j = 1,m
                    itmp = iwork(iindbl + j - 1)
                    w(j) = w(j) + e(iwork(iinspl + itmp - 1))
                 end do
              end if
              if (tryrac) then
                 ! refine computed eigenvalues so that they are relatively accurate
                 ! with respect to the original matrix t.
                 ibegin = 1
                 wbegin = 1
                 loop_39: do jblk = 1,iwork(iindbl + m - 1)
                    iend = iwork(iinspl + jblk - 1)
                    in = iend - ibegin + 1
                    wend = wbegin - 1
                    ! check if any eigenvalues have to be refined in this block
                    36 continue
                    if (wend < m) then
                       if (iwork(iindbl + wend) == jblk) then
                          wend = wend + 1
                          go to 36
                       end if
                    end if
                    if (wend < wbegin) then
                       ibegin = iend + 1
                       cycle loop_39
                    end if
                    offset = iwork(iindw + wbegin - 1) - 1
                    ifirst = iwork(iindw + wbegin - 1)
                    ilast = iwork(iindw + wend - 1)
                    rtol2 = four*eps
                    call la_slarrj(in,work(indd + ibegin - 1),work(inde2 + ibegin - 1),ifirst, &
                    ilast,rtol2,offset,w(wbegin),work(inderr + wbegin - 1),work(indwrk),iwork( &
                               iindwk),pivmin,tnrm,iinfo)
                    ibegin = iend + 1
                    wbegin = wend + 1
                 end do loop_39
              end if
              ! if matrix was scaled, then rescale eigenvalues appropriately.
              if (scale /= one) then
                 call la_sscal(m,one/scale,w,1)
              end if
           end if
           ! if eigenvalues are not in increasing order, then sort them,
           ! possibly along with eigenvectors.
           if (nsplit > 1 .or. n == 2) then
              if (.not. wantz) then
                 call la_slasrt('I',m,w,iinfo)
                 if (iinfo /= 0) then
                    info = 3
                    return
                 end if
              else
                 do j = 1,m - 1
                    i = 0
                    tmp = w(j)
                    do jj = j + 1,m
                       if (w(jj) < tmp) then
                          i = jj
                          tmp = w(jj)
                       end if
                    end do
                    if (i /= 0) then
                       w(i) = w(j)
                       w(j) = tmp
                       if (wantz) then
                          call la_cswap(n,z(1,i),1,z(1,j),1)
                          itmp = isuppz(2*i - 1)
                          isuppz(2*i - 1) = isuppz(2*j - 1)
                          isuppz(2*j - 1) = itmp
                          itmp = isuppz(2*i)
                          isuppz(2*i) = isuppz(2*j)
                          isuppz(2*j) = itmp
                       end if
                    end if
                 end do
              end if
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_cstemr
     !> ZSTEMR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> Depending on the number of desired eigenvalues, these are computed either
     !> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
     !> computed by the use of various suitable L D L^T factorizations near clusters
     !> of close eigenvalues (referred to as RRRs, Relatively Robust
     !> Representations). An informal sketch of the algorithm follows.
     !> For each unreduced block (submatrix) of T,
     !> (a) Compute T - sigma I  = L D L^T, so that L and D
     !> define all the wanted eigenvalues to high relative accuracy.
     !> This means that small relative changes in the entries of D and L
     !> cause only small relative changes in the eigenvalues and
     !> eigenvectors. The standard (unfactored) representation of the
     !> tridiagonal matrix T does not have this property in general.
     !> (b) Compute the eigenvalues to suitable accuracy.
     !> If the eigenvectors are desired, the algorithm attains full
     !> accuracy of the computed eigenvalues only right before
     !> the corresponding vectors have to be computed, see steps c) and d).
     !> (c) For each cluster of close eigenvalues, select a new
     !> shift close to the cluster, find a new factorization, and refine
     !> the shifted eigenvalues to suitable accuracy.
     !> (d) For each eigenvalue with a large enough relative separation compute
     !> the corresponding eigenvector by forming a rank revealing twisted
     !> factorization. Go back to (c) for any clusters that remain.
     !> For more details, see:
     !> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
     !> to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
     !> Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     !> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
     !> Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
     !> 2004.  Also LAPACK Working Note 154.
     !> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem",
     !> Computer Science Division Technical Report No. UCB/CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Further Details
     !> 1.ZSTEMR works only on machines which follow IEEE-754
     !> floating-point standard in their handling of infinities and NaNs.
     !> This permits the use of efficient inner loops avoiding a check for
     !> zero divisors.
     !> 2. LAPACK routines can be used to reduce a complex Hermitean matrix to
     !> real symmetric tridiagonal form.
     !> (Any complex Hermitean tridiagonal matrix has real values on its diagonal
     !> and potentially complex numbers on its off-diagonals. By applying a
     !> similarity transform with an appropriate diagonal matrix
     !> diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex Hermitean
     !> matrix can be transformed into a real symmetric matrix and complex
     !> arithmetic can be entirely avoided.)
     !> While the eigenvectors of the real symmetric tridiagonal matrix are real,
     !> the eigenvectors of original complex Hermitean matrix have complex entries
     !> in general.
     !> Since LAPACK drivers overwrite the matrix data with the eigenvectors,
     !> ZSTEMR accepts complex workspace to facilitate interoperability
     !> with ZUNMTR or ZUPMTR.

     pure subroutine la_zstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc, &
               isuppz,tryrac,work,lwork,iwork,liwork,info)
        use la_constants_dp,only:zero,one,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           logical(lk),intent(inout) :: tryrac
           integer(ilp),intent(in) :: il,iu,ldz,nzc,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*)
           complex(dp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: minrgp = 1.0e-3_dp

           ! Local Scalars
           logical(lk) :: alleig,indeig,lquery,valeig,wantz,zquery
           integer(ilp) :: i,ibegin,iend,ifirst,iil,iindbl,iindw,iindwk,iinfo,iinspl, &
           iiu,ilast,in,indd,inde2,inderr,indgp,indgrs,indwrk,itmp,itmp2,j,jblk,jj, &
                     liwmin,lwmin,nsplit,nzcmin,offset,wbegin,wend
           real(dp) :: bignum,cs,eps,pivmin,r1,r2,rmax,rmin,rtol1,rtol2,safmin,scale, &
                     smlnum,sn,thresh,tmp,tnrm,wl,wu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           zquery = (nzc == -1)
           ! la_dstemr needs work of size 6*n, iwork of size 3*n.
           ! in addition, la_dlarre needs work of size 6*n, iwork of size 5*n.
           ! furthermore, la_zlarrv needs work of size 12*n, iwork of size 7*n.
           if (wantz) then
              lwmin = 18*n
              liwmin = 10*n
           else
              ! need less workspace if only the eigenvalues are wanted
              lwmin = 12*n
              liwmin = 8*n
           end if
           wl = zero
           wu = zero
           iil = 0
           iiu = 0
           nsplit = 0
           if (valeig) then
              ! we do not reference vl, vu in the cases range = 'i','a'
              ! the interval (wl, wu] contains all the wanted eigenvalues.
              ! it is either given by the user or computed in la_dlarre.
              wl = vl
              wu = vu
           elseif (indeig) then
              ! we do not reference il, iu in the cases range = 'v','a'
              iil = il
              iiu = iu
           end if
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (valeig .and. n > 0 .and. wu <= wl) then
              info = -7
           else if (indeig .and. (iil < 1 .or. iil > n)) then
              info = -8
           else if (indeig .and. (iiu < iil .or. iiu > n)) then
              info = -9
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -13
           else if (lwork < lwmin .and. .not. lquery) then
              info = -17
           else if (liwork < liwmin .and. .not. lquery) then
              info = -19
           end if
           ! get machine constants.
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (wantz .and. alleig) then
                 nzcmin = n
              else if (wantz .and. valeig) then
                 call la_dlarrc('T',n,vl,vu,d,e,safmin,nzcmin,itmp,itmp2,info)

              else if (wantz .and. indeig) then
                 nzcmin = iiu - iil + 1
              else
                 ! wantz == false.
                 nzcmin = 0
              end if
              if (zquery .and. info == 0) then
                 z(1,1) = nzcmin
              else if (nzc < nzcmin .and. .not. zquery) then
                 info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZSTEMR',-info)
              return
           else if (lquery .or. zquery) then
              return
           end if
           ! handle n = 0, 1, and 2 cases immediately
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (wl < d(1) .and. wu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz .and. (.not. zquery)) then
                 z(1,1) = one
                 isuppz(1) = 1
                 isuppz(2) = 1
              end if
              return
           end if
           if (n == 2) then
              if (.not. wantz) then
                 call la_dlae2(d(1),e(1),d(2),r1,r2)
              else if (wantz .and. (.not. zquery)) then
                 call la_dlaev2(d(1),e(1),d(2),r1,r2,cs,sn)
              end if
              if (alleig .or. (valeig .and. (r2 > wl) .and. (r2 <= wu)) .or. (indeig .and. (iil == 1))) &
                        then
                 m = m + 1
                 w(m) = r2
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = -sn
                    z(2,m) = cs
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
              if (alleig .or. (valeig .and. (r1 > wl) .and. (r1 <= wu)) .or. (indeig .and. (iiu == 2))) &
                        then
                 m = m + 1
                 w(m) = r1
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = cs
                    z(2,m) = sn
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
           else
              ! continue with general n
              indgrs = 1
              inderr = 2*n + 1
              indgp = 3*n + 1
              indd = 4*n + 1
              inde2 = 5*n + 1
              indwrk = 6*n + 1
              iinspl = 1
              iindbl = n + 1
              iindw = 2*n + 1
              iindwk = 3*n + 1
              ! scale matrix to allowable range, if necessary.
              ! the allowable range is related to the pivmin parameter; see the
              ! comments in la_dlarrd.  the preference for scaling small values
              ! up is heuristic; we expect users' matrices not to be close to the
              ! rmax threshold.
              scale = one
              tnrm = la_dlanst('M',n,d,e)
              if (tnrm > zero .and. tnrm < rmin) then
                 scale = rmin/tnrm
              else if (tnrm > rmax) then
                 scale = rmax/tnrm
              end if
              if (scale /= one) then
                 call la_dscal(n,scale,d,1)
                 call la_dscal(n - 1,scale,e,1)
                 tnrm = tnrm*scale
                 if (valeig) then
                    ! if eigenvalues in interval have to be found,
                    ! scale (wl, wu] accordingly
                    wl = wl*scale
                    wu = wu*scale
                 end if
              end if
              ! compute the desired eigenvalues of the tridiagonal after splitting
              ! into smaller subblocks if the corresponding off-diagonal elements
              ! are small
              ! thresh is the splitting parameter for la_dlarre
              ! a negative thresh forces the old splitting criterion based on the
              ! size of the off-diagonal. a positive thresh switches to splitting
              ! which preserves relative accuracy.
              if (tryrac) then
                 ! test whether the matrix warrants the more expensive relative approach.
                 call la_dlarrr(n,d,e,iinfo)
              else
                 ! the user does not care about relative accurately eigenvalues
                 iinfo = -1
              end if
              ! set the splitting criterion
              if (iinfo == 0) then
                 thresh = eps
              else
                 thresh = -eps
                 ! relative accuracy is desired but t does not guarantee it
                 tryrac = .false.
              end if
              if (tryrac) then
                 ! copy original diagonal, needed to guarantee relative accuracy
                 call la_dcopy(n,d,1,work(indd),1)
              end if
              ! store the squares of the offdiagonal values of t
              do j = 1,n - 1
                 work(inde2 + j - 1) = e(j)**2
              end do
              ! set the tolerance parameters for bisection
              if (.not. wantz) then
                 ! la_dlarre computes the eigenvalues to full precision.
                 rtol1 = four*eps
                 rtol2 = four*eps
              else
                 ! la_dlarre computes the eigenvalues to less than full precision.
                 ! la_zlarrv will refine the eigenvalue approximations, and we only
                 ! need less accurate initial bisection in la_dlarre.
                 ! note: these settings do only affect the subset case and la_dlarre
                 rtol1 = sqrt(eps)
                 rtol2 = max(sqrt(eps)*5.0e-3_dp,four*eps)
              end if
              call la_dlarre(range,n,wl,wu,iil,iiu,d,e,work(inde2),rtol1,rtol2, &
              thresh,nsplit,iwork(iinspl),m,w,work(inderr),work(indgp),iwork(iindbl), &
              iwork(iindw),work(indgrs),pivmin,work(indwrk),iwork(iindwk),iinfo)

              if (iinfo /= 0) then
                 info = 10 + abs(iinfo)
                 return
              end if
              ! note that if range /= 'v', la_dlarre computes bounds on the desired
              ! part of the spectrum. all desired eigenvalues are contained in
              ! (wl,wu]
              if (wantz) then
                 ! compute the desired eigenvectors corresponding to the computed
                 ! eigenvalues
                 call la_zlarrv(n,wl,wu,d,e,pivmin,iwork(iinspl),m,1,m,minrgp, &
                 rtol1,rtol2,w,work(inderr),work(indgp),iwork(iindbl),iwork(iindw), &
                           work(indgrs),z,ldz,isuppz,work(indwrk),iwork(iindwk),iinfo)
                 if (iinfo /= 0) then
                    info = 20 + abs(iinfo)
                    return
                 end if
              else
                 ! la_dlarre computes eigenvalues of the (shifted) root representation
                 ! la_zlarrv returns the eigenvalues of the unshifted matrix.
                 ! however, if the eigenvectors are not desired by the user, we need
                 ! to apply the corresponding shifts from la_dlarre to obtain the
                 ! eigenvalues of the original matrix.
                 do j = 1,m
                    itmp = iwork(iindbl + j - 1)
                    w(j) = w(j) + e(iwork(iinspl + itmp - 1))
                 end do
              end if
              if (tryrac) then
                 ! refine computed eigenvalues so that they are relatively accurate
                 ! with respect to the original matrix t.
                 ibegin = 1
                 wbegin = 1
                 loop_39: do jblk = 1,iwork(iindbl + m - 1)
                    iend = iwork(iinspl + jblk - 1)
                    in = iend - ibegin + 1
                    wend = wbegin - 1
                    ! check if any eigenvalues have to be refined in this block
                    36 continue
                    if (wend < m) then
                       if (iwork(iindbl + wend) == jblk) then
                          wend = wend + 1
                          go to 36
                       end if
                    end if
                    if (wend < wbegin) then
                       ibegin = iend + 1
                       cycle loop_39
                    end if
                    offset = iwork(iindw + wbegin - 1) - 1
                    ifirst = iwork(iindw + wbegin - 1)
                    ilast = iwork(iindw + wend - 1)
                    rtol2 = four*eps
                    call la_dlarrj(in,work(indd + ibegin - 1),work(inde2 + ibegin - 1),ifirst, &
                    ilast,rtol2,offset,w(wbegin),work(inderr + wbegin - 1),work(indwrk),iwork( &
                               iindwk),pivmin,tnrm,iinfo)
                    ibegin = iend + 1
                    wbegin = wend + 1
                 end do loop_39
              end if
              ! if matrix was scaled, then rescale eigenvalues appropriately.
              if (scale /= one) then
                 call la_dscal(m,one/scale,w,1)
              end if
           end if
           ! if eigenvalues are not in increasing order, then sort them,
           ! possibly along with eigenvectors.
           if (nsplit > 1 .or. n == 2) then
              if (.not. wantz) then
                 call la_dlasrt('I',m,w,iinfo)
                 if (iinfo /= 0) then
                    info = 3
                    return
                 end if
              else
                 do j = 1,m - 1
                    i = 0
                    tmp = w(j)
                    do jj = j + 1,m
                       if (w(jj) < tmp) then
                          i = jj
                          tmp = w(jj)
                       end if
                    end do
                    if (i /= 0) then
                       w(i) = w(j)
                       w(j) = tmp
                       if (wantz) then
                          call la_zswap(n,z(1,i),1,z(1,j),1)
                          itmp = isuppz(2*i - 1)
                          isuppz(2*i - 1) = isuppz(2*j - 1)
                          isuppz(2*j - 1) = itmp
                          itmp = isuppz(2*i)
                          isuppz(2*i) = isuppz(2*j)
                          isuppz(2*j) = itmp
                       end if
                    end if
                 end do
              end if
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_zstemr
     !> WSTEMR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> Depending on the number of desired eigenvalues, these are computed either
     !> by bisection or the dqds algorithm. Numerically orthogonal eigenvectors are
     !> computed by the use of various suitable L D L^T factorizations near clusters
     !> of close eigenvalues (referred to as RRRs, Relatively Robust
     !> Representations). An informal sketch of the algorithm follows.
     !> For each unreduced block (submatrix) of T,
     !> (a) Compute T - sigma I  = L D L^T, so that L and D
     !> define all the wanted eigenvalues to high relative accuracy.
     !> This means that small relative changes in the entries of D and L
     !> cause only small relative changes in the eigenvalues and
     !> eigenvectors. The standard (unfactored) representation of the
     !> tridiagonal matrix T does not have this property in general.
     !> (b) Compute the eigenvalues to suitable accuracy.
     !> If the eigenvectors are desired, the algorithm attains full
     !> accuracy of the computed eigenvalues only right before
     !> the corresponding vectors have to be computed, see steps c) and d).
     !> (c) For each cluster of close eigenvalues, select a new
     !> shift close to the cluster, find a new factorization, and refine
     !> the shifted eigenvalues to suitable accuracy.
     !> (d) For each eigenvalue with a large enough relative separation compute
     !> the corresponding eigenvector by forming a rank revealing twisted
     !> factorization. Go back to (c) for any clusters that remain.
     !> For more details, see:
     !> - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
     !> to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
     !> Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
     !> - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
     !> Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
     !> 2004.  Also LAPACK Working Note 154.
     !> - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
     !> tridiagonal eigenvalue/eigenvector problem",
     !> Computer Science Division Technical Report No. UCB/CSD-97-971,
     !> UC Berkeley, May 1997.
     !> Further Details
     !> 1.WSTEMR works only on machines which follow IEEE-754
     !> floating-point standard in their handling of infinities and NaNs.
     !> This permits the use of efficient inner loops avoiding a check for
     !> zero divisors.
     !> 2. LAPACK routines can be used to reduce a complex Hermitean matrix to
     !> real symmetric tridiagonal form.
     !> (Any complex Hermitean tridiagonal matrix has real values on its diagonal
     !> and potentially complex numbers on its off-diagonals. By applying a
     !> similarity transform with an appropriate diagonal matrix
     !> diag(1,e^{i \phy_1}, ... , e^{i \phy_{n-1}}), the complex Hermitean
     !> matrix can be transformed into a real symmetric matrix and complex
     !> arithmetic can be entirely avoided.)
     !> While the eigenvectors of the real symmetric tridiagonal matrix are real,
     !> the eigenvectors of original complex Hermitean matrix have complex entries
     !> in general.
     !> Since LAPACK drivers overwrite the matrix data with the eigenvectors,
     !> WSTEMR accepts complex workspace to facilitate interoperability
     !> with WUNMTR or WUPMTR.

     pure subroutine la_wstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,nzc, &
               isuppz,tryrac,work,lwork,iwork,liwork,info)
        use la_constants_qp,only:zero,one,four
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           logical(lk),intent(inout) :: tryrac
           integer(ilp),intent(in) :: il,iu,ldz,nzc,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*)
           complex(qp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: minrgp = 1.0e-3_qp

           ! Local Scalars
           logical(lk) :: alleig,indeig,lquery,valeig,wantz,zquery
           integer(ilp) :: i,ibegin,iend,ifirst,iil,iindbl,iindw,iindwk,iinfo,iinspl, &
           iiu,ilast,in,indd,inde2,inderr,indgp,indgrs,indwrk,itmp,itmp2,j,jblk,jj, &
                     liwmin,lwmin,nsplit,nzcmin,offset,wbegin,wend
           real(qp) :: bignum,cs,eps,pivmin,r1,r2,rmax,rmin,rtol1,rtol2,safmin,scale, &
                     smlnum,sn,thresh,tmp,tnrm,wl,wu
           ! Intrinsic Functions
           intrinsic :: max,min,sqrt
           ! Executable Statements
           ! test the input parameters.
           wantz = la_lsame(jobz,'V')
           alleig = la_lsame(range,'A')
           valeig = la_lsame(range,'V')
           indeig = la_lsame(range,'I')
           lquery = ((lwork == -1) .or. (liwork == -1))
           zquery = (nzc == -1)
           ! la_qstemr needs work of size 6*n, iwork of size 3*n.
           ! in addition, la_qlarre needs work of size 6*n, iwork of size 5*n.
           ! furthermore, la_wlarrv needs work of size 12*n, iwork of size 7*n.
           if (wantz) then
              lwmin = 18*n
              liwmin = 10*n
           else
              ! need less workspace if only the eigenvalues are wanted
              lwmin = 12*n
              liwmin = 8*n
           end if
           wl = zero
           wu = zero
           iil = 0
           iiu = 0
           nsplit = 0
           if (valeig) then
              ! we do not reference vl, vu in the cases range = 'i','a'
              ! the interval (wl, wu] contains all the wanted eigenvalues.
              ! it is either given by the user or computed in la_qlarre.
              wl = vl
              wu = vu
           elseif (indeig) then
              ! we do not reference il, iu in the cases range = 'v','a'
              iil = il
              iiu = iu
           end if
           info = 0
           if (.not. (wantz .or. la_lsame(jobz,'N'))) then
              info = -1
           else if (.not. (alleig .or. valeig .or. indeig)) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (valeig .and. n > 0 .and. wu <= wl) then
              info = -7
           else if (indeig .and. (iil < 1 .or. iil > n)) then
              info = -8
           else if (indeig .and. (iiu < iil .or. iiu > n)) then
              info = -9
           else if (ldz < 1 .or. (wantz .and. ldz < n)) then
              info = -13
           else if (lwork < lwmin .and. .not. lquery) then
              info = -17
           else if (liwork < liwmin .and. .not. lquery) then
              info = -19
           end if
           ! get machine constants.
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           bignum = one/smlnum
           rmin = sqrt(smlnum)
           rmax = min(sqrt(bignum),one/sqrt(sqrt(safmin)))
           if (info == 0) then
              work(1) = lwmin
              iwork(1) = liwmin
              if (wantz .and. alleig) then
                 nzcmin = n
              else if (wantz .and. valeig) then
                 call la_qlarrc('T',n,vl,vu,d,e,safmin,nzcmin,itmp,itmp2,info)

              else if (wantz .and. indeig) then
                 nzcmin = iiu - iil + 1
              else
                 ! wantz == false.
                 nzcmin = 0
              end if
              if (zquery .and. info == 0) then
                 z(1,1) = nzcmin
              else if (nzc < nzcmin .and. .not. zquery) then
                 info = -14
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WSTEMR',-info)
              return
           else if (lquery .or. zquery) then
              return
           end if
           ! handle n = 0, 1, and 2 cases immediately
           m = 0
           if (n == 0) return
           if (n == 1) then
              if (alleig .or. indeig) then
                 m = 1
                 w(1) = d(1)
              else
                 if (wl < d(1) .and. wu >= d(1)) then
                    m = 1
                    w(1) = d(1)
                 end if
              end if
              if (wantz .and. (.not. zquery)) then
                 z(1,1) = one
                 isuppz(1) = 1
                 isuppz(2) = 1
              end if
              return
           end if
           if (n == 2) then
              if (.not. wantz) then
                 call la_qlae2(d(1),e(1),d(2),r1,r2)
              else if (wantz .and. (.not. zquery)) then
                 call la_qlaev2(d(1),e(1),d(2),r1,r2,cs,sn)
              end if
              if (alleig .or. (valeig .and. (r2 > wl) .and. (r2 <= wu)) .or. (indeig .and. (iil == 1))) &
                        then
                 m = m + 1
                 w(m) = r2
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = -sn
                    z(2,m) = cs
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
              if (alleig .or. (valeig .and. (r1 > wl) .and. (r1 <= wu)) .or. (indeig .and. (iiu == 2))) &
                        then
                 m = m + 1
                 w(m) = r1
                 if (wantz .and. (.not. zquery)) then
                    z(1,m) = cs
                    z(2,m) = sn
                    ! note: at most one of sn and cs can be zero.
                    if (sn /= zero) then
                       if (cs /= zero) then
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 2
                       else
                          isuppz(2*m - 1) = 1
                          isuppz(2*m) = 1
                       end if
                    else
                       isuppz(2*m - 1) = 2
                       isuppz(2*m) = 2
                    end if
                 end if
              end if
           else
              ! continue with general n
              indgrs = 1
              inderr = 2*n + 1
              indgp = 3*n + 1
              indd = 4*n + 1
              inde2 = 5*n + 1
              indwrk = 6*n + 1
              iinspl = 1
              iindbl = n + 1
              iindw = 2*n + 1
              iindwk = 3*n + 1
              ! scale matrix to allowable range, if necessary.
              ! the allowable range is related to the pivmin parameter; see the
              ! comments in la_qlarrd.  the preference for scaling small values
              ! up is heuristic; we expect users' matrices not to be close to the
              ! rmax threshold.
              scale = one
              tnrm = la_qlanst('M',n,d,e)
              if (tnrm > zero .and. tnrm < rmin) then
                 scale = rmin/tnrm
              else if (tnrm > rmax) then
                 scale = rmax/tnrm
              end if
              if (scale /= one) then
                 call la_qscal(n,scale,d,1)
                 call la_qscal(n - 1,scale,e,1)
                 tnrm = tnrm*scale
                 if (valeig) then
                    ! if eigenvalues in interval have to be found,
                    ! scale (wl, wu] accordingly
                    wl = wl*scale
                    wu = wu*scale
                 end if
              end if
              ! compute the desired eigenvalues of the tridiagonal after splitting
              ! into smaller subblocks if the corresponding off-diagonal elements
              ! are small
              ! thresh is the splitting parameter for la_qlarre
              ! a negative thresh forces the old splitting criterion based on the
              ! size of the off-diagonal. a positive thresh switches to splitting
              ! which preserves relative accuracy.
              if (tryrac) then
                 ! test whether the matrix warrants the more expensive relative approach.
                 call la_qlarrr(n,d,e,iinfo)
              else
                 ! the user does not care about relative accurately eigenvalues
                 iinfo = -1
              end if
              ! set the splitting criterion
              if (iinfo == 0) then
                 thresh = eps
              else
                 thresh = -eps
                 ! relative accuracy is desired but t does not guarantee it
                 tryrac = .false.
              end if
              if (tryrac) then
                 ! copy original diagonal, needed to guarantee relative accuracy
                 call la_qcopy(n,d,1,work(indd),1)
              end if
              ! store the squares of the offdiagonal values of t
              do j = 1,n - 1
                 work(inde2 + j - 1) = e(j)**2
              end do
              ! set the tolerance parameters for bisection
              if (.not. wantz) then
                 ! la_qlarre computes the eigenvalues to full precision.
                 rtol1 = four*eps
                 rtol2 = four*eps
              else
                 ! la_qlarre computes the eigenvalues to less than full precision.
                 ! la_wlarrv will refine the eigenvalue approximations, and we only
                 ! need less accurate initial bisection in la_qlarre.
                 ! note: these settings do only affect the subset case and la_qlarre
                 rtol1 = sqrt(eps)
                 rtol2 = max(sqrt(eps)*5.0e-3_qp,four*eps)
              end if
              call la_qlarre(range,n,wl,wu,iil,iiu,d,e,work(inde2),rtol1,rtol2, &
              thresh,nsplit,iwork(iinspl),m,w,work(inderr),work(indgp),iwork(iindbl), &
              iwork(iindw),work(indgrs),pivmin,work(indwrk),iwork(iindwk),iinfo)

              if (iinfo /= 0) then
                 info = 10 + abs(iinfo)
                 return
              end if
              ! note that if range /= 'v', la_qlarre computes bounds on the desired
              ! part of the spectrum. all desired eigenvalues are contained in
              ! (wl,wu]
              if (wantz) then
                 ! compute the desired eigenvectors corresponding to the computed
                 ! eigenvalues
                 call la_wlarrv(n,wl,wu,d,e,pivmin,iwork(iinspl),m,1,m,minrgp, &
                 rtol1,rtol2,w,work(inderr),work(indgp),iwork(iindbl),iwork(iindw), &
                           work(indgrs),z,ldz,isuppz,work(indwrk),iwork(iindwk),iinfo)
                 if (iinfo /= 0) then
                    info = 20 + abs(iinfo)
                    return
                 end if
              else
                 ! la_qlarre computes eigenvalues of the (shifted) root representation
                 ! la_wlarrv returns the eigenvalues of the unshifted matrix.
                 ! however, if the eigenvectors are not desired by the user, we need
                 ! to apply the corresponding shifts from la_qlarre to obtain the
                 ! eigenvalues of the original matrix.
                 do j = 1,m
                    itmp = iwork(iindbl + j - 1)
                    w(j) = w(j) + e(iwork(iinspl + itmp - 1))
                 end do
              end if
              if (tryrac) then
                 ! refine computed eigenvalues so that they are relatively accurate
                 ! with respect to the original matrix t.
                 ibegin = 1
                 wbegin = 1
                 loop_39: do jblk = 1,iwork(iindbl + m - 1)
                    iend = iwork(iinspl + jblk - 1)
                    in = iend - ibegin + 1
                    wend = wbegin - 1
                    ! check if any eigenvalues have to be refined in this block
                    36 continue
                    if (wend < m) then
                       if (iwork(iindbl + wend) == jblk) then
                          wend = wend + 1
                          go to 36
                       end if
                    end if
                    if (wend < wbegin) then
                       ibegin = iend + 1
                       cycle loop_39
                    end if
                    offset = iwork(iindw + wbegin - 1) - 1
                    ifirst = iwork(iindw + wbegin - 1)
                    ilast = iwork(iindw + wend - 1)
                    rtol2 = four*eps
                    call la_qlarrj(in,work(indd + ibegin - 1),work(inde2 + ibegin - 1),ifirst, &
                    ilast,rtol2,offset,w(wbegin),work(inderr + wbegin - 1),work(indwrk),iwork( &
                               iindwk),pivmin,tnrm,iinfo)
                    ibegin = iend + 1
                    wbegin = wend + 1
                 end do loop_39
              end if
              ! if matrix was scaled, then rescale eigenvalues appropriately.
              if (scale /= one) then
                 call la_qscal(m,one/scale,w,1)
              end if
           end if
           ! if eigenvalues are not in increasing order, then sort them,
           ! possibly along with eigenvectors.
           if (nsplit > 1 .or. n == 2) then
              if (.not. wantz) then
                 call la_qlasrt('I',m,w,iinfo)
                 if (iinfo /= 0) then
                    info = 3
                    return
                 end if
              else
                 do j = 1,m - 1
                    i = 0
                    tmp = w(j)
                    do jj = j + 1,m
                       if (w(jj) < tmp) then
                          i = jj
                          tmp = w(jj)
                       end if
                    end do
                    if (i /= 0) then
                       w(i) = w(j)
                       w(j) = tmp
                       if (wantz) then
                          call la_wswap(n,z(1,i),1,z(1,j),1)
                          itmp = isuppz(2*i - 1)
                          isuppz(2*i - 1) = isuppz(2*j - 1)
                          isuppz(2*j - 1) = itmp
                          itmp = isuppz(2*i)
                          isuppz(2*i) = isuppz(2*j)
                          isuppz(2*j) = itmp
                       end if
                    end if
                 end do
              end if
           end if
           work(1) = lwmin
           iwork(1) = liwmin
           return
     end subroutine la_wstemr

     !> CSTEDC: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if CHETRD or CHPTRD or CHBTRD has been used to reduce this
     !> matrix to tridiagonal form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See SLAED3 for details.

     pure subroutine la_cstedc(compz,n,d,e,z,ldz,work,lwork,rwork,lrwork,iwork, &
               liwork,info)
        use la_constants_sp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lrwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: rwork(*)
           complex(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: finish,i,icompz,ii,j,k,lgn,liwmin,ll,lrwmin,lwmin,m,smlsiz, &
                      start
           real(sp) :: eps,orgnrm,p,tiny
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,mod,real,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1 .or. lrwork == -1 .or. liwork == -1)
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
           if (info == 0) then
              ! compute the workspace requirements
              smlsiz = la_ilaenv(9,'CSTEDC',' ',0,0,0,0)
              if (n <= 1 .or. icompz == 0) then
                 lwmin = 1
                 liwmin = 1
                 lrwmin = 1
              else if (n <= smlsiz) then
                 lwmin = 1
                 liwmin = 1
                 lrwmin = 2*(n - 1)
              else if (icompz == 1) then
                 lgn = int(log(real(n,KIND=sp))/log(two),KIND=ilp)
                 if (2**lgn < n) lgn = lgn + 1
                 if (2**lgn < n) lgn = lgn + 1
                 lwmin = n*n
                 lrwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                 liwmin = 6 + 6*n + 5*n*lgn
              else if (icompz == 2) then
                 lwmin = 1
                 lrwmin = 1 + 4*n + 2*n**2
                 liwmin = 3 + 5*n
              end if
              work(1) = lwmin
              rwork(1) = lrwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (lrwork < lrwmin .and. .not. lquery) then
                 info = -10
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('CSTEDC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz /= 0) z(1,1) = one
              return
           end if
           ! if the following conditional clause is removed, then the routine
           ! will use the divide and conquer routine to compute only the
           ! eigenvalues, which requires (3n + 3n**2) real workspace and
           ! (2 + 5n + 2n lg(n)) integer workspace.
           ! since on many architectures la_ssterf is much faster than any other
           ! algorithm for finding eigenvalues only, it is used here
           ! as the default. if the conditional clause is removed, then
           ! information on the size of workspace needs to be changed.
           ! if compz = 'n', use la_ssterf to compute the eigenvalues.
           if (icompz == 0) then
              call la_ssterf(n,d,e,info)
              go to 70
           end if
           ! if n is smaller than the minimum divide size (smlsiz+1), then
           ! solve the problem with another solver.
           if (n <= smlsiz) then
              call la_csteqr(compz,n,d,e,z,ldz,rwork,info)
           else
              ! if compz = 'i', we simply call la_sstedc instead.
              if (icompz == 2) then
                 call la_slaset('FULL',n,n,zero,one,rwork,n)
                 ll = n*n + 1
                 call la_sstedc('I',n,d,e,rwork,n,rwork(ll),lrwork - ll + 1,iwork, &
                           liwork,info)
                 do j = 1,n
                    do i = 1,n
                       z(i,j) = rwork((j - 1)*n + i)
                    end do
                 end do
                 go to 70
              end if
              ! from now on, only option left to be handled is compz = 'v',
              ! i.e. icompz = 1.
              ! scale.
              orgnrm = la_slanst('M',n,d,e)
              if (orgnrm == zero) go to 70
              eps = la_slamch('EPSILON')
              start = 1
              ! while ( start <= n )
              30 continue
              if (start <= n) then
                 ! let finish be the position of the next subdiagonal entry
                 ! such that e( finish ) <= tiny or finish = n if no such
                 ! subdiagonal exists.  the matrix identified by the elements
                 ! between start and finish constitutes an independent
                 ! sub-problem.
                 finish = start
                 40 continue
                 if (finish < n) then
                    tiny = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                    if (abs(e(finish)) > tiny) then
                       finish = finish + 1
                       go to 40
                    end if
                 end if
                 ! (sub) problem determined.  compute its size and solve it.
                 m = finish - start + 1
                 if (m > smlsiz) then
                    ! scale.
                    orgnrm = la_slanst('M',m,d(start),e(start))
                    call la_slascl('G',0,0,orgnrm,one,m,1,d(start),m,info)
                    call la_slascl('G',0,0,orgnrm,one,m - 1,1,e(start),m - 1,info)

                    call la_claed0(n,m,d(start),e(start),z(1,start),ldz,work,n, &
                              rwork,iwork,info)
                    if (info > 0) then
                       info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, (m + 1)) + start - &
                                 1
                       go to 70
                    end if
                    ! scale back.
                    call la_slascl('G',0,0,one,orgnrm,m,1,d(start),m,info)
                 else
                    call la_ssteqr('I',m,d(start),e(start),rwork,m,rwork(m*m + 1), &
                              info)
                    call la_clacrm(n,m,z(1,start),ldz,rwork,m,work,n,rwork(m*m + 1) &
                               )
                    call la_clacpy('A',n,m,work,n,z(1,start),ldz)
                    if (info > 0) then
                       info = start*(n + 1) + finish
                       go to 70
                    end if
                 end if
                 start = finish + 1
                 go to 30
              end if
              ! endwhile
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
           70 continue
           work(1) = lwmin
           rwork(1) = lrwmin
           iwork(1) = liwmin
           return
     end subroutine la_cstedc
     !> ZSTEDC: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
     !> matrix to tridiagonal form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See DLAED3 for details.

     pure subroutine la_zstedc(compz,n,d,e,z,ldz,work,lwork,rwork,lrwork,iwork, &
               liwork,info)
        use la_constants_dp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lrwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: rwork(*)
           complex(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: finish,i,icompz,ii,j,k,lgn,liwmin,ll,lrwmin,lwmin,m,smlsiz, &
                      start
           real(dp) :: eps,orgnrm,p,tiny
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max,mod,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1 .or. lrwork == -1 .or. liwork == -1)
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
           if (info == 0) then
              ! compute the workspace requirements
              smlsiz = la_ilaenv(9,'ZSTEDC',' ',0,0,0,0)
              if (n <= 1 .or. icompz == 0) then
                 lwmin = 1
                 liwmin = 1
                 lrwmin = 1
              else if (n <= smlsiz) then
                 lwmin = 1
                 liwmin = 1
                 lrwmin = 2*(n - 1)
              else if (icompz == 1) then
                 lgn = int(log(real(n,KIND=dp))/log(two),KIND=ilp)
                 if (2**lgn < n) lgn = lgn + 1
                 if (2**lgn < n) lgn = lgn + 1
                 lwmin = n*n
                 lrwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                 liwmin = 6 + 6*n + 5*n*lgn
              else if (icompz == 2) then
                 lwmin = 1
                 lrwmin = 1 + 4*n + 2*n**2
                 liwmin = 3 + 5*n
              end if
              work(1) = lwmin
              rwork(1) = lrwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (lrwork < lrwmin .and. .not. lquery) then
                 info = -10
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('ZSTEDC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz /= 0) z(1,1) = one
              return
           end if
           ! if the following conditional clause is removed, then the routine
           ! will use the divide and conquer routine to compute only the
           ! eigenvalues, which requires (3n + 3n**2) real workspace and
           ! (2 + 5n + 2n lg(n)) integer workspace.
           ! since on many architectures la_dsterf is much faster than any other
           ! algorithm for finding eigenvalues only, it is used here
           ! as the default. if the conditional clause is removed, then
           ! information on the size of workspace needs to be changed.
           ! if compz = 'n', use la_dsterf to compute the eigenvalues.
           if (icompz == 0) then
              call la_dsterf(n,d,e,info)
              go to 70
           end if
           ! if n is smaller than the minimum divide size (smlsiz+1), then
           ! solve the problem with another solver.
           if (n <= smlsiz) then
              call la_zsteqr(compz,n,d,e,z,ldz,rwork,info)
           else
              ! if compz = 'i', we simply call la_dstedc instead.
              if (icompz == 2) then
                 call la_dlaset('FULL',n,n,zero,one,rwork,n)
                 ll = n*n + 1
                 call la_dstedc('I',n,d,e,rwork,n,rwork(ll),lrwork - ll + 1,iwork, &
                           liwork,info)
                 do j = 1,n
                    do i = 1,n
                       z(i,j) = rwork((j - 1)*n + i)
                    end do
                 end do
                 go to 70
              end if
              ! from now on, only option left to be handled is compz = 'v',
              ! i.e. icompz = 1.
              ! scale.
              orgnrm = la_dlanst('M',n,d,e)
              if (orgnrm == zero) go to 70
              eps = la_dlamch('EPSILON')
              start = 1
              ! while ( start <= n )
              30 continue
              if (start <= n) then
                 ! let finish be the position of the next subdiagonal entry
                 ! such that e( finish ) <= tiny or finish = n if no such
                 ! subdiagonal exists.  the matrix identified by the elements
                 ! between start and finish constitutes an independent
                 ! sub-problem.
                 finish = start
                 40 continue
                 if (finish < n) then
                    tiny = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                    if (abs(e(finish)) > tiny) then
                       finish = finish + 1
                       go to 40
                    end if
                 end if
                 ! (sub) problem determined.  compute its size and solve it.
                 m = finish - start + 1
                 if (m > smlsiz) then
                    ! scale.
                    orgnrm = la_dlanst('M',m,d(start),e(start))
                    call la_dlascl('G',0,0,orgnrm,one,m,1,d(start),m,info)
                    call la_dlascl('G',0,0,orgnrm,one,m - 1,1,e(start),m - 1,info)

                    call la_zlaed0(n,m,d(start),e(start),z(1,start),ldz,work,n, &
                              rwork,iwork,info)
                    if (info > 0) then
                       info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, (m + 1)) + start - &
                                 1
                       go to 70
                    end if
                    ! scale back.
                    call la_dlascl('G',0,0,one,orgnrm,m,1,d(start),m,info)
                 else
                    call la_dsteqr('I',m,d(start),e(start),rwork,m,rwork(m*m + 1), &
                              info)
                    call la_zlacrm(n,m,z(1,start),ldz,rwork,m,work,n,rwork(m*m + 1) &
                               )
                    call la_zlacpy('A',n,m,work,n,z(1,start),ldz)
                    if (info > 0) then
                       info = start*(n + 1) + finish
                       go to 70
                    end if
                 end if
                 start = finish + 1
                 go to 30
              end if
              ! endwhile
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
           70 continue
           work(1) = lwmin
           rwork(1) = lrwmin
           iwork(1) = liwmin
           return
     end subroutine la_zstedc
     !> WSTEDC: computes all eigenvalues and, optionally, eigenvectors of a
     !> symmetric tridiagonal matrix using the divide and conquer method.
     !> The eigenvectors of a full or band complex Hermitian matrix can also
     !> be found if WHETRD or WHPTRD or WHBTRD has been used to reduce this
     !> matrix to tridiagonal form.
     !> This code makes very mild assumptions about floating point
     !> arithmetic. It will work on machines with a guard digit in
     !> add/subtract, or on those binary machines without guard digits
     !> which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
     !> It could conceivably fail on hexadecimal or decimal machines
     !> without guard digits, but we know of none.  See QLAED3 for details.

     pure subroutine la_wstedc(compz,n,d,e,z,ldz,work,lwork,rwork,lrwork,iwork, &
               liwork,info)
        use la_constants_qp,only:zero,one,two
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: compz
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: ldz,liwork,lrwork,lwork,n
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: rwork(*)
           complex(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: z(ldz,*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: lquery
           integer(ilp) :: finish,i,icompz,ii,j,k,lgn,liwmin,ll,lrwmin,lwmin,m,smlsiz, &
                      start
           real(qp) :: eps,orgnrm,p,tiny
           ! Intrinsic Functions
           intrinsic :: abs,real,int,log,max,mod,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           lquery = (lwork == -1 .or. lrwork == -1 .or. liwork == -1)
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
           if (info == 0) then
              ! compute the workspace requirements
              smlsiz = la_ilaenv(9,'WSTEDC',' ',0,0,0,0)
              if (n <= 1 .or. icompz == 0) then
                 lwmin = 1
                 liwmin = 1
                 lrwmin = 1
              else if (n <= smlsiz) then
                 lwmin = 1
                 liwmin = 1
                 lrwmin = 2*(n - 1)
              else if (icompz == 1) then
                 lgn = int(log(real(n,KIND=qp))/log(two),KIND=ilp)
                 if (2**lgn < n) lgn = lgn + 1
                 if (2**lgn < n) lgn = lgn + 1
                 lwmin = n*n
                 lrwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                 liwmin = 6 + 6*n + 5*n*lgn
              else if (icompz == 2) then
                 lwmin = 1
                 lrwmin = 1 + 4*n + 2*n**2
                 liwmin = 3 + 5*n
              end if
              work(1) = lwmin
              rwork(1) = lrwmin
              iwork(1) = liwmin
              if (lwork < lwmin .and. .not. lquery) then
                 info = -8
              else if (lrwork < lrwmin .and. .not. lquery) then
                 info = -10
              else if (liwork < liwmin .and. .not. lquery) then
                 info = -12
              end if
           end if
           if (info /= 0) then
              call la_xerbla('WSTEDC',-info)
              return
           else if (lquery) then
              return
           end if
           ! quick return if possible
           if (n == 0) return
           if (n == 1) then
              if (icompz /= 0) z(1,1) = one
              return
           end if
           ! if the following conditional clause is removed, then the routine
           ! will use the divide and conquer routine to compute only the
           ! eigenvalues, which requires (3n + 3n**2) real workspace and
           ! (2 + 5n + 2n lg(n)) integer workspace.
           ! since on many architectures la_qsterf is much faster than any other
           ! algorithm for finding eigenvalues only, it is used here
           ! as the default. if the conditional clause is removed, then
           ! information on the size of workspace needs to be changed.
           ! if compz = 'n', use la_qsterf to compute the eigenvalues.
           if (icompz == 0) then
              call la_qsterf(n,d,e,info)
              go to 70
           end if
           ! if n is smaller than the minimum divide size (smlsiz+1), then
           ! solve the problem with another solver.
           if (n <= smlsiz) then
              call la_wsteqr(compz,n,d,e,z,ldz,rwork,info)
           else
              ! if compz = 'i', we simply call la_qstedc instead.
              if (icompz == 2) then
                 call la_qlaset('FULL',n,n,zero,one,rwork,n)
                 ll = n*n + 1
                 call la_qstedc('I',n,d,e,rwork,n,rwork(ll),lrwork - ll + 1,iwork, &
                           liwork,info)
                 do j = 1,n
                    do i = 1,n
                       z(i,j) = rwork((j - 1)*n + i)
                    end do
                 end do
                 go to 70
              end if
              ! from now on, only option left to be handled is compz = 'v',
              ! i.e. icompz = 1.
              ! scale.
              orgnrm = la_qlanst('M',n,d,e)
              if (orgnrm == zero) go to 70
              eps = la_qlamch('EPSILON')
              start = 1
              ! while ( start <= n )
              30 continue
              if (start <= n) then
                 ! let finish be the position of the next subdiagonal entry
                 ! such that e( finish ) <= tiny or finish = n if no such
                 ! subdiagonal exists.  the matrix identified by the elements
                 ! between start and finish constitutes an independent
                 ! sub-problem.
                 finish = start
                 40 continue
                 if (finish < n) then
                    tiny = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                    if (abs(e(finish)) > tiny) then
                       finish = finish + 1
                       go to 40
                    end if
                 end if
                 ! (sub) problem determined.  compute its size and solve it.
                 m = finish - start + 1
                 if (m > smlsiz) then
                    ! scale.
                    orgnrm = la_qlanst('M',m,d(start),e(start))
                    call la_qlascl('G',0,0,orgnrm,one,m,1,d(start),m,info)
                    call la_qlascl('G',0,0,orgnrm,one,m - 1,1,e(start),m - 1,info)

                    call la_wlaed0(n,m,d(start),e(start),z(1,start),ldz,work,n, &
                              rwork,iwork,info)
                    if (info > 0) then
                       info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, (m + 1)) + start - &
                                 1
                       go to 70
                    end if
                    ! scale back.
                    call la_qlascl('G',0,0,one,orgnrm,m,1,d(start),m,info)
                 else
                    call la_qsteqr('I',m,d(start),e(start),rwork,m,rwork(m*m + 1), &
                              info)
                    call la_wlacrm(n,m,z(1,start),ldz,rwork,m,work,n,rwork(m*m + 1) &
                               )
                    call la_wlacpy('A',n,m,work,n,z(1,start),ldz)
                    if (info > 0) then
                       info = start*(n + 1) + finish
                       go to 70
                    end if
                 end if
                 start = finish + 1
                 go to 30
              end if
              ! endwhile
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
           70 continue
           work(1) = lwmin
           rwork(1) = lrwmin
           iwork(1) = liwmin
           return
     end subroutine la_wstedc

     !> CSTEGR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> CSTEGR is a compatibility wrapper around the improved CSTEMR routine.
     !> See SSTEMR for further details.
     !> One important change is that the ABSTOL parameter no longer provides any
     !> benefit and hence is no longer used.
     !> Note : CSTEGR and CSTEMR work only on machines which follow
     !> IEEE-754 floating-point standard in their handling of infinities and
     !> NaNs.  Normal execution may create these exceptiona values and hence
     !> may abort due to a floating point exception in environments which
     !> do not conform to the IEEE-754 standard.

     pure subroutine la_cstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),e(*)
           real(sp),intent(out) :: w(*),work(*)
           complex(sp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: tryrac
           ! Executable Statements
           info = 0
           tryrac = .false.
           call la_cstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,n,isuppz, &
                     tryrac,work,lwork,iwork,liwork,info)
     end subroutine la_cstegr
     !> ZSTEGR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> ZSTEGR is a compatibility wrapper around the improved ZSTEMR routine.
     !> See ZSTEMR for further details.
     !> One important change is that the ABSTOL parameter no longer provides any
     !> benefit and hence is no longer used.
     !> Note : ZSTEGR and ZSTEMR work only on machines which follow
     !> IEEE-754 floating-point standard in their handling of infinities and
     !> NaNs.  Normal execution may create these exceptiona values and hence
     !> may abort due to a floating point exception in environments which
     !> do not conform to the IEEE-754 standard.

     pure subroutine la_zstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),e(*)
           real(dp),intent(out) :: w(*),work(*)
           complex(dp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: tryrac
           ! Executable Statements
           info = 0
           tryrac = .false.
           call la_zstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,n,isuppz, &
                     tryrac,work,lwork,iwork,liwork,info)
     end subroutine la_zstegr
     !> WSTEGR: computes selected eigenvalues and, optionally, eigenvectors
     !> of a real symmetric tridiagonal matrix T. Any such unreduced matrix has
     !> a well defined set of pairwise different real eigenvalues, the corresponding
     !> real eigenvectors are pairwise orthogonal.
     !> The spectrum may be computed either completely or partially by specifying
     !> either an interval (VL,VU] or a range of indices IL:IU for the desired
     !> eigenvalues.
     !> WSTEGR is a compatibility wrapper around the improved WSTEMR routine.
     !> See WSTEMR for further details.
     !> One important change is that the ABSTOL parameter no longer provides any
     !> benefit and hence is no longer used.
     !> Note : WSTEGR and WSTEMR work only on machines which follow
     !> IEEE-754 floating-point standard in their handling of infinities and
     !> NaNs.  Normal execution may create these exceptiona values and hence
     !> may abort due to a floating point exception in environments which
     !> do not conform to the IEEE-754 standard.

     pure subroutine la_wstegr(jobz,range,n,d,e,vl,vu,il,iu,abstol,m,w,z,ldz, &
               isuppz,work,lwork,iwork,liwork,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobz,range
           integer(ilp),intent(in) :: il,iu,ldz,liwork,lwork,n
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: abstol,vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),e(*)
           real(qp),intent(out) :: w(*),work(*)
           complex(qp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Local Scalars
           logical(lk) :: tryrac
           ! Executable Statements
           info = 0
           tryrac = .false.
           call la_wstemr(jobz,range,n,d,e,vl,vu,il,iu,m,w,z,ldz,n,isuppz, &
                     tryrac,work,lwork,iwork,liwork,info)
     end subroutine la_wstegr

end module la_lapack_eigv_tridiag3
