!> Symmetric tridiagonal eigenvalues: MRRR representation tree, bisection, eigenvector generation
module la_lapack_eigv_tridiag2
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_blas_level2_gen
     use la_lapack_aux
     use la_lapack_auxiliary
     use la_lapack_blas_like_base
     use la_lapack_blas_like_scalar
     use la_lapack_svd_bidiag_qr
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slaebz
     public :: la_slaeda
     public :: la_slamrg
     public :: la_slarra
     public :: la_slarrc
     public :: la_slarrd
     public :: la_slarrj
     public :: la_slarrk
     public :: la_slarrr
     public :: la_slaneg
     public :: la_slar1v
     public :: la_slarrb
     public :: la_slarrf
     public :: la_slarrv
     public :: la_slarre
     public :: la_dlaebz
     public :: la_dlaeda
     public :: la_dlamrg
     public :: la_dlarra
     public :: la_dlarrc
     public :: la_dlarrd
     public :: la_dlarrj
     public :: la_dlarrk
     public :: la_dlarrr
     public :: la_dlaneg
     public :: la_dlar1v
     public :: la_dlarrb
     public :: la_dlarrf
     public :: la_dlarrv
     public :: la_dlarre
     public :: la_qlaebz
     public :: la_qlaeda
     public :: la_qlamrg
     public :: la_qlarra
     public :: la_qlarrc
     public :: la_qlarrd
     public :: la_qlarrj
     public :: la_qlarrk
     public :: la_qlarrr
     public :: la_qlaneg
     public :: la_qlar1v
     public :: la_qlarrb
     public :: la_qlarrf
     public :: la_qlarrv
     public :: la_qlarre
     public :: la_clar1v
     public :: la_clarrv
     public :: la_zlar1v
     public :: la_zlarrv
     public :: la_wlar1v
     public :: la_wlarrv

     contains

     !> SLAEBZ: contains the iteration loops which compute and use the
     !> function N(w), which is the count of eigenvalues of a symmetric
     !> tridiagonal matrix T less than or equal to its argument  w.  It
     !> performs a choice of two types of loops:
     !> IJOB=1, followed by
     !> IJOB=2: It takes as input a list of intervals and returns a list of
     !> sufficiently small intervals whose union contains the same
     !> eigenvalues as the union of the original intervals.
     !> The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
     !> The output interval (AB(j,1),AB(j,2)] will contain
     !> eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
     !> IJOB=3: It performs a binary search in each input interval
     !> (AB(j,1),AB(j,2)] for a point  w(j)  such that
     !> N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
     !> the search.  If such a w(j) is found, then on output
     !> AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
     !> (AB(j,1),AB(j,2)] will be a small interval containing the
     !> point where N(w) jumps through NVAL(j), unless that point
     !> lies outside the initial interval.
     !> Note that the intervals are in all cases half-open intervals,
     !> i.e., of the form  (a,b] , which includes  b  but not  a .
     !> To avoid underflow, the matrix should be scaled so that its largest
     !> element is no greater than  overflow**(1/2) * underflow**(1/4)
     !> in absolute value.  To assure the most accurate computation
     !> of small eigenvalues, the matrix should be scaled to be
     !> not much smaller than that, either.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966
     !> Note: the arguments are, in general, *not* checked for unreasonable
     !> values.

     pure subroutine la_slaebz(ijob,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d, &
               e,e2,nval,ab,c,mout,nab,work,iwork,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ijob,minp,mmax,n,nbmin,nitmax
           integer(ilp),intent(out) :: info,mout
           real(sp),intent(in) :: abstol,pivmin,reltol
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           integer(ilp),intent(inout) :: nab(mmax,*),nval(*)
           real(sp),intent(inout) :: ab(mmax,*),c(*)
           real(sp),intent(in) :: d(*),e(*),e2(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: itmp1,itmp2,j,ji,jit,jp,kf,kfnew,kl,klnew
           real(sp) :: tmp1,tmp2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! check for errors
           info = 0
           if (ijob < 1 .or. ijob > 3) then
              info = -1
              return
           end if
           ! initialize nab
           if (ijob == 1) then
              ! compute the number of eigenvalues in the initial intervals.
              mout = 0
              do ji = 1,minp
                 do jp = 1,2
                    tmp1 = d(1) - ab(ji,jp)
                    if (abs(tmp1) < pivmin) tmp1 = -pivmin
                    nab(ji,jp) = 0
                    if (tmp1 <= zero) nab(ji,jp) = 1
                    do j = 2,n
                       tmp1 = d(j) - e2(j - 1)/tmp1 - ab(ji,jp)
                       if (abs(tmp1) < pivmin) tmp1 = -pivmin
                       if (tmp1 <= zero) nab(ji,jp) = nab(ji,jp) + 1
                    end do
                 end do
                 mout = mout + nab(ji,2) - nab(ji,1)
              end do
              return
           end if
           ! initialize for loop
           ! kf and kl have the following meaning:
              ! intervals 1,...,kf-1 have converged.
              ! intervals kf,...,kl  still need to be refined.
           kf = 1
           kl = minp
           ! if ijob=2, initialize c.
           ! if ijob=3, use the user-supplied starting point.
           if (ijob == 2) then
              do ji = 1,minp
                 c(ji) = half*(ab(ji,1) + ab(ji,2))
              end do
           end if
           ! iteration loop
           loop_130: do jit = 1,nitmax
              ! loop over intervals
              if (kl - kf + 1 >= nbmin .and. nbmin > 0) then
                 ! begin of parallel version of the loop
                 do ji = kf,kl
                    ! compute n(c), the number of eigenvalues less than c
                    work(ji) = d(1) - c(ji)
                    iwork(ji) = 0
                    if (work(ji) <= pivmin) then
                       iwork(ji) = 1
                       work(ji) = min(work(ji),-pivmin)
                    end if
                    do j = 2,n
                       work(ji) = d(j) - e2(j - 1)/work(ji) - c(ji)
                       if (work(ji) <= pivmin) then
                          iwork(ji) = iwork(ji) + 1
                          work(ji) = min(work(ji),-pivmin)
                       end if
                    end do
                 end do
                 if (ijob <= 2) then
                    ! ijob=2: choose all intervals containing eigenvalues.
                    klnew = kl
                    loop_70: do ji = kf,kl
                       ! insure that n(w) is monotone
                       iwork(ji) = min(nab(ji,2),max(nab(ji,1),iwork(ji)))
                       ! update the queue -- add intervals if both halves
                       ! contain eigenvalues.
                       if (iwork(ji) == nab(ji,2)) then
                          ! no eigenvalue in the upper interval:
                          ! just use the lower interval.
                          ab(ji,2) = c(ji)
                       else if (iwork(ji) == nab(ji,1)) then
                          ! no eigenvalue in the lower interval:
                          ! just use the upper interval.
                          ab(ji,1) = c(ji)
                       else
                          klnew = klnew + 1
                          if (klnew <= mmax) then
                             ! eigenvalue in both intervals -- add upper to
                             ! queue.
                             ab(klnew,2) = ab(ji,2)
                             nab(klnew,2) = nab(ji,2)
                             ab(klnew,1) = c(ji)
                             nab(klnew,1) = iwork(ji)
                             ab(ji,2) = c(ji)
                             nab(ji,2) = iwork(ji)
                          else
                             info = mmax + 1
                          end if
                       end if
                    end do loop_70
                    if (info /= 0) return
                    kl = klnew
                 else
                    ! ijob=3: binary search.  keep only the interval containing
                            ! w   s.t. n(w) = nval
                    do ji = kf,kl
                       if (iwork(ji) <= nval(ji)) then
                          ab(ji,1) = c(ji)
                          nab(ji,1) = iwork(ji)
                       end if
                       if (iwork(ji) >= nval(ji)) then
                          ab(ji,2) = c(ji)
                          nab(ji,2) = iwork(ji)
                       end if
                    end do
                 end if
              else
                 ! end of parallel version of the loop
                 ! begin of serial version of the loop
                 klnew = kl
                 loop_100: do ji = kf,kl
                    ! compute n(w), the number of eigenvalues less than w
                    tmp1 = c(ji)
                    tmp2 = d(1) - tmp1
                    itmp1 = 0
                    if (tmp2 <= pivmin) then
                       itmp1 = 1
                       tmp2 = min(tmp2,-pivmin)
                    end if
                    do j = 2,n
                       tmp2 = d(j) - e2(j - 1)/tmp2 - tmp1
                       if (tmp2 <= pivmin) then
                          itmp1 = itmp1 + 1
                          tmp2 = min(tmp2,-pivmin)
                       end if
                    end do
                    if (ijob <= 2) then
                       ! ijob=2: choose all intervals containing eigenvalues.
                       ! insure that n(w) is monotone
                       itmp1 = min(nab(ji,2),max(nab(ji,1),itmp1))
                       ! update the queue -- add intervals if both halves
                       ! contain eigenvalues.
                       if (itmp1 == nab(ji,2)) then
                          ! no eigenvalue in the upper interval:
                          ! just use the lower interval.
                          ab(ji,2) = tmp1
                       else if (itmp1 == nab(ji,1)) then
                          ! no eigenvalue in the lower interval:
                          ! just use the upper interval.
                          ab(ji,1) = tmp1
                       else if (klnew < mmax) then
                          ! eigenvalue in both intervals -- add upper to queue.
                          klnew = klnew + 1
                          ab(klnew,2) = ab(ji,2)
                          nab(klnew,2) = nab(ji,2)
                          ab(klnew,1) = tmp1
                          nab(klnew,1) = itmp1
                          ab(ji,2) = tmp1
                          nab(ji,2) = itmp1
                       else
                          info = mmax + 1
                          return
                       end if
                    else
                       ! ijob=3: binary search.  keep only the interval
                               ! containing  w  s.t. n(w) = nval
                       if (itmp1 <= nval(ji)) then
                          ab(ji,1) = tmp1
                          nab(ji,1) = itmp1
                       end if
                       if (itmp1 >= nval(ji)) then
                          ab(ji,2) = tmp1
                          nab(ji,2) = itmp1
                       end if
                    end if
                 end do loop_100
                 kl = klnew
              end if
              ! check for convergence
              kfnew = kf
              loop_110: do ji = kf,kl
                 tmp1 = abs(ab(ji,2) - ab(ji,1))
                 tmp2 = max(abs(ab(ji,2)),abs(ab(ji,1)))
                 if (tmp1 < max(abstol,pivmin,reltol*tmp2) .or. nab(ji,1) >= nab(ji,2)) &
                           then
                    ! converged -- swap with position kfnew,
                                 ! then increment kfnew
                    if (ji > kfnew) then
                       tmp1 = ab(ji,1)
                       tmp2 = ab(ji,2)
                       itmp1 = nab(ji,1)
                       itmp2 = nab(ji,2)
                       ab(ji,1) = ab(kfnew,1)
                       ab(ji,2) = ab(kfnew,2)
                       nab(ji,1) = nab(kfnew,1)
                       nab(ji,2) = nab(kfnew,2)
                       ab(kfnew,1) = tmp1
                       ab(kfnew,2) = tmp2
                       nab(kfnew,1) = itmp1
                       nab(kfnew,2) = itmp2
                       if (ijob == 3) then
                          itmp1 = nval(ji)
                          nval(ji) = nval(kfnew)
                          nval(kfnew) = itmp1
                       end if
                    end if
                    kfnew = kfnew + 1
                 end if
              end do loop_110
              kf = kfnew
              ! choose midpoints
              do ji = kf,kl
                 c(ji) = half*(ab(ji,1) + ab(ji,2))
              end do
              ! if no more intervals to refine, quit.
              if (kf > kl) go to 140
           end do loop_130
           ! converged
           140 continue
           info = max(kl + 1 - kf,0)
           mout = kl
           return
     end subroutine la_slaebz
     !> DLAEBZ: contains the iteration loops which compute and use the
     !> function N(w), which is the count of eigenvalues of a symmetric
     !> tridiagonal matrix T less than or equal to its argument  w.  It
     !> performs a choice of two types of loops:
     !> IJOB=1, followed by
     !> IJOB=2: It takes as input a list of intervals and returns a list of
     !> sufficiently small intervals whose union contains the same
     !> eigenvalues as the union of the original intervals.
     !> The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
     !> The output interval (AB(j,1),AB(j,2)] will contain
     !> eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
     !> IJOB=3: It performs a binary search in each input interval
     !> (AB(j,1),AB(j,2)] for a point  w(j)  such that
     !> N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
     !> the search.  If such a w(j) is found, then on output
     !> AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
     !> (AB(j,1),AB(j,2)] will be a small interval containing the
     !> point where N(w) jumps through NVAL(j), unless that point
     !> lies outside the initial interval.
     !> Note that the intervals are in all cases half-open intervals,
     !> i.e., of the form  (a,b] , which includes  b  but not  a .
     !> To avoid underflow, the matrix should be scaled so that its largest
     !> element is no greater than  overflow**(1/2) * underflow**(1/4)
     !> in absolute value.  To assure the most accurate computation
     !> of small eigenvalues, the matrix should be scaled to be
     !> not much smaller than that, either.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966
     !> Note: the arguments are, in general, *not* checked for unreasonable
     !> values.

     pure subroutine la_dlaebz(ijob,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d, &
               e,e2,nval,ab,c,mout,nab,work,iwork,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ijob,minp,mmax,n,nbmin,nitmax
           integer(ilp),intent(out) :: info,mout
           real(dp),intent(in) :: abstol,pivmin,reltol
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           integer(ilp),intent(inout) :: nab(mmax,*),nval(*)
           real(dp),intent(inout) :: ab(mmax,*),c(*)
           real(dp),intent(in) :: d(*),e(*),e2(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: itmp1,itmp2,j,ji,jit,jp,kf,kfnew,kl,klnew
           real(dp) :: tmp1,tmp2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! check for errors
           info = 0
           if (ijob < 1 .or. ijob > 3) then
              info = -1
              return
           end if
           ! initialize nab
           if (ijob == 1) then
              ! compute the number of eigenvalues in the initial intervals.
              mout = 0
              do ji = 1,minp
                 do jp = 1,2
                    tmp1 = d(1) - ab(ji,jp)
                    if (abs(tmp1) < pivmin) tmp1 = -pivmin
                    nab(ji,jp) = 0
                    if (tmp1 <= zero) nab(ji,jp) = 1
                    do j = 2,n
                       tmp1 = d(j) - e2(j - 1)/tmp1 - ab(ji,jp)
                       if (abs(tmp1) < pivmin) tmp1 = -pivmin
                       if (tmp1 <= zero) nab(ji,jp) = nab(ji,jp) + 1
                    end do
                 end do
                 mout = mout + nab(ji,2) - nab(ji,1)
              end do
              return
           end if
           ! initialize for loop
           ! kf and kl have the following meaning:
              ! intervals 1,...,kf-1 have converged.
              ! intervals kf,...,kl  still need to be refined.
           kf = 1
           kl = minp
           ! if ijob=2, initialize c.
           ! if ijob=3, use the user-supplied starting point.
           if (ijob == 2) then
              do ji = 1,minp
                 c(ji) = half*(ab(ji,1) + ab(ji,2))
              end do
           end if
           ! iteration loop
           loop_130: do jit = 1,nitmax
              ! loop over intervals
              if (kl - kf + 1 >= nbmin .and. nbmin > 0) then
                 ! begin of parallel version of the loop
                 do ji = kf,kl
                    ! compute n(c), the number of eigenvalues less than c
                    work(ji) = d(1) - c(ji)
                    iwork(ji) = 0
                    if (work(ji) <= pivmin) then
                       iwork(ji) = 1
                       work(ji) = min(work(ji),-pivmin)
                    end if
                    do j = 2,n
                       work(ji) = d(j) - e2(j - 1)/work(ji) - c(ji)
                       if (work(ji) <= pivmin) then
                          iwork(ji) = iwork(ji) + 1
                          work(ji) = min(work(ji),-pivmin)
                       end if
                    end do
                 end do
                 if (ijob <= 2) then
                    ! ijob=2: choose all intervals containing eigenvalues.
                    klnew = kl
                    loop_70: do ji = kf,kl
                       ! insure that n(w) is monotone
                       iwork(ji) = min(nab(ji,2),max(nab(ji,1),iwork(ji)))
                       ! update the queue -- add intervals if both halves
                       ! contain eigenvalues.
                       if (iwork(ji) == nab(ji,2)) then
                          ! no eigenvalue in the upper interval:
                          ! just use the lower interval.
                          ab(ji,2) = c(ji)
                       else if (iwork(ji) == nab(ji,1)) then
                          ! no eigenvalue in the lower interval:
                          ! just use the upper interval.
                          ab(ji,1) = c(ji)
                       else
                          klnew = klnew + 1
                          if (klnew <= mmax) then
                             ! eigenvalue in both intervals -- add upper to
                             ! queue.
                             ab(klnew,2) = ab(ji,2)
                             nab(klnew,2) = nab(ji,2)
                             ab(klnew,1) = c(ji)
                             nab(klnew,1) = iwork(ji)
                             ab(ji,2) = c(ji)
                             nab(ji,2) = iwork(ji)
                          else
                             info = mmax + 1
                          end if
                       end if
                    end do loop_70
                    if (info /= 0) return
                    kl = klnew
                 else
                    ! ijob=3: binary search.  keep only the interval containing
                            ! w   s.t. n(w) = nval
                    do ji = kf,kl
                       if (iwork(ji) <= nval(ji)) then
                          ab(ji,1) = c(ji)
                          nab(ji,1) = iwork(ji)
                       end if
                       if (iwork(ji) >= nval(ji)) then
                          ab(ji,2) = c(ji)
                          nab(ji,2) = iwork(ji)
                       end if
                    end do
                 end if
              else
                 ! end of parallel version of the loop
                 ! begin of serial version of the loop
                 klnew = kl
                 loop_100: do ji = kf,kl
                    ! compute n(w), the number of eigenvalues less than w
                    tmp1 = c(ji)
                    tmp2 = d(1) - tmp1
                    itmp1 = 0
                    if (tmp2 <= pivmin) then
                       itmp1 = 1
                       tmp2 = min(tmp2,-pivmin)
                    end if
                    do j = 2,n
                       tmp2 = d(j) - e2(j - 1)/tmp2 - tmp1
                       if (tmp2 <= pivmin) then
                          itmp1 = itmp1 + 1
                          tmp2 = min(tmp2,-pivmin)
                       end if
                    end do
                    if (ijob <= 2) then
                       ! ijob=2: choose all intervals containing eigenvalues.
                       ! insure that n(w) is monotone
                       itmp1 = min(nab(ji,2),max(nab(ji,1),itmp1))
                       ! update the queue -- add intervals if both halves
                       ! contain eigenvalues.
                       if (itmp1 == nab(ji,2)) then
                          ! no eigenvalue in the upper interval:
                          ! just use the lower interval.
                          ab(ji,2) = tmp1
                       else if (itmp1 == nab(ji,1)) then
                          ! no eigenvalue in the lower interval:
                          ! just use the upper interval.
                          ab(ji,1) = tmp1
                       else if (klnew < mmax) then
                          ! eigenvalue in both intervals -- add upper to queue.
                          klnew = klnew + 1
                          ab(klnew,2) = ab(ji,2)
                          nab(klnew,2) = nab(ji,2)
                          ab(klnew,1) = tmp1
                          nab(klnew,1) = itmp1
                          ab(ji,2) = tmp1
                          nab(ji,2) = itmp1
                       else
                          info = mmax + 1
                          return
                       end if
                    else
                       ! ijob=3: binary search.  keep only the interval
                               ! containing  w  s.t. n(w) = nval
                       if (itmp1 <= nval(ji)) then
                          ab(ji,1) = tmp1
                          nab(ji,1) = itmp1
                       end if
                       if (itmp1 >= nval(ji)) then
                          ab(ji,2) = tmp1
                          nab(ji,2) = itmp1
                       end if
                    end if
                 end do loop_100
                 kl = klnew
              end if
              ! check for convergence
              kfnew = kf
              loop_110: do ji = kf,kl
                 tmp1 = abs(ab(ji,2) - ab(ji,1))
                 tmp2 = max(abs(ab(ji,2)),abs(ab(ji,1)))
                 if (tmp1 < max(abstol,pivmin,reltol*tmp2) .or. nab(ji,1) >= nab(ji,2)) &
                           then
                    ! converged -- swap with position kfnew,
                                 ! then increment kfnew
                    if (ji > kfnew) then
                       tmp1 = ab(ji,1)
                       tmp2 = ab(ji,2)
                       itmp1 = nab(ji,1)
                       itmp2 = nab(ji,2)
                       ab(ji,1) = ab(kfnew,1)
                       ab(ji,2) = ab(kfnew,2)
                       nab(ji,1) = nab(kfnew,1)
                       nab(ji,2) = nab(kfnew,2)
                       ab(kfnew,1) = tmp1
                       ab(kfnew,2) = tmp2
                       nab(kfnew,1) = itmp1
                       nab(kfnew,2) = itmp2
                       if (ijob == 3) then
                          itmp1 = nval(ji)
                          nval(ji) = nval(kfnew)
                          nval(kfnew) = itmp1
                       end if
                    end if
                    kfnew = kfnew + 1
                 end if
              end do loop_110
              kf = kfnew
              ! choose midpoints
              do ji = kf,kl
                 c(ji) = half*(ab(ji,1) + ab(ji,2))
              end do
              ! if no more intervals to refine, quit.
              if (kf > kl) go to 140
           end do loop_130
           ! converged
           140 continue
           info = max(kl + 1 - kf,0)
           mout = kl
           return
     end subroutine la_dlaebz
     !> QLAEBZ: contains the iteration loops which compute and use the
     !> function N(w), which is the count of eigenvalues of a symmetric
     !> tridiagonal matrix T less than or equal to its argument  w.  It
     !> performs a choice of two types of loops:
     !> IJOB=1, followed by
     !> IJOB=2: It takes as input a list of intervals and returns a list of
     !> sufficiently small intervals whose union contains the same
     !> eigenvalues as the union of the original intervals.
     !> The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.
     !> The output interval (AB(j,1),AB(j,2)] will contain
     !> eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.
     !> IJOB=3: It performs a binary search in each input interval
     !> (AB(j,1),AB(j,2)] for a point  w(j)  such that
     !> N(w(j))=NVAL(j), and uses  C(j)  as the starting point of
     !> the search.  If such a w(j) is found, then on output
     !> AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on output
     !> (AB(j,1),AB(j,2)] will be a small interval containing the
     !> point where N(w) jumps through NVAL(j), unless that point
     !> lies outside the initial interval.
     !> Note that the intervals are in all cases half-open intervals,
     !> i.e., of the form  (a,b] , which includes  b  but not  a .
     !> To avoid underflow, the matrix should be scaled so that its largest
     !> element is no greater than  overflow**(1/2) * underflow**(1/4)
     !> in absolute value.  To assure the most accurate computation
     !> of small eigenvalues, the matrix should be scaled to be
     !> not much smaller than that, either.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966
     !> Note: the arguments are, in general, *not* checked for unreasonable
     !> values.

     pure subroutine la_qlaebz(ijob,nitmax,n,mmax,minp,nbmin,abstol,reltol,pivmin,d, &
               e,e2,nval,ab,c,mout,nab,work,iwork,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ijob,minp,mmax,n,nbmin,nitmax
           integer(ilp),intent(out) :: info,mout
           real(qp),intent(in) :: abstol,pivmin,reltol
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           integer(ilp),intent(inout) :: nab(mmax,*),nval(*)
           real(qp),intent(inout) :: ab(mmax,*),c(*)
           real(qp),intent(in) :: d(*),e(*),e2(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: itmp1,itmp2,j,ji,jit,jp,kf,kfnew,kl,klnew
           real(qp) :: tmp1,tmp2
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           ! check for errors
           info = 0
           if (ijob < 1 .or. ijob > 3) then
              info = -1
              return
           end if
           ! initialize nab
           if (ijob == 1) then
              ! compute the number of eigenvalues in the initial intervals.
              mout = 0
              do ji = 1,minp
                 do jp = 1,2
                    tmp1 = d(1) - ab(ji,jp)
                    if (abs(tmp1) < pivmin) tmp1 = -pivmin
                    nab(ji,jp) = 0
                    if (tmp1 <= zero) nab(ji,jp) = 1
                    do j = 2,n
                       tmp1 = d(j) - e2(j - 1)/tmp1 - ab(ji,jp)
                       if (abs(tmp1) < pivmin) tmp1 = -pivmin
                       if (tmp1 <= zero) nab(ji,jp) = nab(ji,jp) + 1
                    end do
                 end do
                 mout = mout + nab(ji,2) - nab(ji,1)
              end do
              return
           end if
           ! initialize for loop
           ! kf and kl have the following meaning:
              ! intervals 1,...,kf-1 have converged.
              ! intervals kf,...,kl  still need to be refined.
           kf = 1
           kl = minp
           ! if ijob=2, initialize c.
           ! if ijob=3, use the user-supplied starting point.
           if (ijob == 2) then
              do ji = 1,minp
                 c(ji) = half*(ab(ji,1) + ab(ji,2))
              end do
           end if
           ! iteration loop
           loop_130: do jit = 1,nitmax
              ! loop over intervals
              if (kl - kf + 1 >= nbmin .and. nbmin > 0) then
                 ! begin of parallel version of the loop
                 do ji = kf,kl
                    ! compute n(c), the number of eigenvalues less than c
                    work(ji) = d(1) - c(ji)
                    iwork(ji) = 0
                    if (work(ji) <= pivmin) then
                       iwork(ji) = 1
                       work(ji) = min(work(ji),-pivmin)
                    end if
                    do j = 2,n
                       work(ji) = d(j) - e2(j - 1)/work(ji) - c(ji)
                       if (work(ji) <= pivmin) then
                          iwork(ji) = iwork(ji) + 1
                          work(ji) = min(work(ji),-pivmin)
                       end if
                    end do
                 end do
                 if (ijob <= 2) then
                    ! ijob=2: choose all intervals containing eigenvalues.
                    klnew = kl
                    loop_70: do ji = kf,kl
                       ! insure that n(w) is monotone
                       iwork(ji) = min(nab(ji,2),max(nab(ji,1),iwork(ji)))
                       ! update the queue -- add intervals if both halves
                       ! contain eigenvalues.
                       if (iwork(ji) == nab(ji,2)) then
                          ! no eigenvalue in the upper interval:
                          ! just use the lower interval.
                          ab(ji,2) = c(ji)
                       else if (iwork(ji) == nab(ji,1)) then
                          ! no eigenvalue in the lower interval:
                          ! just use the upper interval.
                          ab(ji,1) = c(ji)
                       else
                          klnew = klnew + 1
                          if (klnew <= mmax) then
                             ! eigenvalue in both intervals -- add upper to
                             ! queue.
                             ab(klnew,2) = ab(ji,2)
                             nab(klnew,2) = nab(ji,2)
                             ab(klnew,1) = c(ji)
                             nab(klnew,1) = iwork(ji)
                             ab(ji,2) = c(ji)
                             nab(ji,2) = iwork(ji)
                          else
                             info = mmax + 1
                          end if
                       end if
                    end do loop_70
                    if (info /= 0) return
                    kl = klnew
                 else
                    ! ijob=3: binary search.  keep only the interval containing
                            ! w   s.t. n(w) = nval
                    do ji = kf,kl
                       if (iwork(ji) <= nval(ji)) then
                          ab(ji,1) = c(ji)
                          nab(ji,1) = iwork(ji)
                       end if
                       if (iwork(ji) >= nval(ji)) then
                          ab(ji,2) = c(ji)
                          nab(ji,2) = iwork(ji)
                       end if
                    end do
                 end if
              else
                 ! end of parallel version of the loop
                 ! begin of serial version of the loop
                 klnew = kl
                 loop_100: do ji = kf,kl
                    ! compute n(w), the number of eigenvalues less than w
                    tmp1 = c(ji)
                    tmp2 = d(1) - tmp1
                    itmp1 = 0
                    if (tmp2 <= pivmin) then
                       itmp1 = 1
                       tmp2 = min(tmp2,-pivmin)
                    end if
                    do j = 2,n
                       tmp2 = d(j) - e2(j - 1)/tmp2 - tmp1
                       if (tmp2 <= pivmin) then
                          itmp1 = itmp1 + 1
                          tmp2 = min(tmp2,-pivmin)
                       end if
                    end do
                    if (ijob <= 2) then
                       ! ijob=2: choose all intervals containing eigenvalues.
                       ! insure that n(w) is monotone
                       itmp1 = min(nab(ji,2),max(nab(ji,1),itmp1))
                       ! update the queue -- add intervals if both halves
                       ! contain eigenvalues.
                       if (itmp1 == nab(ji,2)) then
                          ! no eigenvalue in the upper interval:
                          ! just use the lower interval.
                          ab(ji,2) = tmp1
                       else if (itmp1 == nab(ji,1)) then
                          ! no eigenvalue in the lower interval:
                          ! just use the upper interval.
                          ab(ji,1) = tmp1
                       else if (klnew < mmax) then
                          ! eigenvalue in both intervals -- add upper to queue.
                          klnew = klnew + 1
                          ab(klnew,2) = ab(ji,2)
                          nab(klnew,2) = nab(ji,2)
                          ab(klnew,1) = tmp1
                          nab(klnew,1) = itmp1
                          ab(ji,2) = tmp1
                          nab(ji,2) = itmp1
                       else
                          info = mmax + 1
                          return
                       end if
                    else
                       ! ijob=3: binary search.  keep only the interval
                               ! containing  w  s.t. n(w) = nval
                       if (itmp1 <= nval(ji)) then
                          ab(ji,1) = tmp1
                          nab(ji,1) = itmp1
                       end if
                       if (itmp1 >= nval(ji)) then
                          ab(ji,2) = tmp1
                          nab(ji,2) = itmp1
                       end if
                    end if
                 end do loop_100
                 kl = klnew
              end if
              ! check for convergence
              kfnew = kf
              loop_110: do ji = kf,kl
                 tmp1 = abs(ab(ji,2) - ab(ji,1))
                 tmp2 = max(abs(ab(ji,2)),abs(ab(ji,1)))
                 if (tmp1 < max(abstol,pivmin,reltol*tmp2) .or. nab(ji,1) >= nab(ji,2)) &
                           then
                    ! converged -- swap with position kfnew,
                                 ! then increment kfnew
                    if (ji > kfnew) then
                       tmp1 = ab(ji,1)
                       tmp2 = ab(ji,2)
                       itmp1 = nab(ji,1)
                       itmp2 = nab(ji,2)
                       ab(ji,1) = ab(kfnew,1)
                       ab(ji,2) = ab(kfnew,2)
                       nab(ji,1) = nab(kfnew,1)
                       nab(ji,2) = nab(kfnew,2)
                       ab(kfnew,1) = tmp1
                       ab(kfnew,2) = tmp2
                       nab(kfnew,1) = itmp1
                       nab(kfnew,2) = itmp2
                       if (ijob == 3) then
                          itmp1 = nval(ji)
                          nval(ji) = nval(kfnew)
                          nval(kfnew) = itmp1
                       end if
                    end if
                    kfnew = kfnew + 1
                 end if
              end do loop_110
              kf = kfnew
              ! choose midpoints
              do ji = kf,kl
                 c(ji) = half*(ab(ji,1) + ab(ji,2))
              end do
              ! if no more intervals to refine, quit.
              if (kf > kl) go to 140
           end do loop_130
           ! converged
           140 continue
           info = max(kl + 1 - kf,0)
           mout = kl
           return
     end subroutine la_qlaebz

     !> SLAEDA: computes the Z vector corresponding to the merge step in the
     !> CURLVLth step of the merge process with TLVLS steps for the CURPBMth
     !> problem.

     pure subroutine la_slaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                q,qptr,z,ztemp,info)
        use la_constants_sp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,n,tlvls
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)
           real(sp),intent(in) :: givnum(2,*),q(*)
           real(sp),intent(out) :: z(*),ztemp(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: bsiz1,bsiz2,curr,i,k,mid,psiz1,psiz2,ptr,zptr1
           ! Intrinsic Functions
           intrinsic :: int,real,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           end if
           if (info /= 0) then
              call la_xerbla('SLAEDA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine location of first number in second half.
           mid = n/2 + 1
           ! gather last/first rows of appropriate eigenblocks into center of z
           ptr = 1
           ! determine location of lowest level subproblem in the full storage
           ! scheme
           curr = ptr + curpbm*2**curlvl + 2**(curlvl - 1) - 1
           ! determine size of these matrices.  we add half to the value of
           ! the sqrt in case the machine underestimates one of these square
           ! roots.
           bsiz1 = int(half + sqrt(real(qptr(curr + 1) - qptr(curr),KIND=sp)),KIND=ilp)
           bsiz2 = int(half + sqrt(real(qptr(curr + 2) - qptr(curr + 1),KIND=sp)),KIND=ilp)

           do k = 1,mid - bsiz1 - 1
              z(k) = zero
           end do
           call la_scopy(bsiz1,q(qptr(curr) + bsiz1 - 1),bsiz1,z(mid - bsiz1),1)
           call la_scopy(bsiz2,q(qptr(curr + 1)),bsiz2,z(mid),1)
           do k = mid + bsiz2,n
              z(k) = zero
           end do
           ! loop through remaining levels 1 -> curlvl applying the givens
           ! rotations and permutation and then multiplying the center matrices
           ! against the current z.
           ptr = 2**tlvls + 1
           loop_70: do k = 1,curlvl - 1
              curr = ptr + curpbm*2**(curlvl - k) + 2**(curlvl - k - 1) - 1
              psiz1 = prmptr(curr + 1) - prmptr(curr)
              psiz2 = prmptr(curr + 2) - prmptr(curr + 1)
              zptr1 = mid - psiz1
             ! apply givens at curr and curr+1
              do i = givptr(curr),givptr(curr + 1) - 1
                 call la_srot(1,z(zptr1 + givcol(1,i) - 1),1,z(zptr1 + givcol(2,i) - 1), &
                           1,givnum(1,i),givnum(2,i))
              end do
              do i = givptr(curr + 1),givptr(curr + 2) - 1
                 call la_srot(1,z(mid - 1 + givcol(1,i)),1,z(mid - 1 + givcol(2,i)),1, &
                           givnum(1,i),givnum(2,i))
              end do
              psiz1 = prmptr(curr + 1) - prmptr(curr)
              psiz2 = prmptr(curr + 2) - prmptr(curr + 1)
              do i = 0,psiz1 - 1
                 ztemp(i + 1) = z(zptr1 + perm(prmptr(curr) + i) - 1)
              end do
              do i = 0,psiz2 - 1
                 ztemp(psiz1 + i + 1) = z(mid + perm(prmptr(curr + 1) + i) - 1)
              end do
              ! multiply blocks at curr and curr+1
              ! determine size of these matrices.  we add half to the value of
              ! the sqrt in case the machine underestimates one of these
              ! square roots.
              bsiz1 = int(half + sqrt(real(qptr(curr + 1) - qptr(curr),KIND=sp)),KIND=ilp)

              bsiz2 = int(half + sqrt(real(qptr(curr + 2) - qptr(curr + 1),KIND=sp)),KIND=ilp)

              if (bsiz1 > 0) then
                 call la_sgemv('T',bsiz1,bsiz1,one,q(qptr(curr)),bsiz1,ztemp(1), &
                           1,zero,z(zptr1),1)
              end if
              call la_scopy(psiz1 - bsiz1,ztemp(bsiz1 + 1),1,z(zptr1 + bsiz1),1)
              if (bsiz2 > 0) then
                 call la_sgemv('T',bsiz2,bsiz2,one,q(qptr(curr + 1)),bsiz2,ztemp( &
                           psiz1 + 1),1,zero,z(mid),1)
              end if
              call la_scopy(psiz2 - bsiz2,ztemp(psiz1 + bsiz2 + 1),1,z(mid + bsiz2),1)

              ptr = ptr + 2**(tlvls - k)
           end do loop_70
           return
     end subroutine la_slaeda
     !> DLAEDA: computes the Z vector corresponding to the merge step in the
     !> CURLVLth step of the merge process with TLVLS steps for the CURPBMth
     !> problem.

     pure subroutine la_dlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                q,qptr,z,ztemp,info)
        use la_constants_dp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,n,tlvls
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)
           real(dp),intent(in) :: givnum(2,*),q(*)
           real(dp),intent(out) :: z(*),ztemp(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: bsiz1,bsiz2,curr,i,k,mid,psiz1,psiz2,ptr,zptr1
           ! Intrinsic Functions
           intrinsic :: real,int,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           end if
           if (info /= 0) then
              call la_xerbla('DLAEDA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine location of first number in second half.
           mid = n/2 + 1
           ! gather last/first rows of appropriate eigenblocks into center of z
           ptr = 1
           ! determine location of lowest level subproblem in the full storage
           ! scheme
           curr = ptr + curpbm*2**curlvl + 2**(curlvl - 1) - 1
           ! determine size of these matrices.  we add half to the value of
           ! the sqrt in case the machine underestimates one of these square
           ! roots.
           bsiz1 = int(half + sqrt(real(qptr(curr + 1) - qptr(curr),KIND=dp)),KIND=ilp)
           bsiz2 = int(half + sqrt(real(qptr(curr + 2) - qptr(curr + 1),KIND=dp)),KIND=ilp)

           do k = 1,mid - bsiz1 - 1
              z(k) = zero
           end do
           call la_dcopy(bsiz1,q(qptr(curr) + bsiz1 - 1),bsiz1,z(mid - bsiz1),1)
           call la_dcopy(bsiz2,q(qptr(curr + 1)),bsiz2,z(mid),1)
           do k = mid + bsiz2,n
              z(k) = zero
           end do
           ! loop through remaining levels 1 -> curlvl applying the givens
           ! rotations and permutation and then multiplying the center matrices
           ! against the current z.
           ptr = 2**tlvls + 1
           loop_70: do k = 1,curlvl - 1
              curr = ptr + curpbm*2**(curlvl - k) + 2**(curlvl - k - 1) - 1
              psiz1 = prmptr(curr + 1) - prmptr(curr)
              psiz2 = prmptr(curr + 2) - prmptr(curr + 1)
              zptr1 = mid - psiz1
             ! apply givens at curr and curr+1
              do i = givptr(curr),givptr(curr + 1) - 1
                 call la_drot(1,z(zptr1 + givcol(1,i) - 1),1,z(zptr1 + givcol(2,i) - 1), &
                           1,givnum(1,i),givnum(2,i))
              end do
              do i = givptr(curr + 1),givptr(curr + 2) - 1
                 call la_drot(1,z(mid - 1 + givcol(1,i)),1,z(mid - 1 + givcol(2,i)),1, &
                           givnum(1,i),givnum(2,i))
              end do
              psiz1 = prmptr(curr + 1) - prmptr(curr)
              psiz2 = prmptr(curr + 2) - prmptr(curr + 1)
              do i = 0,psiz1 - 1
                 ztemp(i + 1) = z(zptr1 + perm(prmptr(curr) + i) - 1)
              end do
              do i = 0,psiz2 - 1
                 ztemp(psiz1 + i + 1) = z(mid + perm(prmptr(curr + 1) + i) - 1)
              end do
              ! multiply blocks at curr and curr+1
              ! determine size of these matrices.  we add half to the value of
              ! the sqrt in case the machine underestimates one of these
              ! square roots.
              bsiz1 = int(half + sqrt(real(qptr(curr + 1) - qptr(curr),KIND=dp)),KIND=ilp)

              bsiz2 = int(half + sqrt(real(qptr(curr + 2) - qptr(curr + 1),KIND=dp)),KIND=ilp)

              if (bsiz1 > 0) then
                 call la_dgemv('T',bsiz1,bsiz1,one,q(qptr(curr)),bsiz1,ztemp(1), &
                           1,zero,z(zptr1),1)
              end if
              call la_dcopy(psiz1 - bsiz1,ztemp(bsiz1 + 1),1,z(zptr1 + bsiz1),1)
              if (bsiz2 > 0) then
                 call la_dgemv('T',bsiz2,bsiz2,one,q(qptr(curr + 1)),bsiz2,ztemp( &
                           psiz1 + 1),1,zero,z(mid),1)
              end if
              call la_dcopy(psiz2 - bsiz2,ztemp(psiz1 + bsiz2 + 1),1,z(mid + bsiz2),1)

              ptr = ptr + 2**(tlvls - k)
           end do loop_70
           return
     end subroutine la_dlaeda
     !> QLAEDA: computes the Z vector corresponding to the merge step in the
     !> CURLVLth step of the merge process with TLVLS steps for the CURPBMth
     !> problem.

     pure subroutine la_qlaeda(n,tlvls,curlvl,curpbm,prmptr,perm,givptr,givcol,givnum, &
                q,qptr,z,ztemp,info)
        use la_constants_qp
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: curlvl,curpbm,n,tlvls
           integer(ilp),intent(out) :: info
           ! Array Arguments
           integer(ilp),intent(in) :: givcol(2,*),givptr(*),perm(*),prmptr(*),qptr(*)
           real(qp),intent(in) :: givnum(2,*),q(*)
           real(qp),intent(out) :: z(*),ztemp(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: bsiz1,bsiz2,curr,i,k,mid,psiz1,psiz2,ptr,zptr1
           ! Intrinsic Functions
           intrinsic :: real,int,sqrt
           ! Executable Statements
           ! test the input parameters.
           info = 0
           if (n < 0) then
              info = -1
           end if
           if (info /= 0) then
              call la_xerbla('QLAEDA',-info)
              return
           end if
           ! quick return if possible
           if (n == 0) return
           ! determine location of first number in second half.
           mid = n/2 + 1
           ! gather last/first rows of appropriate eigenblocks into center of z
           ptr = 1
           ! determine location of lowest level subproblem in the full storage
           ! scheme
           curr = ptr + curpbm*2**curlvl + 2**(curlvl - 1) - 1
           ! determine size of these matrices.  we add half to the value of
           ! the sqrt in case the machine underestimates one of these square
           ! roots.
           bsiz1 = int(half + sqrt(real(qptr(curr + 1) - qptr(curr),KIND=qp)),KIND=ilp)
           bsiz2 = int(half + sqrt(real(qptr(curr + 2) - qptr(curr + 1),KIND=qp)),KIND=ilp)

           do k = 1,mid - bsiz1 - 1
              z(k) = zero
           end do
           call la_qcopy(bsiz1,q(qptr(curr) + bsiz1 - 1),bsiz1,z(mid - bsiz1),1)
           call la_qcopy(bsiz2,q(qptr(curr + 1)),bsiz2,z(mid),1)
           do k = mid + bsiz2,n
              z(k) = zero
           end do
           ! loop through remaining levels 1 -> curlvl applying the givens
           ! rotations and permutation and then multiplying the center matrices
           ! against the current z.
           ptr = 2**tlvls + 1
           loop_70: do k = 1,curlvl - 1
              curr = ptr + curpbm*2**(curlvl - k) + 2**(curlvl - k - 1) - 1
              psiz1 = prmptr(curr + 1) - prmptr(curr)
              psiz2 = prmptr(curr + 2) - prmptr(curr + 1)
              zptr1 = mid - psiz1
             ! apply givens at curr and curr+1
              do i = givptr(curr),givptr(curr + 1) - 1
                 call la_qrot(1,z(zptr1 + givcol(1,i) - 1),1,z(zptr1 + givcol(2,i) - 1), &
                           1,givnum(1,i),givnum(2,i))
              end do
              do i = givptr(curr + 1),givptr(curr + 2) - 1
                 call la_qrot(1,z(mid - 1 + givcol(1,i)),1,z(mid - 1 + givcol(2,i)),1, &
                           givnum(1,i),givnum(2,i))
              end do
              psiz1 = prmptr(curr + 1) - prmptr(curr)
              psiz2 = prmptr(curr + 2) - prmptr(curr + 1)
              do i = 0,psiz1 - 1
                 ztemp(i + 1) = z(zptr1 + perm(prmptr(curr) + i) - 1)
              end do
              do i = 0,psiz2 - 1
                 ztemp(psiz1 + i + 1) = z(mid + perm(prmptr(curr + 1) + i) - 1)
              end do
              ! multiply blocks at curr and curr+1
              ! determine size of these matrices.  we add half to the value of
              ! the sqrt in case the machine underestimates one of these
              ! square roots.
              bsiz1 = int(half + sqrt(real(qptr(curr + 1) - qptr(curr),KIND=qp)),KIND=ilp)

              bsiz2 = int(half + sqrt(real(qptr(curr + 2) - qptr(curr + 1),KIND=qp)),KIND=ilp)

              if (bsiz1 > 0) then
                 call la_qgemv('T',bsiz1,bsiz1,one,q(qptr(curr)),bsiz1,ztemp(1), &
                           1,zero,z(zptr1),1)
              end if
              call la_qcopy(psiz1 - bsiz1,ztemp(bsiz1 + 1),1,z(zptr1 + bsiz1),1)
              if (bsiz2 > 0) then
                 call la_qgemv('T',bsiz2,bsiz2,one,q(qptr(curr + 1)),bsiz2,ztemp( &
                           psiz1 + 1),1,zero,z(mid),1)
              end if
              call la_qcopy(psiz2 - bsiz2,ztemp(psiz1 + bsiz2 + 1),1,z(mid + bsiz2),1)

              ptr = ptr + 2**(tlvls - k)
           end do loop_70
           return
     end subroutine la_qlaeda

     !> SLAMRG: will create a permutation list which will merge the elements
     !> of A (which is composed of two independently sorted sets) into a
     !> single set which is sorted in ascending order.

     pure subroutine la_slamrg(n1,n2,a,strd1,strd2,index)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n1,n2,strd1,strd2
           ! Array Arguments
           integer(ilp),intent(out) :: index(*)
           real(sp),intent(in) :: a(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ind1,ind2,n1sv,n2sv
           ! Executable Statements
           n1sv = n1
           n2sv = n2
           if (strd1 > 0) then
              ind1 = 1
           else
              ind1 = n1
           end if
           if (strd2 > 0) then
              ind2 = 1 + n1
           else
              ind2 = n1 + n2
           end if
           i = 1
           ! while ( (n1sv > 0)
           10 continue
           if (n1sv > 0 .and. n2sv > 0) then
              if (a(ind1) <= a(ind2)) then
                 index(i) = ind1
                 i = i + 1
                 ind1 = ind1 + strd1
                 n1sv = n1sv - 1
              else
                 index(i) = ind2
                 i = i + 1
                 ind2 = ind2 + strd2
                 n2sv = n2sv - 1
              end if
              go to 10
           end if
           ! end while
           if (n1sv == 0) then
              do n1sv = 1,n2sv
                 index(i) = ind2
                 i = i + 1
                 ind2 = ind2 + strd2
              end do
           else
           ! n2sv == 0
              do n2sv = 1,n1sv
                 index(i) = ind1
                 i = i + 1
                 ind1 = ind1 + strd1
              end do
           end if
           return
     end subroutine la_slamrg
     !> DLAMRG: will create a permutation list which will merge the elements
     !> of A (which is composed of two independently sorted sets) into a
     !> single set which is sorted in ascending order.

     pure subroutine la_dlamrg(n1,n2,a,dtrd1,dtrd2,index)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dtrd1,dtrd2,n1,n2
           ! Array Arguments
           integer(ilp),intent(out) :: index(*)
           real(dp),intent(in) :: a(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ind1,ind2,n1sv,n2sv
           ! Executable Statements
           n1sv = n1
           n2sv = n2
           if (dtrd1 > 0) then
              ind1 = 1
           else
              ind1 = n1
           end if
           if (dtrd2 > 0) then
              ind2 = 1 + n1
           else
              ind2 = n1 + n2
           end if
           i = 1
           ! while ( (n1sv > 0)
           10 continue
           if (n1sv > 0 .and. n2sv > 0) then
              if (a(ind1) <= a(ind2)) then
                 index(i) = ind1
                 i = i + 1
                 ind1 = ind1 + dtrd1
                 n1sv = n1sv - 1
              else
                 index(i) = ind2
                 i = i + 1
                 ind2 = ind2 + dtrd2
                 n2sv = n2sv - 1
              end if
              go to 10
           end if
           ! end while
           if (n1sv == 0) then
              do n1sv = 1,n2sv
                 index(i) = ind2
                 i = i + 1
                 ind2 = ind2 + dtrd2
              end do
           else
           ! n2sv == 0
              do n2sv = 1,n1sv
                 index(i) = ind1
                 i = i + 1
                 ind1 = ind1 + dtrd1
              end do
           end if
           return
     end subroutine la_dlamrg
     !> QLAMRG: will create a permutation list which will merge the elements
     !> of A (which is composed of two independently sorted sets) into a
     !> single set which is sorted in ascending order.

     pure subroutine la_qlamrg(n1,n2,a,qtrd1,qtrd2,index)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: qtrd1,qtrd2,n1,n2
           ! Array Arguments
           integer(ilp),intent(out) :: index(*)
           real(qp),intent(in) :: a(*)
        ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ind1,ind2,n1sv,n2sv
           ! Executable Statements
           n1sv = n1
           n2sv = n2
           if (qtrd1 > 0) then
              ind1 = 1
           else
              ind1 = n1
           end if
           if (qtrd2 > 0) then
              ind2 = 1 + n1
           else
              ind2 = n1 + n2
           end if
           i = 1
           ! while ( (n1sv > 0)
           10 continue
           if (n1sv > 0 .and. n2sv > 0) then
              if (a(ind1) <= a(ind2)) then
                 index(i) = ind1
                 i = i + 1
                 ind1 = ind1 + qtrd1
                 n1sv = n1sv - 1
              else
                 index(i) = ind2
                 i = i + 1
                 ind2 = ind2 + qtrd2
                 n2sv = n2sv - 1
              end if
              go to 10
           end if
           ! end while
           if (n1sv == 0) then
              do n1sv = 1,n2sv
                 index(i) = ind2
                 i = i + 1
                 ind2 = ind2 + qtrd2
              end do
           else
           ! n2sv == 0
              do n2sv = 1,n1sv
                 index(i) = ind1
                 i = i + 1
                 ind1 = ind1 + qtrd1
              end do
           end if
           return
     end subroutine la_qlamrg

     !> Compute the splitting points with threshold SPLTOL.
     !> SLARRA: sets any "small" off-diagonal elements to zero.

     pure subroutine la_slarra(n,d,e,e2,spltol,tnrm,nsplit,isplit,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,nsplit
           integer(ilp),intent(in) :: n
           real(sp),intent(in) :: spltol,tnrm
           ! Array Arguments
           integer(ilp),intent(out) :: isplit(*)
           real(sp),intent(in) :: d(*)
           real(sp),intent(inout) :: e(*),e2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(sp) :: eabs,tmp1
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! compute splitting points
           nsplit = 1
           if (spltol < zero) then
              ! criterion based on absolute off-diagonal value
              tmp1 = abs(spltol)*tnrm
              do i = 1,n - 1
                 eabs = abs(e(i))
                 if (eabs <= tmp1) then
                    e(i) = zero
                    e2(i) = zero
                    isplit(nsplit) = i
                    nsplit = nsplit + 1
                 end if
              end do
           else
              ! criterion that guarantees relative accuracy
              do i = 1,n - 1
                 eabs = abs(e(i))
                 if (eabs <= spltol*sqrt(abs(d(i)))*sqrt(abs(d(i + 1)))) then
                    e(i) = zero
                    e2(i) = zero
                    isplit(nsplit) = i
                    nsplit = nsplit + 1
                 end if
              end do
           end if
           isplit(nsplit) = n
           return
     end subroutine la_slarra
     !> Compute the splitting points with threshold SPLTOL.
     !> DLARRA: sets any "small" off-diagonal elements to zero.

     pure subroutine la_dlarra(n,d,e,e2,spltol,tnrm,nsplit,isplit,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,nsplit
           integer(ilp),intent(in) :: n
           real(dp),intent(in) :: spltol,tnrm
           ! Array Arguments
           integer(ilp),intent(out) :: isplit(*)
           real(dp),intent(in) :: d(*)
           real(dp),intent(inout) :: e(*),e2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(dp) :: eabs,tmp1
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! compute splitting points
           nsplit = 1
           if (spltol < zero) then
              ! criterion based on absolute off-diagonal value
              tmp1 = abs(spltol)*tnrm
              do i = 1,n - 1
                 eabs = abs(e(i))
                 if (eabs <= tmp1) then
                    e(i) = zero
                    e2(i) = zero
                    isplit(nsplit) = i
                    nsplit = nsplit + 1
                 end if
              end do
           else
              ! criterion that guarantees relative accuracy
              do i = 1,n - 1
                 eabs = abs(e(i))
                 if (eabs <= spltol*sqrt(abs(d(i)))*sqrt(abs(d(i + 1)))) then
                    e(i) = zero
                    e2(i) = zero
                    isplit(nsplit) = i
                    nsplit = nsplit + 1
                 end if
              end do
           end if
           isplit(nsplit) = n
           return
     end subroutine la_dlarra
     !> Compute the splitting points with threshold SPLTOL.
     !> QLARRA: sets any "small" off-diagonal elements to zero.

     pure subroutine la_qlarra(n,d,e,e2,spltol,tnrm,nsplit,isplit,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info,nsplit
           integer(ilp),intent(in) :: n
           real(qp),intent(in) :: spltol,tnrm
           ! Array Arguments
           integer(ilp),intent(out) :: isplit(*)
           real(qp),intent(in) :: d(*)
           real(qp),intent(inout) :: e(*),e2(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           real(qp) :: eabs,tmp1
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! compute splitting points
           nsplit = 1
           if (spltol < zero) then
              ! criterion based on absolute off-diagonal value
              tmp1 = abs(spltol)*tnrm
              do i = 1,n - 1
                 eabs = abs(e(i))
                 if (eabs <= tmp1) then
                    e(i) = zero
                    e2(i) = zero
                    isplit(nsplit) = i
                    nsplit = nsplit + 1
                 end if
              end do
           else
              ! criterion that guarantees relative accuracy
              do i = 1,n - 1
                 eabs = abs(e(i))
                 if (eabs <= spltol*sqrt(abs(d(i)))*sqrt(abs(d(i + 1)))) then
                    e(i) = zero
                    e2(i) = zero
                    isplit(nsplit) = i
                    nsplit = nsplit + 1
                 end if
              end do
           end if
           isplit(nsplit) = n
           return
     end subroutine la_qlarra

     !> Find the number of eigenvalues of the symmetric tridiagonal matrix T
     !> that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T
     !> if JOBT = 'L'.

     pure subroutine la_slarrc(jobt,n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info)
        use la_constants_sp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobt
           integer(ilp),intent(out) :: eigcnt,info,lcnt,rcnt
           integer(ilp),intent(in) :: n
           real(sp),intent(in) :: pivmin,vl,vu
           ! Array Arguments
           real(sp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           logical(lk) :: matt
           real(sp) :: lpivot,rpivot,sl,su,tmp,tmp2
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           lcnt = 0
           rcnt = 0
           eigcnt = 0
           matt = la_lsame(jobt,'T')
           if (matt) then
              ! sturm sequence count on t
              lpivot = d(1) - vl
              rpivot = d(1) - vu
              if (lpivot <= zero) then
                 lcnt = lcnt + 1
              end if
              if (rpivot <= zero) then
                 rcnt = rcnt + 1
              end if
              do i = 1,n - 1
                 tmp = e(i)**2
                 lpivot = (d(i + 1) - vl) - tmp/lpivot
                 rpivot = (d(i + 1) - vu) - tmp/rpivot
                 if (lpivot <= zero) then
                    lcnt = lcnt + 1
                 end if
                 if (rpivot <= zero) then
                    rcnt = rcnt + 1
                 end if
              end do
           else
              ! sturm sequence count on l d l^t
              sl = -vl
              su = -vu
              do i = 1,n - 1
                 lpivot = d(i) + sl
                 rpivot = d(i) + su
                 if (lpivot <= zero) then
                    lcnt = lcnt + 1
                 end if
                 if (rpivot <= zero) then
                    rcnt = rcnt + 1
                 end if
                 tmp = e(i)*d(i)*e(i)
                 tmp2 = tmp/lpivot
                 if (tmp2 == zero) then
                    sl = tmp - vl
                 else
                    sl = sl*tmp2 - vl
                 end if
                 tmp2 = tmp/rpivot
                 if (tmp2 == zero) then
                    su = tmp - vu
                 else
                    su = su*tmp2 - vu
                 end if
              end do
              lpivot = d(n) + sl
              rpivot = d(n) + su
              if (lpivot <= zero) then
                 lcnt = lcnt + 1
              end if
              if (rpivot <= zero) then
                 rcnt = rcnt + 1
              end if
           end if
           eigcnt = rcnt - lcnt
           return
     end subroutine la_slarrc
     !> Find the number of eigenvalues of the symmetric tridiagonal matrix T
     !> that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T
     !> if JOBT = 'L'.

     pure subroutine la_dlarrc(jobt,n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info)
        use la_constants_dp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobt
           integer(ilp),intent(out) :: eigcnt,info,lcnt,rcnt
           integer(ilp),intent(in) :: n
           real(dp),intent(in) :: pivmin,vl,vu
           ! Array Arguments
           real(dp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           logical(lk) :: matt
           real(dp) :: lpivot,rpivot,sl,su,tmp,tmp2
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           lcnt = 0
           rcnt = 0
           eigcnt = 0
           matt = la_lsame(jobt,'T')
           if (matt) then
              ! sturm sequence count on t
              lpivot = d(1) - vl
              rpivot = d(1) - vu
              if (lpivot <= zero) then
                 lcnt = lcnt + 1
              end if
              if (rpivot <= zero) then
                 rcnt = rcnt + 1
              end if
              do i = 1,n - 1
                 tmp = e(i)**2
                 lpivot = (d(i + 1) - vl) - tmp/lpivot
                 rpivot = (d(i + 1) - vu) - tmp/rpivot
                 if (lpivot <= zero) then
                    lcnt = lcnt + 1
                 end if
                 if (rpivot <= zero) then
                    rcnt = rcnt + 1
                 end if
              end do
           else
              ! sturm sequence count on l d l^t
              sl = -vl
              su = -vu
              do i = 1,n - 1
                 lpivot = d(i) + sl
                 rpivot = d(i) + su
                 if (lpivot <= zero) then
                    lcnt = lcnt + 1
                 end if
                 if (rpivot <= zero) then
                    rcnt = rcnt + 1
                 end if
                 tmp = e(i)*d(i)*e(i)
                 tmp2 = tmp/lpivot
                 if (tmp2 == zero) then
                    sl = tmp - vl
                 else
                    sl = sl*tmp2 - vl
                 end if
                 tmp2 = tmp/rpivot
                 if (tmp2 == zero) then
                    su = tmp - vu
                 else
                    su = su*tmp2 - vu
                 end if
              end do
              lpivot = d(n) + sl
              rpivot = d(n) + su
              if (lpivot <= zero) then
                 lcnt = lcnt + 1
              end if
              if (rpivot <= zero) then
                 rcnt = rcnt + 1
              end if
           end if
           eigcnt = rcnt - lcnt
           return
     end subroutine la_dlarrc
     !> Find the number of eigenvalues of the symmetric tridiagonal matrix T
     !> that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T
     !> if JOBT = 'L'.

     pure subroutine la_qlarrc(jobt,n,vl,vu,d,e,pivmin,eigcnt,lcnt,rcnt,info)
        use la_constants_qp

        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: jobt
           integer(ilp),intent(out) :: eigcnt,info,lcnt,rcnt
           integer(ilp),intent(in) :: n
           real(qp),intent(in) :: pivmin,vl,vu
           ! Array Arguments
           real(qp),intent(in) :: d(*),e(*)
        ! =====================================================================

           ! Local Scalars
           integer(ilp) :: i
           logical(lk) :: matt
           real(qp) :: lpivot,rpivot,sl,su,tmp,tmp2
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           lcnt = 0
           rcnt = 0
           eigcnt = 0
           matt = la_lsame(jobt,'T')
           if (matt) then
              ! sturm sequence count on t
              lpivot = d(1) - vl
              rpivot = d(1) - vu
              if (lpivot <= zero) then
                 lcnt = lcnt + 1
              end if
              if (rpivot <= zero) then
                 rcnt = rcnt + 1
              end if
              do i = 1,n - 1
                 tmp = e(i)**2
                 lpivot = (d(i + 1) - vl) - tmp/lpivot
                 rpivot = (d(i + 1) - vu) - tmp/rpivot
                 if (lpivot <= zero) then
                    lcnt = lcnt + 1
                 end if
                 if (rpivot <= zero) then
                    rcnt = rcnt + 1
                 end if
              end do
           else
              ! sturm sequence count on l d l^t
              sl = -vl
              su = -vu
              do i = 1,n - 1
                 lpivot = d(i) + sl
                 rpivot = d(i) + su
                 if (lpivot <= zero) then
                    lcnt = lcnt + 1
                 end if
                 if (rpivot <= zero) then
                    rcnt = rcnt + 1
                 end if
                 tmp = e(i)*d(i)*e(i)
                 tmp2 = tmp/lpivot
                 if (tmp2 == zero) then
                    sl = tmp - vl
                 else
                    sl = sl*tmp2 - vl
                 end if
                 tmp2 = tmp/rpivot
                 if (tmp2 == zero) then
                    su = tmp - vu
                 else
                    su = su*tmp2 - vu
                 end if
              end do
              lpivot = d(n) + sl
              rpivot = d(n) + su
              if (lpivot <= zero) then
                 lcnt = lcnt + 1
              end if
              if (rpivot <= zero) then
                 rcnt = rcnt + 1
              end if
           end if
           eigcnt = rcnt - lcnt
           return
     end subroutine la_qlarrc

     !> SLARRD: computes the eigenvalues of a symmetric tridiagonal
     !> matrix T to suitable accuracy. This is an auxiliary code to be
     !> called from SSTEMR.
     !> The user may ask for all eigenvalues, all eigenvalues
     !> in the half-open interval (VL, VU], or the IL-th through IU-th
     !> eigenvalues.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_slarrd(range,order,n,vl,vu,il,iu,gers,reltol,d,e,e2, &
               pivmin,nsplit,isplit,m,w,werr,wl,wu,iblock,indexw,work,iwork,info)
        use la_constants_sp,only:zero,half,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: order,range
           integer(ilp),intent(in) :: il,iu,n,nsplit
           integer(ilp),intent(out) :: info,m
           real(sp),intent(in) :: pivmin,reltol,vl,vu
           real(sp),intent(out) :: wl,wu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),indexw(*),iwork(*)
           integer(ilp),intent(in) :: isplit(*)
           real(sp),intent(in) :: d(*),e(*),e2(*),gers(*)
           real(sp),intent(out) :: w(*),werr(*),work(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: fudge = two
           integer(ilp),parameter :: allrng = 1
           integer(ilp),parameter :: valrng = 2
           integer(ilp),parameter :: indrng = 3

           ! Local Scalars
           logical(lk) :: ncnvrg,toofew
           integer(ilp) :: i,ib,ibegin,idiscl,idiscu,ie,iend,iinfo,im,in,ioff,iout, &
                     irange,itmax,itmp1,itmp2,iw,iwoff,j,jblk,jdisc,je,jee,nb,nwl,nwu
           real(sp) :: atoli,eps,gl,gu,rtoli,tmp1,tmp2,tnorm,uflow,wkill,wlu, &
                     wul
           ! Local Arrays
           integer(ilp) :: idumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! decode range
           if (la_lsame(range,'A')) then
              irange = allrng
           else if (la_lsame(range,'V')) then
              irange = valrng
           else if (la_lsame(range,'I')) then
              irange = indrng
           else
              irange = 0
           end if
           ! check for errors
           if (irange <= 0) then
              info = -1
           else if (.not. (la_lsame(order,'B') .or. la_lsame(order,'E'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (irange == valrng) then
              if (vl >= vu) info = -5
           else if (irange == indrng .and. (il < 1 .or. il > max(1,n))) then
              info = -6
           else if (irange == indrng .and. (iu < min(n,il) .or. iu > n)) then
              info = -7
           end if
           if (info /= 0) then
              return
           end if
           ! initialize error flags
           info = 0
           ncnvrg = .false.
           toofew = .false.
           ! quick return if possible
           m = 0
           if (n == 0) return
           ! simplification:
           if (irange == indrng .and. il == 1 .and. iu == n) irange = 1
           ! get machine constants
           eps = la_slamch('P')
           uflow = la_slamch('U')
           ! special case when n=1
           ! treat case of 1x1 matrix for quick return
           if (n == 1) then
              if ((irange == allrng) .or. ((irange == valrng) .and. (d(1) > vl) .and. (d(1) <= vu)) .or. (( &
                        irange == indrng) .and. (il == 1) .and. (iu == 1))) then
                 m = 1
                 w(1) = d(1)
                 ! the computation error of the eigenvalue is zero
                 werr(1) = zero
                 iblock(1) = 1
                 indexw(1) = 1
              end if
              return
           end if
           ! nb is the minimum vector length for vector bisection, or 0
           ! if only scalar is to be done.
           nb = la_ilaenv(1,'SSTEBZ',' ',n,-1,-1,-1)
           if (nb <= 1) nb = 0
           ! find global spectral radius
           gl = d(1)
           gu = d(1)
           do i = 1,n
              gl = min(gl,gers(2*i - 1))
              gu = max(gu,gers(2*i))
           end do
           ! compute global gerschgorin bounds and spectral diameter
           tnorm = max(abs(gl),abs(gu))
           gl = gl - fudge*tnorm*eps*n - fudge*two*pivmin
           gu = gu + fudge*tnorm*eps*n + fudge*two*pivmin
           ! [jan/28/2009] remove the line below since spdiam variable not use
           ! spdiam = gu - gl
           ! input arguments for la_slaebz:
           ! the relative tolerance.  an interval (a,b] lies within
           ! "relative tolerance" if  b-a < reltol*max(|a|,|b|),
           rtoli = reltol
           ! set the absolute tolerance for interval convergence to zero to force
           ! interval convergence based on relative size of the interval.
           ! this is dangerous because intervals might not converge when reltol is
           ! small. but at least a very small number should be selected so that for
           ! strongly graded matrices, the code can get relatively accurate
           ! eigenvalues.
           atoli = fudge*two*uflow + fudge*two*pivmin
           if (irange == indrng) then
              ! range='i': compute an interval containing eigenvalues
              ! il through iu. the initial interval [gl,gu] from the global
              ! gerschgorin bounds gl and gu is refined by la_slaebz.
              itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
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
              call la_slaebz(3,itmax,n,2,2,nb,atoli,rtoli,pivmin,d,e,e2,iwork(5) &
                        ,work(n + 1),work(n + 5),iout,iwork,w,iblock,iinfo)
              if (iinfo /= 0) then
                 info = iinfo
                 return
              end if
              ! on exit, output intervals may not be ordered by ascending negcount
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
              ! on exit, the interval [wl, wlu] contains a value with negcount nwl,
              ! and [wul, wu] contains a value with negcount nwu.
              if (nwl < 0 .or. nwl >= n .or. nwu < 1 .or. nwu > n) then
                 info = 4
                 return
              end if
           elseif (irange == valrng) then
              wl = vl
              wu = vu
           elseif (irange == allrng) then
              wl = gl
              wu = gu
           end if
           ! find eigenvalues -- loop over blocks and recompute nwl and nwu.
           ! nwl accumulates the number of eigenvalues .le. wl,
           ! nwu accumulates the number of eigenvalues .le. wu
           m = 0
           iend = 0
           info = 0
           nwl = 0
           nwu = 0
           loop_70: do jblk = 1,nsplit
              ioff = iend
              ibegin = ioff + 1
              iend = isplit(jblk)
              in = iend - ioff
              if (in == 1) then
                 ! 1x1 block
                 if (wl >= d(ibegin) - pivmin) nwl = nwl + 1
                 if (wu >= d(ibegin) - pivmin) nwu = nwu + 1
                 if (irange == allrng .or. (wl < d(ibegin) - pivmin .and. wu >= d(ibegin) - pivmin)) &
                           then
                    m = m + 1
                    w(m) = d(ibegin)
                    werr(m) = zero
                    ! the gap for a single block doesn't matter for the later
                    ! algorithm and is assigned an arbitrary large value
                    iblock(m) = jblk
                    indexw(m) = 1
                 end if
              ! disabled 2x2 case because of a failure on the following matrix
              ! range = 'i', il = iu = 4
                ! original tridiagonal, d = [
                 ! -0.150102010615740e+00_sp
                 ! -0.849897989384260e+00_sp
                 ! -0.128208148052635e-15_sp
                  ! 0.128257718286320e-15_sp
                ! ];
                ! e = [
                 ! -0.357171383266986e+00_sp
                 ! -0.180411241501588e-15_sp
                 ! -0.175152352710251e-15_sp
                ! ];
               ! else if( in==2 ) then
      ! *           2x2 block
                  ! disc = sqrt( (half*(d(ibegin)-d(iend)))**2 + e(ibegin)**2 )
                  ! tmp1 = half*(d(ibegin)+d(iend))
                  ! l1 = tmp1 - disc
                  ! if( wl>= l1-pivmin )
           ! $         nwl = nwl + 1
                  ! if( wu>= l1-pivmin )
           ! $         nwu = nwu + 1
                  ! if( irange==allrng .or. ( wl<l1-pivmin .and. wu>=
           ! $          l1-pivmin ) ) then
                     ! m = m + 1
                     ! w( m ) = l1
      ! *              the uncertainty of eigenvalues of a 2x2 matrix is very small
                     ! werr( m ) = eps * abs( w( m ) ) * two
                     ! iblock( m ) = jblk
                     ! indexw( m ) = 1
                  ! endif
                  ! l2 = tmp1 + disc
                  ! if( wl>= l2-pivmin )
           ! $         nwl = nwl + 1
                  ! if( wu>= l2-pivmin )
           ! $         nwu = nwu + 1
                  ! if( irange==allrng .or. ( wl<l2-pivmin .and. wu>=
           ! $          l2-pivmin ) ) then
                     ! m = m + 1
                     ! w( m ) = l2
      ! *              the uncertainty of eigenvalues of a 2x2 matrix is very small
                     ! werr( m ) = eps * abs( w( m ) ) * two
                     ! iblock( m ) = jblk
                     ! indexw( m ) = 2
                  ! endif
              else
                 ! general case - block of size in >= 2
                 ! compute local gerschgorin interval and use it as the initial
                 ! interval for la_slaebz
                 gu = d(ibegin)
                 gl = d(ibegin)
                 tmp1 = zero
                 do j = ibegin,iend
                    gl = min(gl,gers(2*j - 1))
                    gu = max(gu,gers(2*j))
                 end do
                 ! [jan/28/2009]
                 ! change spdiam by tnorm in lines 2 and 3 thereafter
                 ! line 1: remove computation of spdiam (not useful anymore)
                 ! spdiam = gu - gl
                 ! gl = gl - fudge*spdiam*eps*in - fudge*pivmin
                 ! gu = gu + fudge*spdiam*eps*in + fudge*pivmin
                 gl = gl - fudge*tnorm*eps*in - fudge*pivmin
                 gu = gu + fudge*tnorm*eps*in + fudge*pivmin
                 if (irange > 1) then
                    if (gu < wl) then
                       ! the local block contains none of the wanted eigenvalues
                       nwl = nwl + in
                       nwu = nwu + in
                       cycle loop_70
                    end if
                    ! refine search interval if possible, only range (wl,wu] matters
                    gl = max(gl,wl)
                    gu = min(gu,wu)
                    if (gl >= gu) cycle loop_70
                 end if
                 ! find negcount of initial interval boundaries gl and gu
                 work(n + 1) = gl
                 work(n + in + 1) = gu
                 call la_slaebz(1,0,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                 ibegin),e2(ibegin),idumma,work(n + 1),work(n + 2*in + 1),im,iwork,w(m + 1), &
                            iblock(m + 1),iinfo)
                 if (iinfo /= 0) then
                    info = iinfo
                    return
                 end if
                 nwl = nwl + iwork(1)
                 nwu = nwu + iwork(in + 1)
                 iwoff = m - iwork(1)
                 ! compute eigenvalues
                 itmax = int((log(gu - gl + pivmin) - log(pivmin))/log(two),KIND=ilp) + &
                           2
                 call la_slaebz(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                  ibegin),e2(ibegin),idumma,work(n + 1),work(n + 2*in + 1),iout,iwork,w(m + &
                            1),iblock(m + 1),iinfo)
                 if (iinfo /= 0) then
                    info = iinfo
                    return
                 end if
                 ! copy eigenvalues into w and iblock
                 ! use -jblk for block number for unconverged eigenvalues.
                 ! loop over the number of output intervals from la_slaebz
                 do j = 1,iout
                    ! eigenvalue approximation is middle point of interval
                    tmp1 = half*(work(j + n) + work(j + in + n))
                    ! semi length of error interval
                    tmp2 = half*abs(work(j + n) - work(j + in + n))
                    if (j > iout - iinfo) then
                       ! flag non-convergence.
                       ncnvrg = .true.
                       ib = -jblk
                    else
                       ib = jblk
                    end if
                    do je = iwork(j) + 1 + iwoff,iwork(j + in) + iwoff
                       w(je) = tmp1
                       werr(je) = tmp2
                       indexw(je) = je - iwoff
                       iblock(je) = ib
                    end do
                 end do
                 m = m + im
              end if
           end do loop_70
           ! if range='i', then (wl,wu) contains eigenvalues nwl+1,...,nwu
           ! if nwl+1 < il or nwu > iu, discard extra eigenvalues.
           if (irange == indrng) then
              idiscl = il - 1 - nwl
              idiscu = nwu - iu
              if (idiscl > 0) then
                 im = 0
                 do je = 1,m
                    ! remove some of the smallest eigenvalues from the left so that
                    ! at the end idiscl =0. move all eigenvalues up to the left.
                    if (w(je) <= wlu .and. idiscl > 0) then
                       idiscl = idiscl - 1
                    else
                       im = im + 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscu > 0) then
                 ! remove some of the largest eigenvalues from the right so that
                 ! at the end idiscu =0. move all eigenvalues up to the left.
                 im = m + 1
                 do je = m,1,-1
                    if (w(je) >= wul .and. idiscu > 0) then
                       idiscu = idiscu - 1
                    else
                       im = im - 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 jee = 0
                 do je = im,m
                    jee = jee + 1
                    w(jee) = w(je)
                    werr(jee) = werr(je)
                    indexw(jee) = indexw(je)
                    iblock(jee) = iblock(je)
                 end do
                 m = m - im + 1
              end if
              if (idiscl > 0 .or. idiscu > 0) then
                 ! code to deal with effects of bad arithmetic. (if n(w) is
                 ! monotone non-decreasing, this should never happen.)
                 ! some low eigenvalues to be discarded are not in (wl,wlu],
                 ! or high eigenvalues to be discarded are not in (wul,wu]
                 ! so just kill off the smallest idiscl/largest idiscu
                 ! eigenvalues, by marking the corresponding iblock = 0
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
                          if (iblock(je) /= 0 .and. (w(je) >= wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 ! now erase all eigenvalues with iblock set to zero
                 im = 0
                 do je = 1,m
                    if (iblock(je) /= 0) then
                       im = im + 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl < 0 .or. idiscu < 0) then
                 toofew = .true.
              end if
           end if
           if ((irange == allrng .and. m /= n) .or. (irange == indrng .and. m /= iu - il + 1)) then
              toofew = .true.
           end if
           ! if order='b', do nothing the eigenvalues are already sorted by
              ! block.
           ! if order='e', sort the eigenvalues from smallest to largest
           if (la_lsame(order,'E') .and. nsplit > 1) then
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
                    tmp2 = werr(ie)
                    itmp1 = iblock(ie)
                    itmp2 = indexw(ie)
                    w(ie) = w(je)
                    werr(ie) = werr(je)
                    iblock(ie) = iblock(je)
                    indexw(ie) = indexw(je)
                    w(je) = tmp1
                    werr(je) = tmp2
                    iblock(je) = itmp1
                    indexw(je) = itmp2
                 end if
              end do
           end if
           info = 0
           if (ncnvrg) info = info + 1
           if (toofew) info = info + 2
           return
     end subroutine la_slarrd
     !> DLARRD: computes the eigenvalues of a symmetric tridiagonal
     !> matrix T to suitable accuracy. This is an auxiliary code to be
     !> called from DSTEMR.
     !> The user may ask for all eigenvalues, all eigenvalues
     !> in the half-open interval (VL, VU], or the IL-th through IU-th
     !> eigenvalues.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_dlarrd(range,order,n,vl,vu,il,iu,gers,reltol,d,e,e2, &
               pivmin,nsplit,isplit,m,w,werr,wl,wu,iblock,indexw,work,iwork,info)
        use la_constants_dp,only:zero,half,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: order,range
           integer(ilp),intent(in) :: il,iu,n,nsplit
           integer(ilp),intent(out) :: info,m
           real(dp),intent(in) :: pivmin,reltol,vl,vu
           real(dp),intent(out) :: wl,wu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),indexw(*),iwork(*)
           integer(ilp),intent(in) :: isplit(*)
           real(dp),intent(in) :: d(*),e(*),e2(*),gers(*)
           real(dp),intent(out) :: w(*),werr(*),work(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: fudge = two
           integer(ilp),parameter :: allrng = 1
           integer(ilp),parameter :: valrng = 2
           integer(ilp),parameter :: indrng = 3

           ! Local Scalars
           logical(lk) :: ncnvrg,toofew
           integer(ilp) :: i,ib,ibegin,idiscl,idiscu,ie,iend,iinfo,im,in,ioff,iout, &
                     irange,itmax,itmp1,itmp2,iw,iwoff,j,jblk,jdisc,je,jee,nb,nwl,nwu
           real(dp) :: atoli,eps,gl,gu,rtoli,tmp1,tmp2,tnorm,uflow,wkill,wlu, &
                     wul
           ! Local Arrays
           integer(ilp) :: idumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! decode range
           if (la_lsame(range,'A')) then
              irange = allrng
           else if (la_lsame(range,'V')) then
              irange = valrng
           else if (la_lsame(range,'I')) then
              irange = indrng
           else
              irange = 0
           end if
           ! check for errors
           if (irange <= 0) then
              info = -1
           else if (.not. (la_lsame(order,'B') .or. la_lsame(order,'E'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (irange == valrng) then
              if (vl >= vu) info = -5
           else if (irange == indrng .and. (il < 1 .or. il > max(1,n))) then
              info = -6
           else if (irange == indrng .and. (iu < min(n,il) .or. iu > n)) then
              info = -7
           end if
           if (info /= 0) then
              return
           end if
           ! initialize error flags
           info = 0
           ncnvrg = .false.
           toofew = .false.
           ! quick return if possible
           m = 0
           if (n == 0) return
           ! simplification:
           if (irange == indrng .and. il == 1 .and. iu == n) irange = 1
           ! get machine constants
           eps = la_dlamch('P')
           uflow = la_dlamch('U')
           ! special case when n=1
           ! treat case of 1x1 matrix for quick return
           if (n == 1) then
              if ((irange == allrng) .or. ((irange == valrng) .and. (d(1) > vl) .and. (d(1) <= vu)) .or. (( &
                        irange == indrng) .and. (il == 1) .and. (iu == 1))) then
                 m = 1
                 w(1) = d(1)
                 ! the computation error of the eigenvalue is zero
                 werr(1) = zero
                 iblock(1) = 1
                 indexw(1) = 1
              end if
              return
           end if
           ! nb is the minimum vector length for vector bisection, or 0
           ! if only scalar is to be done.
           nb = la_ilaenv(1,'DSTEBZ',' ',n,-1,-1,-1)
           if (nb <= 1) nb = 0
           ! find global spectral radius
           gl = d(1)
           gu = d(1)
           do i = 1,n
              gl = min(gl,gers(2*i - 1))
              gu = max(gu,gers(2*i))
           end do
           ! compute global gerschgorin bounds and spectral diameter
           tnorm = max(abs(gl),abs(gu))
           gl = gl - fudge*tnorm*eps*n - fudge*two*pivmin
           gu = gu + fudge*tnorm*eps*n + fudge*two*pivmin
           ! [jan/28/2009] remove the line below since spdiam variable not use
           ! spdiam = gu - gl
           ! input arguments for la_dlaebz:
           ! the relative tolerance.  an interval (a,b] lies within
           ! "relative tolerance" if  b-a < reltol*max(|a|,|b|),
           rtoli = reltol
           ! set the absolute tolerance for interval convergence to zero to force
           ! interval convergence based on relative size of the interval.
           ! this is dangerous because intervals might not converge when reltol is
           ! small. but at least a very small number should be selected so that for
           ! strongly graded matrices, the code can get relatively accurate
           ! eigenvalues.
           atoli = fudge*two*uflow + fudge*two*pivmin
           if (irange == indrng) then
              ! range='i': compute an interval containing eigenvalues
              ! il through iu. the initial interval [gl,gu] from the global
              ! gerschgorin bounds gl and gu is refined by la_dlaebz.
              itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
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
              call la_dlaebz(3,itmax,n,2,2,nb,atoli,rtoli,pivmin,d,e,e2,iwork(5) &
                        ,work(n + 1),work(n + 5),iout,iwork,w,iblock,iinfo)
              if (iinfo /= 0) then
                 info = iinfo
                 return
              end if
              ! on exit, output intervals may not be ordered by ascending negcount
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
              ! on exit, the interval [wl, wlu] contains a value with negcount nwl,
              ! and [wul, wu] contains a value with negcount nwu.
              if (nwl < 0 .or. nwl >= n .or. nwu < 1 .or. nwu > n) then
                 info = 4
                 return
              end if
           elseif (irange == valrng) then
              wl = vl
              wu = vu
           elseif (irange == allrng) then
              wl = gl
              wu = gu
           end if
           ! find eigenvalues -- loop over blocks and recompute nwl and nwu.
           ! nwl accumulates the number of eigenvalues .le. wl,
           ! nwu accumulates the number of eigenvalues .le. wu
           m = 0
           iend = 0
           info = 0
           nwl = 0
           nwu = 0
           loop_70: do jblk = 1,nsplit
              ioff = iend
              ibegin = ioff + 1
              iend = isplit(jblk)
              in = iend - ioff
              if (in == 1) then
                 ! 1x1 block
                 if (wl >= d(ibegin) - pivmin) nwl = nwl + 1
                 if (wu >= d(ibegin) - pivmin) nwu = nwu + 1
                 if (irange == allrng .or. (wl < d(ibegin) - pivmin .and. wu >= d(ibegin) - pivmin)) &
                           then
                    m = m + 1
                    w(m) = d(ibegin)
                    werr(m) = zero
                    ! the gap for a single block doesn't matter for the later
                    ! algorithm and is assigned an arbitrary large value
                    iblock(m) = jblk
                    indexw(m) = 1
                 end if
              ! disabled 2x2 case because of a failure on the following matrix
              ! range = 'i', il = iu = 4
                ! original tridiagonal, d = [
                 ! -0.150102010615740e+00_dp
                 ! -0.849897989384260e+00_dp
                 ! -0.128208148052635e-15_dp
                  ! 0.128257718286320e-15_dp
                ! ];
                ! e = [
                 ! -0.357171383266986e+00_dp
                 ! -0.180411241501588e-15_dp
                 ! -0.175152352710251e-15_dp
                ! ];
               ! else if( in==2 ) then
      ! *           2x2 block
                  ! disc = sqrt( (half*(d(ibegin)-d(iend)))**2 + e(ibegin)**2 )
                  ! tmp1 = half*(d(ibegin)+d(iend))
                  ! l1 = tmp1 - disc
                  ! if( wl>= l1-pivmin )
           ! $         nwl = nwl + 1
                  ! if( wu>= l1-pivmin )
           ! $         nwu = nwu + 1
                  ! if( irange==allrng .or. ( wl<l1-pivmin .and. wu>=
           ! $          l1-pivmin ) ) then
                     ! m = m + 1
                     ! w( m ) = l1
      ! *              the uncertainty of eigenvalues of a 2x2 matrix is very small
                     ! werr( m ) = eps * abs( w( m ) ) * two
                     ! iblock( m ) = jblk
                     ! indexw( m ) = 1
                  ! endif
                  ! l2 = tmp1 + disc
                  ! if( wl>= l2-pivmin )
           ! $         nwl = nwl + 1
                  ! if( wu>= l2-pivmin )
           ! $         nwu = nwu + 1
                  ! if( irange==allrng .or. ( wl<l2-pivmin .and. wu>=
           ! $          l2-pivmin ) ) then
                     ! m = m + 1
                     ! w( m ) = l2
      ! *              the uncertainty of eigenvalues of a 2x2 matrix is very small
                     ! werr( m ) = eps * abs( w( m ) ) * two
                     ! iblock( m ) = jblk
                     ! indexw( m ) = 2
                  ! endif
              else
                 ! general case - block of size in >= 2
                 ! compute local gerschgorin interval and use it as the initial
                 ! interval for la_dlaebz
                 gu = d(ibegin)
                 gl = d(ibegin)
                 tmp1 = zero
                 do j = ibegin,iend
                    gl = min(gl,gers(2*j - 1))
                    gu = max(gu,gers(2*j))
                 end do
                 ! [jan/28/2009]
                 ! change spdiam by tnorm in lines 2 and 3 thereafter
                 ! line 1: remove computation of spdiam (not useful anymore)
                 ! spdiam = gu - gl
                 ! gl = gl - fudge*spdiam*eps*in - fudge*pivmin
                 ! gu = gu + fudge*spdiam*eps*in + fudge*pivmin
                 gl = gl - fudge*tnorm*eps*in - fudge*pivmin
                 gu = gu + fudge*tnorm*eps*in + fudge*pivmin
                 if (irange > 1) then
                    if (gu < wl) then
                       ! the local block contains none of the wanted eigenvalues
                       nwl = nwl + in
                       nwu = nwu + in
                       cycle loop_70
                    end if
                    ! refine search interval if possible, only range (wl,wu] matters
                    gl = max(gl,wl)
                    gu = min(gu,wu)
                    if (gl >= gu) cycle loop_70
                 end if
                 ! find negcount of initial interval boundaries gl and gu
                 work(n + 1) = gl
                 work(n + in + 1) = gu
                 call la_dlaebz(1,0,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                 ibegin),e2(ibegin),idumma,work(n + 1),work(n + 2*in + 1),im,iwork,w(m + 1), &
                            iblock(m + 1),iinfo)
                 if (iinfo /= 0) then
                    info = iinfo
                    return
                 end if
                 nwl = nwl + iwork(1)
                 nwu = nwu + iwork(in + 1)
                 iwoff = m - iwork(1)
                 ! compute eigenvalues
                 itmax = int((log(gu - gl + pivmin) - log(pivmin))/log(two),KIND=ilp) + &
                           2
                 call la_dlaebz(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                  ibegin),e2(ibegin),idumma,work(n + 1),work(n + 2*in + 1),iout,iwork,w(m + &
                            1),iblock(m + 1),iinfo)
                 if (iinfo /= 0) then
                    info = iinfo
                    return
                 end if
                 ! copy eigenvalues into w and iblock
                 ! use -jblk for block number for unconverged eigenvalues.
                 ! loop over the number of output intervals from la_dlaebz
                 do j = 1,iout
                    ! eigenvalue approximation is middle point of interval
                    tmp1 = half*(work(j + n) + work(j + in + n))
                    ! semi length of error interval
                    tmp2 = half*abs(work(j + n) - work(j + in + n))
                    if (j > iout - iinfo) then
                       ! flag non-convergence.
                       ncnvrg = .true.
                       ib = -jblk
                    else
                       ib = jblk
                    end if
                    do je = iwork(j) + 1 + iwoff,iwork(j + in) + iwoff
                       w(je) = tmp1
                       werr(je) = tmp2
                       indexw(je) = je - iwoff
                       iblock(je) = ib
                    end do
                 end do
                 m = m + im
              end if
           end do loop_70
           ! if range='i', then (wl,wu) contains eigenvalues nwl+1,...,nwu
           ! if nwl+1 < il or nwu > iu, discard extra eigenvalues.
           if (irange == indrng) then
              idiscl = il - 1 - nwl
              idiscu = nwu - iu
              if (idiscl > 0) then
                 im = 0
                 do je = 1,m
                    ! remove some of the smallest eigenvalues from the left so that
                    ! at the end idiscl =0. move all eigenvalues up to the left.
                    if (w(je) <= wlu .and. idiscl > 0) then
                       idiscl = idiscl - 1
                    else
                       im = im + 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscu > 0) then
                 ! remove some of the largest eigenvalues from the right so that
                 ! at the end idiscu =0. move all eigenvalues up to the left.
                 im = m + 1
                 do je = m,1,-1
                    if (w(je) >= wul .and. idiscu > 0) then
                       idiscu = idiscu - 1
                    else
                       im = im - 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 jee = 0
                 do je = im,m
                    jee = jee + 1
                    w(jee) = w(je)
                    werr(jee) = werr(je)
                    indexw(jee) = indexw(je)
                    iblock(jee) = iblock(je)
                 end do
                 m = m - im + 1
              end if
              if (idiscl > 0 .or. idiscu > 0) then
                 ! code to deal with effects of bad arithmetic. (if n(w) is
                 ! monotone non-decreasing, this should never happen.)
                 ! some low eigenvalues to be discarded are not in (wl,wlu],
                 ! or high eigenvalues to be discarded are not in (wul,wu]
                 ! so just kill off the smallest idiscl/largest idiscu
                 ! eigenvalues, by marking the corresponding iblock = 0
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
                          if (iblock(je) /= 0 .and. (w(je) >= wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 ! now erase all eigenvalues with iblock set to zero
                 im = 0
                 do je = 1,m
                    if (iblock(je) /= 0) then
                       im = im + 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl < 0 .or. idiscu < 0) then
                 toofew = .true.
              end if
           end if
           if ((irange == allrng .and. m /= n) .or. (irange == indrng .and. m /= iu - il + 1)) then
              toofew = .true.
           end if
           ! if order='b', do nothing the eigenvalues are already sorted by
              ! block.
           ! if order='e', sort the eigenvalues from smallest to largest
           if (la_lsame(order,'E') .and. nsplit > 1) then
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
                    tmp2 = werr(ie)
                    itmp1 = iblock(ie)
                    itmp2 = indexw(ie)
                    w(ie) = w(je)
                    werr(ie) = werr(je)
                    iblock(ie) = iblock(je)
                    indexw(ie) = indexw(je)
                    w(je) = tmp1
                    werr(je) = tmp2
                    iblock(je) = itmp1
                    indexw(je) = itmp2
                 end if
              end do
           end if
           info = 0
           if (ncnvrg) info = info + 1
           if (toofew) info = info + 2
           return
     end subroutine la_dlarrd
     !> QLARRD: computes the eigenvalues of a symmetric tridiagonal
     !> matrix T to suitable accuracy. This is an auxiliary code to be
     !> called from QSTEMR.
     !> The user may ask for all eigenvalues, all eigenvalues
     !> in the half-open interval (VL, VU], or the IL-th through IU-th
     !> eigenvalues.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_qlarrd(range,order,n,vl,vu,il,iu,gers,reltol,d,e,e2, &
               pivmin,nsplit,isplit,m,w,werr,wl,wu,iblock,indexw,work,iwork,info)
        use la_constants_qp,only:zero,half,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: order,range
           integer(ilp),intent(in) :: il,iu,n,nsplit
           integer(ilp),intent(out) :: info,m
           real(qp),intent(in) :: pivmin,reltol,vl,vu
           real(qp),intent(out) :: wl,wu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),indexw(*),iwork(*)
           integer(ilp),intent(in) :: isplit(*)
           real(qp),intent(in) :: d(*),e(*),e2(*),gers(*)
           real(qp),intent(out) :: w(*),werr(*),work(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: fudge = two
           integer(ilp),parameter :: allrng = 1
           integer(ilp),parameter :: valrng = 2
           integer(ilp),parameter :: indrng = 3

           ! Local Scalars
           logical(lk) :: ncnvrg,toofew
           integer(ilp) :: i,ib,ibegin,idiscl,idiscu,ie,iend,iinfo,im,in,ioff,iout, &
                     irange,itmax,itmp1,itmp2,iw,iwoff,j,jblk,jdisc,je,jee,nb,nwl,nwu
           real(qp) :: atoli,eps,gl,gu,rtoli,tmp1,tmp2,tnorm,uflow,wkill,wlu, &
                     wul
           ! Local Arrays
           integer(ilp) :: idumma(1)
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! decode range
           if (la_lsame(range,'A')) then
              irange = allrng
           else if (la_lsame(range,'V')) then
              irange = valrng
           else if (la_lsame(range,'I')) then
              irange = indrng
           else
              irange = 0
           end if
           ! check for errors
           if (irange <= 0) then
              info = -1
           else if (.not. (la_lsame(order,'B') .or. la_lsame(order,'E'))) then
              info = -2
           else if (n < 0) then
              info = -3
           else if (irange == valrng) then
              if (vl >= vu) info = -5
           else if (irange == indrng .and. (il < 1 .or. il > max(1,n))) then
              info = -6
           else if (irange == indrng .and. (iu < min(n,il) .or. iu > n)) then
              info = -7
           end if
           if (info /= 0) then
              return
           end if
           ! initialize error flags
           info = 0
           ncnvrg = .false.
           toofew = .false.
           ! quick return if possible
           m = 0
           if (n == 0) return
           ! simplification:
           if (irange == indrng .and. il == 1 .and. iu == n) irange = 1
           ! get machine constants
           eps = la_qlamch('P')
           uflow = la_qlamch('U')
           ! special case when n=1
           ! treat case of 1x1 matrix for quick return
           if (n == 1) then
              if ((irange == allrng) .or. ((irange == valrng) .and. (d(1) > vl) .and. (d(1) <= vu)) .or. (( &
                        irange == indrng) .and. (il == 1) .and. (iu == 1))) then
                 m = 1
                 w(1) = d(1)
                 ! the computation error of the eigenvalue is zero
                 werr(1) = zero
                 iblock(1) = 1
                 indexw(1) = 1
              end if
              return
           end if
           ! nb is the minimum vector length for vector bisection, or 0
           ! if only scalar is to be done.
           nb = la_ilaenv(1,'QSTEBZ',' ',n,-1,-1,-1)
           if (nb <= 1) nb = 0
           ! find global spectral radius
           gl = d(1)
           gu = d(1)
           do i = 1,n
              gl = min(gl,gers(2*i - 1))
              gu = max(gu,gers(2*i))
           end do
           ! compute global gerschgorin bounds and spectral diameter
           tnorm = max(abs(gl),abs(gu))
           gl = gl - fudge*tnorm*eps*n - fudge*two*pivmin
           gu = gu + fudge*tnorm*eps*n + fudge*two*pivmin
           ! [jan/28/2009] remove the line below since spdiam variable not use
           ! spdiam = gu - gl
           ! input arguments for la_qlaebz:
           ! the relative tolerance.  an interval (a,b] lies within
           ! "relative tolerance" if  b-a < reltol*max(|a|,|b|),
           rtoli = reltol
           ! set the absolute tolerance for interval convergence to zero to force
           ! interval convergence based on relative size of the interval.
           ! this is dangerous because intervals might not converge when reltol is
           ! small. but at least a very small number should be selected so that for
           ! strongly graded matrices, the code can get relatively accurate
           ! eigenvalues.
           atoli = fudge*two*uflow + fudge*two*pivmin
           if (irange == indrng) then
              ! range='i': compute an interval containing eigenvalues
              ! il through iu. the initial interval [gl,gu] from the global
              ! gerschgorin bounds gl and gu is refined by la_qlaebz.
              itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
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
              call la_qlaebz(3,itmax,n,2,2,nb,atoli,rtoli,pivmin,d,e,e2,iwork(5) &
                        ,work(n + 1),work(n + 5),iout,iwork,w,iblock,iinfo)
              if (iinfo /= 0) then
                 info = iinfo
                 return
              end if
              ! on exit, output intervals may not be ordered by ascending negcount
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
              ! on exit, the interval [wl, wlu] contains a value with negcount nwl,
              ! and [wul, wu] contains a value with negcount nwu.
              if (nwl < 0 .or. nwl >= n .or. nwu < 1 .or. nwu > n) then
                 info = 4
                 return
              end if
           elseif (irange == valrng) then
              wl = vl
              wu = vu
           elseif (irange == allrng) then
              wl = gl
              wu = gu
           end if
           ! find eigenvalues -- loop over blocks and recompute nwl and nwu.
           ! nwl accumulates the number of eigenvalues .le. wl,
           ! nwu accumulates the number of eigenvalues .le. wu
           m = 0
           iend = 0
           info = 0
           nwl = 0
           nwu = 0
           loop_70: do jblk = 1,nsplit
              ioff = iend
              ibegin = ioff + 1
              iend = isplit(jblk)
              in = iend - ioff
              if (in == 1) then
                 ! 1x1 block
                 if (wl >= d(ibegin) - pivmin) nwl = nwl + 1
                 if (wu >= d(ibegin) - pivmin) nwu = nwu + 1
                 if (irange == allrng .or. (wl < d(ibegin) - pivmin .and. wu >= d(ibegin) - pivmin)) &
                           then
                    m = m + 1
                    w(m) = d(ibegin)
                    werr(m) = zero
                    ! the gap for a single block doesn't matter for the later
                    ! algorithm and is assigned an arbitrary large value
                    iblock(m) = jblk
                    indexw(m) = 1
                 end if
              ! disabled 2x2 case because of a failure on the following matrix
              ! range = 'i', il = iu = 4
                ! original tridiagonal, d = [
                 ! -0.150102010615740e+00_qp
                 ! -0.849897989384260e+00_qp
                 ! -0.128208148052635e-15_qp
                  ! 0.128257718286320e-15_qp
                ! ];
                ! e = [
                 ! -0.357171383266986e+00_qp
                 ! -0.180411241501588e-15_qp
                 ! -0.175152352710251e-15_qp
                ! ];
               ! else if( in==2 ) then
      ! *           2x2 block
                  ! disc = sqrt( (half*(d(ibegin)-d(iend)))**2 + e(ibegin)**2 )
                  ! tmp1 = half*(d(ibegin)+d(iend))
                  ! l1 = tmp1 - disc
                  ! if( wl>= l1-pivmin )
           ! $         nwl = nwl + 1
                  ! if( wu>= l1-pivmin )
           ! $         nwu = nwu + 1
                  ! if( irange==allrng .or. ( wl<l1-pivmin .and. wu>=
           ! $          l1-pivmin ) ) then
                     ! m = m + 1
                     ! w( m ) = l1
      ! *              the uncertainty of eigenvalues of a 2x2 matrix is very small
                     ! werr( m ) = eps * abs( w( m ) ) * two
                     ! iblock( m ) = jblk
                     ! indexw( m ) = 1
                  ! endif
                  ! l2 = tmp1 + disc
                  ! if( wl>= l2-pivmin )
           ! $         nwl = nwl + 1
                  ! if( wu>= l2-pivmin )
           ! $         nwu = nwu + 1
                  ! if( irange==allrng .or. ( wl<l2-pivmin .and. wu>=
           ! $          l2-pivmin ) ) then
                     ! m = m + 1
                     ! w( m ) = l2
      ! *              the uncertainty of eigenvalues of a 2x2 matrix is very small
                     ! werr( m ) = eps * abs( w( m ) ) * two
                     ! iblock( m ) = jblk
                     ! indexw( m ) = 2
                  ! endif
              else
                 ! general case - block of size in >= 2
                 ! compute local gerschgorin interval and use it as the initial
                 ! interval for la_qlaebz
                 gu = d(ibegin)
                 gl = d(ibegin)
                 tmp1 = zero
                 do j = ibegin,iend
                    gl = min(gl,gers(2*j - 1))
                    gu = max(gu,gers(2*j))
                 end do
                 ! [jan/28/2009]
                 ! change spdiam by tnorm in lines 2 and 3 thereafter
                 ! line 1: remove computation of spdiam (not useful anymore)
                 ! spdiam = gu - gl
                 ! gl = gl - fudge*spdiam*eps*in - fudge*pivmin
                 ! gu = gu + fudge*spdiam*eps*in + fudge*pivmin
                 gl = gl - fudge*tnorm*eps*in - fudge*pivmin
                 gu = gu + fudge*tnorm*eps*in + fudge*pivmin
                 if (irange > 1) then
                    if (gu < wl) then
                       ! the local block contains none of the wanted eigenvalues
                       nwl = nwl + in
                       nwu = nwu + in
                       cycle loop_70
                    end if
                    ! refine search interval if possible, only range (wl,wu] matters
                    gl = max(gl,wl)
                    gu = min(gu,wu)
                    if (gl >= gu) cycle loop_70
                 end if
                 ! find negcount of initial interval boundaries gl and gu
                 work(n + 1) = gl
                 work(n + in + 1) = gu
                 call la_qlaebz(1,0,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                 ibegin),e2(ibegin),idumma,work(n + 1),work(n + 2*in + 1),im,iwork,w(m + 1), &
                            iblock(m + 1),iinfo)
                 if (iinfo /= 0) then
                    info = iinfo
                    return
                 end if
                 nwl = nwl + iwork(1)
                 nwu = nwu + iwork(in + 1)
                 iwoff = m - iwork(1)
                 ! compute eigenvalues
                 itmax = int((log(gu - gl + pivmin) - log(pivmin))/log(two),KIND=ilp) + &
                           2
                 call la_qlaebz(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,d(ibegin),e( &
                  ibegin),e2(ibegin),idumma,work(n + 1),work(n + 2*in + 1),iout,iwork,w(m + &
                            1),iblock(m + 1),iinfo)
                 if (iinfo /= 0) then
                    info = iinfo
                    return
                 end if
                 ! copy eigenvalues into w and iblock
                 ! use -jblk for block number for unconverged eigenvalues.
                 ! loop over the number of output intervals from la_qlaebz
                 do j = 1,iout
                    ! eigenvalue approximation is middle point of interval
                    tmp1 = half*(work(j + n) + work(j + in + n))
                    ! semi length of error interval
                    tmp2 = half*abs(work(j + n) - work(j + in + n))
                    if (j > iout - iinfo) then
                       ! flag non-convergence.
                       ncnvrg = .true.
                       ib = -jblk
                    else
                       ib = jblk
                    end if
                    do je = iwork(j) + 1 + iwoff,iwork(j + in) + iwoff
                       w(je) = tmp1
                       werr(je) = tmp2
                       indexw(je) = je - iwoff
                       iblock(je) = ib
                    end do
                 end do
                 m = m + im
              end if
           end do loop_70
           ! if range='i', then (wl,wu) contains eigenvalues nwl+1,...,nwu
           ! if nwl+1 < il or nwu > iu, discard extra eigenvalues.
           if (irange == indrng) then
              idiscl = il - 1 - nwl
              idiscu = nwu - iu
              if (idiscl > 0) then
                 im = 0
                 do je = 1,m
                    ! remove some of the smallest eigenvalues from the left so that
                    ! at the end idiscl =0. move all eigenvalues up to the left.
                    if (w(je) <= wlu .and. idiscl > 0) then
                       idiscl = idiscl - 1
                    else
                       im = im + 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscu > 0) then
                 ! remove some of the largest eigenvalues from the right so that
                 ! at the end idiscu =0. move all eigenvalues up to the left.
                 im = m + 1
                 do je = m,1,-1
                    if (w(je) >= wul .and. idiscu > 0) then
                       idiscu = idiscu - 1
                    else
                       im = im - 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 jee = 0
                 do je = im,m
                    jee = jee + 1
                    w(jee) = w(je)
                    werr(jee) = werr(je)
                    indexw(jee) = indexw(je)
                    iblock(jee) = iblock(je)
                 end do
                 m = m - im + 1
              end if
              if (idiscl > 0 .or. idiscu > 0) then
                 ! code to deal with effects of bad arithmetic. (if n(w) is
                 ! monotone non-decreasing, this should never happen.)
                 ! some low eigenvalues to be discarded are not in (wl,wlu],
                 ! or high eigenvalues to be discarded are not in (wul,wu]
                 ! so just kill off the smallest idiscl/largest idiscu
                 ! eigenvalues, by marking the corresponding iblock = 0
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
                          if (iblock(je) /= 0 .and. (w(je) >= wkill .or. iw == 0)) then
                             iw = je
                             wkill = w(je)
                          end if
                       end do
                       iblock(iw) = 0
                    end do
                 end if
                 ! now erase all eigenvalues with iblock set to zero
                 im = 0
                 do je = 1,m
                    if (iblock(je) /= 0) then
                       im = im + 1
                       w(im) = w(je)
                       werr(im) = werr(je)
                       indexw(im) = indexw(je)
                       iblock(im) = iblock(je)
                    end if
                 end do
                 m = im
              end if
              if (idiscl < 0 .or. idiscu < 0) then
                 toofew = .true.
              end if
           end if
           if ((irange == allrng .and. m /= n) .or. (irange == indrng .and. m /= iu - il + 1)) then
              toofew = .true.
           end if
           ! if order='b', do nothing the eigenvalues are already sorted by
              ! block.
           ! if order='e', sort the eigenvalues from smallest to largest
           if (la_lsame(order,'E') .and. nsplit > 1) then
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
                    tmp2 = werr(ie)
                    itmp1 = iblock(ie)
                    itmp2 = indexw(ie)
                    w(ie) = w(je)
                    werr(ie) = werr(je)
                    iblock(ie) = iblock(je)
                    indexw(ie) = indexw(je)
                    w(je) = tmp1
                    werr(je) = tmp2
                    iblock(je) = itmp1
                    indexw(je) = itmp2
                 end if
              end do
           end if
           info = 0
           if (ncnvrg) info = info + 1
           if (toofew) info = info + 2
           return
     end subroutine la_qlarrd

     !> Given the initial eigenvalue approximations of T, SLARRJ:
     !> does  bisection to refine the eigenvalues of T,
     !> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
     !> guesses for these eigenvalues are input in W, the corresponding estimate
     !> of the error in these guesses in WERR. During bisection, intervals
     !> [left, right] are maintained by storing their mid-points and
     !> semi-widths in the arrays W and WERR respectively.

     pure subroutine la_slarrj(n,d,e2,ifirst,ilast,rtol,offset,w,werr,work,iwork, &
               pivmin,spdiam,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ifirst,ilast,n,offset
           integer(ilp),intent(out) :: info
           real(sp),intent(in) :: pivmin,rtol,spdiam
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: d(*),e2(*)
           real(sp),intent(inout) :: w(*),werr(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           integer(ilp) :: maxitr
           ! Local Scalars
           integer(ilp) :: cnt,i,i1,i2,ii,iter,j,k,next,nint,olnint,p,prev, &
                     savi1
           real(sp) :: dplus,fac,left,mid,right,s,tmp,width
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           maxitr = int((log(spdiam + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           ! initialize unconverged intervals in [ work(2*i-1), work(2*i) ].
           ! the sturm count, count( work(2*i-1) ) is arranged to be i-1, while
           ! count( work(2*i) ) is stored in iwork( 2*i ). the integer iwork( 2*i-1 )
           ! for an unconverged interval is set to the index of the next unconverged
           ! interval, and is -1 or 0 for a converged interval. thus a linked
           ! list of unconverged intervals is set up.
           i1 = ifirst
           i2 = ilast
           ! the number of unconverged intervals
           nint = 0
           ! the last unconverged interval found
           prev = 0
           loop_75: do i = i1,i2
              k = 2*i
              ii = i - offset
              left = w(ii) - werr(ii)
              mid = w(ii)
              right = w(ii) + werr(ii)
              width = right - mid
              tmp = max(abs(left),abs(right))
              ! the following test prevents the test of converged intervals
              if (width < rtol*tmp) then
                 ! this interval has already converged and does not need refinement.
                 ! (note that the gaps might change through refining the
                  ! eigenvalues, however, they can only get bigger.)
                 ! remove it from the list.
                 iwork(k - 1) = -1
                 ! make sure that i1 always points to the first unconverged interval
                 if ((i == i1) .and. (i < i2)) i1 = i + 1
                 if ((prev >= i1) .and. (i <= i2)) iwork(2*prev - 1) = i + 1
              else
                 ! unconverged interval found
                 prev = i
                 ! make sure that [left,right] contains the desired eigenvalue
                 ! do while( cnt(left)>i-1 )
                 fac = one
                 20 continue
                 cnt = 0
                 s = left
                 dplus = d(1) - s
                 if (dplus < zero) cnt = cnt + 1
                 do j = 2,n
                    dplus = d(j) - s - e2(j - 1)/dplus
                    if (dplus < zero) cnt = cnt + 1
                 end do
                 if (cnt > i - 1) then
                    left = left - werr(ii)*fac
                    fac = two*fac
                    go to 20
                 end if
                 ! do while( cnt(right)<i )
                 fac = one
                 50 continue
                 cnt = 0
                 s = right
                 dplus = d(1) - s
                 if (dplus < zero) cnt = cnt + 1
                 do j = 2,n
                    dplus = d(j) - s - e2(j - 1)/dplus
                    if (dplus < zero) cnt = cnt + 1
                 end do
                 if (cnt < i) then
                    right = right + werr(ii)*fac
                    fac = two*fac
                    go to 50
                 end if
                 nint = nint + 1
                 iwork(k - 1) = i + 1
                 iwork(k) = cnt
              end if
              work(k - 1) = left
              work(k) = right
           end do loop_75
           savi1 = i1
           ! do while( nint>0 ), i.e. there are still unconverged intervals
           ! and while (iter<maxitr)
           iter = 0
           80 continue
           prev = i1 - 1
           i = i1
           olnint = nint
           loop_100: do p = 1,olnint
              k = 2*i
              ii = i - offset
              next = iwork(k - 1)
              left = work(k - 1)
              right = work(k)
              mid = half*(left + right)
              ! semiwidth of interval
              width = right - mid
              tmp = max(abs(left),abs(right))
              if ((width < rtol*tmp) .or. (iter == maxitr)) then
                 ! reduce number of unconverged intervals
                 nint = nint - 1
                 ! mark interval as converged.
                 iwork(k - 1) = 0
                 if (i1 == i) then
                    i1 = next
                 else
                    ! prev holds the last unconverged interval previously examined
                    if (prev >= i1) iwork(2*prev - 1) = next
                 end if
                 i = next
                 cycle loop_100
              end if
              prev = i
              ! perform one bisection step
              cnt = 0
              s = mid
              dplus = d(1) - s
              if (dplus < zero) cnt = cnt + 1
              do j = 2,n
                 dplus = d(j) - s - e2(j - 1)/dplus
                 if (dplus < zero) cnt = cnt + 1
              end do
              if (cnt <= i - 1) then
                 work(k - 1) = mid
              else
                 work(k) = mid
              end if
              i = next
           end do loop_100
           iter = iter + 1
           ! do another loop if there are still unconverged intervals
           ! however, in the last iteration, all intervals are accepted
           ! since this is the best we can do.
           if ((nint > 0) .and. (iter <= maxitr)) go to 80
           ! at this point, all the intervals have converged
           do i = savi1,ilast
              k = 2*i
              ii = i - offset
              ! all intervals marked by '0' have been refined.
              if (iwork(k - 1) == 0) then
                 w(ii) = half*(work(k - 1) + work(k))
                 werr(ii) = work(k) - w(ii)
              end if
           end do
           return
     end subroutine la_slarrj
     !> Given the initial eigenvalue approximations of T, DLARRJ:
     !> does  bisection to refine the eigenvalues of T,
     !> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
     !> guesses for these eigenvalues are input in W, the corresponding estimate
     !> of the error in these guesses in WERR. During bisection, intervals
     !> [left, right] are maintained by storing their mid-points and
     !> semi-widths in the arrays W and WERR respectively.

     pure subroutine la_dlarrj(n,d,e2,ifirst,ilast,rtol,offset,w,werr,work,iwork, &
               pivmin,spdiam,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ifirst,ilast,n,offset
           integer(ilp),intent(out) :: info
           real(dp),intent(in) :: pivmin,rtol,spdiam
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: d(*),e2(*)
           real(dp),intent(inout) :: w(*),werr(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           integer(ilp) :: maxitr
           ! Local Scalars
           integer(ilp) :: cnt,i,i1,i2,ii,iter,j,k,next,nint,olnint,p,prev, &
                     savi1
           real(dp) :: dplus,fac,left,mid,right,s,tmp,width
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           maxitr = int((log(spdiam + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           ! initialize unconverged intervals in [ work(2*i-1), work(2*i) ].
           ! the sturm count, count( work(2*i-1) ) is arranged to be i-1, while
           ! count( work(2*i) ) is stored in iwork( 2*i ). the integer iwork( 2*i-1 )
           ! for an unconverged interval is set to the index of the next unconverged
           ! interval, and is -1 or 0 for a converged interval. thus a linked
           ! list of unconverged intervals is set up.
           i1 = ifirst
           i2 = ilast
           ! the number of unconverged intervals
           nint = 0
           ! the last unconverged interval found
           prev = 0
           loop_75: do i = i1,i2
              k = 2*i
              ii = i - offset
              left = w(ii) - werr(ii)
              mid = w(ii)
              right = w(ii) + werr(ii)
              width = right - mid
              tmp = max(abs(left),abs(right))
              ! the following test prevents the test of converged intervals
              if (width < rtol*tmp) then
                 ! this interval has already converged and does not need refinement.
                 ! (note that the gaps might change through refining the
                  ! eigenvalues, however, they can only get bigger.)
                 ! remove it from the list.
                 iwork(k - 1) = -1
                 ! make sure that i1 always points to the first unconverged interval
                 if ((i == i1) .and. (i < i2)) i1 = i + 1
                 if ((prev >= i1) .and. (i <= i2)) iwork(2*prev - 1) = i + 1
              else
                 ! unconverged interval found
                 prev = i
                 ! make sure that [left,right] contains the desired eigenvalue
                 ! do while( cnt(left)>i-1 )
                 fac = one
                 20 continue
                 cnt = 0
                 s = left
                 dplus = d(1) - s
                 if (dplus < zero) cnt = cnt + 1
                 do j = 2,n
                    dplus = d(j) - s - e2(j - 1)/dplus
                    if (dplus < zero) cnt = cnt + 1
                 end do
                 if (cnt > i - 1) then
                    left = left - werr(ii)*fac
                    fac = two*fac
                    go to 20
                 end if
                 ! do while( cnt(right)<i )
                 fac = one
                 50 continue
                 cnt = 0
                 s = right
                 dplus = d(1) - s
                 if (dplus < zero) cnt = cnt + 1
                 do j = 2,n
                    dplus = d(j) - s - e2(j - 1)/dplus
                    if (dplus < zero) cnt = cnt + 1
                 end do
                 if (cnt < i) then
                    right = right + werr(ii)*fac
                    fac = two*fac
                    go to 50
                 end if
                 nint = nint + 1
                 iwork(k - 1) = i + 1
                 iwork(k) = cnt
              end if
              work(k - 1) = left
              work(k) = right
           end do loop_75
           savi1 = i1
           ! do while( nint>0 ), i.e. there are still unconverged intervals
           ! and while (iter<maxitr)
           iter = 0
           80 continue
           prev = i1 - 1
           i = i1
           olnint = nint
           loop_100: do p = 1,olnint
              k = 2*i
              ii = i - offset
              next = iwork(k - 1)
              left = work(k - 1)
              right = work(k)
              mid = half*(left + right)
              ! semiwidth of interval
              width = right - mid
              tmp = max(abs(left),abs(right))
              if ((width < rtol*tmp) .or. (iter == maxitr)) then
                 ! reduce number of unconverged intervals
                 nint = nint - 1
                 ! mark interval as converged.
                 iwork(k - 1) = 0
                 if (i1 == i) then
                    i1 = next
                 else
                    ! prev holds the last unconverged interval previously examined
                    if (prev >= i1) iwork(2*prev - 1) = next
                 end if
                 i = next
                 cycle loop_100
              end if
              prev = i
              ! perform one bisection step
              cnt = 0
              s = mid
              dplus = d(1) - s
              if (dplus < zero) cnt = cnt + 1
              do j = 2,n
                 dplus = d(j) - s - e2(j - 1)/dplus
                 if (dplus < zero) cnt = cnt + 1
              end do
              if (cnt <= i - 1) then
                 work(k - 1) = mid
              else
                 work(k) = mid
              end if
              i = next
           end do loop_100
           iter = iter + 1
           ! do another loop if there are still unconverged intervals
           ! however, in the last iteration, all intervals are accepted
           ! since this is the best we can do.
           if ((nint > 0) .and. (iter <= maxitr)) go to 80
           ! at this point, all the intervals have converged
           do i = savi1,ilast
              k = 2*i
              ii = i - offset
              ! all intervals marked by '0' have been refined.
              if (iwork(k - 1) == 0) then
                 w(ii) = half*(work(k - 1) + work(k))
                 werr(ii) = work(k) - w(ii)
              end if
           end do
           return
     end subroutine la_dlarrj
     !> Given the initial eigenvalue approximations of T, QLARRJ:
     !> does  bisection to refine the eigenvalues of T,
     !> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
     !> guesses for these eigenvalues are input in W, the corresponding estimate
     !> of the error in these guesses in WERR. During bisection, intervals
     !> [left, right] are maintained by storing their mid-points and
     !> semi-widths in the arrays W and WERR respectively.

     pure subroutine la_qlarrj(n,d,e2,ifirst,ilast,rtol,offset,w,werr,work,iwork, &
               pivmin,spdiam,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ifirst,ilast,n,offset
           integer(ilp),intent(out) :: info
           real(qp),intent(in) :: pivmin,rtol,spdiam
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: d(*),e2(*)
           real(qp),intent(inout) :: w(*),werr(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           integer(ilp) :: maxitr
           ! Local Scalars
           integer(ilp) :: cnt,i,i1,i2,ii,iter,j,k,next,nint,olnint,p,prev, &
                     savi1
           real(qp) :: dplus,fac,left,mid,right,s,tmp,width
           ! Intrinsic Functions
           intrinsic :: abs,max
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           maxitr = int((log(spdiam + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           ! initialize unconverged intervals in [ work(2*i-1), work(2*i) ].
           ! the sturm count, count( work(2*i-1) ) is arranged to be i-1, while
           ! count( work(2*i) ) is stored in iwork( 2*i ). the integer iwork( 2*i-1 )
           ! for an unconverged interval is set to the index of the next unconverged
           ! interval, and is -1 or 0 for a converged interval. thus a linked
           ! list of unconverged intervals is set up.
           i1 = ifirst
           i2 = ilast
           ! the number of unconverged intervals
           nint = 0
           ! the last unconverged interval found
           prev = 0
           loop_75: do i = i1,i2
              k = 2*i
              ii = i - offset
              left = w(ii) - werr(ii)
              mid = w(ii)
              right = w(ii) + werr(ii)
              width = right - mid
              tmp = max(abs(left),abs(right))
              ! the following test prevents the test of converged intervals
              if (width < rtol*tmp) then
                 ! this interval has already converged and does not need refinement.
                 ! (note that the gaps might change through refining the
                  ! eigenvalues, however, they can only get bigger.)
                 ! remove it from the list.
                 iwork(k - 1) = -1
                 ! make sure that i1 always points to the first unconverged interval
                 if ((i == i1) .and. (i < i2)) i1 = i + 1
                 if ((prev >= i1) .and. (i <= i2)) iwork(2*prev - 1) = i + 1
              else
                 ! unconverged interval found
                 prev = i
                 ! make sure that [left,right] contains the desired eigenvalue
                 ! do while( cnt(left)>i-1 )
                 fac = one
                 20 continue
                 cnt = 0
                 s = left
                 dplus = d(1) - s
                 if (dplus < zero) cnt = cnt + 1
                 do j = 2,n
                    dplus = d(j) - s - e2(j - 1)/dplus
                    if (dplus < zero) cnt = cnt + 1
                 end do
                 if (cnt > i - 1) then
                    left = left - werr(ii)*fac
                    fac = two*fac
                    go to 20
                 end if
                 ! do while( cnt(right)<i )
                 fac = one
                 50 continue
                 cnt = 0
                 s = right
                 dplus = d(1) - s
                 if (dplus < zero) cnt = cnt + 1
                 do j = 2,n
                    dplus = d(j) - s - e2(j - 1)/dplus
                    if (dplus < zero) cnt = cnt + 1
                 end do
                 if (cnt < i) then
                    right = right + werr(ii)*fac
                    fac = two*fac
                    go to 50
                 end if
                 nint = nint + 1
                 iwork(k - 1) = i + 1
                 iwork(k) = cnt
              end if
              work(k - 1) = left
              work(k) = right
           end do loop_75
           savi1 = i1
           ! do while( nint>0 ), i.e. there are still unconverged intervals
           ! and while (iter<maxitr)
           iter = 0
           80 continue
           prev = i1 - 1
           i = i1
           olnint = nint
           loop_100: do p = 1,olnint
              k = 2*i
              ii = i - offset
              next = iwork(k - 1)
              left = work(k - 1)
              right = work(k)
              mid = half*(left + right)
              ! semiwidth of interval
              width = right - mid
              tmp = max(abs(left),abs(right))
              if ((width < rtol*tmp) .or. (iter == maxitr)) then
                 ! reduce number of unconverged intervals
                 nint = nint - 1
                 ! mark interval as converged.
                 iwork(k - 1) = 0
                 if (i1 == i) then
                    i1 = next
                 else
                    ! prev holds the last unconverged interval previously examined
                    if (prev >= i1) iwork(2*prev - 1) = next
                 end if
                 i = next
                 cycle loop_100
              end if
              prev = i
              ! perform one bisection step
              cnt = 0
              s = mid
              dplus = d(1) - s
              if (dplus < zero) cnt = cnt + 1
              do j = 2,n
                 dplus = d(j) - s - e2(j - 1)/dplus
                 if (dplus < zero) cnt = cnt + 1
              end do
              if (cnt <= i - 1) then
                 work(k - 1) = mid
              else
                 work(k) = mid
              end if
              i = next
           end do loop_100
           iter = iter + 1
           ! do another loop if there are still unconverged intervals
           ! however, in the last iteration, all intervals are accepted
           ! since this is the best we can do.
           if ((nint > 0) .and. (iter <= maxitr)) go to 80
           ! at this point, all the intervals have converged
           do i = savi1,ilast
              k = 2*i
              ii = i - offset
              ! all intervals marked by '0' have been refined.
              if (iwork(k - 1) == 0) then
                 w(ii) = half*(work(k - 1) + work(k))
                 werr(ii) = work(k) - w(ii)
              end if
           end do
           return
     end subroutine la_qlarrj

     !> SLARRK: computes one eigenvalue of a symmetric tridiagonal
     !> matrix T to suitable accuracy. This is an auxiliary code to be
     !> called from SSTEMR.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_slarrk(n,iw,gl,gu,d,e2,pivmin,reltol,w,werr,info)
        use la_constants_sp,only:zero,half,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: iw,n
           real(sp),intent(in) :: pivmin,reltol,gl,gu
           real(sp),intent(out) :: w,werr
           ! Array Arguments
           real(sp),intent(in) :: d(*),e2(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: fudge = two

           ! Local Scalars
           integer(ilp) :: i,it,itmax,negcnt
           real(sp) :: atoli,eps,left,mid,right,rtoli,tmp1,tmp2,tnorm
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              info = 0
              return
           end if
           ! get machine constants
           eps = la_slamch('P')
           tnorm = max(abs(gl),abs(gu))
           rtoli = reltol
           atoli = fudge*two*pivmin
           itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           info = -1
           left = gl - fudge*tnorm*eps*n - fudge*two*pivmin
           right = gu + fudge*tnorm*eps*n + fudge*two*pivmin
           it = 0
           10 continue
           ! check if interval converged or maximum number of iterations reached
           tmp1 = abs(right - left)
           tmp2 = max(abs(right),abs(left))
           if (tmp1 < max(atoli,pivmin,rtoli*tmp2)) then
              info = 0
              goto 30
           end if
           if (it > itmax) goto 30
           ! count number of negative pivots for mid-point
           it = it + 1
           mid = half*(left + right)
           negcnt = 0
           tmp1 = d(1) - mid
           if (abs(tmp1) < pivmin) tmp1 = -pivmin
           if (tmp1 <= zero) negcnt = negcnt + 1
           do i = 2,n
              tmp1 = d(i) - e2(i - 1)/tmp1 - mid
              if (abs(tmp1) < pivmin) tmp1 = -pivmin
              if (tmp1 <= zero) negcnt = negcnt + 1
           end do
           if (negcnt >= iw) then
              right = mid
           else
              left = mid
           end if
           goto 10
           30 continue
           ! converged or maximum number of iterations reached
           w = half*(left + right)
           werr = half*abs(right - left)
           return
     end subroutine la_slarrk
     !> DLARRK: computes one eigenvalue of a symmetric tridiagonal
     !> matrix T to suitable accuracy. This is an auxiliary code to be
     !> called from DSTEMR.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_dlarrk(n,iw,gl,gu,d,e2,pivmin,reltol,w,werr,info)
        use la_constants_dp,only:zero,half,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: iw,n
           real(dp),intent(in) :: pivmin,reltol,gl,gu
           real(dp),intent(out) :: w,werr
           ! Array Arguments
           real(dp),intent(in) :: d(*),e2(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: fudge = two

           ! Local Scalars
           integer(ilp) :: i,it,itmax,negcnt
           real(dp) :: atoli,eps,left,mid,right,rtoli,tmp1,tmp2,tnorm
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              info = 0
              return
           end if
           ! get machine constants
           eps = la_dlamch('P')
           tnorm = max(abs(gl),abs(gu))
           rtoli = reltol
           atoli = fudge*two*pivmin
           itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           info = -1
           left = gl - fudge*tnorm*eps*n - fudge*two*pivmin
           right = gu + fudge*tnorm*eps*n + fudge*two*pivmin
           it = 0
           10 continue
           ! check if interval converged or maximum number of iterations reached
           tmp1 = abs(right - left)
           tmp2 = max(abs(right),abs(left))
           if (tmp1 < max(atoli,pivmin,rtoli*tmp2)) then
              info = 0
              goto 30
           end if
           if (it > itmax) goto 30
           ! count number of negative pivots for mid-point
           it = it + 1
           mid = half*(left + right)
           negcnt = 0
           tmp1 = d(1) - mid
           if (abs(tmp1) < pivmin) tmp1 = -pivmin
           if (tmp1 <= zero) negcnt = negcnt + 1
           do i = 2,n
              tmp1 = d(i) - e2(i - 1)/tmp1 - mid
              if (abs(tmp1) < pivmin) tmp1 = -pivmin
              if (tmp1 <= zero) negcnt = negcnt + 1
           end do
           if (negcnt >= iw) then
              right = mid
           else
              left = mid
           end if
           goto 10
           30 continue
           ! converged or maximum number of iterations reached
           w = half*(left + right)
           werr = half*abs(right - left)
           return
     end subroutine la_dlarrk
     !> QLARRK: computes one eigenvalue of a symmetric tridiagonal
     !> matrix T to suitable accuracy. This is an auxiliary code to be
     !> called from QSTEMR.
     !> To avoid overflow, the matrix must be scaled so that its
     !> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
     !> accuracy, it should not be much smaller than that.
     !> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
     !> Matrix", Report CS41, Computer Science Dept., Stanford
     !> University, July 21, 1966.

     pure subroutine la_qlarrk(n,iw,gl,gu,d,e2,pivmin,reltol,w,werr,info)
        use la_constants_qp,only:zero,half,two
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: iw,n
           real(qp),intent(in) :: pivmin,reltol,gl,gu
           real(qp),intent(out) :: w,werr
           ! Array Arguments
           real(qp),intent(in) :: d(*),e2(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: fudge = two

           ! Local Scalars
           integer(ilp) :: i,it,itmax,negcnt
           real(qp) :: atoli,eps,left,mid,right,rtoli,tmp1,tmp2,tnorm
           ! Intrinsic Functions
           intrinsic :: abs,int,log,max
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              info = 0
              return
           end if
           ! get machine constants
           eps = la_qlamch('P')
           tnorm = max(abs(gl),abs(gu))
           rtoli = reltol
           atoli = fudge*two*pivmin
           itmax = int((log(tnorm + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           info = -1
           left = gl - fudge*tnorm*eps*n - fudge*two*pivmin
           right = gu + fudge*tnorm*eps*n + fudge*two*pivmin
           it = 0
           10 continue
           ! check if interval converged or maximum number of iterations reached
           tmp1 = abs(right - left)
           tmp2 = max(abs(right),abs(left))
           if (tmp1 < max(atoli,pivmin,rtoli*tmp2)) then
              info = 0
              goto 30
           end if
           if (it > itmax) goto 30
           ! count number of negative pivots for mid-point
           it = it + 1
           mid = half*(left + right)
           negcnt = 0
           tmp1 = d(1) - mid
           if (abs(tmp1) < pivmin) tmp1 = -pivmin
           if (tmp1 <= zero) negcnt = negcnt + 1
           do i = 2,n
              tmp1 = d(i) - e2(i - 1)/tmp1 - mid
              if (abs(tmp1) < pivmin) tmp1 = -pivmin
              if (tmp1 <= zero) negcnt = negcnt + 1
           end do
           if (negcnt >= iw) then
              right = mid
           else
              left = mid
           end if
           goto 10
           30 continue
           ! converged or maximum number of iterations reached
           w = half*(left + right)
           werr = half*abs(right - left)
           return
     end subroutine la_qlarrk

     !> Perform tests to decide whether the symmetric tridiagonal matrix T
     !> warrants expensive computations which guarantee high relative accuracy
     !> in the eigenvalues.

     pure subroutine la_slarrr(n,d,e,info)
        use la_constants_sp,only:zero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(sp),intent(in) :: d(*)
           real(sp),intent(inout) :: e(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: relcond = 0.999_sp

           ! Local Scalars
           integer(ilp) :: i
           logical(lk) :: yesrel
           real(sp) :: eps,safmin,smlnum,rmin,tmp,tmp2,offdig,offdig2
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              info = 0
              return
           end if
           ! as a default, do not go for relative-accuracy preserving computations.
           info = 1
           safmin = la_slamch('SAFE MINIMUM')
           eps = la_slamch('PRECISION')
           smlnum = safmin/eps
           rmin = sqrt(smlnum)
           ! tests for relative accuracy
           ! test for scaled diagonal dominance
           ! scale the diagonal entries to one and check whether the sum of the
           ! off-diagonals is less than one
           ! the sdd relative error bounds have a 1/(1- 2*x) factor in them,
           ! x = max(offdig + offdig2), so when x is close to 1/2, no relative
           ! accuracy is promised.  in the notation of the code fragment below,
           ! 1/(1 - (offdig + offdig2)) is the condition number.
           ! we don't think it is worth going into "sdd mode" unless the relative
           ! condition number is reasonable, not 1/macheps.
           ! the threshold should be compatible with other thresholds used in the
           ! code. we set  offdig + offdig2 <= .999_sp =: relcond, it corresponds
           ! to losing at most 3 decimal digits: 1 / (1 - (offdig + offdig2)) <= 1000
           ! instead of the current offdig + offdig2 < 1
           yesrel = .true.
           offdig = zero
           tmp = sqrt(abs(d(1)))
           if (tmp < rmin) yesrel = .false.
           if (.not. yesrel) goto 11
           do i = 2,n
              tmp2 = sqrt(abs(d(i)))
              if (tmp2 < rmin) yesrel = .false.
              if (.not. yesrel) goto 11
              offdig2 = abs(e(i - 1))/(tmp*tmp2)
              if (offdig + offdig2 >= relcond) yesrel = .false.
              if (.not. yesrel) goto 11
              tmp = tmp2
              offdig = offdig2
           end do
           11 continue
           if (yesrel) then
              info = 0
              return
           else
           end if
           ! *** more to be implemented ***
           ! test if the lower bidiagonal matrix l from t = l d l^t
           ! (zero shift facto) is well conditioned
           ! test if the upper bidiagonal matrix u from t = u d u^t
           ! (zero shift facto) is well conditioned.
           ! in this case, the matrix needs to be flipped and, at the end
           ! of the eigenvector computation, the flip needs to be applied
           ! to the computed eigenvectors (and the support)
           return
     end subroutine la_slarrr
     !> Perform tests to decide whether the symmetric tridiagonal matrix T
     !> warrants expensive computations which guarantee high relative accuracy
     !> in the eigenvalues.

     pure subroutine la_dlarrr(n,d,e,info)
        use la_constants_dp,only:zero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(dp),intent(in) :: d(*)
           real(dp),intent(inout) :: e(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: relcond = 0.999_dp

           ! Local Scalars
           integer(ilp) :: i
           logical(lk) :: yesrel
           real(dp) :: eps,safmin,smlnum,rmin,tmp,tmp2,offdig,offdig2
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              info = 0
              return
           end if
           ! as a default, do not go for relative-accuracy preserving computations.
           info = 1
           safmin = la_dlamch('SAFE MINIMUM')
           eps = la_dlamch('PRECISION')
           smlnum = safmin/eps
           rmin = sqrt(smlnum)
           ! tests for relative accuracy
           ! test for scaled diagonal dominance
           ! scale the diagonal entries to one and check whether the sum of the
           ! off-diagonals is less than one
           ! the sdd relative error bounds have a 1/(1- 2*x) factor in them,
           ! x = max(offdig + offdig2), so when x is close to 1/2, no relative
           ! accuracy is promised.  in the notation of the code fragment below,
           ! 1/(1 - (offdig + offdig2)) is the condition number.
           ! we don't think it is worth going into "sdd mode" unless the relative
           ! condition number is reasonable, not 1/macheps.
           ! the threshold should be compatible with other thresholds used in the
           ! code. we set  offdig + offdig2 <= .999_dp =: relcond, it corresponds
           ! to losing at most 3 decimal digits: 1 / (1 - (offdig + offdig2)) <= 1000
           ! instead of the current offdig + offdig2 < 1
           yesrel = .true.
           offdig = zero
           tmp = sqrt(abs(d(1)))
           if (tmp < rmin) yesrel = .false.
           if (.not. yesrel) goto 11
           do i = 2,n
              tmp2 = sqrt(abs(d(i)))
              if (tmp2 < rmin) yesrel = .false.
              if (.not. yesrel) goto 11
              offdig2 = abs(e(i - 1))/(tmp*tmp2)
              if (offdig + offdig2 >= relcond) yesrel = .false.
              if (.not. yesrel) goto 11
              tmp = tmp2
              offdig = offdig2
           end do
           11 continue
           if (yesrel) then
              info = 0
              return
           else
           end if
           ! *** more to be implemented ***
           ! test if the lower bidiagonal matrix l from t = l d l^t
           ! (zero shift facto) is well conditioned
           ! test if the upper bidiagonal matrix u from t = u d u^t
           ! (zero shift facto) is well conditioned.
           ! in this case, the matrix needs to be flipped and, at the end
           ! of the eigenvector computation, the flip needs to be applied
           ! to the computed eigenvectors (and the support)
           return
     end subroutine la_dlarrr
     !> Perform tests to decide whether the symmetric tridiagonal matrix T
     !> warrants expensive computations which guarantee high relative accuracy
     !> in the eigenvalues.

     pure subroutine la_qlarrr(n,d,e,info)
        use la_constants_qp,only:zero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n
           integer(ilp),intent(out) :: info
           ! Array Arguments
           real(qp),intent(in) :: d(*)
           real(qp),intent(inout) :: e(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: relcond = 0.999_qp

           ! Local Scalars
           integer(ilp) :: i
           logical(lk) :: yesrel
           real(qp) :: eps,safmin,smlnum,rmin,tmp,tmp2,offdig,offdig2
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) then
              info = 0
              return
           end if
           ! as a default, do not go for relative-accuracy preserving computations.
           info = 1
           safmin = la_qlamch('SAFE MINIMUM')
           eps = la_qlamch('PRECISION')
           smlnum = safmin/eps
           rmin = sqrt(smlnum)
           ! tests for relative accuracy
           ! test for scaled diagonal dominance
           ! scale the diagonal entries to one and check whether the sum of the
           ! off-diagonals is less than one
           ! the sdd relative error bounds have a 1/(1- 2*x) factor in them,
           ! x = max(offdig + offdig2), so when x is close to 1/2, no relative
           ! accuracy is promised.  in the notation of the code fragment below,
           ! 1/(1 - (offdig + offdig2)) is the condition number.
           ! we don't think it is worth going into "sdd mode" unless the relative
           ! condition number is reasonable, not 1/macheps.
           ! the threshold should be compatible with other thresholds used in the
           ! code. we set  offdig + offdig2 <= .999_qp =: relcond, it corresponds
           ! to losing at most 3 decimal digits: 1 / (1 - (offdig + offdig2)) <= 1000
           ! instead of the current offdig + offdig2 < 1
           yesrel = .true.
           offdig = zero
           tmp = sqrt(abs(d(1)))
           if (tmp < rmin) yesrel = .false.
           if (.not. yesrel) goto 11
           do i = 2,n
              tmp2 = sqrt(abs(d(i)))
              if (tmp2 < rmin) yesrel = .false.
              if (.not. yesrel) goto 11
              offdig2 = abs(e(i - 1))/(tmp*tmp2)
              if (offdig + offdig2 >= relcond) yesrel = .false.
              if (.not. yesrel) goto 11
              tmp = tmp2
              offdig = offdig2
           end do
           11 continue
           if (yesrel) then
              info = 0
              return
           else
           end if
           ! *** more to be implemented ***
           ! test if the lower bidiagonal matrix l from t = l d l^t
           ! (zero shift facto) is well conditioned
           ! test if the upper bidiagonal matrix u from t = u d u^t
           ! (zero shift facto) is well conditioned.
           ! in this case, the matrix needs to be flipped and, at the end
           ! of the eigenvector computation, the flip needs to be applied
           ! to the computed eigenvectors (and the support)
           return
     end subroutine la_qlarrr

     !> SLANEG: computes the Sturm count, the number of negative pivots
     !> encountered while factoring tridiagonal T - sigma I = L D L^T.
     !> This implementation works directly on the factors without forming
     !> the tridiagonal matrix T.  The Sturm count is also the number of
     !> eigenvalues of T less than sigma.
     !> This routine is called from SLARRB.
     !> The current routine does not use the PIVMIN parameter but rather
     !> requires IEEE-754 propagation of Infinities and NaNs.  This
     !> routine also has no input range restrictions but does require
     !> default exception handling such that x/0 produces Inf when x is
     !> non-zero, and Inf/Inf produces NaN.  For more information, see:
     !> Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in
     !> Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on
     !> Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624
     !> (Tech report version in LAWN 172 with the same title.)

     pure integer(ilp) function la_slaneg(n,d,lld,sigma,pivmin,r)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,r
           real(sp),intent(in) :: pivmin,sigma
           ! Array Arguments
           real(sp),intent(in) :: d(*),lld(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blklen = 128

           ! some architectures propagate infinities and nans very slowly, so
           ! the code computes counts in blklen chunks.  then a nan can
           ! propagate at most blklen columns before being detected.  this is
           ! not a general tuning parameter; it needs only to be just large
           ! enough that the overhead is tiny in common cases.

           ! Local Scalars
           integer(ilp) :: bj,j,neg1,neg2,negcnt
           real(sp) :: bsav,dminus,dplus,gamma,p,t,tmp
           logical(lk) :: sawnan
           ! Intrinsic Functions
           intrinsic :: min,max
           ! Executable Statements
           negcnt = 0
           ! i) upper part: l d l^t - sigma i = l+ d+ l+^t
           t = -sigma
           loop_210: do bj = 1,r - 1,blklen
              neg1 = 0
              bsav = t
              do j = bj,min(bj + blklen - 1,r - 1)
                 dplus = d(j) + t
                 if (dplus < zero) neg1 = neg1 + 1
                 tmp = t/dplus
                 t = tmp*lld(j) - sigma
              end do
              sawnan = la_sisnan(t)
           ! run a slower version of the above loop if a nan is detected.
           ! a nan should occur only with a zero pivot after an infinite
           ! pivot.  in that case, substituting 1 for t/dplus is the
           ! correct limit.
              if (sawnan) then
                 neg1 = 0
                 t = bsav
                 do j = bj,min(bj + blklen - 1,r - 1)
                    dplus = d(j) + t
                    if (dplus < zero) neg1 = neg1 + 1
                    tmp = t/dplus
                    if (la_sisnan(tmp)) tmp = one
                    t = tmp*lld(j) - sigma
                 end do
              end if
              negcnt = negcnt + neg1
           end do loop_210
           ! ii) lower part: l d l^t - sigma i = u- d- u-^t
           p = d(n) - sigma
           do bj = n - 1,r,-blklen
              neg2 = 0
              bsav = p
              do j = bj,max(bj - blklen + 1,r),-1
                 dminus = lld(j) + p
                 if (dminus < zero) neg2 = neg2 + 1
                 tmp = p/dminus
                 p = tmp*d(j) - sigma
              end do
              sawnan = la_sisnan(p)
           ! as above, run a slower version that substitutes 1 for inf/inf.
              if (sawnan) then
                 neg2 = 0
                 p = bsav
                 do j = bj,max(bj - blklen + 1,r),-1
                    dminus = lld(j) + p
                    if (dminus < zero) neg2 = neg2 + 1
                    tmp = p/dminus
                    if (la_sisnan(tmp)) tmp = one
                    p = tmp*d(j) - sigma
                 end do
              end if
              negcnt = negcnt + neg2
           end do
           ! iii) twist index
             ! t was shifted by sigma initially.
           gamma = (t + sigma) + p
           if (gamma < zero) negcnt = negcnt + 1
           la_slaneg = negcnt
     end function la_slaneg
     !> DLANEG: computes the Sturm count, the number of negative pivots
     !> encountered while factoring tridiagonal T - sigma I = L D L^T.
     !> This implementation works directly on the factors without forming
     !> the tridiagonal matrix T.  The Sturm count is also the number of
     !> eigenvalues of T less than sigma.
     !> This routine is called from DLARRB.
     !> The current routine does not use the PIVMIN parameter but rather
     !> requires IEEE-754 propagation of Infinities and NaNs.  This
     !> routine also has no input range restrictions but does require
     !> default exception handling such that x/0 produces Inf when x is
     !> non-zero, and Inf/Inf produces NaN.  For more information, see:
     !> Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in
     !> Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on
     !> Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624
     !> (Tech report version in LAWN 172 with the same title.)

     pure integer(ilp) function la_dlaneg(n,d,lld,sigma,pivmin,r)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,r
           real(dp),intent(in) :: pivmin,sigma
           ! Array Arguments
           real(dp),intent(in) :: d(*),lld(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blklen = 128

           ! some architectures propagate infinities and nans very slowly, so
           ! the code computes counts in blklen chunks.  then a nan can
           ! propagate at most blklen columns before being detected.  this is
           ! not a general tuning parameter; it needs only to be just large
           ! enough that the overhead is tiny in common cases.

           ! Local Scalars
           integer(ilp) :: bj,j,neg1,neg2,negcnt
           real(dp) :: bsav,dminus,dplus,gamma,p,t,tmp
           logical(lk) :: sawnan
           ! Intrinsic Functions
           intrinsic :: min,max
           ! Executable Statements
           negcnt = 0
           ! i) upper part: l d l^t - sigma i = l+ d+ l+^t
           t = -sigma
           loop_210: do bj = 1,r - 1,blklen
              neg1 = 0
              bsav = t
              do j = bj,min(bj + blklen - 1,r - 1)
                 dplus = d(j) + t
                 if (dplus < zero) neg1 = neg1 + 1
                 tmp = t/dplus
                 t = tmp*lld(j) - sigma
              end do
              sawnan = la_disnan(t)
           ! run a slower version of the above loop if a nan is detected.
           ! a nan should occur only with a zero pivot after an infinite
           ! pivot.  in that case, substituting 1 for t/dplus is the
           ! correct limit.
              if (sawnan) then
                 neg1 = 0
                 t = bsav
                 do j = bj,min(bj + blklen - 1,r - 1)
                    dplus = d(j) + t
                    if (dplus < zero) neg1 = neg1 + 1
                    tmp = t/dplus
                    if (la_disnan(tmp)) tmp = one
                    t = tmp*lld(j) - sigma
                 end do
              end if
              negcnt = negcnt + neg1
           end do loop_210
           ! ii) lower part: l d l^t - sigma i = u- d- u-^t
           p = d(n) - sigma
           do bj = n - 1,r,-blklen
              neg2 = 0
              bsav = p
              do j = bj,max(bj - blklen + 1,r),-1
                 dminus = lld(j) + p
                 if (dminus < zero) neg2 = neg2 + 1
                 tmp = p/dminus
                 p = tmp*d(j) - sigma
              end do
              sawnan = la_disnan(p)
           ! as above, run a slower version that substitutes 1 for inf/inf.
              if (sawnan) then
                 neg2 = 0
                 p = bsav
                 do j = bj,max(bj - blklen + 1,r),-1
                    dminus = lld(j) + p
                    if (dminus < zero) neg2 = neg2 + 1
                    tmp = p/dminus
                    if (la_disnan(tmp)) tmp = one
                    p = tmp*d(j) - sigma
                 end do
              end if
              negcnt = negcnt + neg2
           end do
           ! iii) twist index
             ! t was shifted by sigma initially.
           gamma = (t + sigma) + p
           if (gamma < zero) negcnt = negcnt + 1
           la_dlaneg = negcnt
     end function la_dlaneg
     !> QLANEG: computes the Sturm count, the number of negative pivots
     !> encountered while factoring tridiagonal T - sigma I = L D L^T.
     !> This implementation works directly on the factors without forming
     !> the tridiagonal matrix T.  The Sturm count is also the number of
     !> eigenvalues of T less than sigma.
     !> This routine is called from QLARRB.
     !> The current routine does not use the PIVMIN parameter but rather
     !> requires IEEE-754 propagation of Infinities and NaNs.  This
     !> routine also has no input range restrictions but does require
     !> default exception handling such that x/0 produces Inf when x is
     !> non-zero, and Inf/Inf produces NaN.  For more information, see:
     !> Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in
     !> Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on
     !> Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624
     !> (Tech report version in LAWN 172 with the same title.)

     pure integer(ilp) function la_qlaneg(n,d,lld,sigma,pivmin,r)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: n,r
           real(qp),intent(in) :: pivmin,sigma
           ! Array Arguments
           real(qp),intent(in) :: d(*),lld(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: blklen = 128

           ! some architectures propagate infinities and nans very slowly, so
           ! the code computes counts in blklen chunks.  then a nan can
           ! propagate at most blklen columns before being detected.  this is
           ! not a general tuning parameter; it needs only to be just large
           ! enough that the overhead is tiny in common cases.

           ! Local Scalars
           integer(ilp) :: bj,j,neg1,neg2,negcnt
           real(qp) :: bsav,dminus,dplus,gamma,p,t,tmp
           logical(lk) :: sawnan
           ! Intrinsic Functions
           intrinsic :: min,max
           ! Executable Statements
           negcnt = 0
           ! i) upper part: l d l^t - sigma i = l+ d+ l+^t
           t = -sigma
           loop_210: do bj = 1,r - 1,blklen
              neg1 = 0
              bsav = t
              do j = bj,min(bj + blklen - 1,r - 1)
                 dplus = d(j) + t
                 if (dplus < zero) neg1 = neg1 + 1
                 tmp = t/dplus
                 t = tmp*lld(j) - sigma
              end do
              sawnan = la_qisnan(t)
           ! run a slower version of the above loop if a nan is detected.
           ! a nan should occur only with a zero pivot after an infinite
           ! pivot.  in that case, substituting 1 for t/dplus is the
           ! correct limit.
              if (sawnan) then
                 neg1 = 0
                 t = bsav
                 do j = bj,min(bj + blklen - 1,r - 1)
                    dplus = d(j) + t
                    if (dplus < zero) neg1 = neg1 + 1
                    tmp = t/dplus
                    if (la_qisnan(tmp)) tmp = one
                    t = tmp*lld(j) - sigma
                 end do
              end if
              negcnt = negcnt + neg1
           end do loop_210
           ! ii) lower part: l d l^t - sigma i = u- d- u-^t
           p = d(n) - sigma
           do bj = n - 1,r,-blklen
              neg2 = 0
              bsav = p
              do j = bj,max(bj - blklen + 1,r),-1
                 dminus = lld(j) + p
                 if (dminus < zero) neg2 = neg2 + 1
                 tmp = p/dminus
                 p = tmp*d(j) - sigma
              end do
              sawnan = la_qisnan(p)
           ! as above, run a slower version that substitutes 1 for inf/inf.
              if (sawnan) then
                 neg2 = 0
                 p = bsav
                 do j = bj,max(bj - blklen + 1,r),-1
                    dminus = lld(j) + p
                    if (dminus < zero) neg2 = neg2 + 1
                    tmp = p/dminus
                    if (la_qisnan(tmp)) tmp = one
                    p = tmp*d(j) - sigma
                 end do
              end if
              negcnt = negcnt + neg2
           end do
           ! iii) twist index
             ! t was shifted by sigma initially.
           gamma = (t + sigma) + p
           if (gamma < zero) negcnt = negcnt + 1
           la_qlaneg = negcnt
     end function la_qlaneg

     !> SLAR1V: computes the (scaled) r-th column of the inverse of
     !> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
     !> L D L**T - sigma I. When sigma is close to an eigenvalue, the
     !> computed vector is an accurate eigenvector. Usually, r corresponds
     !> to the index where the eigenvector is largest in magnitude.
     !> The following steps accomplish this computation :
     !> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
     !> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
     !> (c) Computation of the diagonal elements of the inverse of
     !> L D L**T - sigma I by combining the above transforms, and choosing
     !> r as the index where the diagonal of the inverse is (one of the)
     !> largest in magnitude.
     !> (d) Computation of the (scaled) r-th column of the inverse using the
     !> twisted factorization obtained by combining the top part of the
     !> the stationary and the bottom part of the progressive transform.

     pure subroutine la_slar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc, &
               negcnt,ztz,mingma,r,isuppz,nrminv,resid,rqcorr,work)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantnc
           integer(ilp),intent(in) :: b1,bn,n
           integer(ilp),intent(out) :: negcnt
           integer(ilp),intent(inout) :: r
           real(sp),intent(in) :: gaptol,lambda,pivmin
           real(sp),intent(out) :: mingma,nrminv,resid,rqcorr,ztz
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*)
           real(sp),intent(in) :: d(*),l(*),ld(*),lld(*)
           real(sp),intent(out) :: work(*)
           real(sp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: sawnan1,sawnan2
           integer(ilp) :: i,indlpl,indp,inds,indumn,neg1,neg2,r1,r2
           real(sp) :: dminus,dplus,eps,s,tmp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           eps = la_slamch('PRECISION')
           if (r == 0) then
              r1 = b1
              r2 = bn
           else
              r1 = r
              r2 = r
           end if
           ! storage for lplus
           indlpl = 0
           ! storage for uminus
           indumn = n
           inds = 2*n + 1
           indp = 3*n + 1
           if (b1 == 1) then
              work(inds) = zero
           else
              work(inds + b1 - 1) = lld(b1 - 1)
           end if
           ! compute the stationary transform (using the differential form)
           ! until the index r2.
           sawnan1 = .false.
           neg1 = 0
           s = work(inds + b1 - 1) - lambda
           do i = b1,r1 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              if (dplus < zero) neg1 = neg1 + 1
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_sisnan(s)
           if (sawnan1) goto 60
           do i = r1,r2 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_sisnan(s)
           60 continue
           if (sawnan1) then
              ! runs a slower version of the above loop if a nan is detected
              neg1 = 0
              s = work(inds + b1 - 1) - lambda
              do i = b1,r1 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 if (dplus < zero) neg1 = neg1 + 1
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
              do i = r1,r2 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
           end if
           ! compute the progressive transform (using the differential form)
           ! until the index r1
           sawnan2 = .false.
           neg2 = 0
           work(indp + bn - 1) = d(bn) - lambda
           do i = bn - 1,r1,-1
              dminus = lld(i) + work(indp + i)
              tmp = d(i)/dminus
              if (dminus < zero) neg2 = neg2 + 1
              work(indumn + i) = l(i)*tmp
              work(indp + i - 1) = work(indp + i)*tmp - lambda
           end do
           tmp = work(indp + r1 - 1)
           sawnan2 = la_sisnan(tmp)
           if (sawnan2) then
              ! runs a slower version of the above loop if a nan is detected
              neg2 = 0
              do i = bn - 1,r1,-1
                 dminus = lld(i) + work(indp + i)
                 if (abs(dminus) < pivmin) dminus = -pivmin
                 tmp = d(i)/dminus
                 if (dminus < zero) neg2 = neg2 + 1
                 work(indumn + i) = l(i)*tmp
                 work(indp + i - 1) = work(indp + i)*tmp - lambda
                 if (tmp == zero) work(indp + i - 1) = d(i) - lambda
              end do
           end if
           ! find the index (from r1 to r2) of the largest (in magnitude)
           ! diagonal element of the inverse
           mingma = work(inds + r1 - 1) + work(indp + r1 - 1)
           if (mingma < zero) neg1 = neg1 + 1
           if (wantnc) then
              negcnt = neg1 + neg2
           else
              negcnt = -1
           end if
           if (abs(mingma) == zero) mingma = eps*work(inds + r1 - 1)
           r = r1
           do i = r1,r2 - 1
              tmp = work(inds + i) + work(indp + i)
              if (tmp == zero) tmp = eps*work(inds + i)
              if (abs(tmp) <= abs(mingma)) then
                 mingma = tmp
                 r = i + 1
              end if
           end do
           ! compute the fp vector: solve n^t v = e_r
           isuppz(1) = b1
           isuppz(2) = bn
           z(r) = one
           ztz = one
           ! compute the fp vector upwards from r
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r - 1,b1,-1
                 z(i) = -(work(indlpl + i)*z(i + 1))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    goto 220
                 end if
                 ztz = ztz + z(i)*z(i)
              end do
              220 continue
           else
              ! run slower loop if nan occurred.
              do i = r - 1,b1,-1
                 if (z(i + 1) == zero) then
                    z(i) = -(ld(i + 1)/ld(i))*z(i + 2)
                 else
                    z(i) = -(work(indlpl + i)*z(i + 1))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    go to 240
                 end if
                 ztz = ztz + z(i)*z(i)
              end do
              240 continue
           end if
           ! compute the fp vector downwards from r in blocks of size blksiz
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r,bn - 1
                 z(i + 1) = -(work(indumn + i)*z(i))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 260
                 end if
                 ztz = ztz + z(i + 1)*z(i + 1)
              end do
              260 continue
           else
              ! run slower loop if nan occurred.
              do i = r,bn - 1
                 if (z(i) == zero) then
                    z(i + 1) = -(ld(i - 1)/ld(i))*z(i - 1)
                 else
                    z(i + 1) = -(work(indumn + i)*z(i))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 280
                 end if
                 ztz = ztz + z(i + 1)*z(i + 1)
              end do
              280 continue
           end if
           ! compute quantities for convergence test
           tmp = one/ztz
           nrminv = sqrt(tmp)
           resid = abs(mingma)*nrminv
           rqcorr = mingma*tmp
           return
     end subroutine la_slar1v
     !> DLAR1V: computes the (scaled) r-th column of the inverse of
     !> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
     !> L D L**T - sigma I. When sigma is close to an eigenvalue, the
     !> computed vector is an accurate eigenvector. Usually, r corresponds
     !> to the index where the eigenvector is largest in magnitude.
     !> The following steps accomplish this computation :
     !> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
     !> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
     !> (c) Computation of the diagonal elements of the inverse of
     !> L D L**T - sigma I by combining the above transforms, and choosing
     !> r as the index where the diagonal of the inverse is (one of the)
     !> largest in magnitude.
     !> (d) Computation of the (scaled) r-th column of the inverse using the
     !> twisted factorization obtained by combining the top part of the
     !> the stationary and the bottom part of the progressive transform.

     pure subroutine la_dlar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc, &
               negcnt,ztz,mingma,r,isuppz,nrminv,resid,rqcorr,work)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantnc
           integer(ilp),intent(in) :: b1,bn,n
           integer(ilp),intent(out) :: negcnt
           integer(ilp),intent(inout) :: r
           real(dp),intent(in) :: gaptol,lambda,pivmin
           real(dp),intent(out) :: mingma,nrminv,resid,rqcorr,ztz
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*)
           real(dp),intent(in) :: d(*),l(*),ld(*),lld(*)
           real(dp),intent(out) :: work(*)
           real(dp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: sawnan1,sawnan2
           integer(ilp) :: i,indlpl,indp,inds,indumn,neg1,neg2,r1,r2
           real(dp) :: dminus,dplus,eps,s,tmp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           eps = la_dlamch('PRECISION')
           if (r == 0) then
              r1 = b1
              r2 = bn
           else
              r1 = r
              r2 = r
           end if
           ! storage for lplus
           indlpl = 0
           ! storage for uminus
           indumn = n
           inds = 2*n + 1
           indp = 3*n + 1
           if (b1 == 1) then
              work(inds) = zero
           else
              work(inds + b1 - 1) = lld(b1 - 1)
           end if
           ! compute the stationary transform (using the differential form)
           ! until the index r2.
           sawnan1 = .false.
           neg1 = 0
           s = work(inds + b1 - 1) - lambda
           do i = b1,r1 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              if (dplus < zero) neg1 = neg1 + 1
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_disnan(s)
           if (sawnan1) goto 60
           do i = r1,r2 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_disnan(s)
           60 continue
           if (sawnan1) then
              ! runs a slower version of the above loop if a nan is detected
              neg1 = 0
              s = work(inds + b1 - 1) - lambda
              do i = b1,r1 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 if (dplus < zero) neg1 = neg1 + 1
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
              do i = r1,r2 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
           end if
           ! compute the progressive transform (using the differential form)
           ! until the index r1
           sawnan2 = .false.
           neg2 = 0
           work(indp + bn - 1) = d(bn) - lambda
           do i = bn - 1,r1,-1
              dminus = lld(i) + work(indp + i)
              tmp = d(i)/dminus
              if (dminus < zero) neg2 = neg2 + 1
              work(indumn + i) = l(i)*tmp
              work(indp + i - 1) = work(indp + i)*tmp - lambda
           end do
           tmp = work(indp + r1 - 1)
           sawnan2 = la_disnan(tmp)
           if (sawnan2) then
              ! runs a slower version of the above loop if a nan is detected
              neg2 = 0
              do i = bn - 1,r1,-1
                 dminus = lld(i) + work(indp + i)
                 if (abs(dminus) < pivmin) dminus = -pivmin
                 tmp = d(i)/dminus
                 if (dminus < zero) neg2 = neg2 + 1
                 work(indumn + i) = l(i)*tmp
                 work(indp + i - 1) = work(indp + i)*tmp - lambda
                 if (tmp == zero) work(indp + i - 1) = d(i) - lambda
              end do
           end if
           ! find the index (from r1 to r2) of the largest (in magnitude)
           ! diagonal element of the inverse
           mingma = work(inds + r1 - 1) + work(indp + r1 - 1)
           if (mingma < zero) neg1 = neg1 + 1
           if (wantnc) then
              negcnt = neg1 + neg2
           else
              negcnt = -1
           end if
           if (abs(mingma) == zero) mingma = eps*work(inds + r1 - 1)
           r = r1
           do i = r1,r2 - 1
              tmp = work(inds + i) + work(indp + i)
              if (tmp == zero) tmp = eps*work(inds + i)
              if (abs(tmp) <= abs(mingma)) then
                 mingma = tmp
                 r = i + 1
              end if
           end do
           ! compute the fp vector: solve n^t v = e_r
           isuppz(1) = b1
           isuppz(2) = bn
           z(r) = one
           ztz = one
           ! compute the fp vector upwards from r
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r - 1,b1,-1
                 z(i) = -(work(indlpl + i)*z(i + 1))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    goto 220
                 end if
                 ztz = ztz + z(i)*z(i)
              end do
              220 continue
           else
              ! run slower loop if nan occurred.
              do i = r - 1,b1,-1
                 if (z(i + 1) == zero) then
                    z(i) = -(ld(i + 1)/ld(i))*z(i + 2)
                 else
                    z(i) = -(work(indlpl + i)*z(i + 1))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    go to 240
                 end if
                 ztz = ztz + z(i)*z(i)
              end do
              240 continue
           end if
           ! compute the fp vector downwards from r in blocks of size blksiz
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r,bn - 1
                 z(i + 1) = -(work(indumn + i)*z(i))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 260
                 end if
                 ztz = ztz + z(i + 1)*z(i + 1)
              end do
              260 continue
           else
              ! run slower loop if nan occurred.
              do i = r,bn - 1
                 if (z(i) == zero) then
                    z(i + 1) = -(ld(i - 1)/ld(i))*z(i - 1)
                 else
                    z(i + 1) = -(work(indumn + i)*z(i))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 280
                 end if
                 ztz = ztz + z(i + 1)*z(i + 1)
              end do
              280 continue
           end if
           ! compute quantities for convergence test
           tmp = one/ztz
           nrminv = sqrt(tmp)
           resid = abs(mingma)*nrminv
           rqcorr = mingma*tmp
           return
     end subroutine la_dlar1v
     !> QLAR1V: computes the (scaled) r-th column of the inverse of
     !> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
     !> L D L**T - sigma I. When sigma is close to an eigenvalue, the
     !> computed vector is an accurate eigenvector. Usually, r corresponds
     !> to the index where the eigenvector is largest in magnitude.
     !> The following steps accomplish this computation :
     !> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
     !> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
     !> (c) Computation of the diagonal elements of the inverse of
     !> L D L**T - sigma I by combining the above transforms, and choosing
     !> r as the index where the diagonal of the inverse is (one of the)
     !> largest in magnitude.
     !> (d) Computation of the (scaled) r-th column of the inverse using the
     !> twisted factorization obtained by combining the top part of the
     !> the stationary and the bottom part of the progressive transform.

     pure subroutine la_qlar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc, &
               negcnt,ztz,mingma,r,isuppz,nrminv,resid,rqcorr,work)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantnc
           integer(ilp),intent(in) :: b1,bn,n
           integer(ilp),intent(out) :: negcnt
           integer(ilp),intent(inout) :: r
           real(qp),intent(in) :: gaptol,lambda,pivmin
           real(qp),intent(out) :: mingma,nrminv,resid,rqcorr,ztz
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*)
           real(qp),intent(in) :: d(*),l(*),ld(*),lld(*)
           real(qp),intent(out) :: work(*)
           real(qp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: sawnan1,sawnan2
           integer(ilp) :: i,indlpl,indp,inds,indumn,neg1,neg2,r1,r2
           real(qp) :: dminus,dplus,eps,s,tmp
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           eps = la_qlamch('PRECISION')
           if (r == 0) then
              r1 = b1
              r2 = bn
           else
              r1 = r
              r2 = r
           end if
           ! storage for lplus
           indlpl = 0
           ! storage for uminus
           indumn = n
           inds = 2*n + 1
           indp = 3*n + 1
           if (b1 == 1) then
              work(inds) = zero
           else
              work(inds + b1 - 1) = lld(b1 - 1)
           end if
           ! compute the stationary transform (using the differential form)
           ! until the index r2.
           sawnan1 = .false.
           neg1 = 0
           s = work(inds + b1 - 1) - lambda
           do i = b1,r1 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              if (dplus < zero) neg1 = neg1 + 1
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_qisnan(s)
           if (sawnan1) goto 60
           do i = r1,r2 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_qisnan(s)
           60 continue
           if (sawnan1) then
              ! runs a slower version of the above loop if a nan is detected
              neg1 = 0
              s = work(inds + b1 - 1) - lambda
              do i = b1,r1 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 if (dplus < zero) neg1 = neg1 + 1
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
              do i = r1,r2 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
           end if
           ! compute the progressive transform (using the differential form)
           ! until the index r1
           sawnan2 = .false.
           neg2 = 0
           work(indp + bn - 1) = d(bn) - lambda
           do i = bn - 1,r1,-1
              dminus = lld(i) + work(indp + i)
              tmp = d(i)/dminus
              if (dminus < zero) neg2 = neg2 + 1
              work(indumn + i) = l(i)*tmp
              work(indp + i - 1) = work(indp + i)*tmp - lambda
           end do
           tmp = work(indp + r1 - 1)
           sawnan2 = la_qisnan(tmp)
           if (sawnan2) then
              ! runs a slower version of the above loop if a nan is detected
              neg2 = 0
              do i = bn - 1,r1,-1
                 dminus = lld(i) + work(indp + i)
                 if (abs(dminus) < pivmin) dminus = -pivmin
                 tmp = d(i)/dminus
                 if (dminus < zero) neg2 = neg2 + 1
                 work(indumn + i) = l(i)*tmp
                 work(indp + i - 1) = work(indp + i)*tmp - lambda
                 if (tmp == zero) work(indp + i - 1) = d(i) - lambda
              end do
           end if
           ! find the index (from r1 to r2) of the largest (in magnitude)
           ! diagonal element of the inverse
           mingma = work(inds + r1 - 1) + work(indp + r1 - 1)
           if (mingma < zero) neg1 = neg1 + 1
           if (wantnc) then
              negcnt = neg1 + neg2
           else
              negcnt = -1
           end if
           if (abs(mingma) == zero) mingma = eps*work(inds + r1 - 1)
           r = r1
           do i = r1,r2 - 1
              tmp = work(inds + i) + work(indp + i)
              if (tmp == zero) tmp = eps*work(inds + i)
              if (abs(tmp) <= abs(mingma)) then
                 mingma = tmp
                 r = i + 1
              end if
           end do
           ! compute the fp vector: solve n^t v = e_r
           isuppz(1) = b1
           isuppz(2) = bn
           z(r) = one
           ztz = one
           ! compute the fp vector upwards from r
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r - 1,b1,-1
                 z(i) = -(work(indlpl + i)*z(i + 1))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    goto 220
                 end if
                 ztz = ztz + z(i)*z(i)
              end do
              220 continue
           else
              ! run slower loop if nan occurred.
              do i = r - 1,b1,-1
                 if (z(i + 1) == zero) then
                    z(i) = -(ld(i + 1)/ld(i))*z(i + 2)
                 else
                    z(i) = -(work(indlpl + i)*z(i + 1))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    go to 240
                 end if
                 ztz = ztz + z(i)*z(i)
              end do
              240 continue
           end if
           ! compute the fp vector downwards from r in blocks of size blksiz
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r,bn - 1
                 z(i + 1) = -(work(indumn + i)*z(i))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 260
                 end if
                 ztz = ztz + z(i + 1)*z(i + 1)
              end do
              260 continue
           else
              ! run slower loop if nan occurred.
              do i = r,bn - 1
                 if (z(i) == zero) then
                    z(i + 1) = -(ld(i - 1)/ld(i))*z(i - 1)
                 else
                    z(i + 1) = -(work(indumn + i)*z(i))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 280
                 end if
                 ztz = ztz + z(i + 1)*z(i + 1)
              end do
              280 continue
           end if
           ! compute quantities for convergence test
           tmp = one/ztz
           nrminv = sqrt(tmp)
           resid = abs(mingma)*nrminv
           rqcorr = mingma*tmp
           return
     end subroutine la_qlar1v

     !> Given the relatively robust representation(RRR) L D L^T, SLARRB:
     !> does "limited" bisection to refine the eigenvalues of L D L^T,
     !> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
     !> guesses for these eigenvalues are input in W, the corresponding estimate
     !> of the error in these guesses and their gaps are input in WERR
     !> and WGAP, respectively. During bisection, intervals
     !> [left, right] are maintained by storing their mid-points and
     !> semi-widths in the arrays W and WERR respectively.

     pure subroutine la_slarrb(n,d,lld,ifirst,ilast,rtol1,rtol2,offset,w,wgap,werr, &
               work,iwork,pivmin,spdiam,twist,info)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ifirst,ilast,n,offset,twist
           integer(ilp),intent(out) :: info
           real(sp),intent(in) :: pivmin,rtol1,rtol2,spdiam
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(sp),intent(in) :: d(*),lld(*)
           real(sp),intent(inout) :: w(*),werr(*),wgap(*)
           real(sp),intent(out) :: work(*)
        ! =====================================================================

           integer(ilp) :: maxitr
           ! Local Scalars
           integer(ilp) :: i,i1,ii,ip,iter,k,negcnt,next,nint,olnint,prev,r
           real(sp) :: back,cvrgd,gap,left,lgap,mid,mnwdth,rgap,right,tmp,width
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           maxitr = int((log(spdiam + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           mnwdth = two*pivmin
           r = twist
           if ((r < 1) .or. (r > n)) r = n
           ! initialize unconverged intervals in [ work(2*i-1), work(2*i) ].
           ! the sturm count, count( work(2*i-1) ) is arranged to be i-1, while
           ! count( work(2*i) ) is stored in iwork( 2*i ). the integer iwork( 2*i-1 )
           ! for an unconverged interval is set to the index of the next unconverged
           ! interval, and is -1 or 0 for a converged interval. thus a linked
           ! list of unconverged intervals is set up.
           i1 = ifirst
           ! the number of unconverged intervals
           nint = 0
           ! the last unconverged interval found
           prev = 0
           rgap = wgap(i1 - offset)
           loop_75: do i = i1,ilast
              k = 2*i
              ii = i - offset
              left = w(ii) - werr(ii)
              right = w(ii) + werr(ii)
              lgap = rgap
              rgap = wgap(ii)
              gap = min(lgap,rgap)
              ! make sure that [left,right] contains the desired eigenvalue
              ! compute negcount from dstqds facto l+d+l+^t = l d l^t - left
              ! do while( negcnt(left)>i-1 )
              back = werr(ii)
              20 continue
              negcnt = la_slaneg(n,d,lld,left,pivmin,r)
              if (negcnt > i - 1) then
                 left = left - back
                 back = two*back
                 go to 20
              end if
              ! do while( negcnt(right)<i )
              ! compute negcount from dstqds facto l+d+l+^t = l d l^t - right
              back = werr(ii)
              50 continue
              negcnt = la_slaneg(n,d,lld,right,pivmin,r)
               if (negcnt < i) then
                  right = right + back
                  back = two*back
                  go to 50
               end if
              width = half*abs(left - right)
              tmp = max(abs(left),abs(right))
              cvrgd = max(rtol1*gap,rtol2*tmp)
              if (width <= cvrgd .or. width <= mnwdth) then
                 ! this interval has already converged and does not need refinement.
                 ! (note that the gaps might change through refining the
                  ! eigenvalues, however, they can only get bigger.)
                 ! remove it from the list.
                 iwork(k - 1) = -1
                 ! make sure that i1 always points to the first unconverged interval
                 if ((i == i1) .and. (i < ilast)) i1 = i + 1
                 if ((prev >= i1) .and. (i <= ilast)) iwork(2*prev - 1) = i + 1
              else
                 ! unconverged interval found
                 prev = i
                 nint = nint + 1
                 iwork(k - 1) = i + 1
                 iwork(k) = negcnt
              end if
              work(k - 1) = left
              work(k) = right
           end do loop_75
           ! do while( nint>0 ), i.e. there are still unconverged intervals
           ! and while (iter<maxitr)
           iter = 0
           80 continue
           prev = i1 - 1
           i = i1
           olnint = nint
           loop_100: do ip = 1,olnint
              k = 2*i
              ii = i - offset
              rgap = wgap(ii)
              lgap = rgap
              if (ii > 1) lgap = wgap(ii - 1)
              gap = min(lgap,rgap)
              next = iwork(k - 1)
              left = work(k - 1)
              right = work(k)
              mid = half*(left + right)
              ! semiwidth of interval
              width = right - mid
              tmp = max(abs(left),abs(right))
              cvrgd = max(rtol1*gap,rtol2*tmp)
              if ((width <= cvrgd) .or. (width <= mnwdth) .or. (iter == maxitr)) then
                 ! reduce number of unconverged intervals
                 nint = nint - 1
                 ! mark interval as converged.
                 iwork(k - 1) = 0
                 if (i1 == i) then
                    i1 = next
                 else
                    ! prev holds the last unconverged interval previously examined
                    if (prev >= i1) iwork(2*prev - 1) = next
                 end if
                 i = next
                 cycle loop_100
              end if
              prev = i
              ! perform one bisection step
              negcnt = la_slaneg(n,d,lld,mid,pivmin,r)
              if (negcnt <= i - 1) then
                 work(k - 1) = mid
              else
                 work(k) = mid
              end if
              i = next
           end do loop_100
           iter = iter + 1
           ! do another loop if there are still unconverged intervals
           ! however, in the last iteration, all intervals are accepted
           ! since this is the best we can do.
           if ((nint > 0) .and. (iter <= maxitr)) go to 80
           ! at this point, all the intervals have converged
           do i = ifirst,ilast
              k = 2*i
              ii = i - offset
              ! all intervals marked by '0' have been refined.
              if (iwork(k - 1) == 0) then
                 w(ii) = half*(work(k - 1) + work(k))
                 werr(ii) = work(k) - w(ii)
              end if
           end do
           do i = ifirst + 1,ilast
              k = 2*i
              ii = i - offset
              wgap(ii - 1) = max(zero,w(ii) - werr(ii) - w(ii - 1) - werr(ii - 1))
           end do
           return
     end subroutine la_slarrb
     !> Given the relatively robust representation(RRR) L D L^T, DLARRB:
     !> does "limited" bisection to refine the eigenvalues of L D L^T,
     !> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
     !> guesses for these eigenvalues are input in W, the corresponding estimate
     !> of the error in these guesses and their gaps are input in WERR
     !> and WGAP, respectively. During bisection, intervals
     !> [left, right] are maintained by storing their mid-points and
     !> semi-widths in the arrays W and WERR respectively.

     pure subroutine la_dlarrb(n,d,lld,ifirst,ilast,rtol1,rtol2,offset,w,wgap,werr, &
               work,iwork,pivmin,spdiam,twist,info)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ifirst,ilast,n,offset,twist
           integer(ilp),intent(out) :: info
           real(dp),intent(in) :: pivmin,rtol1,rtol2,spdiam
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(dp),intent(in) :: d(*),lld(*)
           real(dp),intent(inout) :: w(*),werr(*),wgap(*)
           real(dp),intent(out) :: work(*)
        ! =====================================================================

           integer(ilp) :: maxitr
           ! Local Scalars
           integer(ilp) :: i,i1,ii,ip,iter,k,negcnt,next,nint,olnint,prev,r
           real(dp) :: back,cvrgd,gap,left,lgap,mid,mnwdth,rgap,right,tmp,width
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           maxitr = int((log(spdiam + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           mnwdth = two*pivmin
           r = twist
           if ((r < 1) .or. (r > n)) r = n
           ! initialize unconverged intervals in [ work(2*i-1), work(2*i) ].
           ! the sturm count, count( work(2*i-1) ) is arranged to be i-1, while
           ! count( work(2*i) ) is stored in iwork( 2*i ). the integer iwork( 2*i-1 )
           ! for an unconverged interval is set to the index of the next unconverged
           ! interval, and is -1 or 0 for a converged interval. thus a linked
           ! list of unconverged intervals is set up.
           i1 = ifirst
           ! the number of unconverged intervals
           nint = 0
           ! the last unconverged interval found
           prev = 0
           rgap = wgap(i1 - offset)
           loop_75: do i = i1,ilast
              k = 2*i
              ii = i - offset
              left = w(ii) - werr(ii)
              right = w(ii) + werr(ii)
              lgap = rgap
              rgap = wgap(ii)
              gap = min(lgap,rgap)
              ! make sure that [left,right] contains the desired eigenvalue
              ! compute negcount from dstqds facto l+d+l+^t = l d l^t - left
              ! do while( negcnt(left)>i-1 )
              back = werr(ii)
              20 continue
              negcnt = la_dlaneg(n,d,lld,left,pivmin,r)
              if (negcnt > i - 1) then
                 left = left - back
                 back = two*back
                 go to 20
              end if
              ! do while( negcnt(right)<i )
              ! compute negcount from dstqds facto l+d+l+^t = l d l^t - right
              back = werr(ii)
              50 continue
              negcnt = la_dlaneg(n,d,lld,right,pivmin,r)
               if (negcnt < i) then
                  right = right + back
                  back = two*back
                  go to 50
               end if
              width = half*abs(left - right)
              tmp = max(abs(left),abs(right))
              cvrgd = max(rtol1*gap,rtol2*tmp)
              if (width <= cvrgd .or. width <= mnwdth) then
                 ! this interval has already converged and does not need refinement.
                 ! (note that the gaps might change through refining the
                  ! eigenvalues, however, they can only get bigger.)
                 ! remove it from the list.
                 iwork(k - 1) = -1
                 ! make sure that i1 always points to the first unconverged interval
                 if ((i == i1) .and. (i < ilast)) i1 = i + 1
                 if ((prev >= i1) .and. (i <= ilast)) iwork(2*prev - 1) = i + 1
              else
                 ! unconverged interval found
                 prev = i
                 nint = nint + 1
                 iwork(k - 1) = i + 1
                 iwork(k) = negcnt
              end if
              work(k - 1) = left
              work(k) = right
           end do loop_75
           ! do while( nint>0 ), i.e. there are still unconverged intervals
           ! and while (iter<maxitr)
           iter = 0
           80 continue
           prev = i1 - 1
           i = i1
           olnint = nint
           loop_100: do ip = 1,olnint
              k = 2*i
              ii = i - offset
              rgap = wgap(ii)
              lgap = rgap
              if (ii > 1) lgap = wgap(ii - 1)
              gap = min(lgap,rgap)
              next = iwork(k - 1)
              left = work(k - 1)
              right = work(k)
              mid = half*(left + right)
              ! semiwidth of interval
              width = right - mid
              tmp = max(abs(left),abs(right))
              cvrgd = max(rtol1*gap,rtol2*tmp)
              if ((width <= cvrgd) .or. (width <= mnwdth) .or. (iter == maxitr)) then
                 ! reduce number of unconverged intervals
                 nint = nint - 1
                 ! mark interval as converged.
                 iwork(k - 1) = 0
                 if (i1 == i) then
                    i1 = next
                 else
                    ! prev holds the last unconverged interval previously examined
                    if (prev >= i1) iwork(2*prev - 1) = next
                 end if
                 i = next
                 cycle loop_100
              end if
              prev = i
              ! perform one bisection step
              negcnt = la_dlaneg(n,d,lld,mid,pivmin,r)
              if (negcnt <= i - 1) then
                 work(k - 1) = mid
              else
                 work(k) = mid
              end if
              i = next
           end do loop_100
           iter = iter + 1
           ! do another loop if there are still unconverged intervals
           ! however, in the last iteration, all intervals are accepted
           ! since this is the best we can do.
           if ((nint > 0) .and. (iter <= maxitr)) go to 80
           ! at this point, all the intervals have converged
           do i = ifirst,ilast
              k = 2*i
              ii = i - offset
              ! all intervals marked by '0' have been refined.
              if (iwork(k - 1) == 0) then
                 w(ii) = half*(work(k - 1) + work(k))
                 werr(ii) = work(k) - w(ii)
              end if
           end do
           do i = ifirst + 1,ilast
              k = 2*i
              ii = i - offset
              wgap(ii - 1) = max(zero,w(ii) - werr(ii) - w(ii - 1) - werr(ii - 1))
           end do
           return
     end subroutine la_dlarrb
     !> Given the relatively robust representation(RRR) L D L^T, QLARRB:
     !> does "limited" bisection to refine the eigenvalues of L D L^T,
     !> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
     !> guesses for these eigenvalues are input in W, the corresponding estimate
     !> of the error in these guesses and their gaps are input in WERR
     !> and WGAP, respectively. During bisection, intervals
     !> [left, right] are maintained by storing their mid-points and
     !> semi-widths in the arrays W and WERR respectively.

     pure subroutine la_qlarrb(n,d,lld,ifirst,ilast,rtol1,rtol2,offset,w,wgap,werr, &
               work,iwork,pivmin,spdiam,twist,info)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: ifirst,ilast,n,offset,twist
           integer(ilp),intent(out) :: info
           real(qp),intent(in) :: pivmin,rtol1,rtol2,spdiam
           ! Array Arguments
           integer(ilp),intent(out) :: iwork(*)
           real(qp),intent(in) :: d(*),lld(*)
           real(qp),intent(inout) :: w(*),werr(*),wgap(*)
           real(qp),intent(out) :: work(*)
        ! =====================================================================

           integer(ilp) :: maxitr
           ! Local Scalars
           integer(ilp) :: i,i1,ii,ip,iter,k,negcnt,next,nint,olnint,prev,r
           real(qp) :: back,cvrgd,gap,left,lgap,mid,mnwdth,rgap,right,tmp,width
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           maxitr = int((log(spdiam + pivmin) - log(pivmin))/log(two),KIND=ilp) + 2
           mnwdth = two*pivmin
           r = twist
           if ((r < 1) .or. (r > n)) r = n
           ! initialize unconverged intervals in [ work(2*i-1), work(2*i) ].
           ! the sturm count, count( work(2*i-1) ) is arranged to be i-1, while
           ! count( work(2*i) ) is stored in iwork( 2*i ). the integer iwork( 2*i-1 )
           ! for an unconverged interval is set to the index of the next unconverged
           ! interval, and is -1 or 0 for a converged interval. thus a linked
           ! list of unconverged intervals is set up.
           i1 = ifirst
           ! the number of unconverged intervals
           nint = 0
           ! the last unconverged interval found
           prev = 0
           rgap = wgap(i1 - offset)
           loop_75: do i = i1,ilast
              k = 2*i
              ii = i - offset
              left = w(ii) - werr(ii)
              right = w(ii) + werr(ii)
              lgap = rgap
              rgap = wgap(ii)
              gap = min(lgap,rgap)
              ! make sure that [left,right] contains the desired eigenvalue
              ! compute negcount from dstqds facto l+d+l+^t = l d l^t - left
              ! do while( negcnt(left)>i-1 )
              back = werr(ii)
              20 continue
              negcnt = la_qlaneg(n,d,lld,left,pivmin,r)
              if (negcnt > i - 1) then
                 left = left - back
                 back = two*back
                 go to 20
              end if
              ! do while( negcnt(right)<i )
              ! compute negcount from dstqds facto l+d+l+^t = l d l^t - right
              back = werr(ii)
              50 continue
              negcnt = la_qlaneg(n,d,lld,right,pivmin,r)
               if (negcnt < i) then
                  right = right + back
                  back = two*back
                  go to 50
               end if
              width = half*abs(left - right)
              tmp = max(abs(left),abs(right))
              cvrgd = max(rtol1*gap,rtol2*tmp)
              if (width <= cvrgd .or. width <= mnwdth) then
                 ! this interval has already converged and does not need refinement.
                 ! (note that the gaps might change through refining the
                  ! eigenvalues, however, they can only get bigger.)
                 ! remove it from the list.
                 iwork(k - 1) = -1
                 ! make sure that i1 always points to the first unconverged interval
                 if ((i == i1) .and. (i < ilast)) i1 = i + 1
                 if ((prev >= i1) .and. (i <= ilast)) iwork(2*prev - 1) = i + 1
              else
                 ! unconverged interval found
                 prev = i
                 nint = nint + 1
                 iwork(k - 1) = i + 1
                 iwork(k) = negcnt
              end if
              work(k - 1) = left
              work(k) = right
           end do loop_75
           ! do while( nint>0 ), i.e. there are still unconverged intervals
           ! and while (iter<maxitr)
           iter = 0
           80 continue
           prev = i1 - 1
           i = i1
           olnint = nint
           loop_100: do ip = 1,olnint
              k = 2*i
              ii = i - offset
              rgap = wgap(ii)
              lgap = rgap
              if (ii > 1) lgap = wgap(ii - 1)
              gap = min(lgap,rgap)
              next = iwork(k - 1)
              left = work(k - 1)
              right = work(k)
              mid = half*(left + right)
              ! semiwidth of interval
              width = right - mid
              tmp = max(abs(left),abs(right))
              cvrgd = max(rtol1*gap,rtol2*tmp)
              if ((width <= cvrgd) .or. (width <= mnwdth) .or. (iter == maxitr)) then
                 ! reduce number of unconverged intervals
                 nint = nint - 1
                 ! mark interval as converged.
                 iwork(k - 1) = 0
                 if (i1 == i) then
                    i1 = next
                 else
                    ! prev holds the last unconverged interval previously examined
                    if (prev >= i1) iwork(2*prev - 1) = next
                 end if
                 i = next
                 cycle loop_100
              end if
              prev = i
              ! perform one bisection step
              negcnt = la_qlaneg(n,d,lld,mid,pivmin,r)
              if (negcnt <= i - 1) then
                 work(k - 1) = mid
              else
                 work(k) = mid
              end if
              i = next
           end do loop_100
           iter = iter + 1
           ! do another loop if there are still unconverged intervals
           ! however, in the last iteration, all intervals are accepted
           ! since this is the best we can do.
           if ((nint > 0) .and. (iter <= maxitr)) go to 80
           ! at this point, all the intervals have converged
           do i = ifirst,ilast
              k = 2*i
              ii = i - offset
              ! all intervals marked by '0' have been refined.
              if (iwork(k - 1) == 0) then
                 w(ii) = half*(work(k - 1) + work(k))
                 werr(ii) = work(k) - w(ii)
              end if
           end do
           do i = ifirst + 1,ilast
              k = 2*i
              ii = i - offset
              wgap(ii - 1) = max(zero,w(ii) - werr(ii) - w(ii - 1) - werr(ii - 1))
           end do
           return
     end subroutine la_qlarrb

     !> Given the initial representation L D L^T and its cluster of close
     !> eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...
     !> W( CLEND ), SLARRF: finds a new relatively robust representation
     !> L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the
     !> eigenvalues of L(+) D(+) L(+)^T is relatively isolated.

     pure subroutine la_slarrf(n,d,l,ld,clstrt,clend,w,wgap,werr,spdiam,clgapl, &
               clgapr,pivmin,sigma,dplus,lplus,work,info)
        use la_constants_sp,only:one,two,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: clstrt,clend,n
           integer(ilp),intent(out) :: info
           real(sp),intent(in) :: clgapl,clgapr,pivmin,spdiam
           real(sp),intent(out) :: sigma
           ! Array Arguments
           real(sp),intent(in) :: d(*),l(*),ld(*),w(*),werr(*)
           real(sp),intent(out) :: dplus(*),lplus(*),work(*)
           real(sp),intent(inout) :: wgap(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: quart = 0.25_sp
           real(sp),parameter :: maxgrowth1 = 8._sp
           real(sp),parameter :: maxgrowth2 = 8._sp
           integer(ilp),parameter :: ktrymax = 1
           integer(ilp),parameter :: sleft = 1
           integer(ilp),parameter :: sright = 2

           ! Local Scalars
           logical(lk) :: dorrr1,forcer,nofail,sawnan1,sawnan2,tryrrr1
           integer(ilp) :: i,indx,ktry,shift
           real(sp) :: avgap,bestshift,clwdth,eps,fact,fail,fail2,growthbound,ldelta, &
           ldmax,lsigma,max1,max2,mingap,oldp,prod,rdelta,rdmax,rrr1,rrr2,rsigma,s, &
                     smlgrowth,tmp,znm2
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           fact = real(2**ktrymax,KIND=sp)
           eps = la_slamch('PRECISION')
           shift = 0
           forcer = .false.
           ! note that we cannot guarantee that for any of the shifts tried,
           ! the factorization has a small or even moderate element growth.
           ! there could be ritz values at both ends of the cluster and despite
           ! backing off, there are examples where all factorizations tried
           ! (in ieee mode, allowing zero pivots
           ! element growth.
           ! for this reason, we should use pivmin in this subroutine so that at
           ! least the l d l^t factorization exists. it can be checked afterwards
           ! whether the element growth caused bad residuals/orthogonality.
           ! decide whether the code should accept the best among all
           ! representations despite large element growth or signal info=1
           ! setting nofail to .false. for quick fix for bug 113
           nofail = .false.
           ! compute the average gap length of the cluster
           clwdth = abs(w(clend) - w(clstrt)) + werr(clend) + werr(clstrt)
           avgap = clwdth/real(clend - clstrt,KIND=sp)
           mingap = min(clgapl,clgapr)
           ! initial values for shifts to both ends of cluster
           lsigma = min(w(clstrt),w(clend)) - werr(clstrt)
           rsigma = max(w(clstrt),w(clend)) + werr(clend)
           ! use a small fudge to make sure that we really shift to the outside
           lsigma = lsigma - abs(lsigma)*two*eps
           rsigma = rsigma + abs(rsigma)*two*eps
           ! compute upper bounds for how much to back off the initial shifts
           ldmax = quart*mingap + two*pivmin
           rdmax = quart*mingap + two*pivmin
           ldelta = max(avgap,wgap(clstrt))/fact
           rdelta = max(avgap,wgap(clend - 1))/fact
           ! initialize the record of the best representation found
           s = la_slamch('S')
           smlgrowth = one/s
           fail = real(n - 1,KIND=sp)*mingap/(spdiam*eps)
           fail2 = real(n - 1,KIND=sp)*mingap/(spdiam*sqrt(eps))
           bestshift = lsigma
           ! while (ktry <= ktrymax)
           ktry = 0
           growthbound = maxgrowth1*spdiam
           5 continue
           sawnan1 = .false.
           sawnan2 = .false.
           ! ensure that we do not back off too much of the initial shifts
           ldelta = min(ldmax,ldelta)
           rdelta = min(rdmax,rdelta)
           ! compute the element growth when shifting to both ends of the cluster
           ! accept the shift if there is no element growth at one of the two ends
           ! left end
           s = -lsigma
           dplus(1) = d(1) + s
           if (abs(dplus(1)) < pivmin) then
              dplus(1) = -pivmin
              ! need to set sawnan1 because refined rrr test should not be used
              ! in this case
              sawnan1 = .true.
           end if
           max1 = abs(dplus(1))
           do i = 1,n - 1
              lplus(i) = ld(i)/dplus(i)
              s = s*lplus(i)*l(i) - lsigma
              dplus(i + 1) = d(i + 1) + s
              if (abs(dplus(i + 1)) < pivmin) then
                 dplus(i + 1) = -pivmin
                 ! need to set sawnan1 because refined rrr test should not be used
                 ! in this case
                 sawnan1 = .true.
              end if
              max1 = max(max1,abs(dplus(i + 1)))
           end do
           sawnan1 = sawnan1 .or. la_sisnan(max1)
           if (forcer .or. (max1 <= growthbound .and. .not. sawnan1)) then
              sigma = lsigma
              shift = sleft
              goto 100
           end if
           ! right end
           s = -rsigma
           work(1) = d(1) + s
           if (abs(work(1)) < pivmin) then
              work(1) = -pivmin
              ! need to set sawnan2 because refined rrr test should not be used
              ! in this case
              sawnan2 = .true.
           end if
           max2 = abs(work(1))
           do i = 1,n - 1
              work(n + i) = ld(i)/work(i)
              s = s*work(n + i)*l(i) - rsigma
              work(i + 1) = d(i + 1) + s
              if (abs(work(i + 1)) < pivmin) then
                 work(i + 1) = -pivmin
                 ! need to set sawnan2 because refined rrr test should not be used
                 ! in this case
                 sawnan2 = .true.
              end if
              max2 = max(max2,abs(work(i + 1)))
           end do
           sawnan2 = sawnan2 .or. la_sisnan(max2)
           if (forcer .or. (max2 <= growthbound .and. .not. sawnan2)) then
              sigma = rsigma
              shift = sright
              goto 100
           end if
           ! if we are at this point, both shifts led to too much element growth
           ! record the better of the two shifts (provided it didn't lead to nan)
           if (sawnan1 .and. sawnan2) then
              ! both max1 and max2 are nan
              goto 50
           else
              if (.not. sawnan1) then
                 indx = 1
                 if (max1 <= smlgrowth) then
                    smlgrowth = max1
                    bestshift = lsigma
                 end if
              end if
              if (.not. sawnan2) then
                 if (sawnan1 .or. max2 <= max1) indx = 2
                 if (max2 <= smlgrowth) then
                    smlgrowth = max2
                    bestshift = rsigma
                 end if
              end if
           end if
           ! if we are here, both the left and the right shift led to
           ! element growth. if the element growth is moderate, then
           ! we may still accept the representation, if it passes a
           ! refined test for rrr. this test supposes that no nan occurred.
           ! moreover, we use the refined rrr test only for isolated clusters.
           if ((clwdth < mingap/real(128,KIND=sp)) .and. (min(max1,max2) < fail2) .and. (.not. sawnan1) &
                     .and. (.not. sawnan2)) then
              dorrr1 = .true.
           else
              dorrr1 = .false.
           end if
           tryrrr1 = .true.
           if (tryrrr1 .and. dorrr1) then
           if (indx == 1) then
              tmp = abs(dplus(n))
              znm2 = one
              prod = one
              oldp = one
              do i = n - 1,1,-1
                 if (prod <= eps) then
                    prod = ((dplus(i + 1)*work(n + i + 1))/(dplus(i)*work(n + i)))*oldp
                 else
                    prod = prod*abs(work(n + i))
                 end if
                 oldp = prod
                 znm2 = znm2 + prod**2
                 tmp = max(tmp,abs(dplus(i)*prod))
              end do
              rrr1 = tmp/(spdiam*sqrt(znm2))
              if (rrr1 <= maxgrowth2) then
                 sigma = lsigma
                 shift = sleft
                 goto 100
              end if
           else if (indx == 2) then
              tmp = abs(work(n))
              znm2 = one
              prod = one
              oldp = one
              do i = n - 1,1,-1
                 if (prod <= eps) then
                    prod = ((work(i + 1)*lplus(i + 1))/(work(i)*lplus(i)))*oldp
                 else
                    prod = prod*abs(lplus(i))
                 end if
                 oldp = prod
                 znm2 = znm2 + prod**2
                 tmp = max(tmp,abs(work(i)*prod))
              end do
              rrr2 = tmp/(spdiam*sqrt(znm2))
              if (rrr2 <= maxgrowth2) then
                 sigma = rsigma
                 shift = sright
                 goto 100
              end if
           end if
           end if
           50 continue
           if (ktry < ktrymax) then
              ! if we are here, both shifts failed also the rrr test.
              ! back off to the outside
              lsigma = max(lsigma - ldelta,lsigma - ldmax)
              rsigma = min(rsigma + rdelta,rsigma + rdmax)
              ldelta = two*ldelta
              rdelta = two*rdelta
              ktry = ktry + 1
              goto 5
           else
              ! none of the representations investigated satisfied our
              ! criteria. take the best one we found.
              if ((smlgrowth < fail) .or. nofail) then
                 lsigma = bestshift
                 rsigma = bestshift
                 forcer = .true.
                 goto 5
              else
                 info = 1
                 return
              end if
           end if
           100 continue
           if (shift == sleft) then
           elseif (shift == sright) then
              ! store new l and d back into dplus, lplus
              call la_scopy(n,work,1,dplus,1)
              call la_scopy(n - 1,work(n + 1),1,lplus,1)
           end if
           return
     end subroutine la_slarrf
     !> Given the initial representation L D L^T and its cluster of close
     !> eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...
     !> W( CLEND ), DLARRF: finds a new relatively robust representation
     !> L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the
     !> eigenvalues of L(+) D(+) L(+)^T is relatively isolated.

     pure subroutine la_dlarrf(n,d,l,ld,clstrt,clend,w,wgap,werr,spdiam,clgapl, &
               clgapr,pivmin,sigma,dplus,lplus,work,info)
        use la_constants_dp,only:one,two,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: clstrt,clend,n
           integer(ilp),intent(out) :: info
           real(dp),intent(in) :: clgapl,clgapr,pivmin,spdiam
           real(dp),intent(out) :: sigma
           ! Array Arguments
           real(dp),intent(in) :: d(*),l(*),ld(*),w(*),werr(*)
           real(dp),intent(out) :: dplus(*),lplus(*),work(*)
           real(dp),intent(inout) :: wgap(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: quart = 0.25_dp
           real(dp),parameter :: maxgrowth1 = 8._dp
           real(dp),parameter :: maxgrowth2 = 8._dp
           integer(ilp),parameter :: ktrymax = 1
           integer(ilp),parameter :: sleft = 1
           integer(ilp),parameter :: sright = 2

           ! Local Scalars
           logical(lk) :: dorrr1,forcer,nofail,sawnan1,sawnan2,tryrrr1
           integer(ilp) :: i,indx,ktry,shift
           real(dp) :: avgap,bestshift,clwdth,eps,fact,fail,fail2,growthbound,ldelta, &
           ldmax,lsigma,max1,max2,mingap,oldp,prod,rdelta,rdmax,rrr1,rrr2,rsigma,s, &
                     smlgrowth,tmp,znm2
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           fact = real(2**ktrymax,KIND=dp)
           eps = la_dlamch('PRECISION')
           shift = 0
           forcer = .false.
           ! note that we cannot guarantee that for any of the shifts tried,
           ! the factorization has a small or even moderate element growth.
           ! there could be ritz values at both ends of the cluster and despite
           ! backing off, there are examples where all factorizations tried
           ! (in ieee mode, allowing zero pivots
           ! element growth.
           ! for this reason, we should use pivmin in this subroutine so that at
           ! least the l d l^t factorization exists. it can be checked afterwards
           ! whether the element growth caused bad residuals/orthogonality.
           ! decide whether the code should accept the best among all
           ! representations despite large element growth or signal info=1
           ! setting nofail to .false. for quick fix for bug 113
           nofail = .false.
           ! compute the average gap length of the cluster
           clwdth = abs(w(clend) - w(clstrt)) + werr(clend) + werr(clstrt)
           avgap = clwdth/real(clend - clstrt,KIND=dp)
           mingap = min(clgapl,clgapr)
           ! initial values for shifts to both ends of cluster
           lsigma = min(w(clstrt),w(clend)) - werr(clstrt)
           rsigma = max(w(clstrt),w(clend)) + werr(clend)
           ! use a small fudge to make sure that we really shift to the outside
           lsigma = lsigma - abs(lsigma)*four*eps
           rsigma = rsigma + abs(rsigma)*four*eps
           ! compute upper bounds for how much to back off the initial shifts
           ldmax = quart*mingap + two*pivmin
           rdmax = quart*mingap + two*pivmin
           ldelta = max(avgap,wgap(clstrt))/fact
           rdelta = max(avgap,wgap(clend - 1))/fact
           ! initialize the record of the best representation found
           s = la_dlamch('S')
           smlgrowth = one/s
           fail = real(n - 1,KIND=dp)*mingap/(spdiam*eps)
           fail2 = real(n - 1,KIND=dp)*mingap/(spdiam*sqrt(eps))
           bestshift = lsigma
           ! while (ktry <= ktrymax)
           ktry = 0
           growthbound = maxgrowth1*spdiam
           5 continue
           sawnan1 = .false.
           sawnan2 = .false.
           ! ensure that we do not back off too much of the initial shifts
           ldelta = min(ldmax,ldelta)
           rdelta = min(rdmax,rdelta)
           ! compute the element growth when shifting to both ends of the cluster
           ! accept the shift if there is no element growth at one of the two ends
           ! left end
           s = -lsigma
           dplus(1) = d(1) + s
           if (abs(dplus(1)) < pivmin) then
              dplus(1) = -pivmin
              ! need to set sawnan1 because refined rrr test should not be used
              ! in this case
              sawnan1 = .true.
           end if
           max1 = abs(dplus(1))
           do i = 1,n - 1
              lplus(i) = ld(i)/dplus(i)
              s = s*lplus(i)*l(i) - lsigma
              dplus(i + 1) = d(i + 1) + s
              if (abs(dplus(i + 1)) < pivmin) then
                 dplus(i + 1) = -pivmin
                 ! need to set sawnan1 because refined rrr test should not be used
                 ! in this case
                 sawnan1 = .true.
              end if
              max1 = max(max1,abs(dplus(i + 1)))
           end do
           sawnan1 = sawnan1 .or. la_disnan(max1)
           if (forcer .or. (max1 <= growthbound .and. .not. sawnan1)) then
              sigma = lsigma
              shift = sleft
              goto 100
           end if
           ! right end
           s = -rsigma
           work(1) = d(1) + s
           if (abs(work(1)) < pivmin) then
              work(1) = -pivmin
              ! need to set sawnan2 because refined rrr test should not be used
              ! in this case
              sawnan2 = .true.
           end if
           max2 = abs(work(1))
           do i = 1,n - 1
              work(n + i) = ld(i)/work(i)
              s = s*work(n + i)*l(i) - rsigma
              work(i + 1) = d(i + 1) + s
              if (abs(work(i + 1)) < pivmin) then
                 work(i + 1) = -pivmin
                 ! need to set sawnan2 because refined rrr test should not be used
                 ! in this case
                 sawnan2 = .true.
              end if
              max2 = max(max2,abs(work(i + 1)))
           end do
           sawnan2 = sawnan2 .or. la_disnan(max2)
           if (forcer .or. (max2 <= growthbound .and. .not. sawnan2)) then
              sigma = rsigma
              shift = sright
              goto 100
           end if
           ! if we are at this point, both shifts led to too much element growth
           ! record the better of the two shifts (provided it didn't lead to nan)
           if (sawnan1 .and. sawnan2) then
              ! both max1 and max2 are nan
              goto 50
           else
              if (.not. sawnan1) then
                 indx = 1
                 if (max1 <= smlgrowth) then
                    smlgrowth = max1
                    bestshift = lsigma
                 end if
              end if
              if (.not. sawnan2) then
                 if (sawnan1 .or. max2 <= max1) indx = 2
                 if (max2 <= smlgrowth) then
                    smlgrowth = max2
                    bestshift = rsigma
                 end if
              end if
           end if
           ! if we are here, both the left and the right shift led to
           ! element growth. if the element growth is moderate, then
           ! we may still accept the representation, if it passes a
           ! refined test for rrr. this test supposes that no nan occurred.
           ! moreover, we use the refined rrr test only for isolated clusters.
           if ((clwdth < mingap/real(128,KIND=dp)) .and. (min(max1,max2) < fail2) .and. (.not. sawnan1) &
                     .and. (.not. sawnan2)) then
              dorrr1 = .true.
           else
              dorrr1 = .false.
           end if
           tryrrr1 = .true.
           if (tryrrr1 .and. dorrr1) then
           if (indx == 1) then
              tmp = abs(dplus(n))
              znm2 = one
              prod = one
              oldp = one
              do i = n - 1,1,-1
                 if (prod <= eps) then
                    prod = ((dplus(i + 1)*work(n + i + 1))/(dplus(i)*work(n + i)))*oldp
                 else
                    prod = prod*abs(work(n + i))
                 end if
                 oldp = prod
                 znm2 = znm2 + prod**2
                 tmp = max(tmp,abs(dplus(i)*prod))
              end do
              rrr1 = tmp/(spdiam*sqrt(znm2))
              if (rrr1 <= maxgrowth2) then
                 sigma = lsigma
                 shift = sleft
                 goto 100
              end if
           else if (indx == 2) then
              tmp = abs(work(n))
              znm2 = one
              prod = one
              oldp = one
              do i = n - 1,1,-1
                 if (prod <= eps) then
                    prod = ((work(i + 1)*lplus(i + 1))/(work(i)*lplus(i)))*oldp
                 else
                    prod = prod*abs(lplus(i))
                 end if
                 oldp = prod
                 znm2 = znm2 + prod**2
                 tmp = max(tmp,abs(work(i)*prod))
              end do
              rrr2 = tmp/(spdiam*sqrt(znm2))
              if (rrr2 <= maxgrowth2) then
                 sigma = rsigma
                 shift = sright
                 goto 100
              end if
           end if
           end if
           50 continue
           if (ktry < ktrymax) then
              ! if we are here, both shifts failed also the rrr test.
              ! back off to the outside
              lsigma = max(lsigma - ldelta,lsigma - ldmax)
              rsigma = min(rsigma + rdelta,rsigma + rdmax)
              ldelta = two*ldelta
              rdelta = two*rdelta
              ktry = ktry + 1
              goto 5
           else
              ! none of the representations investigated satisfied our
              ! criteria. take the best one we found.
              if ((smlgrowth < fail) .or. nofail) then
                 lsigma = bestshift
                 rsigma = bestshift
                 forcer = .true.
                 goto 5
              else
                 info = 1
                 return
              end if
           end if
           100 continue
           if (shift == sleft) then
           elseif (shift == sright) then
              ! store new l and d back into dplus, lplus
              call la_dcopy(n,work,1,dplus,1)
              call la_dcopy(n - 1,work(n + 1),1,lplus,1)
           end if
           return
     end subroutine la_dlarrf
     !> Given the initial representation L D L^T and its cluster of close
     !> eigenvalues (in a relative measure), W( CLSTRT ), W( CLSTRT+1 ), ...
     !> W( CLEND ), QLARRF: finds a new relatively robust representation
     !> L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the
     !> eigenvalues of L(+) D(+) L(+)^T is relatively isolated.

     pure subroutine la_qlarrf(n,d,l,ld,clstrt,clend,w,wgap,werr,spdiam,clgapl, &
               clgapr,pivmin,sigma,dplus,lplus,work,info)
        use la_constants_qp,only:one,two,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: clstrt,clend,n
           integer(ilp),intent(out) :: info
           real(qp),intent(in) :: clgapl,clgapr,pivmin,spdiam
           real(qp),intent(out) :: sigma
           ! Array Arguments
           real(qp),intent(in) :: d(*),l(*),ld(*),w(*),werr(*)
           real(qp),intent(out) :: dplus(*),lplus(*),work(*)
           real(qp),intent(inout) :: wgap(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: quart = 0.25_qp
           real(qp),parameter :: maxgrowth1 = 8._qp
           real(qp),parameter :: maxgrowth2 = 8._qp
           integer(ilp),parameter :: ktrymax = 1
           integer(ilp),parameter :: sleft = 1
           integer(ilp),parameter :: sright = 2

           ! Local Scalars
           logical(lk) :: dorrr1,forcer,nofail,sawnan1,sawnan2,tryrrr1
           integer(ilp) :: i,indx,ktry,shift
           real(qp) :: avgap,bestshift,clwdth,eps,fact,fail,fail2,growthbound,ldelta, &
           ldmax,lsigma,max1,max2,mingap,oldp,prod,rdelta,rdmax,rrr1,rrr2,rsigma,s, &
                     smlgrowth,tmp,znm2
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           fact = real(2**ktrymax,KIND=qp)
           eps = la_qlamch('PRECISION')
           shift = 0
           forcer = .false.
           ! note that we cannot guarantee that for any of the shifts tried,
           ! the factorization has a small or even moderate element growth.
           ! there could be ritz values at both ends of the cluster and despite
           ! backing off, there are examples where all factorizations tried
           ! (in ieee mode, allowing zero pivots
           ! element growth.
           ! for this reason, we should use pivmin in this subroutine so that at
           ! least the l d l^t factorization exists. it can be checked afterwards
           ! whether the element growth caused bad residuals/orthogonality.
           ! decide whether the code should accept the best among all
           ! representations despite large element growth or signal info=1
           ! setting nofail to .false. for quick fix for bug 113
           nofail = .false.
           ! compute the average gap length of the cluster
           clwdth = abs(w(clend) - w(clstrt)) + werr(clend) + werr(clstrt)
           avgap = clwdth/real(clend - clstrt,KIND=qp)
           mingap = min(clgapl,clgapr)
           ! initial values for shifts to both ends of cluster
           lsigma = min(w(clstrt),w(clend)) - werr(clstrt)
           rsigma = max(w(clstrt),w(clend)) + werr(clend)
           ! use a small fudge to make sure that we really shift to the outside
           lsigma = lsigma - abs(lsigma)*four*eps
           rsigma = rsigma + abs(rsigma)*four*eps
           ! compute upper bounds for how much to back off the initial shifts
           ldmax = quart*mingap + two*pivmin
           rdmax = quart*mingap + two*pivmin
           ldelta = max(avgap,wgap(clstrt))/fact
           rdelta = max(avgap,wgap(clend - 1))/fact
           ! initialize the record of the best representation found
           s = la_qlamch('S')
           smlgrowth = one/s
           fail = real(n - 1,KIND=qp)*mingap/(spdiam*eps)
           fail2 = real(n - 1,KIND=qp)*mingap/(spdiam*sqrt(eps))
           bestshift = lsigma
           ! while (ktry <= ktrymax)
           ktry = 0
           growthbound = maxgrowth1*spdiam
           5 continue
           sawnan1 = .false.
           sawnan2 = .false.
           ! ensure that we do not back off too much of the initial shifts
           ldelta = min(ldmax,ldelta)
           rdelta = min(rdmax,rdelta)
           ! compute the element growth when shifting to both ends of the cluster
           ! accept the shift if there is no element growth at one of the two ends
           ! left end
           s = -lsigma
           dplus(1) = d(1) + s
           if (abs(dplus(1)) < pivmin) then
              dplus(1) = -pivmin
              ! need to set sawnan1 because refined rrr test should not be used
              ! in this case
              sawnan1 = .true.
           end if
           max1 = abs(dplus(1))
           do i = 1,n - 1
              lplus(i) = ld(i)/dplus(i)
              s = s*lplus(i)*l(i) - lsigma
              dplus(i + 1) = d(i + 1) + s
              if (abs(dplus(i + 1)) < pivmin) then
                 dplus(i + 1) = -pivmin
                 ! need to set sawnan1 because refined rrr test should not be used
                 ! in this case
                 sawnan1 = .true.
              end if
              max1 = max(max1,abs(dplus(i + 1)))
           end do
           sawnan1 = sawnan1 .or. la_qisnan(max1)
           if (forcer .or. (max1 <= growthbound .and. .not. sawnan1)) then
              sigma = lsigma
              shift = sleft
              goto 100
           end if
           ! right end
           s = -rsigma
           work(1) = d(1) + s
           if (abs(work(1)) < pivmin) then
              work(1) = -pivmin
              ! need to set sawnan2 because refined rrr test should not be used
              ! in this case
              sawnan2 = .true.
           end if
           max2 = abs(work(1))
           do i = 1,n - 1
              work(n + i) = ld(i)/work(i)
              s = s*work(n + i)*l(i) - rsigma
              work(i + 1) = d(i + 1) + s
              if (abs(work(i + 1)) < pivmin) then
                 work(i + 1) = -pivmin
                 ! need to set sawnan2 because refined rrr test should not be used
                 ! in this case
                 sawnan2 = .true.
              end if
              max2 = max(max2,abs(work(i + 1)))
           end do
           sawnan2 = sawnan2 .or. la_qisnan(max2)
           if (forcer .or. (max2 <= growthbound .and. .not. sawnan2)) then
              sigma = rsigma
              shift = sright
              goto 100
           end if
           ! if we are at this point, both shifts led to too much element growth
           ! record the better of the two shifts (provided it didn't lead to nan)
           if (sawnan1 .and. sawnan2) then
              ! both max1 and max2 are nan
              goto 50
           else
              if (.not. sawnan1) then
                 indx = 1
                 if (max1 <= smlgrowth) then
                    smlgrowth = max1
                    bestshift = lsigma
                 end if
              end if
              if (.not. sawnan2) then
                 if (sawnan1 .or. max2 <= max1) indx = 2
                 if (max2 <= smlgrowth) then
                    smlgrowth = max2
                    bestshift = rsigma
                 end if
              end if
           end if
           ! if we are here, both the left and the right shift led to
           ! element growth. if the element growth is moderate, then
           ! we may still accept the representation, if it passes a
           ! refined test for rrr. this test supposes that no nan occurred.
           ! moreover, we use the refined rrr test only for isolated clusters.
           if ((clwdth < mingap/real(128,KIND=qp)) .and. (min(max1,max2) < fail2) .and. (.not. sawnan1) &
                     .and. (.not. sawnan2)) then
              dorrr1 = .true.
           else
              dorrr1 = .false.
           end if
           tryrrr1 = .true.
           if (tryrrr1 .and. dorrr1) then
           if (indx == 1) then
              tmp = abs(dplus(n))
              znm2 = one
              prod = one
              oldp = one
              do i = n - 1,1,-1
                 if (prod <= eps) then
                    prod = ((dplus(i + 1)*work(n + i + 1))/(dplus(i)*work(n + i)))*oldp
                 else
                    prod = prod*abs(work(n + i))
                 end if
                 oldp = prod
                 znm2 = znm2 + prod**2
                 tmp = max(tmp,abs(dplus(i)*prod))
              end do
              rrr1 = tmp/(spdiam*sqrt(znm2))
              if (rrr1 <= maxgrowth2) then
                 sigma = lsigma
                 shift = sleft
                 goto 100
              end if
           else if (indx == 2) then
              tmp = abs(work(n))
              znm2 = one
              prod = one
              oldp = one
              do i = n - 1,1,-1
                 if (prod <= eps) then
                    prod = ((work(i + 1)*lplus(i + 1))/(work(i)*lplus(i)))*oldp
                 else
                    prod = prod*abs(lplus(i))
                 end if
                 oldp = prod
                 znm2 = znm2 + prod**2
                 tmp = max(tmp,abs(work(i)*prod))
              end do
              rrr2 = tmp/(spdiam*sqrt(znm2))
              if (rrr2 <= maxgrowth2) then
                 sigma = rsigma
                 shift = sright
                 goto 100
              end if
           end if
           end if
           50 continue
           if (ktry < ktrymax) then
              ! if we are here, both shifts failed also the rrr test.
              ! back off to the outside
              lsigma = max(lsigma - ldelta,lsigma - ldmax)
              rsigma = min(rsigma + rdelta,rsigma + rdmax)
              ldelta = two*ldelta
              rdelta = two*rdelta
              ktry = ktry + 1
              goto 5
           else
              ! none of the representations investigated satisfied our
              ! criteria. take the best one we found.
              if ((smlgrowth < fail) .or. nofail) then
                 lsigma = bestshift
                 rsigma = bestshift
                 forcer = .true.
                 goto 5
              else
                 info = 1
                 return
              end if
           end if
           100 continue
           if (shift == sleft) then
           elseif (shift == sright) then
              ! store new l and d back into dplus, lplus
              call la_qcopy(n,work,1,dplus,1)
              call la_qcopy(n - 1,work(n + 1),1,lplus,1)
           end if
           return
     end subroutine la_qlarrf

     !> SLARRV: computes the eigenvectors of the tridiagonal matrix
     !> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
     !> The input eigenvalues should have been computed by SLARRE.

     pure subroutine la_slarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1, &
               rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz,isuppz,work,iwork,info)
        use la_constants_sp,only:zero,half,one,two,three,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dol,dou,ldz,m,n
           integer(ilp),intent(out) :: info
           real(sp),intent(in) :: minrgp,pivmin,vl,vu
           real(sp),intent(inout) :: rtol1,rtol2
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),indexw(*),isplit(*)
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),l(*),w(*),werr(*),wgap(*)
           real(sp),intent(in) :: gers(*)
           real(sp),intent(out) :: work(*)
           real(sp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 10

           ! Local Scalars
           logical(lk) :: eskip,needbs,stp2ii,tryrqc,usedbs,usedrq
           integer(ilp) :: done,i,ibegin,idone,iend,ii,iindc1,iindc2,iindr,iindwk,iinfo, &
            im,in,indeig,indld,indlld,indwrk,isupmn,isupmx,iter,itmp1,j,jblk,k, &
            miniwsize,minwsize,nclus,ndepth,negcnt,newcls,newfst,newftt,newlst,newsiz, &
            offset,oldcls,oldfst,oldien,oldlst,oldncl,p,parity,q,wbegin,wend,windex, &
                      windmn,windpl,zfrom,zto,zusedl,zusedu,zusedw
           real(sp) :: bstres,bstw,eps,fudge,gap,gaptol,gl,gu,lambda,left,lgap,mingma, &
           nrminv,resid,rgap,right,rqcorr,rqtol,savgap,sgndef,sigma,spdiam,ssigma,tau, &
                     tmp,tol,ztz
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if ((n <= 0) .or. (m <= 0)) then
              return
           end if
           ! the first n entries of work are reserved for the eigenvalues
           indld = n + 1
           indlld = 2*n + 1
           indwrk = 3*n + 1
           minwsize = 12*n
           do i = 1,minwsize
              work(i) = zero
           end do
           ! iwork(iindr+1:iindr+n) hold the twist indices r for the
           ! factorization used to compute the fp vector
           iindr = 0
           ! iwork(iindc1+1:iinc2+n) are used to store the clusters of the current
           ! layer and the one above.
           iindc1 = n
           iindc2 = 2*n
           iindwk = 3*n + 1
           miniwsize = 7*n
           do i = 1,miniwsize
              iwork(i) = 0
           end do
           zusedl = 1
           if (dol > 1) then
              ! set lower bound for use of z
              zusedl = dol - 1
           end if
           zusedu = m
           if (dou < m) then
              ! set lower bound for use of z
              zusedu = dou + 1
           end if
           ! the width of the part of z that is used
           zusedw = zusedu - zusedl + 1
           call la_slaset('FULL',n,zusedw,zero,zero,z(1,zusedl),ldz)
           eps = la_slamch('PRECISION')
           rqtol = two*eps
           ! set expert flags for standard code.
           tryrqc = .true.
           if ((dol == 1) .and. (dou == m)) then
           else
              ! only selected eigenpairs are computed. since the other evalues
              ! are not refined by rq iteration, bisection has to compute to full
              ! accuracy.
              rtol1 = four*eps
              rtol2 = four*eps
           end if
           ! the entries wbegin:wend in w, werr, wgap correspond to the
           ! desired eigenvalues. the support of the nonzero eigenvector
           ! entries is contained in the interval ibegin:iend.
           ! remark that if k eigenpairs are desired, then the eigenvectors
           ! are stored in k contiguous columns of z.
           ! done is the number of eigenvectors already computed
           done = 0
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,iblock(m)
              iend = isplit(jblk)
              sigma = l(iend)
              ! find the eigenvectors of the submatrix indexed ibegin
              ! through iend.
              wend = wbegin - 1
              15 continue
              if (wend < m) then
                 if (iblock(wend + 1) == jblk) then
                    wend = wend + 1
                    go to 15
                 end if
              end if
              if (wend < wbegin) then
                 ibegin = iend + 1
                 cycle loop_170
              elseif ((wend < dol) .or. (wbegin > dou)) then
                 ibegin = iend + 1
                 wbegin = wend + 1
                 cycle loop_170
              end if
              ! find local spectral diameter of the block
              gl = gers(2*ibegin - 1)
              gu = gers(2*ibegin)
              do i = ibegin + 1,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              ! oldien is the last index of the previous block
              oldien = ibegin - 1
              ! calculate the size of the current block
              in = iend - ibegin + 1
              ! the number of eigenvalues in the current block
              im = wend - wbegin + 1
              ! this is for a 1x1 block
              if (ibegin == iend) then
                 done = done + 1
                 z(ibegin,wbegin) = one
                 isuppz(2*wbegin - 1) = ibegin
                 isuppz(2*wbegin) = ibegin
                 w(wbegin) = w(wbegin) + sigma
                 work(wbegin) = w(wbegin)
                 ibegin = iend + 1
                 wbegin = wbegin + 1
                 cycle loop_170
              end if
              ! the desired (shifted) eigenvalues are stored in w(wbegin:wend)
              ! note that these can be approximations, in this case, the corresp.
              ! entries of werr give the size of the uncertainty interval.
              ! the eigenvalue approximations will be refined when necessary as
              ! high relative accuracy is required for the computation of the
              ! corresponding eigenvectors.
              call la_scopy(im,w(wbegin),1,work(wbegin),1)
              ! we store in w the eigenvalue approximations w.r.t. the original
              ! matrix t.
              do i = 1,im
                 w(wbegin + i - 1) = w(wbegin + i - 1) + sigma
              end do
              ! ndepth is the current depth of the representation tree
              ndepth = 0
              ! parity is either 1 or 0
              parity = 1
              ! nclus is the number of clusters for the next level of the
              ! representation tree, we start with nclus = 1 for the root
              nclus = 1
              iwork(iindc1 + 1) = 1
              iwork(iindc1 + 2) = im
              ! idone is the number of eigenvectors already computed in the current
              ! block
              idone = 0
              ! loop while( idone<im )
              ! generate the representation tree for the current block and
              ! compute the eigenvectors
              40 continue
              if (idone < im) then
                 ! this is a crude protection against infinitely deep trees
                 if (ndepth > m) then
                    info = -2
                    return
                 end if
                 ! breadth first processing of the current level of the representation
                 ! tree: oldncl = number of clusters on current level
                 oldncl = nclus
                 ! reset nclus to count the number of child clusters
                 nclus = 0
                 parity = 1 - parity
                 if (parity == 0) then
                    oldcls = iindc1
                    newcls = iindc2
                 else
                    oldcls = iindc2
                    newcls = iindc1
                 end if
                 ! process the clusters on the current level
                 loop_150: do i = 1,oldncl
                    j = oldcls + 2*i
                    ! oldfst, oldlst = first, last index of current cluster.
                                     ! cluster indices start with 1 and are relative
                                     ! to wbegin when accessing w, wgap, werr, z
                    oldfst = iwork(j - 1)
                    oldlst = iwork(j)
                    if (ndepth > 0) then
                       ! retrieve relatively robust representation (rrr) of cluster
                       ! that has been computed at the previous level
                       ! the rrr is stored in z and overwritten once the eigenvectors
                       ! have been computed or when the cluster is refined
                       if ((dol == 1) .and. (dou == m)) then
                          ! get representation from location of the leftmost evalue
                          ! of the cluster
                          j = wbegin + oldfst - 1
                       else
                          if (wbegin + oldfst - 1 < dol) then
                             ! get representation from the left end of z array
                             j = dol - 1
                          elseif (wbegin + oldfst - 1 > dou) then
                             ! get representation from the right end of z array
                             j = dou
                          else
                             j = wbegin + oldfst - 1
                          end if
                       end if
                       call la_scopy(in,z(ibegin,j),1,d(ibegin),1)
                       call la_scopy(in - 1,z(ibegin,j + 1),1,l(ibegin),1)
                       sigma = z(iend,j + 1)
                       ! set the corresponding entries in z to zero
                       call la_slaset('FULL',in,2,zero,zero,z(ibegin,j),ldz)
                    end if
                    ! compute dl and dll of current rrr
                    do j = ibegin,iend - 1
                       tmp = d(j)*l(j)
                       work(indld - 1 + j) = tmp
                       work(indlld - 1 + j) = tmp*l(j)
                    end do
                    if (ndepth > 0) then
                       ! p and q are index of the first and last eigenvalue to compute
                       ! within the current block
                       p = indexw(wbegin - 1 + oldfst)
                       q = indexw(wbegin - 1 + oldlst)
                       ! offset for the arrays work, wgap and werr, i.e., the p-offset
                       ! through the q-offset elements of these arrays are to be used.
                        ! offset = p-oldfst
                       offset = indexw(wbegin) - 1
                       ! perform limited bisection (if necessary) to get approximate
                       ! eigenvalues to the precision needed.
                       call la_slarrb(in,d(ibegin),work(indlld + ibegin - 1),p,q,rtol1, &
                       rtol2,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk),iwork( &
                                  iindwk),pivmin,spdiam,in,iinfo)
                       if (iinfo /= 0) then
                          info = -1
                          return
                       end if
                       ! we also recompute the extremal gaps. w holds all eigenvalues
                       ! of the unshifted matrix and must be used for computation
                       ! of wgap, the entries of work might stem from rrrs with
                       ! different shifts. the gaps from wbegin-1+oldfst to
                       ! wbegin-1+oldlst are correctly computed in la_slarrb.
                       ! however, we only allow the gaps to become greater since
                       ! this is what should happen when we decrease werr
                       if (oldfst > 1) then
                          wgap(wbegin + oldfst - 2) = max(wgap(wbegin + oldfst - 2),w(wbegin + oldfst - 1) - &
                          werr(wbegin + oldfst - 1) - w(wbegin + oldfst - 2) - werr(wbegin + oldfst - 2))

                       end if
                       if (wbegin + oldlst - 1 < wend) then
                          wgap(wbegin + oldlst - 1) = max(wgap(wbegin + oldlst - 1),w(wbegin + oldlst) - &
                                    werr(wbegin + oldlst) - w(wbegin + oldlst - 1) - werr(wbegin + oldlst - 1))
                       end if
                       ! each time the eigenvalues in work get refined, we store
                       ! the newly found approximation with all shifts applied in w
                       do j = oldfst,oldlst
                          w(wbegin + j - 1) = work(wbegin + j - 1) + sigma
                       end do
                    end if
                    ! process the current node.
                    newfst = oldfst
                    loop_140: do j = oldfst,oldlst
                       if (j == oldlst) then
                          ! we are at the right end of the cluster, this is also the
                          ! boundary of the child cluster
                          newlst = j
                       else if (wgap(wbegin + j - 1) >= minrgp*abs(work(wbegin + j - 1))) &
                                 then
                          ! the right relative gap is big enough, the child cluster
                          ! (newfst,..,newlst) is well separated from the following
                          newlst = j
                        else
                          ! inside a child cluster, the relative gap is not
                          ! big enough.
                          cycle loop_140
                       end if
                       ! compute size of child cluster found
                       newsiz = newlst - newfst + 1
                       ! newftt is the place in z where the new rrr or the computed
                       ! eigenvector is to be stored
                       if ((dol == 1) .and. (dou == m)) then
                          ! store representation at location of the leftmost evalue
                          ! of the cluster
                          newftt = wbegin + newfst - 1
                       else
                          if (wbegin + newfst - 1 < dol) then
                             ! store representation at the left end of z array
                             newftt = dol - 1
                          elseif (wbegin + newfst - 1 > dou) then
                             ! store representation at the right end of z array
                             newftt = dou
                          else
                             newftt = wbegin + newfst - 1
                          end if
                       end if
                       if (newsiz > 1) then
                          ! current child is not a singleton but a cluster.
                          ! compute and store new representation of child.
                          ! compute left and right cluster gap.
                          ! lgap and rgap are not computed from work because
                          ! the eigenvalue approximations may stem from rrrs
                          ! different shifts. however, w hold all eigenvalues
                          ! of the unshifted matrix. still, the entries in wgap
                          ! have to be computed from work since the entries
                          ! in w might be of the same order so that gaps are not
                          ! exhibited correctly for very close eigenvalues.
                          if (newfst == 1) then
                             lgap = max(zero,w(wbegin) - werr(wbegin) - vl)
                         else
                             lgap = wgap(wbegin + newfst - 2)
                          end if
                          rgap = wgap(wbegin + newlst - 1)
                          ! compute left- and rightmost eigenvalue of child
                          ! to high precision in order to shift as close
                          ! as possible and obtain as large relative gaps
                          ! as possible
                          do k = 1,2
                             if (k == 1) then
                                p = indexw(wbegin - 1 + newfst)
                             else
                                p = indexw(wbegin - 1 + newlst)
                             end if
                             offset = indexw(wbegin) - 1
                             call la_slarrb(in,d(ibegin),work(indlld + ibegin - 1),p,p,rqtol, &
                             rqtol,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk), &
                                       iwork(iindwk),pivmin,spdiam,in,iinfo)
                          end do
                          if ((wbegin + newlst - 1 < dol) .or. (wbegin + newfst - 1 > dou)) then
                             ! if the cluster contains no desired eigenvalues
                             ! skip the computation of that branch of the rep. tree
                             ! we could skip before the refinement of the extremal
                             ! eigenvalues of the child, but then the representation
                             ! tree could be different from the one when nothing is
                             ! skipped. for this reason we skip at this place.
                             idone = idone + newlst - newfst + 1
                             goto 139
                          end if
                          ! compute rrr of child cluster.
                          ! note that the new rrr is stored in z
                          ! la_slarrf needs lwork = 2*n
                          call la_slarrf(in,d(ibegin),l(ibegin),work(indld + ibegin - 1), &
                          newfst,newlst,work(wbegin),wgap(wbegin),werr(wbegin),spdiam,lgap, &
                          rgap,pivmin,tau,z(ibegin,newftt),z(ibegin,newftt + 1),work(indwrk), &
                                    iinfo)
                          if (iinfo == 0) then
                             ! a new rrr for the cluster was found by la_slarrf
                             ! update shift and store it
                             ssigma = sigma + tau
                             z(iend,newftt + 1) = ssigma
                             ! work() are the midpoints and werr() the semi-width
                             ! note that the entries in w are unchanged.
                             do k = newfst,newlst
                                fudge = three*eps*abs(work(wbegin + k - 1))
                                work(wbegin + k - 1) = work(wbegin + k - 1) - tau
                                fudge = fudge + four*eps*abs(work(wbegin + k - 1))
                                ! fudge errors
                                werr(wbegin + k - 1) = werr(wbegin + k - 1) + fudge
                                ! gaps are not fudged. provided that werr is small
                                ! when eigenvalues are close, a zero gap indicates
                                ! that a new representation is needed for resolving
                                ! the cluster. a fudge could lead to a wrong decision
                                ! of judging eigenvalues 'separated' which in
                                ! reality are not. this could have a negative impact
                                ! on the orthogonality of the computed eigenvectors.
                             end do
                             nclus = nclus + 1
                             k = newcls + 2*nclus
                             iwork(k - 1) = newfst
                             iwork(k) = newlst
                          else
                             info = -2
                             return
                          end if
                       else
                          ! compute eigenvector of singleton
                          iter = 0
                          tol = four*log(real(in,KIND=sp))*eps
                          k = newfst
                          windex = wbegin + k - 1
                          windmn = max(windex - 1,1)
                          windpl = min(windex + 1,m)
                          lambda = work(windex)
                          done = done + 1
                          ! check if eigenvector computation is to be skipped
                          if ((windex < dol) .or. (windex > dou)) then
                             eskip = .true.
                             goto 125
                          else
                             eskip = .false.
                          end if
                          left = work(windex) - werr(windex)
                          right = work(windex) + werr(windex)
                          indeig = indexw(windex)
                          ! note that since we compute the eigenpairs for a child,
                          ! all eigenvalue approximations are w.r.t the same shift.
                          ! in this case, the entries in work should be used for
                          ! computing the gaps since they exhibit even very small
                          ! differences in the eigenvalues, as opposed to the
                          ! entries in w which might "look" the same.
                          if (k == 1) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vl, the formula
                             ! lgap = max( zero, (sigma - vl) + lambda )
                             ! can lead to an overestimation of the left gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small left gap.
                             lgap = eps*max(abs(left),abs(right))
                          else
                             lgap = wgap(windmn)
                          end if
                          if (k == im) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vu, the formula
                             ! can lead to an overestimation of the right gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small right gap.
                             rgap = eps*max(abs(left),abs(right))
                          else
                             rgap = wgap(windex)
                          end if
                          gap = min(lgap,rgap)
                          if ((k == 1) .or. (k == im)) then
                             ! the eigenvector support can become wrong
                             ! because significant entries could be cut off due to a
                             ! large gaptol parameter in lar1v. prevent this.
                             gaptol = zero
                          else
                             gaptol = gap*eps
                          end if
                          isupmn = in
                          isupmx = 1
                          ! update wgap so that it holds the minimum gap
                          ! to the left or the right. this is crucial in the
                          ! case where bisection is used to ensure that the
                          ! eigenvalue is refined up to the required precision.
                          ! the correct value is restored afterwards.
                          savgap = wgap(windex)
                          wgap(windex) = gap
                          ! we want to use the rayleigh quotient correction
                          ! as often as possible since it converges quadratically
                          ! when we are close enough to the desired eigenvalue.
                          ! however, the rayleigh quotient can have the wrong sign
                          ! and lead us away from the desired eigenvalue. in this
                          ! case, the best we can do is to use bisection.
                          usedbs = .false.
                          usedrq = .false.
                          ! bisection is initially turned off unless it is forced
                          needbs = .not. tryrqc
                          120 continue
                          ! check if bisection should be used to refine eigenvalue
                          if (needbs) then
                             ! take the bisection as new iterate
                             usedbs = .true.
                             itmp1 = iwork(iindr + windex)
                             offset = indexw(wbegin) - 1
                             call la_slarrb(in,d(ibegin),work(indlld + ibegin - 1),indeig, &
                             indeig,zero,two*eps,offset,work(wbegin),wgap(wbegin),werr(wbegin), &
                                       work(indwrk),iwork(iindwk),pivmin,spdiam,itmp1,iinfo)
                             if (iinfo /= 0) then
                                info = -3
                                return
                             end if
                             lambda = work(windex)
                             ! reset twist index from inaccurate lambda to
                             ! force computation of true mingma
                             iwork(iindr + windex) = 0
                          end if
                          ! given lambda, compute the eigenvector.
                          call la_slar1v(in,1,in,lambda,d(ibegin),l(ibegin),work( &
                          indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z(ibegin,windex &
                          ),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + windex),isuppz( &
                                    2*windex - 1),nrminv,resid,rqcorr,work(indwrk))
                          if (iter == 0) then
                             bstres = resid
                             bstw = lambda
                          elseif (resid < bstres) then
                             bstres = resid
                             bstw = lambda
                          end if
                          isupmn = min(isupmn,isuppz(2*windex - 1))
                          isupmx = max(isupmx,isuppz(2*windex))
                          iter = iter + 1
                          ! sin alpha <= |resid|/gap
                          ! note that both the residual and the gap are
                          ! proportional to the matrix, so ||t|| doesn't play
                          ! a role in the quotient
                          ! convergence test for rayleigh-quotient iteration
                          ! (omitted when bisection has been used)
                          if (resid > tol*gap .and. abs(rqcorr) > rqtol*abs(lambda) .and. .not. &
                                    usedbs) then
                             ! we need to check that the rqcorr update doesn't
                             ! move the eigenvalue away from the desired one and
                             ! towards a neighbor. -> protection with bisection
                             if (indeig <= negcnt) then
                                ! the wanted eigenvalue lies to the left
                                sgndef = -one
                             else
                                ! the wanted eigenvalue lies to the right
                                sgndef = one
                             end if
                             ! we only use the rqcorr if it improves the
                             ! the iterate reasonably.
                             if ((rqcorr*sgndef >= zero) .and. (lambda + rqcorr <= right) .and. ( &
                                       lambda + rqcorr >= left)) then
                                usedrq = .true.
                                ! store new midpoint of bisection interval in work
                                if (sgndef == one) then
                                   ! the current lambda is on the left of the true
                                   ! eigenvalue
                                   left = lambda
                                   ! we prefer to assume that the error estimate
                                   ! is correct. we could make the interval not
                                   ! as a bracket but to be modified if the rqcorr
                                   ! chooses to. in this case, the right side should
                                   ! be modified as follows:
                                    ! right = max(right, lambda + rqcorr)
                                else
                                   ! the current lambda is on the right of the true
                                   ! eigenvalue
                                   right = lambda
                                   ! see comment about assuming the error estimate is
                                   ! correct above.
                                    ! left = min(left, lambda + rqcorr)
                                end if
                                work(windex) = half*(right + left)
                                ! take rqcorr since it has the correct sign and
                                ! improves the iterate reasonably
                                lambda = lambda + rqcorr
                                ! update width of error interval
                                werr(windex) = half*(right - left)
                             else
                                needbs = .true.
                             end if
                             if (right - left < rqtol*abs(lambda)) then
                                   ! the eigenvalue is computed to bisection accuracy
                                   ! compute eigenvector and stop
                                usedbs = .true.
                                goto 120
                             elseif (iter < maxitr) then
                                goto 120
                             elseif (iter == maxitr) then
                                needbs = .true.
                                goto 120
                             else
                                info = 5
                                return
                             end if
                          else
                             stp2ii = .false.
             if (usedrq .and. usedbs .and. bstres <= resid) then
                                lambda = bstw
                                stp2ii = .true.
                             end if
                             if (stp2ii) then
                                ! improve error angle by second step
                                call la_slar1v(in,1,in,lambda,d(ibegin),l(ibegin), &
                                work(indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z( &
                                ibegin,windex),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + &
                                windex),isuppz(2*windex - 1),nrminv,resid,rqcorr,work(indwrk &
                                          ))
                             end if
                             work(windex) = lambda
                          end if
                          ! compute fp-vector support w.r.t. whole matrix
                          isuppz(2*windex - 1) = isuppz(2*windex - 1) + oldien
                          isuppz(2*windex) = isuppz(2*windex) + oldien
                          zfrom = isuppz(2*windex - 1)
                          zto = isuppz(2*windex)
                          isupmn = isupmn + oldien
                          isupmx = isupmx + oldien
                          ! ensure vector is ok if support in the rqi has changed
                          if (isupmn < zfrom) then
                             do ii = isupmn,zfrom - 1
                                z(ii,windex) = zero
                             end do
                          end if
                          if (isupmx > zto) then
                             do ii = zto + 1,isupmx
                                z(ii,windex) = zero
                             end do
                          end if
                          call la_sscal(zto - zfrom + 1,nrminv,z(zfrom,windex),1)
                          125 continue
                          ! update w
                          w(windex) = lambda + sigma
                          ! recompute the gaps on the left and right
                          ! but only allow them to become larger and not
                          ! smaller (which can only happen through "bad"
                          ! cancellation and doesn't reflect the theory
                          ! where the initial gaps are underestimated due
                          ! to werr being too crude.)
                          if (.not. eskip) then
                             if (k > 1) then
                                wgap(windmn) = max(wgap(windmn),w(windex) - werr(windex) - w( &
                                          windmn) - werr(windmn))
                             end if
                             if (windex < wend) then
                                wgap(windex) = max(savgap,w(windpl) - werr(windpl) - w( &
                                          windex) - werr(windex))
                             end if
                          end if
                          idone = idone + 1
                       end if
                       ! here ends the code for the current child
                       139 continue
                       ! proceed to any remaining child nodes
                       newfst = j + 1
                    end do loop_140
                 end do loop_150
                 ndepth = ndepth + 1
                 go to 40
              end if
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_slarrv
     !> DLARRV: computes the eigenvectors of the tridiagonal matrix
     !> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
     !> The input eigenvalues should have been computed by DLARRE.

     pure subroutine la_dlarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1, &
               rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz,isuppz,work,iwork,info)
        use la_constants_dp,only:zero,half,one,two,three,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dol,dou,ldz,m,n
           integer(ilp),intent(out) :: info
           real(dp),intent(in) :: minrgp,pivmin,vl,vu
           real(dp),intent(inout) :: rtol1,rtol2
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),indexw(*),isplit(*)
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),l(*),w(*),werr(*),wgap(*)
           real(dp),intent(in) :: gers(*)
           real(dp),intent(out) :: work(*)
           real(dp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 10

           ! Local Scalars
           logical(lk) :: eskip,needbs,stp2ii,tryrqc,usedbs,usedrq
           integer(ilp) :: done,i,ibegin,idone,iend,ii,iindc1,iindc2,iindr,iindwk,iinfo, &
            im,in,indeig,indld,indlld,indwrk,isupmn,isupmx,iter,itmp1,j,jblk,k, &
            miniwsize,minwsize,nclus,ndepth,negcnt,newcls,newfst,newftt,newlst,newsiz, &
            offset,oldcls,oldfst,oldien,oldlst,oldncl,p,parity,q,wbegin,wend,windex, &
                      windmn,windpl,zfrom,zto,zusedl,zusedu,zusedw
           real(dp) :: bstres,bstw,eps,fudge,gap,gaptol,gl,gu,lambda,left,lgap,mingma, &
           nrminv,resid,rgap,right,rqcorr,rqtol,savgap,sgndef,sigma,spdiam,ssigma,tau, &
                     tmp,tol,ztz
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if ((n <= 0) .or. (m <= 0)) then
              return
           end if
           ! the first n entries of work are reserved for the eigenvalues
           indld = n + 1
           indlld = 2*n + 1
           indwrk = 3*n + 1
           minwsize = 12*n
           do i = 1,minwsize
              work(i) = zero
           end do
           ! iwork(iindr+1:iindr+n) hold the twist indices r for the
           ! factorization used to compute the fp vector
           iindr = 0
           ! iwork(iindc1+1:iinc2+n) are used to store the clusters of the current
           ! layer and the one above.
           iindc1 = n
           iindc2 = 2*n
           iindwk = 3*n + 1
           miniwsize = 7*n
           do i = 1,miniwsize
              iwork(i) = 0
           end do
           zusedl = 1
           if (dol > 1) then
              ! set lower bound for use of z
              zusedl = dol - 1
           end if
           zusedu = m
           if (dou < m) then
              ! set lower bound for use of z
              zusedu = dou + 1
           end if
           ! the width of the part of z that is used
           zusedw = zusedu - zusedl + 1
           call la_dlaset('FULL',n,zusedw,zero,zero,z(1,zusedl),ldz)
           eps = la_dlamch('PRECISION')
           rqtol = two*eps
           ! set expert flags for standard code.
           tryrqc = .true.
           if ((dol == 1) .and. (dou == m)) then
           else
              ! only selected eigenpairs are computed. since the other evalues
              ! are not refined by rq iteration, bisection has to compute to full
              ! accuracy.
              rtol1 = four*eps
              rtol2 = four*eps
           end if
           ! the entries wbegin:wend in w, werr, wgap correspond to the
           ! desired eigenvalues. the support of the nonzero eigenvector
           ! entries is contained in the interval ibegin:iend.
           ! remark that if k eigenpairs are desired, then the eigenvectors
           ! are stored in k contiguous columns of z.
           ! done is the number of eigenvectors already computed
           done = 0
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,iblock(m)
              iend = isplit(jblk)
              sigma = l(iend)
              ! find the eigenvectors of the submatrix indexed ibegin
              ! through iend.
              wend = wbegin - 1
              15 continue
              if (wend < m) then
                 if (iblock(wend + 1) == jblk) then
                    wend = wend + 1
                    go to 15
                 end if
              end if
              if (wend < wbegin) then
                 ibegin = iend + 1
                 cycle loop_170
              elseif ((wend < dol) .or. (wbegin > dou)) then
                 ibegin = iend + 1
                 wbegin = wend + 1
                 cycle loop_170
              end if
              ! find local spectral diameter of the block
              gl = gers(2*ibegin - 1)
              gu = gers(2*ibegin)
              do i = ibegin + 1,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              ! oldien is the last index of the previous block
              oldien = ibegin - 1
              ! calculate the size of the current block
              in = iend - ibegin + 1
              ! the number of eigenvalues in the current block
              im = wend - wbegin + 1
              ! this is for a 1x1 block
              if (ibegin == iend) then
                 done = done + 1
                 z(ibegin,wbegin) = one
                 isuppz(2*wbegin - 1) = ibegin
                 isuppz(2*wbegin) = ibegin
                 w(wbegin) = w(wbegin) + sigma
                 work(wbegin) = w(wbegin)
                 ibegin = iend + 1
                 wbegin = wbegin + 1
                 cycle loop_170
              end if
              ! the desired (shifted) eigenvalues are stored in w(wbegin:wend)
              ! note that these can be approximations, in this case, the corresp.
              ! entries of werr give the size of the uncertainty interval.
              ! the eigenvalue approximations will be refined when necessary as
              ! high relative accuracy is required for the computation of the
              ! corresponding eigenvectors.
              call la_dcopy(im,w(wbegin),1,work(wbegin),1)
              ! we store in w the eigenvalue approximations w.r.t. the original
              ! matrix t.
              do i = 1,im
                 w(wbegin + i - 1) = w(wbegin + i - 1) + sigma
              end do
              ! ndepth is the current depth of the representation tree
              ndepth = 0
              ! parity is either 1 or 0
              parity = 1
              ! nclus is the number of clusters for the next level of the
              ! representation tree, we start with nclus = 1 for the root
              nclus = 1
              iwork(iindc1 + 1) = 1
              iwork(iindc1 + 2) = im
              ! idone is the number of eigenvectors already computed in the current
              ! block
              idone = 0
              ! loop while( idone<im )
              ! generate the representation tree for the current block and
              ! compute the eigenvectors
              40 continue
              if (idone < im) then
                 ! this is a crude protection against infinitely deep trees
                 if (ndepth > m) then
                    info = -2
                    return
                 end if
                 ! breadth first processing of the current level of the representation
                 ! tree: oldncl = number of clusters on current level
                 oldncl = nclus
                 ! reset nclus to count the number of child clusters
                 nclus = 0
                 parity = 1 - parity
                 if (parity == 0) then
                    oldcls = iindc1
                    newcls = iindc2
                 else
                    oldcls = iindc2
                    newcls = iindc1
                 end if
                 ! process the clusters on the current level
                 loop_150: do i = 1,oldncl
                    j = oldcls + 2*i
                    ! oldfst, oldlst = first, last index of current cluster.
                                     ! cluster indices start with 1 and are relative
                                     ! to wbegin when accessing w, wgap, werr, z
                    oldfst = iwork(j - 1)
                    oldlst = iwork(j)
                    if (ndepth > 0) then
                       ! retrieve relatively robust representation (rrr) of cluster
                       ! that has been computed at the previous level
                       ! the rrr is stored in z and overwritten once the eigenvectors
                       ! have been computed or when the cluster is refined
                       if ((dol == 1) .and. (dou == m)) then
                          ! get representation from location of the leftmost evalue
                          ! of the cluster
                          j = wbegin + oldfst - 1
                       else
                          if (wbegin + oldfst - 1 < dol) then
                             ! get representation from the left end of z array
                             j = dol - 1
                          elseif (wbegin + oldfst - 1 > dou) then
                             ! get representation from the right end of z array
                             j = dou
                          else
                             j = wbegin + oldfst - 1
                          end if
                       end if
                       call la_dcopy(in,z(ibegin,j),1,d(ibegin),1)
                       call la_dcopy(in - 1,z(ibegin,j + 1),1,l(ibegin),1)
                       sigma = z(iend,j + 1)
                       ! set the corresponding entries in z to zero
                       call la_dlaset('FULL',in,2,zero,zero,z(ibegin,j),ldz)
                    end if
                    ! compute dl and dll of current rrr
                    do j = ibegin,iend - 1
                       tmp = d(j)*l(j)
                       work(indld - 1 + j) = tmp
                       work(indlld - 1 + j) = tmp*l(j)
                    end do
                    if (ndepth > 0) then
                       ! p and q are index of the first and last eigenvalue to compute
                       ! within the current block
                       p = indexw(wbegin - 1 + oldfst)
                       q = indexw(wbegin - 1 + oldlst)
                       ! offset for the arrays work, wgap and werr, i.e., the p-offset
                       ! through the q-offset elements of these arrays are to be used.
                        ! offset = p-oldfst
                       offset = indexw(wbegin) - 1
                       ! perform limited bisection (if necessary) to get approximate
                       ! eigenvalues to the precision needed.
                       call la_dlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,q,rtol1, &
                       rtol2,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk),iwork( &
                                  iindwk),pivmin,spdiam,in,iinfo)
                       if (iinfo /= 0) then
                          info = -1
                          return
                       end if
                       ! we also recompute the extremal gaps. w holds all eigenvalues
                       ! of the unshifted matrix and must be used for computation
                       ! of wgap, the entries of work might stem from rrrs with
                       ! different shifts. the gaps from wbegin-1+oldfst to
                       ! wbegin-1+oldlst are correctly computed in la_dlarrb.
                       ! however, we only allow the gaps to become greater since
                       ! this is what should happen when we decrease werr
                       if (oldfst > 1) then
                          wgap(wbegin + oldfst - 2) = max(wgap(wbegin + oldfst - 2),w(wbegin + oldfst - 1) - &
                          werr(wbegin + oldfst - 1) - w(wbegin + oldfst - 2) - werr(wbegin + oldfst - 2))

                       end if
                       if (wbegin + oldlst - 1 < wend) then
                          wgap(wbegin + oldlst - 1) = max(wgap(wbegin + oldlst - 1),w(wbegin + oldlst) - &
                                    werr(wbegin + oldlst) - w(wbegin + oldlst - 1) - werr(wbegin + oldlst - 1))
                       end if
                       ! each time the eigenvalues in work get refined, we store
                       ! the newly found approximation with all shifts applied in w
                       do j = oldfst,oldlst
                          w(wbegin + j - 1) = work(wbegin + j - 1) + sigma
                       end do
                    end if
                    ! process the current node.
                    newfst = oldfst
                    loop_140: do j = oldfst,oldlst
                       if (j == oldlst) then
                          ! we are at the right end of the cluster, this is also the
                          ! boundary of the child cluster
                          newlst = j
                       else if (wgap(wbegin + j - 1) >= minrgp*abs(work(wbegin + j - 1))) &
                                 then
                          ! the right relative gap is big enough, the child cluster
                          ! (newfst,..,newlst) is well separated from the following
                          newlst = j
                        else
                          ! inside a child cluster, the relative gap is not
                          ! big enough.
                          cycle loop_140
                       end if
                       ! compute size of child cluster found
                       newsiz = newlst - newfst + 1
                       ! newftt is the place in z where the new rrr or the computed
                       ! eigenvector is to be stored
                       if ((dol == 1) .and. (dou == m)) then
                          ! store representation at location of the leftmost evalue
                          ! of the cluster
                          newftt = wbegin + newfst - 1
                       else
                          if (wbegin + newfst - 1 < dol) then
                             ! store representation at the left end of z array
                             newftt = dol - 1
                          elseif (wbegin + newfst - 1 > dou) then
                             ! store representation at the right end of z array
                             newftt = dou
                          else
                             newftt = wbegin + newfst - 1
                          end if
                       end if
                       if (newsiz > 1) then
                          ! current child is not a singleton but a cluster.
                          ! compute and store new representation of child.
                          ! compute left and right cluster gap.
                          ! lgap and rgap are not computed from work because
                          ! the eigenvalue approximations may stem from rrrs
                          ! different shifts. however, w hold all eigenvalues
                          ! of the unshifted matrix. still, the entries in wgap
                          ! have to be computed from work since the entries
                          ! in w might be of the same order so that gaps are not
                          ! exhibited correctly for very close eigenvalues.
                          if (newfst == 1) then
                             lgap = max(zero,w(wbegin) - werr(wbegin) - vl)
                         else
                             lgap = wgap(wbegin + newfst - 2)
                          end if
                          rgap = wgap(wbegin + newlst - 1)
                          ! compute left- and rightmost eigenvalue of child
                          ! to high precision in order to shift as close
                          ! as possible and obtain as large relative gaps
                          ! as possible
                          do k = 1,2
                             if (k == 1) then
                                p = indexw(wbegin - 1 + newfst)
                             else
                                p = indexw(wbegin - 1 + newlst)
                             end if
                             offset = indexw(wbegin) - 1
                             call la_dlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,p,rqtol, &
                             rqtol,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk), &
                                       iwork(iindwk),pivmin,spdiam,in,iinfo)
                          end do
                          if ((wbegin + newlst - 1 < dol) .or. (wbegin + newfst - 1 > dou)) then
                             ! if the cluster contains no desired eigenvalues
                             ! skip the computation of that branch of the rep. tree
                             ! we could skip before the refinement of the extremal
                             ! eigenvalues of the child, but then the representation
                             ! tree could be different from the one when nothing is
                             ! skipped. for this reason we skip at this place.
                             idone = idone + newlst - newfst + 1
                             goto 139
                          end if
                          ! compute rrr of child cluster.
                          ! note that the new rrr is stored in z
                          ! la_dlarrf needs lwork = 2*n
                          call la_dlarrf(in,d(ibegin),l(ibegin),work(indld + ibegin - 1), &
                          newfst,newlst,work(wbegin),wgap(wbegin),werr(wbegin),spdiam,lgap, &
                          rgap,pivmin,tau,z(ibegin,newftt),z(ibegin,newftt + 1),work(indwrk), &
                                    iinfo)
                          if (iinfo == 0) then
                             ! a new rrr for the cluster was found by la_dlarrf
                             ! update shift and store it
                             ssigma = sigma + tau
                             z(iend,newftt + 1) = ssigma
                             ! work() are the midpoints and werr() the semi-width
                             ! note that the entries in w are unchanged.
                             do k = newfst,newlst
                                fudge = three*eps*abs(work(wbegin + k - 1))
                                work(wbegin + k - 1) = work(wbegin + k - 1) - tau
                                fudge = fudge + four*eps*abs(work(wbegin + k - 1))
                                ! fudge errors
                                werr(wbegin + k - 1) = werr(wbegin + k - 1) + fudge
                                ! gaps are not fudged. provided that werr is small
                                ! when eigenvalues are close, a zero gap indicates
                                ! that a new representation is needed for resolving
                                ! the cluster. a fudge could lead to a wrong decision
                                ! of judging eigenvalues 'separated' which in
                                ! reality are not. this could have a negative impact
                                ! on the orthogonality of the computed eigenvectors.
                             end do
                             nclus = nclus + 1
                             k = newcls + 2*nclus
                             iwork(k - 1) = newfst
                             iwork(k) = newlst
                          else
                             info = -2
                             return
                          end if
                       else
                          ! compute eigenvector of singleton
                          iter = 0
                          tol = four*log(real(in,KIND=dp))*eps
                          k = newfst
                          windex = wbegin + k - 1
                          windmn = max(windex - 1,1)
                          windpl = min(windex + 1,m)
                          lambda = work(windex)
                          done = done + 1
                          ! check if eigenvector computation is to be skipped
                          if ((windex < dol) .or. (windex > dou)) then
                             eskip = .true.
                             goto 125
                          else
                             eskip = .false.
                          end if
                          left = work(windex) - werr(windex)
                          right = work(windex) + werr(windex)
                          indeig = indexw(windex)
                          ! note that since we compute the eigenpairs for a child,
                          ! all eigenvalue approximations are w.r.t the same shift.
                          ! in this case, the entries in work should be used for
                          ! computing the gaps since they exhibit even very small
                          ! differences in the eigenvalues, as opposed to the
                          ! entries in w which might "look" the same.
                          if (k == 1) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vl, the formula
                             ! lgap = max( zero, (sigma - vl) + lambda )
                             ! can lead to an overestimation of the left gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small left gap.
                             lgap = eps*max(abs(left),abs(right))
                          else
                             lgap = wgap(windmn)
                          end if
                          if (k == im) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vu, the formula
                             ! can lead to an overestimation of the right gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small right gap.
                             rgap = eps*max(abs(left),abs(right))
                          else
                             rgap = wgap(windex)
                          end if
                          gap = min(lgap,rgap)
                          if ((k == 1) .or. (k == im)) then
                             ! the eigenvector support can become wrong
                             ! because significant entries could be cut off due to a
                             ! large gaptol parameter in lar1v. prevent this.
                             gaptol = zero
                          else
                             gaptol = gap*eps
                          end if
                          isupmn = in
                          isupmx = 1
                          ! update wgap so that it holds the minimum gap
                          ! to the left or the right. this is crucial in the
                          ! case where bisection is used to ensure that the
                          ! eigenvalue is refined up to the required precision.
                          ! the correct value is restored afterwards.
                          savgap = wgap(windex)
                          wgap(windex) = gap
                          ! we want to use the rayleigh quotient correction
                          ! as often as possible since it converges quadratically
                          ! when we are close enough to the desired eigenvalue.
                          ! however, the rayleigh quotient can have the wrong sign
                          ! and lead us away from the desired eigenvalue. in this
                          ! case, the best we can do is to use bisection.
                          usedbs = .false.
                          usedrq = .false.
                          ! bisection is initially turned off unless it is forced
                          needbs = .not. tryrqc
                          120 continue
                          ! check if bisection should be used to refine eigenvalue
                          if (needbs) then
                             ! take the bisection as new iterate
                             usedbs = .true.
                             itmp1 = iwork(iindr + windex)
                             offset = indexw(wbegin) - 1
                             call la_dlarrb(in,d(ibegin),work(indlld + ibegin - 1),indeig, &
                             indeig,zero,two*eps,offset,work(wbegin),wgap(wbegin),werr(wbegin), &
                                       work(indwrk),iwork(iindwk),pivmin,spdiam,itmp1,iinfo)
                             if (iinfo /= 0) then
                                info = -3
                                return
                             end if
                             lambda = work(windex)
                             ! reset twist index from inaccurate lambda to
                             ! force computation of true mingma
                             iwork(iindr + windex) = 0
                          end if
                          ! given lambda, compute the eigenvector.
                          call la_dlar1v(in,1,in,lambda,d(ibegin),l(ibegin),work( &
                          indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z(ibegin,windex &
                          ),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + windex),isuppz( &
                                    2*windex - 1),nrminv,resid,rqcorr,work(indwrk))
                          if (iter == 0) then
                             bstres = resid
                             bstw = lambda
                          elseif (resid < bstres) then
                             bstres = resid
                             bstw = lambda
                          end if
                          isupmn = min(isupmn,isuppz(2*windex - 1))
                          isupmx = max(isupmx,isuppz(2*windex))
                          iter = iter + 1
                          ! sin alpha <= |resid|/gap
                          ! note that both the residual and the gap are
                          ! proportional to the matrix, so ||t|| doesn't play
                          ! a role in the quotient
                          ! convergence test for rayleigh-quotient iteration
                          ! (omitted when bisection has been used)
                          if (resid > tol*gap .and. abs(rqcorr) > rqtol*abs(lambda) .and. .not. &
                                    usedbs) then
                             ! we need to check that the rqcorr update doesn't
                             ! move the eigenvalue away from the desired one and
                             ! towards a neighbor. -> protection with bisection
                             if (indeig <= negcnt) then
                                ! the wanted eigenvalue lies to the left
                                sgndef = -one
                             else
                                ! the wanted eigenvalue lies to the right
                                sgndef = one
                             end if
                             ! we only use the rqcorr if it improves the
                             ! the iterate reasonably.
                             if ((rqcorr*sgndef >= zero) .and. (lambda + rqcorr <= right) .and. ( &
                                       lambda + rqcorr >= left)) then
                                usedrq = .true.
                                ! store new midpoint of bisection interval in work
                                if (sgndef == one) then
                                   ! the current lambda is on the left of the true
                                   ! eigenvalue
                                   left = lambda
                                   ! we prefer to assume that the error estimate
                                   ! is correct. we could make the interval not
                                   ! as a bracket but to be modified if the rqcorr
                                   ! chooses to. in this case, the right side should
                                   ! be modified as follows:
                                    ! right = max(right, lambda + rqcorr)
                                else
                                   ! the current lambda is on the right of the true
                                   ! eigenvalue
                                   right = lambda
                                   ! see comment about assuming the error estimate is
                                   ! correct above.
                                    ! left = min(left, lambda + rqcorr)
                                end if
                                work(windex) = half*(right + left)
                                ! take rqcorr since it has the correct sign and
                                ! improves the iterate reasonably
                                lambda = lambda + rqcorr
                                ! update width of error interval
                                werr(windex) = half*(right - left)
                             else
                                needbs = .true.
                             end if
                             if (right - left < rqtol*abs(lambda)) then
                                   ! the eigenvalue is computed to bisection accuracy
                                   ! compute eigenvector and stop
                                usedbs = .true.
                                goto 120
                             elseif (iter < maxitr) then
                                goto 120
                             elseif (iter == maxitr) then
                                needbs = .true.
                                goto 120
                             else
                                info = 5
                                return
                             end if
                          else
                             stp2ii = .false.
             if (usedrq .and. usedbs .and. bstres <= resid) then
                                lambda = bstw
                                stp2ii = .true.
                             end if
                             if (stp2ii) then
                                ! improve error angle by second step
                                call la_dlar1v(in,1,in,lambda,d(ibegin),l(ibegin), &
                                work(indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z( &
                                ibegin,windex),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + &
                                windex),isuppz(2*windex - 1),nrminv,resid,rqcorr,work(indwrk &
                                          ))
                             end if
                             work(windex) = lambda
                          end if
                          ! compute fp-vector support w.r.t. whole matrix
                          isuppz(2*windex - 1) = isuppz(2*windex - 1) + oldien
                          isuppz(2*windex) = isuppz(2*windex) + oldien
                          zfrom = isuppz(2*windex - 1)
                          zto = isuppz(2*windex)
                          isupmn = isupmn + oldien
                          isupmx = isupmx + oldien
                          ! ensure vector is ok if support in the rqi has changed
                          if (isupmn < zfrom) then
                             do ii = isupmn,zfrom - 1
                                z(ii,windex) = zero
                             end do
                          end if
                          if (isupmx > zto) then
                             do ii = zto + 1,isupmx
                                z(ii,windex) = zero
                             end do
                          end if
                          call la_dscal(zto - zfrom + 1,nrminv,z(zfrom,windex),1)
                          125 continue
                          ! update w
                          w(windex) = lambda + sigma
                          ! recompute the gaps on the left and right
                          ! but only allow them to become larger and not
                          ! smaller (which can only happen through "bad"
                          ! cancellation and doesn't reflect the theory
                          ! where the initial gaps are underestimated due
                          ! to werr being too crude.)
                          if (.not. eskip) then
                             if (k > 1) then
                                wgap(windmn) = max(wgap(windmn),w(windex) - werr(windex) - w( &
                                          windmn) - werr(windmn))
                             end if
                             if (windex < wend) then
                                wgap(windex) = max(savgap,w(windpl) - werr(windpl) - w( &
                                          windex) - werr(windex))
                             end if
                          end if
                          idone = idone + 1
                       end if
                       ! here ends the code for the current child
                       139 continue
                       ! proceed to any remaining child nodes
                       newfst = j + 1
                    end do loop_140
                 end do loop_150
                 ndepth = ndepth + 1
                 go to 40
              end if
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_dlarrv
     !> QLARRV: computes the eigenvectors of the tridiagonal matrix
     !> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
     !> The input eigenvalues should have been computed by QLARRE.

     pure subroutine la_qlarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1, &
               rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz,isuppz,work,iwork,info)
        use la_constants_qp,only:zero,half,one,two,three,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dol,dou,ldz,m,n
           integer(ilp),intent(out) :: info
           real(qp),intent(in) :: minrgp,pivmin,vl,vu
           real(qp),intent(inout) :: rtol1,rtol2
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),indexw(*),isplit(*)
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),l(*),w(*),werr(*),wgap(*)
           real(qp),intent(in) :: gers(*)
           real(qp),intent(out) :: work(*)
           real(qp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 10

           ! Local Scalars
           logical(lk) :: eskip,needbs,stp2ii,tryrqc,usedbs,usedrq
           integer(ilp) :: done,i,ibegin,idone,iend,ii,iindc1,iindc2,iindr,iindwk,iinfo, &
            im,in,indeig,indld,indlld,indwrk,isupmn,isupmx,iter,itmp1,j,jblk,k, &
            miniwsize,minwsize,nclus,ndepth,negcnt,newcls,newfst,newftt,newlst,newsiz, &
            offset,oldcls,oldfst,oldien,oldlst,oldncl,p,parity,q,wbegin,wend,windex, &
                      windmn,windpl,zfrom,zto,zusedl,zusedu,zusedw
           real(qp) :: bstres,bstw,eps,fudge,gap,gaptol,gl,gu,lambda,left,lgap,mingma, &
           nrminv,resid,rgap,right,rqcorr,rqtol,savgap,sgndef,sigma,spdiam,ssigma,tau, &
                     tmp,tol,ztz
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if ((n <= 0) .or. (m <= 0)) then
              return
           end if
           ! the first n entries of work are reserved for the eigenvalues
           indld = n + 1
           indlld = 2*n + 1
           indwrk = 3*n + 1
           minwsize = 12*n
           do i = 1,minwsize
              work(i) = zero
           end do
           ! iwork(iindr+1:iindr+n) hold the twist indices r for the
           ! factorization used to compute the fp vector
           iindr = 0
           ! iwork(iindc1+1:iinc2+n) are used to store the clusters of the current
           ! layer and the one above.
           iindc1 = n
           iindc2 = 2*n
           iindwk = 3*n + 1
           miniwsize = 7*n
           do i = 1,miniwsize
              iwork(i) = 0
           end do
           zusedl = 1
           if (dol > 1) then
              ! set lower bound for use of z
              zusedl = dol - 1
           end if
           zusedu = m
           if (dou < m) then
              ! set lower bound for use of z
              zusedu = dou + 1
           end if
           ! the width of the part of z that is used
           zusedw = zusedu - zusedl + 1
           call la_qlaset('FULL',n,zusedw,zero,zero,z(1,zusedl),ldz)
           eps = la_qlamch('PRECISION')
           rqtol = two*eps
           ! set expert flags for standard code.
           tryrqc = .true.
           if ((dol == 1) .and. (dou == m)) then
           else
              ! only selected eigenpairs are computed. since the other evalues
              ! are not refined by rq iteration, bisection has to compute to full
              ! accuracy.
              rtol1 = four*eps
              rtol2 = four*eps
           end if
           ! the entries wbegin:wend in w, werr, wgap correspond to the
           ! desired eigenvalues. the support of the nonzero eigenvector
           ! entries is contained in the interval ibegin:iend.
           ! remark that if k eigenpairs are desired, then the eigenvectors
           ! are stored in k contiguous columns of z.
           ! done is the number of eigenvectors already computed
           done = 0
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,iblock(m)
              iend = isplit(jblk)
              sigma = l(iend)
              ! find the eigenvectors of the submatrix indexed ibegin
              ! through iend.
              wend = wbegin - 1
              15 continue
              if (wend < m) then
                 if (iblock(wend + 1) == jblk) then
                    wend = wend + 1
                    go to 15
                 end if
              end if
              if (wend < wbegin) then
                 ibegin = iend + 1
                 cycle loop_170
              elseif ((wend < dol) .or. (wbegin > dou)) then
                 ibegin = iend + 1
                 wbegin = wend + 1
                 cycle loop_170
              end if
              ! find local spectral diameter of the block
              gl = gers(2*ibegin - 1)
              gu = gers(2*ibegin)
              do i = ibegin + 1,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              ! oldien is the last index of the previous block
              oldien = ibegin - 1
              ! calculate the size of the current block
              in = iend - ibegin + 1
              ! the number of eigenvalues in the current block
              im = wend - wbegin + 1
              ! this is for a 1x1 block
              if (ibegin == iend) then
                 done = done + 1
                 z(ibegin,wbegin) = one
                 isuppz(2*wbegin - 1) = ibegin
                 isuppz(2*wbegin) = ibegin
                 w(wbegin) = w(wbegin) + sigma
                 work(wbegin) = w(wbegin)
                 ibegin = iend + 1
                 wbegin = wbegin + 1
                 cycle loop_170
              end if
              ! the desired (shifted) eigenvalues are stored in w(wbegin:wend)
              ! note that these can be approximations, in this case, the corresp.
              ! entries of werr give the size of the uncertainty interval.
              ! the eigenvalue approximations will be refined when necessary as
              ! high relative accuracy is required for the computation of the
              ! corresponding eigenvectors.
              call la_qcopy(im,w(wbegin),1,work(wbegin),1)
              ! we store in w the eigenvalue approximations w.r.t. the original
              ! matrix t.
              do i = 1,im
                 w(wbegin + i - 1) = w(wbegin + i - 1) + sigma
              end do
              ! ndepth is the current depth of the representation tree
              ndepth = 0
              ! parity is either 1 or 0
              parity = 1
              ! nclus is the number of clusters for the next level of the
              ! representation tree, we start with nclus = 1 for the root
              nclus = 1
              iwork(iindc1 + 1) = 1
              iwork(iindc1 + 2) = im
              ! idone is the number of eigenvectors already computed in the current
              ! block
              idone = 0
              ! loop while( idone<im )
              ! generate the representation tree for the current block and
              ! compute the eigenvectors
              40 continue
              if (idone < im) then
                 ! this is a crude protection against infinitely deep trees
                 if (ndepth > m) then
                    info = -2
                    return
                 end if
                 ! breadth first processing of the current level of the representation
                 ! tree: oldncl = number of clusters on current level
                 oldncl = nclus
                 ! reset nclus to count the number of child clusters
                 nclus = 0
                 parity = 1 - parity
                 if (parity == 0) then
                    oldcls = iindc1
                    newcls = iindc2
                 else
                    oldcls = iindc2
                    newcls = iindc1
                 end if
                 ! process the clusters on the current level
                 loop_150: do i = 1,oldncl
                    j = oldcls + 2*i
                    ! oldfst, oldlst = first, last index of current cluster.
                                     ! cluster indices start with 1 and are relative
                                     ! to wbegin when accessing w, wgap, werr, z
                    oldfst = iwork(j - 1)
                    oldlst = iwork(j)
                    if (ndepth > 0) then
                       ! retrieve relatively robust representation (rrr) of cluster
                       ! that has been computed at the previous level
                       ! the rrr is stored in z and overwritten once the eigenvectors
                       ! have been computed or when the cluster is refined
                       if ((dol == 1) .and. (dou == m)) then
                          ! get representation from location of the leftmost evalue
                          ! of the cluster
                          j = wbegin + oldfst - 1
                       else
                          if (wbegin + oldfst - 1 < dol) then
                             ! get representation from the left end of z array
                             j = dol - 1
                          elseif (wbegin + oldfst - 1 > dou) then
                             ! get representation from the right end of z array
                             j = dou
                          else
                             j = wbegin + oldfst - 1
                          end if
                       end if
                       call la_qcopy(in,z(ibegin,j),1,d(ibegin),1)
                       call la_qcopy(in - 1,z(ibegin,j + 1),1,l(ibegin),1)
                       sigma = z(iend,j + 1)
                       ! set the corresponding entries in z to zero
                       call la_qlaset('FULL',in,2,zero,zero,z(ibegin,j),ldz)
                    end if
                    ! compute dl and dll of current rrr
                    do j = ibegin,iend - 1
                       tmp = d(j)*l(j)
                       work(indld - 1 + j) = tmp
                       work(indlld - 1 + j) = tmp*l(j)
                    end do
                    if (ndepth > 0) then
                       ! p and q are index of the first and last eigenvalue to compute
                       ! within the current block
                       p = indexw(wbegin - 1 + oldfst)
                       q = indexw(wbegin - 1 + oldlst)
                       ! offset for the arrays work, wgap and werr, i.e., the p-offset
                       ! through the q-offset elements of these arrays are to be used.
                        ! offset = p-oldfst
                       offset = indexw(wbegin) - 1
                       ! perform limited bisection (if necessary) to get approximate
                       ! eigenvalues to the precision needed.
                       call la_qlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,q,rtol1, &
                       rtol2,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk),iwork( &
                                  iindwk),pivmin,spdiam,in,iinfo)
                       if (iinfo /= 0) then
                          info = -1
                          return
                       end if
                       ! we also recompute the extremal gaps. w holds all eigenvalues
                       ! of the unshifted matrix and must be used for computation
                       ! of wgap, the entries of work might stem from rrrs with
                       ! different shifts. the gaps from wbegin-1+oldfst to
                       ! wbegin-1+oldlst are correctly computed in la_qlarrb.
                       ! however, we only allow the gaps to become greater since
                       ! this is what should happen when we decrease werr
                       if (oldfst > 1) then
                          wgap(wbegin + oldfst - 2) = max(wgap(wbegin + oldfst - 2),w(wbegin + oldfst - 1) - &
                          werr(wbegin + oldfst - 1) - w(wbegin + oldfst - 2) - werr(wbegin + oldfst - 2))

                       end if
                       if (wbegin + oldlst - 1 < wend) then
                          wgap(wbegin + oldlst - 1) = max(wgap(wbegin + oldlst - 1),w(wbegin + oldlst) - &
                                    werr(wbegin + oldlst) - w(wbegin + oldlst - 1) - werr(wbegin + oldlst - 1))
                       end if
                       ! each time the eigenvalues in work get refined, we store
                       ! the newly found approximation with all shifts applied in w
                       do j = oldfst,oldlst
                          w(wbegin + j - 1) = work(wbegin + j - 1) + sigma
                       end do
                    end if
                    ! process the current node.
                    newfst = oldfst
                    loop_140: do j = oldfst,oldlst
                       if (j == oldlst) then
                          ! we are at the right end of the cluster, this is also the
                          ! boundary of the child cluster
                          newlst = j
                       else if (wgap(wbegin + j - 1) >= minrgp*abs(work(wbegin + j - 1))) &
                                 then
                          ! the right relative gap is big enough, the child cluster
                          ! (newfst,..,newlst) is well separated from the following
                          newlst = j
                        else
                          ! inside a child cluster, the relative gap is not
                          ! big enough.
                          cycle loop_140
                       end if
                       ! compute size of child cluster found
                       newsiz = newlst - newfst + 1
                       ! newftt is the place in z where the new rrr or the computed
                       ! eigenvector is to be stored
                       if ((dol == 1) .and. (dou == m)) then
                          ! store representation at location of the leftmost evalue
                          ! of the cluster
                          newftt = wbegin + newfst - 1
                       else
                          if (wbegin + newfst - 1 < dol) then
                             ! store representation at the left end of z array
                             newftt = dol - 1
                          elseif (wbegin + newfst - 1 > dou) then
                             ! store representation at the right end of z array
                             newftt = dou
                          else
                             newftt = wbegin + newfst - 1
                          end if
                       end if
                       if (newsiz > 1) then
                          ! current child is not a singleton but a cluster.
                          ! compute and store new representation of child.
                          ! compute left and right cluster gap.
                          ! lgap and rgap are not computed from work because
                          ! the eigenvalue approximations may stem from rrrs
                          ! different shifts. however, w hold all eigenvalues
                          ! of the unshifted matrix. still, the entries in wgap
                          ! have to be computed from work since the entries
                          ! in w might be of the same order so that gaps are not
                          ! exhibited correctly for very close eigenvalues.
                          if (newfst == 1) then
                             lgap = max(zero,w(wbegin) - werr(wbegin) - vl)
                         else
                             lgap = wgap(wbegin + newfst - 2)
                          end if
                          rgap = wgap(wbegin + newlst - 1)
                          ! compute left- and rightmost eigenvalue of child
                          ! to high precision in order to shift as close
                          ! as possible and obtain as large relative gaps
                          ! as possible
                          do k = 1,2
                             if (k == 1) then
                                p = indexw(wbegin - 1 + newfst)
                             else
                                p = indexw(wbegin - 1 + newlst)
                             end if
                             offset = indexw(wbegin) - 1
                             call la_qlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,p,rqtol, &
                             rqtol,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk), &
                                       iwork(iindwk),pivmin,spdiam,in,iinfo)
                          end do
                          if ((wbegin + newlst - 1 < dol) .or. (wbegin + newfst - 1 > dou)) then
                             ! if the cluster contains no desired eigenvalues
                             ! skip the computation of that branch of the rep. tree
                             ! we could skip before the refinement of the extremal
                             ! eigenvalues of the child, but then the representation
                             ! tree could be different from the one when nothing is
                             ! skipped. for this reason we skip at this place.
                             idone = idone + newlst - newfst + 1
                             goto 139
                          end if
                          ! compute rrr of child cluster.
                          ! note that the new rrr is stored in z
                          ! la_qlarrf needs lwork = 2*n
                          call la_qlarrf(in,d(ibegin),l(ibegin),work(indld + ibegin - 1), &
                          newfst,newlst,work(wbegin),wgap(wbegin),werr(wbegin),spdiam,lgap, &
                          rgap,pivmin,tau,z(ibegin,newftt),z(ibegin,newftt + 1),work(indwrk), &
                                    iinfo)
                          if (iinfo == 0) then
                             ! a new rrr for the cluster was found by la_qlarrf
                             ! update shift and store it
                             ssigma = sigma + tau
                             z(iend,newftt + 1) = ssigma
                             ! work() are the midpoints and werr() the semi-width
                             ! note that the entries in w are unchanged.
                             do k = newfst,newlst
                                fudge = three*eps*abs(work(wbegin + k - 1))
                                work(wbegin + k - 1) = work(wbegin + k - 1) - tau
                                fudge = fudge + four*eps*abs(work(wbegin + k - 1))
                                ! fudge errors
                                werr(wbegin + k - 1) = werr(wbegin + k - 1) + fudge
                                ! gaps are not fudged. provided that werr is small
                                ! when eigenvalues are close, a zero gap indicates
                                ! that a new representation is needed for resolving
                                ! the cluster. a fudge could lead to a wrong decision
                                ! of judging eigenvalues 'separated' which in
                                ! reality are not. this could have a negative impact
                                ! on the orthogonality of the computed eigenvectors.
                             end do
                             nclus = nclus + 1
                             k = newcls + 2*nclus
                             iwork(k - 1) = newfst
                             iwork(k) = newlst
                          else
                             info = -2
                             return
                          end if
                       else
                          ! compute eigenvector of singleton
                          iter = 0
                          tol = four*log(real(in,KIND=qp))*eps
                          k = newfst
                          windex = wbegin + k - 1
                          windmn = max(windex - 1,1)
                          windpl = min(windex + 1,m)
                          lambda = work(windex)
                          done = done + 1
                          ! check if eigenvector computation is to be skipped
                          if ((windex < dol) .or. (windex > dou)) then
                             eskip = .true.
                             goto 125
                          else
                             eskip = .false.
                          end if
                          left = work(windex) - werr(windex)
                          right = work(windex) + werr(windex)
                          indeig = indexw(windex)
                          ! note that since we compute the eigenpairs for a child,
                          ! all eigenvalue approximations are w.r.t the same shift.
                          ! in this case, the entries in work should be used for
                          ! computing the gaps since they exhibit even very small
                          ! differences in the eigenvalues, as opposed to the
                          ! entries in w which might "look" the same.
                          if (k == 1) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vl, the formula
                             ! lgap = max( zero, (sigma - vl) + lambda )
                             ! can lead to an overestimation of the left gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small left gap.
                             lgap = eps*max(abs(left),abs(right))
                          else
                             lgap = wgap(windmn)
                          end if
                          if (k == im) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vu, the formula
                             ! can lead to an overestimation of the right gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small right gap.
                             rgap = eps*max(abs(left),abs(right))
                          else
                             rgap = wgap(windex)
                          end if
                          gap = min(lgap,rgap)
                          if ((k == 1) .or. (k == im)) then
                             ! the eigenvector support can become wrong
                             ! because significant entries could be cut off due to a
                             ! large gaptol parameter in lar1v. prevent this.
                             gaptol = zero
                          else
                             gaptol = gap*eps
                          end if
                          isupmn = in
                          isupmx = 1
                          ! update wgap so that it holds the minimum gap
                          ! to the left or the right. this is crucial in the
                          ! case where bisection is used to ensure that the
                          ! eigenvalue is refined up to the required precision.
                          ! the correct value is restored afterwards.
                          savgap = wgap(windex)
                          wgap(windex) = gap
                          ! we want to use the rayleigh quotient correction
                          ! as often as possible since it converges quadratically
                          ! when we are close enough to the desired eigenvalue.
                          ! however, the rayleigh quotient can have the wrong sign
                          ! and lead us away from the desired eigenvalue. in this
                          ! case, the best we can do is to use bisection.
                          usedbs = .false.
                          usedrq = .false.
                          ! bisection is initially turned off unless it is forced
                          needbs = .not. tryrqc
                          120 continue
                          ! check if bisection should be used to refine eigenvalue
                          if (needbs) then
                             ! take the bisection as new iterate
                             usedbs = .true.
                             itmp1 = iwork(iindr + windex)
                             offset = indexw(wbegin) - 1
                             call la_qlarrb(in,d(ibegin),work(indlld + ibegin - 1),indeig, &
                             indeig,zero,two*eps,offset,work(wbegin),wgap(wbegin),werr(wbegin), &
                                       work(indwrk),iwork(iindwk),pivmin,spdiam,itmp1,iinfo)
                             if (iinfo /= 0) then
                                info = -3
                                return
                             end if
                             lambda = work(windex)
                             ! reset twist index from inaccurate lambda to
                             ! force computation of true mingma
                             iwork(iindr + windex) = 0
                          end if
                          ! given lambda, compute the eigenvector.
                          call la_qlar1v(in,1,in,lambda,d(ibegin),l(ibegin),work( &
                          indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z(ibegin,windex &
                          ),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + windex),isuppz( &
                                    2*windex - 1),nrminv,resid,rqcorr,work(indwrk))
                          if (iter == 0) then
                             bstres = resid
                             bstw = lambda
                          elseif (resid < bstres) then
                             bstres = resid
                             bstw = lambda
                          end if
                          isupmn = min(isupmn,isuppz(2*windex - 1))
                          isupmx = max(isupmx,isuppz(2*windex))
                          iter = iter + 1
                          ! sin alpha <= |resid|/gap
                          ! note that both the residual and the gap are
                          ! proportional to the matrix, so ||t|| doesn't play
                          ! a role in the quotient
                          ! convergence test for rayleigh-quotient iteration
                          ! (omitted when bisection has been used)
                          if (resid > tol*gap .and. abs(rqcorr) > rqtol*abs(lambda) .and. .not. &
                                    usedbs) then
                             ! we need to check that the rqcorr update doesn't
                             ! move the eigenvalue away from the desired one and
                             ! towards a neighbor. -> protection with bisection
                             if (indeig <= negcnt) then
                                ! the wanted eigenvalue lies to the left
                                sgndef = -one
                             else
                                ! the wanted eigenvalue lies to the right
                                sgndef = one
                             end if
                             ! we only use the rqcorr if it improves the
                             ! the iterate reasonably.
                             if ((rqcorr*sgndef >= zero) .and. (lambda + rqcorr <= right) .and. ( &
                                       lambda + rqcorr >= left)) then
                                usedrq = .true.
                                ! store new midpoint of bisection interval in work
                                if (sgndef == one) then
                                   ! the current lambda is on the left of the true
                                   ! eigenvalue
                                   left = lambda
                                   ! we prefer to assume that the error estimate
                                   ! is correct. we could make the interval not
                                   ! as a bracket but to be modified if the rqcorr
                                   ! chooses to. in this case, the right side should
                                   ! be modified as follows:
                                    ! right = max(right, lambda + rqcorr)
                                else
                                   ! the current lambda is on the right of the true
                                   ! eigenvalue
                                   right = lambda
                                   ! see comment about assuming the error estimate is
                                   ! correct above.
                                    ! left = min(left, lambda + rqcorr)
                                end if
                                work(windex) = half*(right + left)
                                ! take rqcorr since it has the correct sign and
                                ! improves the iterate reasonably
                                lambda = lambda + rqcorr
                                ! update width of error interval
                                werr(windex) = half*(right - left)
                             else
                                needbs = .true.
                             end if
                             if (right - left < rqtol*abs(lambda)) then
                                   ! the eigenvalue is computed to bisection accuracy
                                   ! compute eigenvector and stop
                                usedbs = .true.
                                goto 120
                             elseif (iter < maxitr) then
                                goto 120
                             elseif (iter == maxitr) then
                                needbs = .true.
                                goto 120
                             else
                                info = 5
                                return
                             end if
                          else
                             stp2ii = .false.
             if (usedrq .and. usedbs .and. bstres <= resid) then
                                lambda = bstw
                                stp2ii = .true.
                             end if
                             if (stp2ii) then
                                ! improve error angle by second step
                                call la_qlar1v(in,1,in,lambda,d(ibegin),l(ibegin), &
                                work(indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z( &
                                ibegin,windex),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + &
                                windex),isuppz(2*windex - 1),nrminv,resid,rqcorr,work(indwrk &
                                          ))
                             end if
                             work(windex) = lambda
                          end if
                          ! compute fp-vector support w.r.t. whole matrix
                          isuppz(2*windex - 1) = isuppz(2*windex - 1) + oldien
                          isuppz(2*windex) = isuppz(2*windex) + oldien
                          zfrom = isuppz(2*windex - 1)
                          zto = isuppz(2*windex)
                          isupmn = isupmn + oldien
                          isupmx = isupmx + oldien
                          ! ensure vector is ok if support in the rqi has changed
                          if (isupmn < zfrom) then
                             do ii = isupmn,zfrom - 1
                                z(ii,windex) = zero
                             end do
                          end if
                          if (isupmx > zto) then
                             do ii = zto + 1,isupmx
                                z(ii,windex) = zero
                             end do
                          end if
                          call la_qscal(zto - zfrom + 1,nrminv,z(zfrom,windex),1)
                          125 continue
                          ! update w
                          w(windex) = lambda + sigma
                          ! recompute the gaps on the left and right
                          ! but only allow them to become larger and not
                          ! smaller (which can only happen through "bad"
                          ! cancellation and doesn't reflect the theory
                          ! where the initial gaps are underestimated due
                          ! to werr being too crude.)
                          if (.not. eskip) then
                             if (k > 1) then
                                wgap(windmn) = max(wgap(windmn),w(windex) - werr(windex) - w( &
                                          windmn) - werr(windmn))
                             end if
                             if (windex < wend) then
                                wgap(windex) = max(savgap,w(windpl) - werr(windpl) - w( &
                                          windex) - werr(windex))
                             end if
                          end if
                          idone = idone + 1
                       end if
                       ! here ends the code for the current child
                       139 continue
                       ! proceed to any remaining child nodes
                       newfst = j + 1
                    end do loop_140
                 end do loop_150
                 ndepth = ndepth + 1
                 go to 40
              end if
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_qlarrv

     !> To find the desired eigenvalues of a given real symmetric
     !> tridiagonal matrix T, SLARRE: sets any "small" off-diagonal
     !> elements to zero, and for each unreduced block T_i, it finds
     !> (a) a suitable shift at one end of the block's spectrum,
     !> (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and
     !> (c) eigenvalues of each L_i D_i L_i^T.
     !> The representations and eigenvalues found are then used by
     !> SSTEMR to compute the eigenvectors of T.
     !> The accuracy varies depending on whether bisection is used to
     !> find a few eigenvalues or the dqds algorithm (subroutine SLASQ2) to
     !> conpute all and then discard any unwanted one.
     !> As an added benefit, SLARRE also outputs the n
     !> Gerschgorin intervals for the matrices L_i D_i L_i^T.

     pure subroutine la_slarre(range,n,vl,vu,il,iu,d,e,e2,rtol1,rtol2,spltol, &
               nsplit,isplit,m,w,werr,wgap,iblock,indexw,gers,pivmin,work,iwork,info)
        use la_constants_sp,only:zero,half,one,two,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: range
           integer(ilp),intent(in) :: il,iu,n
           integer(ilp),intent(out) :: info,m,nsplit
           real(sp),intent(out) :: pivmin
           real(sp),intent(in) :: rtol1,rtol2,spltol
           real(sp),intent(inout) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),isplit(*),iwork(*),indexw(*)
           real(sp),intent(inout) :: d(*),e(*),e2(*)
           real(sp),intent(out) :: gers(*),w(*),werr(*),wgap(*),work(*)
        ! =====================================================================
           ! Parameters
           real(sp),parameter :: hndrd = 100.0_sp
           real(sp),parameter :: pert = 4.0_sp
           real(sp),parameter :: fourth = one/four
           real(sp),parameter :: fac = half
           real(sp),parameter :: maxgrowth = 64.0_sp
           real(sp),parameter :: fudge = 2.0_sp
           integer(ilp),parameter :: maxtry = 6
           integer(ilp),parameter :: allrng = 1
           integer(ilp),parameter :: indrng = 2
           integer(ilp),parameter :: valrng = 3

           ! Local Scalars
           logical(lk) :: forceb,norep,usedqd
           integer(ilp) :: cnt,cnt1,cnt2,i,ibegin,idum,iend,iinfo,in,indl,indu,irange, &
                     j,jblk,mb,mm,wbegin,wend
           real(sp) :: avgap,bsrtol,clwdth,dmax,dpivot,eabs,emax,eold,eps,gl,gu,isleft, &
                      isrght,rtl,rtol,s1,s2,safmin,sgndef,sigma,spdiam,tau,tmp,tmp1
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! decode range
           if (la_lsame(range,'A')) then
              irange = allrng
           else if (la_lsame(range,'V')) then
              irange = valrng
           else if (la_lsame(range,'I')) then
              irange = indrng
           end if
           m = 0
           ! get machine constants
           safmin = la_slamch('S')
           eps = la_slamch('P')
           ! set parameters
           rtl = hndrd*eps
           ! if one were ever to ask for less initial precision in bsrtol,
           ! one should keep in mind that for the subset case, the extremal
           ! eigenvalues must be at least as accurate as the current setting
           ! (eigenvalues in the middle need not as much accuracy)
           bsrtol = sqrt(eps)*(0.5e-3_sp)
           ! treat case of 1x1 matrix for quick return
           if (n == 1) then
              if ((irange == allrng) .or. ((irange == valrng) .and. (d(1) > vl) .and. (d(1) <= vu)) .or. (( &
                        irange == indrng) .and. (il == 1) .and. (iu == 1))) then
                 m = 1
                 w(1) = d(1)
                 ! the computation error of the eigenvalue is zero
                 werr(1) = zero
                 wgap(1) = zero
                 iblock(1) = 1
                 indexw(1) = 1
                 gers(1) = d(1)
                 gers(2) = d(1)
              end if
              ! store the shift for the initial rrr, which is zero in this case
              e(1) = zero
              return
           end if
           ! general case: tridiagonal matrix of order > 1
           ! init werr, wgap. compute gerschgorin intervals and spectral diameter.
           ! compute maximum off-diagonal entry and pivmin.
           gl = d(1)
           gu = d(1)
           eold = zero
           emax = zero
           e(n) = zero
           do i = 1,n
              werr(i) = zero
              wgap(i) = zero
              eabs = abs(e(i))
              if (eabs >= emax) then
                 emax = eabs
              end if
              tmp1 = eabs + eold
              gers(2*i - 1) = d(i) - tmp1
              gl = min(gl,gers(2*i - 1))
              gers(2*i) = d(i) + tmp1
              gu = max(gu,gers(2*i))
              eold = eabs
           end do
           ! the minimum pivot allowed in the sturm sequence for t
           pivmin = safmin*max(one,emax**2)
           ! compute spectral diameter. the gerschgorin bounds give an
           ! estimate that is wrong by at most a factor of sqrt(2)
           spdiam = gu - gl
           ! compute splitting points
           call la_slarra(n,d,e,e2,spltol,spdiam,nsplit,isplit,iinfo)
           ! can force use of bisection instead of faster dqds.
           ! option left in the code for future multisection work.
           forceb = .false.
           ! initialize usedqd, dqds should be used for allrng unless someone
           ! explicitly wants bisection.
           usedqd = ((irange == allrng) .and. (.not. forceb))
           if ((irange == allrng) .and. (.not. forceb)) then
              ! set interval [vl,vu] that contains all eigenvalues
              vl = gl
              vu = gu
           else
              ! we call la_slarrd to find crude approximations to the eigenvalues
              ! in the desired range. in case irange = indrng, we also obtain the
              ! interval (vl,vu] that contains all the wanted eigenvalues.
              ! an interval [left,right] has converged if
              ! right-left<rtol*max(abs(left),abs(right))
              ! la_slarrd needs a work of size 4*n, iwork of size 3*n
              call la_slarrd(range,'B',n,vl,vu,il,iu,gers,bsrtol,d,e,e2,pivmin, &
                        nsplit,isplit,mm,w,werr,vl,vu,iblock,indexw,work,iwork,iinfo)
              if (iinfo /= 0) then
                 info = -1
                 return
              end if
              ! make sure that the entries m+1 to n in w, werr, iblock, indexw are 0
              do i = mm + 1,n
                 w(i) = zero
                 werr(i) = zero
                 iblock(i) = 0
                 indexw(i) = 0
              end do
           end if
      ! **
           ! loop over unreduced blocks
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,nsplit
              iend = isplit(jblk)
              in = iend - ibegin + 1
              ! 1 x 1 block
              if (in == 1) then
                 if ((irange == allrng) .or. ((irange == valrng) .and. (d(ibegin) > vl) .and. (d( &
                           ibegin) <= vu)) .or. ((irange == indrng) .and. (iblock(wbegin) == jblk))) then
                    m = m + 1
                    w(m) = d(ibegin)
                    werr(m) = zero
                    ! the gap for a single block doesn't matter for the later
                    ! algorithm and is assigned an arbitrary large value
                    wgap(m) = zero
                    iblock(m) = jblk
                    indexw(m) = 1
                    wbegin = wbegin + 1
                 end if
                 ! e( iend ) holds the shift for the initial rrr
                 e(iend) = zero
                 ibegin = iend + 1
                 cycle loop_170
              end if
              ! blocks of size larger than 1x1
              ! e( iend ) will hold the shift for the initial rrr, for now set it =0
              e(iend) = zero
              ! find local outer bounds gl,gu for the block
              gl = d(ibegin)
              gu = d(ibegin)
              do i = ibegin,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              if (.not. ((irange == allrng) .and. (.not. forceb))) then
                 ! count the number of eigenvalues in the current block.
                 mb = 0
                 do i = wbegin,mm
                    if (iblock(i) == jblk) then
                       mb = mb + 1
                    else
                       goto 21
                    end if
                 end do
                 21 continue
                 if (mb == 0) then
                    ! no eigenvalue in the current block lies in the desired range
                    ! e( iend ) holds the shift for the initial rrr
                    e(iend) = zero
                    ibegin = iend + 1
                    cycle loop_170
                 else
                    ! decide whether dqds or bisection is more efficient
                    usedqd = ((mb > fac*in) .and. (.not. forceb))
                    wend = wbegin + mb - 1
                    ! calculate gaps for the current block
                    ! in later stages, when representations for individual
                    ! eigenvalues are different, we use sigma = e( iend ).
                    sigma = zero
                    do i = wbegin,wend - 1
                       wgap(i) = max(zero,w(i + 1) - werr(i + 1) - (w(i) + werr(i)))
                    end do
                    wgap(wend) = max(zero,vu - sigma - (w(wend) + werr(wend)))
                    ! find local index of the first and last desired evalue.
                    indl = indexw(wbegin)
                    indu = indexw(wend)
                 end if
              end if
              if (((irange == allrng) .and. (.not. forceb)) .or. usedqd) then
                 ! case of dqds
                 ! find approximations to the extremal eigenvalues of the block
                 call la_slarrk(in,1,gl,gu,d(ibegin),e2(ibegin),pivmin,rtl,tmp,tmp1, &
                           iinfo)
                 if (iinfo /= 0) then
                    info = -1
                    return
                 end if
                 isleft = max(gl,tmp - tmp1 - hndrd*eps*abs(tmp - tmp1))
                 call la_slarrk(in,in,gl,gu,d(ibegin),e2(ibegin),pivmin,rtl,tmp,tmp1, &
                            iinfo)
                 if (iinfo /= 0) then
                    info = -1
                    return
                 end if
                 isrght = min(gu,tmp + tmp1 + hndrd*eps*abs(tmp + tmp1))
                 ! improve the estimate of the spectral diameter
                 spdiam = isrght - isleft
              else
                 ! case of bisection
                 ! find approximations to the wanted extremal eigenvalues
                 isleft = max(gl,w(wbegin) - werr(wbegin) - hndrd*eps*abs(w(wbegin) - werr( &
                           wbegin)))
                 isrght = min(gu,w(wend) + werr(wend) + hndrd*eps*abs(w(wend) + werr(wend)))

              end if
              ! decide whether the base representation for the current block
              ! l_jblk d_jblk l_jblk^t = t_jblk - sigma_jblk i
              ! should be on the left or the right end of the current block.
              ! the strategy is to shift to the end which is "more populated"
              ! furthermore, decide whether to use dqds for the computation of
              ! dqds is chosen if all eigenvalues are desired or the number of
              ! eigenvalues to be computed is large compared to the blocksize.
              if ((irange == allrng) .and. (.not. forceb)) then
                 ! if all the eigenvalues have to be computed, we use dqd
                 usedqd = .true.
                 ! indl is the local index of the first eigenvalue to compute
                 indl = 1
                 indu = in
                 ! mb =  number of eigenvalues to compute
                 mb = in
                 wend = wbegin + mb - 1
                 ! define 1/4 and 3/4 points of the spectrum
                 s1 = isleft + fourth*spdiam
                 s2 = isrght - fourth*spdiam
              else
                 ! la_slarrd has computed iblock and indexw for each eigenvalue
                 ! approximation.
                 ! choose sigma
                 if (usedqd) then
                    s1 = isleft + fourth*spdiam
                    s2 = isrght - fourth*spdiam
                 else
                    tmp = min(isrght,vu) - max(isleft,vl)
                    s1 = max(isleft,vl) + fourth*tmp
                    s2 = min(isrght,vu) - fourth*tmp
                 end if
              end if
              ! compute the negcount at the 1/4 and 3/4 points
              if (mb > 1) then
                 call la_slarrc('T',in,s1,s2,d(ibegin),e(ibegin),pivmin,cnt,cnt1, &
                           cnt2,iinfo)
              end if
              if (mb == 1) then
                 sigma = gl
                 sgndef = one
              elseif (cnt1 - indl >= indu - cnt2) then
                 if ((irange == allrng) .and. (.not. forceb)) then
                    sigma = max(isleft,gl)
                 elseif (usedqd) then
                    ! use gerschgorin bound as shift to get pos def matrix
                    ! for dqds
                    sigma = isleft
                 else
                    ! use approximation of the first desired eigenvalue of the
                    ! block as shift
                    sigma = max(isleft,vl)
                 end if
                 sgndef = one
              else
                 if ((irange == allrng) .and. (.not. forceb)) then
                    sigma = min(isrght,gu)
                 elseif (usedqd) then
                    ! use gerschgorin bound as shift to get neg def matrix
                    ! for dqds
                    sigma = isrght
                 else
                    ! use approximation of the first desired eigenvalue of the
                    ! block as shift
                    sigma = min(isrght,vu)
                 end if
                 sgndef = -one
              end if
              ! an initial sigma has been chosen that will be used for computing
              ! t - sigma i = l d l^t
              ! define the increment tau of the shift in case the initial shift
              ! needs to be refined to obtain a factorization with not too much
              ! element growth.
              if (usedqd) then
                 ! the initial sigma was to the outer end of the spectrum
                 ! the matrix is definite and we need not retreat.
                 tau = spdiam*eps*n + two*pivmin
                 tau = max(tau,two*eps*abs(sigma))
              else
                 if (mb > 1) then
                    clwdth = w(wend) + werr(wend) - w(wbegin) - werr(wbegin)
                    avgap = abs(clwdth/real(wend - wbegin,KIND=sp))
                    if (sgndef == one) then
                       tau = half*max(wgap(wbegin),avgap)
                       tau = max(tau,werr(wbegin))
                    else
                       tau = half*max(wgap(wend - 1),avgap)
                       tau = max(tau,werr(wend))
                    end if
                 else
                    tau = werr(wbegin)
                 end if
              end if
              loop_80: do idum = 1,maxtry
                 ! compute l d l^t factorization of tridiagonal matrix t - sigma i.
                 ! store d in work(1:in), l in work(in+1:2*in), and reciprocals of
                 ! pivots in work(2*in+1:3*in)
                 dpivot = d(ibegin) - sigma
                 work(1) = dpivot
                 dmax = abs(work(1))
                 j = ibegin
                 do i = 1,in - 1
                    work(2*in + i) = one/work(i)
                    tmp = e(j)*work(2*in + i)
                    work(in + i) = tmp
                    dpivot = (d(j + 1) - sigma) - tmp*e(j)
                    work(i + 1) = dpivot
                    dmax = max(dmax,abs(dpivot))
                    j = j + 1
                 end do
                 ! check for element growth
                 if (dmax > maxgrowth*spdiam) then
                    norep = .true.
                 else
                    norep = .false.
                 end if
                 if (usedqd .and. .not. norep) then
                    ! ensure the definiteness of the representation
                    ! all entries of d (of l d l^t) must have the same sign
                    do i = 1,in
                       tmp = sgndef*work(i)
                       if (tmp < zero) norep = .true.
                    end do
                 end if
                 if (norep) then
                    ! note that in the case of irange=allrng, we use the gerschgorin
                    ! shift which makes the matrix definite. so we should end up
                    ! here really only in the case of irange = valrng or indrng.
                    if (idum == maxtry - 1) then
                       if (sgndef == one) then
                          ! the fudged gerschgorin shift should succeed
                          sigma = gl - fudge*spdiam*eps*n - fudge*two*pivmin
                       else
                          sigma = gu + fudge*spdiam*eps*n + fudge*two*pivmin
                       end if
                    else
                       sigma = sigma - sgndef*tau
                       tau = two*tau
                    end if
                 else
                    ! an initial rrr is found
                    go to 83
                 end if
              end do loop_80
              ! if the program reaches this point, no base representation could be
              ! found in maxtry iterations.
              info = 2
              return
              83 continue
              ! at this point, we have found an initial base representation
              ! t - sigma i = l d l^t with not too much element growth.
              ! store the shift.
              e(iend) = sigma
              ! store d and l.
              call la_scopy(in,work,1,d(ibegin),1)
              call la_scopy(in - 1,work(in + 1),1,e(ibegin),1)
              if (mb > 1) then
                 ! perturb each entry of the base representation by a small
                 ! (but random) relative amount to overcome difficulties with
                 ! glued matrices.
                 do i = 1,4
                    iseed(i) = 1
                 end do
                 call la_slarnv(2,iseed,2*in - 1,work(1))
                 do i = 1,in - 1
                    d(ibegin + i - 1) = d(ibegin + i - 1)*(one + eps*pert*work(i))
                    e(ibegin + i - 1) = e(ibegin + i - 1)*(one + eps*pert*work(in + i))
                 end do
                 d(iend) = d(iend)*(one + eps*four*work(in))
              end if
              ! don't update the gerschgorin intervals because keeping track
              ! of the updates would be too much work in la_slarrv.
              ! we update w instead and use it to locate the proper gerschgorin
              ! intervals.
              ! compute the required eigenvalues of l d l' by bisection or dqds
              if (.not. usedqd) then
                 ! if la_slarrd has been used, shift the eigenvalue approximations
                 ! according to their representation. this is necessary for
                 ! a uniform la_slarrv since dqds computes eigenvalues of the
                 ! shifted representation. in la_slarrv, w will always hold the
                 ! unshifted eigenvalue approximation.
                 do j = wbegin,wend
                    w(j) = w(j) - sigma
                    werr(j) = werr(j) + abs(w(j))*eps
                 end do
                 ! call la_slarrb to reduce eigenvalue error of the approximations
                 ! from la_slarrd
                 do i = ibegin,iend - 1
                    work(i) = d(i)*e(i)**2
                 end do
                 ! use bisection to find ev from indl to indu
                 call la_slarrb(in,d(ibegin),work(ibegin),indl,indu,rtol1,rtol2,indl - 1, &
                 w(wbegin),wgap(wbegin),werr(wbegin),work(2*n + 1),iwork,pivmin,spdiam,in, &
                           iinfo)
                 if (iinfo /= 0) then
                    info = -4
                    return
                 end if
                 ! la_slarrb computes all gaps correctly except for the last one
                 ! record distance to vu/gu
                 wgap(wend) = max(zero, (vu - sigma) - (w(wend) + werr(wend)))
                 do i = indl,indu
                    m = m + 1
                    iblock(m) = jblk
                    indexw(m) = i
                 end do
              else
                 ! call dqds to get all eigs (and then possibly delete unwanted
                 ! eigenvalues).
                 ! note that dqds finds the eigenvalues of the l d l^t representation
                 ! of t to high relative accuracy. high relative accuracy
                 ! might be lost when the shift of the rrr is subtracted to obtain
                 ! the eigenvalues of t. however, t is not guaranteed to define its
                 ! eigenvalues to high relative accuracy anyway.
                 ! set rtol to the order of the tolerance used in la_slasq2
                 ! this is an estimated error, the worst case bound is 4*n*eps
                 ! which is usually too large and requires unnecessary work to be
                 ! done by bisection when computing the eigenvectors
                 rtol = log(real(in,KIND=sp))*four*eps
                 j = ibegin
                 do i = 1,in - 1
                    work(2*i - 1) = abs(d(j))
                    work(2*i) = e(j)*e(j)*work(2*i - 1)
                    j = j + 1
                 end do
                 work(2*in - 1) = abs(d(iend))
                 work(2*in) = zero
                 call la_slasq2(in,work,iinfo)
                 if (iinfo /= 0) then
                    ! if iinfo = -5 then an index is part of a tight cluster
                    ! and should be changed. the index is in iwork(1) and the
                    ! gap is in work(n+1)
                    info = -5
                    return
                 else
                    ! test that all eigenvalues are positive as expected
                    do i = 1,in
                       if (work(i) < zero) then
                          info = -6
                          return
                       end if
                    end do
                 end if
                 if (sgndef > zero) then
                    do i = indl,indu
                       m = m + 1
                       w(m) = work(in - i + 1)
                       iblock(m) = jblk
                       indexw(m) = i
                    end do
                 else
                    do i = indl,indu
                       m = m + 1
                       w(m) = -work(i)
                       iblock(m) = jblk
                       indexw(m) = i
                    end do
                 end if
                 do i = m - mb + 1,m
                    ! the value of rtol below should be the tolerance in la_slasq2
                    werr(i) = rtol*abs(w(i))
                 end do
                 do i = m - mb + 1,m - 1
                    ! compute the right gap between the intervals
                    wgap(i) = max(zero,w(i + 1) - werr(i + 1) - (w(i) + werr(i)))
                 end do
                 wgap(m) = max(zero, (vu - sigma) - (w(m) + werr(m)))
              end if
              ! proceed with next block
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_slarre
     !> To find the desired eigenvalues of a given real symmetric
     !> tridiagonal matrix T, DLARRE: sets any "small" off-diagonal
     !> elements to zero, and for each unreduced block T_i, it finds
     !> (a) a suitable shift at one end of the block's spectrum,
     !> (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and
     !> (c) eigenvalues of each L_i D_i L_i^T.
     !> The representations and eigenvalues found are then used by
     !> DSTEMR to compute the eigenvectors of T.
     !> The accuracy varies depending on whether bisection is used to
     !> find a few eigenvalues or the dqds algorithm (subroutine DLASQ2) to
     !> conpute all and then discard any unwanted one.
     !> As an added benefit, DLARRE also outputs the n
     !> Gerschgorin intervals for the matrices L_i D_i L_i^T.

     pure subroutine la_dlarre(range,n,vl,vu,il,iu,d,e,e2,rtol1,rtol2,spltol, &
               nsplit,isplit,m,w,werr,wgap,iblock,indexw,gers,pivmin,work,iwork,info)
        use la_constants_dp,only:zero,half,one,two,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: range
           integer(ilp),intent(in) :: il,iu,n
           integer(ilp),intent(out) :: info,m,nsplit
           real(dp),intent(out) :: pivmin
           real(dp),intent(in) :: rtol1,rtol2,spltol
           real(dp),intent(inout) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),isplit(*),iwork(*),indexw(*)
           real(dp),intent(inout) :: d(*),e(*),e2(*)
           real(dp),intent(out) :: gers(*),w(*),werr(*),wgap(*),work(*)
        ! =====================================================================
           ! Parameters
           real(dp),parameter :: hndrd = 100.0_dp
           real(dp),parameter :: pert = 8.0_dp
           real(dp),parameter :: fourth = one/four
           real(dp),parameter :: fac = half
           real(dp),parameter :: maxgrowth = 64.0_dp
           real(dp),parameter :: fudge = 2.0_dp
           integer(ilp),parameter :: maxtry = 6
           integer(ilp),parameter :: allrng = 1
           integer(ilp),parameter :: indrng = 2
           integer(ilp),parameter :: valrng = 3

           ! Local Scalars
           logical(lk) :: forceb,norep,usedqd
           integer(ilp) :: cnt,cnt1,cnt2,i,ibegin,idum,iend,iinfo,in,indl,indu,irange, &
                     j,jblk,mb,mm,wbegin,wend
           real(dp) :: avgap,bsrtol,clwdth,dmax,dpivot,eabs,emax,eold,eps,gl,gu,isleft, &
                      isrght,rtl,rtol,s1,s2,safmin,sgndef,sigma,spdiam,tau,tmp,tmp1
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! decode range
           if (la_lsame(range,'A')) then
              irange = allrng
           else if (la_lsame(range,'V')) then
              irange = valrng
           else if (la_lsame(range,'I')) then
              irange = indrng
           end if
           m = 0
           ! get machine constants
           safmin = la_dlamch('S')
           eps = la_dlamch('P')
           ! set parameters
           rtl = sqrt(eps)
           bsrtol = sqrt(eps)
           ! treat case of 1x1 matrix for quick return
           if (n == 1) then
              if ((irange == allrng) .or. ((irange == valrng) .and. (d(1) > vl) .and. (d(1) <= vu)) .or. (( &
                        irange == indrng) .and. (il == 1) .and. (iu == 1))) then
                 m = 1
                 w(1) = d(1)
                 ! the computation error of the eigenvalue is zero
                 werr(1) = zero
                 wgap(1) = zero
                 iblock(1) = 1
                 indexw(1) = 1
                 gers(1) = d(1)
                 gers(2) = d(1)
              end if
              ! store the shift for the initial rrr, which is zero in this case
              e(1) = zero
              return
           end if
           ! general case: tridiagonal matrix of order > 1
           ! init werr, wgap. compute gerschgorin intervals and spectral diameter.
           ! compute maximum off-diagonal entry and pivmin.
           gl = d(1)
           gu = d(1)
           eold = zero
           emax = zero
           e(n) = zero
           do i = 1,n
              werr(i) = zero
              wgap(i) = zero
              eabs = abs(e(i))
              if (eabs >= emax) then
                 emax = eabs
              end if
              tmp1 = eabs + eold
              gers(2*i - 1) = d(i) - tmp1
              gl = min(gl,gers(2*i - 1))
              gers(2*i) = d(i) + tmp1
              gu = max(gu,gers(2*i))
              eold = eabs
           end do
           ! the minimum pivot allowed in the sturm sequence for t
           pivmin = safmin*max(one,emax**2)
           ! compute spectral diameter. the gerschgorin bounds give an
           ! estimate that is wrong by at most a factor of sqrt(2)
           spdiam = gu - gl
           ! compute splitting points
           call la_dlarra(n,d,e,e2,spltol,spdiam,nsplit,isplit,iinfo)
           ! can force use of bisection instead of faster dqds.
           ! option left in the code for future multisection work.
           forceb = .false.
           ! initialize usedqd, dqds should be used for allrng unless someone
           ! explicitly wants bisection.
           usedqd = ((irange == allrng) .and. (.not. forceb))
           if ((irange == allrng) .and. (.not. forceb)) then
              ! set interval [vl,vu] that contains all eigenvalues
              vl = gl
              vu = gu
           else
              ! we call la_dlarrd to find crude approximations to the eigenvalues
              ! in the desired range. in case irange = indrng, we also obtain the
              ! interval (vl,vu] that contains all the wanted eigenvalues.
              ! an interval [left,right] has converged if
              ! right-left<rtol*max(abs(left),abs(right))
              ! la_dlarrd needs a work of size 4*n, iwork of size 3*n
              call la_dlarrd(range,'B',n,vl,vu,il,iu,gers,bsrtol,d,e,e2,pivmin, &
                        nsplit,isplit,mm,w,werr,vl,vu,iblock,indexw,work,iwork,iinfo)
              if (iinfo /= 0) then
                 info = -1
                 return
              end if
              ! make sure that the entries m+1 to n in w, werr, iblock, indexw are 0
              do i = mm + 1,n
                 w(i) = zero
                 werr(i) = zero
                 iblock(i) = 0
                 indexw(i) = 0
              end do
           end if
      ! **
           ! loop over unreduced blocks
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,nsplit
              iend = isplit(jblk)
              in = iend - ibegin + 1
              ! 1 x 1 block
              if (in == 1) then
                 if ((irange == allrng) .or. ((irange == valrng) .and. (d(ibegin) > vl) .and. (d( &
                           ibegin) <= vu)) .or. ((irange == indrng) .and. (iblock(wbegin) == jblk))) then
                    m = m + 1
                    w(m) = d(ibegin)
                    werr(m) = zero
                    ! the gap for a single block doesn't matter for the later
                    ! algorithm and is assigned an arbitrary large value
                    wgap(m) = zero
                    iblock(m) = jblk
                    indexw(m) = 1
                    wbegin = wbegin + 1
                 end if
                 ! e( iend ) holds the shift for the initial rrr
                 e(iend) = zero
                 ibegin = iend + 1
                 cycle loop_170
              end if
              ! blocks of size larger than 1x1
              ! e( iend ) will hold the shift for the initial rrr, for now set it =0
              e(iend) = zero
              ! find local outer bounds gl,gu for the block
              gl = d(ibegin)
              gu = d(ibegin)
              do i = ibegin,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              if (.not. ((irange == allrng) .and. (.not. forceb))) then
                 ! count the number of eigenvalues in the current block.
                 mb = 0
                 do i = wbegin,mm
                    if (iblock(i) == jblk) then
                       mb = mb + 1
                    else
                       goto 21
                    end if
                 end do
                 21 continue
                 if (mb == 0) then
                    ! no eigenvalue in the current block lies in the desired range
                    ! e( iend ) holds the shift for the initial rrr
                    e(iend) = zero
                    ibegin = iend + 1
                    cycle loop_170
                 else
                    ! decide whether dqds or bisection is more efficient
                    usedqd = ((mb > fac*in) .and. (.not. forceb))
                    wend = wbegin + mb - 1
                    ! calculate gaps for the current block
                    ! in later stages, when representations for individual
                    ! eigenvalues are different, we use sigma = e( iend ).
                    sigma = zero
                    do i = wbegin,wend - 1
                       wgap(i) = max(zero,w(i + 1) - werr(i + 1) - (w(i) + werr(i)))
                    end do
                    wgap(wend) = max(zero,vu - sigma - (w(wend) + werr(wend)))
                    ! find local index of the first and last desired evalue.
                    indl = indexw(wbegin)
                    indu = indexw(wend)
                 end if
              end if
              if (((irange == allrng) .and. (.not. forceb)) .or. usedqd) then
                 ! case of dqds
                 ! find approximations to the extremal eigenvalues of the block
                 call la_dlarrk(in,1,gl,gu,d(ibegin),e2(ibegin),pivmin,rtl,tmp,tmp1, &
                           iinfo)
                 if (iinfo /= 0) then
                    info = -1
                    return
                 end if
                 isleft = max(gl,tmp - tmp1 - hndrd*eps*abs(tmp - tmp1))
                 call la_dlarrk(in,in,gl,gu,d(ibegin),e2(ibegin),pivmin,rtl,tmp,tmp1, &
                            iinfo)
                 if (iinfo /= 0) then
                    info = -1
                    return
                 end if
                 isrght = min(gu,tmp + tmp1 + hndrd*eps*abs(tmp + tmp1))
                 ! improve the estimate of the spectral diameter
                 spdiam = isrght - isleft
              else
                 ! case of bisection
                 ! find approximations to the wanted extremal eigenvalues
                 isleft = max(gl,w(wbegin) - werr(wbegin) - hndrd*eps*abs(w(wbegin) - werr( &
                           wbegin)))
                 isrght = min(gu,w(wend) + werr(wend) + hndrd*eps*abs(w(wend) + werr(wend)))

              end if
              ! decide whether the base representation for the current block
              ! l_jblk d_jblk l_jblk^t = t_jblk - sigma_jblk i
              ! should be on the left or the right end of the current block.
              ! the strategy is to shift to the end which is "more populated"
              ! furthermore, decide whether to use dqds for the computation of
              ! dqds is chosen if all eigenvalues are desired or the number of
              ! eigenvalues to be computed is large compared to the blocksize.
              if ((irange == allrng) .and. (.not. forceb)) then
                 ! if all the eigenvalues have to be computed, we use dqd
                 usedqd = .true.
                 ! indl is the local index of the first eigenvalue to compute
                 indl = 1
                 indu = in
                 ! mb =  number of eigenvalues to compute
                 mb = in
                 wend = wbegin + mb - 1
                 ! define 1/4 and 3/4 points of the spectrum
                 s1 = isleft + fourth*spdiam
                 s2 = isrght - fourth*spdiam
              else
                 ! la_dlarrd has computed iblock and indexw for each eigenvalue
                 ! approximation.
                 ! choose sigma
                 if (usedqd) then
                    s1 = isleft + fourth*spdiam
                    s2 = isrght - fourth*spdiam
                 else
                    tmp = min(isrght,vu) - max(isleft,vl)
                    s1 = max(isleft,vl) + fourth*tmp
                    s2 = min(isrght,vu) - fourth*tmp
                 end if
              end if
              ! compute the negcount at the 1/4 and 3/4 points
              if (mb > 1) then
                 call la_dlarrc('T',in,s1,s2,d(ibegin),e(ibegin),pivmin,cnt,cnt1, &
                           cnt2,iinfo)
              end if
              if (mb == 1) then
                 sigma = gl
                 sgndef = one
              elseif (cnt1 - indl >= indu - cnt2) then
                 if ((irange == allrng) .and. (.not. forceb)) then
                    sigma = max(isleft,gl)
                 elseif (usedqd) then
                    ! use gerschgorin bound as shift to get pos def matrix
                    ! for dqds
                    sigma = isleft
                 else
                    ! use approximation of the first desired eigenvalue of the
                    ! block as shift
                    sigma = max(isleft,vl)
                 end if
                 sgndef = one
              else
                 if ((irange == allrng) .and. (.not. forceb)) then
                    sigma = min(isrght,gu)
                 elseif (usedqd) then
                    ! use gerschgorin bound as shift to get neg def matrix
                    ! for dqds
                    sigma = isrght
                 else
                    ! use approximation of the first desired eigenvalue of the
                    ! block as shift
                    sigma = min(isrght,vu)
                 end if
                 sgndef = -one
              end if
              ! an initial sigma has been chosen that will be used for computing
              ! t - sigma i = l d l^t
              ! define the increment tau of the shift in case the initial shift
              ! needs to be refined to obtain a factorization with not too much
              ! element growth.
              if (usedqd) then
                 ! the initial sigma was to the outer end of the spectrum
                 ! the matrix is definite and we need not retreat.
                 tau = spdiam*eps*n + two*pivmin
                 tau = max(tau,two*eps*abs(sigma))
              else
                 if (mb > 1) then
                    clwdth = w(wend) + werr(wend) - w(wbegin) - werr(wbegin)
                    avgap = abs(clwdth/real(wend - wbegin,KIND=dp))
                    if (sgndef == one) then
                       tau = half*max(wgap(wbegin),avgap)
                       tau = max(tau,werr(wbegin))
                    else
                       tau = half*max(wgap(wend - 1),avgap)
                       tau = max(tau,werr(wend))
                    end if
                 else
                    tau = werr(wbegin)
                 end if
              end if
              loop_80: do idum = 1,maxtry
                 ! compute l d l^t factorization of tridiagonal matrix t - sigma i.
                 ! store d in work(1:in), l in work(in+1:2*in), and reciprocals of
                 ! pivots in work(2*in+1:3*in)
                 dpivot = d(ibegin) - sigma
                 work(1) = dpivot
                 dmax = abs(work(1))
                 j = ibegin
                 do i = 1,in - 1
                    work(2*in + i) = one/work(i)
                    tmp = e(j)*work(2*in + i)
                    work(in + i) = tmp
                    dpivot = (d(j + 1) - sigma) - tmp*e(j)
                    work(i + 1) = dpivot
                    dmax = max(dmax,abs(dpivot))
                    j = j + 1
                 end do
                 ! check for element growth
                 if (dmax > maxgrowth*spdiam) then
                    norep = .true.
                 else
                    norep = .false.
                 end if
                 if (usedqd .and. .not. norep) then
                    ! ensure the definiteness of the representation
                    ! all entries of d (of l d l^t) must have the same sign
                    do i = 1,in
                       tmp = sgndef*work(i)
                       if (tmp < zero) norep = .true.
                    end do
                 end if
                 if (norep) then
                    ! note that in the case of irange=allrng, we use the gerschgorin
                    ! shift which makes the matrix definite. so we should end up
                    ! here really only in the case of irange = valrng or indrng.
                    if (idum == maxtry - 1) then
                       if (sgndef == one) then
                          ! the fudged gerschgorin shift should succeed
                          sigma = gl - fudge*spdiam*eps*n - fudge*two*pivmin
                       else
                          sigma = gu + fudge*spdiam*eps*n + fudge*two*pivmin
                       end if
                    else
                       sigma = sigma - sgndef*tau
                       tau = two*tau
                    end if
                 else
                    ! an initial rrr is found
                    go to 83
                 end if
              end do loop_80
              ! if the program reaches this point, no base representation could be
              ! found in maxtry iterations.
              info = 2
              return
              83 continue
              ! at this point, we have found an initial base representation
              ! t - sigma i = l d l^t with not too much element growth.
              ! store the shift.
              e(iend) = sigma
              ! store d and l.
              call la_dcopy(in,work,1,d(ibegin),1)
              call la_dcopy(in - 1,work(in + 1),1,e(ibegin),1)
              if (mb > 1) then
                 ! perturb each entry of the base representation by a small
                 ! (but random) relative amount to overcome difficulties with
                 ! glued matrices.
                 do i = 1,4
                    iseed(i) = 1
                 end do
                 call la_dlarnv(2,iseed,2*in - 1,work(1))
                 do i = 1,in - 1
                    d(ibegin + i - 1) = d(ibegin + i - 1)*(one + eps*pert*work(i))
                    e(ibegin + i - 1) = e(ibegin + i - 1)*(one + eps*pert*work(in + i))
                 end do
                 d(iend) = d(iend)*(one + eps*four*work(in))
              end if
              ! don't update the gerschgorin intervals because keeping track
              ! of the updates would be too much work in la_dlarrv.
              ! we update w instead and use it to locate the proper gerschgorin
              ! intervals.
              ! compute the required eigenvalues of l d l' by bisection or dqds
              if (.not. usedqd) then
                 ! if la_dlarrd has been used, shift the eigenvalue approximations
                 ! according to their representation. this is necessary for
                 ! a uniform la_dlarrv since dqds computes eigenvalues of the
                 ! shifted representation. in la_dlarrv, w will always hold the
                 ! unshifted eigenvalue approximation.
                 do j = wbegin,wend
                    w(j) = w(j) - sigma
                    werr(j) = werr(j) + abs(w(j))*eps
                 end do
                 ! call la_dlarrb to reduce eigenvalue error of the approximations
                 ! from la_dlarrd
                 do i = ibegin,iend - 1
                    work(i) = d(i)*e(i)**2
                 end do
                 ! use bisection to find ev from indl to indu
                 call la_dlarrb(in,d(ibegin),work(ibegin),indl,indu,rtol1,rtol2,indl - 1, &
                 w(wbegin),wgap(wbegin),werr(wbegin),work(2*n + 1),iwork,pivmin,spdiam,in, &
                           iinfo)
                 if (iinfo /= 0) then
                    info = -4
                    return
                 end if
                 ! la_dlarrb computes all gaps correctly except for the last one
                 ! record distance to vu/gu
                 wgap(wend) = max(zero, (vu - sigma) - (w(wend) + werr(wend)))
                 do i = indl,indu
                    m = m + 1
                    iblock(m) = jblk
                    indexw(m) = i
                 end do
              else
                 ! call dqds to get all eigs (and then possibly delete unwanted
                 ! eigenvalues).
                 ! note that dqds finds the eigenvalues of the l d l^t representation
                 ! of t to high relative accuracy. high relative accuracy
                 ! might be lost when the shift of the rrr is subtracted to obtain
                 ! the eigenvalues of t. however, t is not guaranteed to define its
                 ! eigenvalues to high relative accuracy anyway.
                 ! set rtol to the order of the tolerance used in la_dlasq2
                 ! this is an estimated error, the worst case bound is 4*n*eps
                 ! which is usually too large and requires unnecessary work to be
                 ! done by bisection when computing the eigenvectors
                 rtol = log(real(in,KIND=dp))*four*eps
                 j = ibegin
                 do i = 1,in - 1
                    work(2*i - 1) = abs(d(j))
                    work(2*i) = e(j)*e(j)*work(2*i - 1)
                    j = j + 1
                 end do
                 work(2*in - 1) = abs(d(iend))
                 work(2*in) = zero
                 call la_dlasq2(in,work,iinfo)
                 if (iinfo /= 0) then
                    ! if iinfo = -5 then an index is part of a tight cluster
                    ! and should be changed. the index is in iwork(1) and the
                    ! gap is in work(n+1)
                    info = -5
                    return
                 else
                    ! test that all eigenvalues are positive as expected
                    do i = 1,in
                       if (work(i) < zero) then
                          info = -6
                          return
                       end if
                    end do
                 end if
                 if (sgndef > zero) then
                    do i = indl,indu
                       m = m + 1
                       w(m) = work(in - i + 1)
                       iblock(m) = jblk
                       indexw(m) = i
                    end do
                 else
                    do i = indl,indu
                       m = m + 1
                       w(m) = -work(i)
                       iblock(m) = jblk
                       indexw(m) = i
                    end do
                 end if
                 do i = m - mb + 1,m
                    ! the value of rtol below should be the tolerance in la_dlasq2
                    werr(i) = rtol*abs(w(i))
                 end do
                 do i = m - mb + 1,m - 1
                    ! compute the right gap between the intervals
                    wgap(i) = max(zero,w(i + 1) - werr(i + 1) - (w(i) + werr(i)))
                 end do
                 wgap(m) = max(zero, (vu - sigma) - (w(m) + werr(m)))
              end if
              ! proceed with next block
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_dlarre
     !> To find the desired eigenvalues of a given real symmetric
     !> tridiagonal matrix T, QLARRE: sets any "small" off-diagonal
     !> elements to zero, and for each unreduced block T_i, it finds
     !> (a) a suitable shift at one end of the block's spectrum,
     !> (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and
     !> (c) eigenvalues of each L_i D_i L_i^T.
     !> The representations and eigenvalues found are then used by
     !> QSTEMR to compute the eigenvectors of T.
     !> The accuracy varies depending on whether bisection is used to
     !> find a few eigenvalues or the dqds algorithm (subroutine QLASQ2) to
     !> conpute all and then discard any unwanted one.
     !> As an added benefit, QLARRE also outputs the n
     !> Gerschgorin intervals for the matrices L_i D_i L_i^T.

     pure subroutine la_qlarre(range,n,vl,vu,il,iu,d,e,e2,rtol1,rtol2,spltol, &
               nsplit,isplit,m,w,werr,wgap,iblock,indexw,gers,pivmin,work,iwork,info)
        use la_constants_qp,only:zero,half,one,two,four
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: range
           integer(ilp),intent(in) :: il,iu,n
           integer(ilp),intent(out) :: info,m,nsplit
           real(qp),intent(out) :: pivmin
           real(qp),intent(in) :: rtol1,rtol2,spltol
           real(qp),intent(inout) :: vl,vu
           ! Array Arguments
           integer(ilp),intent(out) :: iblock(*),isplit(*),iwork(*),indexw(*)
           real(qp),intent(inout) :: d(*),e(*),e2(*)
           real(qp),intent(out) :: gers(*),w(*),werr(*),wgap(*),work(*)
        ! =====================================================================
           ! Parameters
           real(qp),parameter :: hndrd = 100.0_qp
           real(qp),parameter :: pert = 8.0_qp
           real(qp),parameter :: fourth = one/four
           real(qp),parameter :: fac = half
           real(qp),parameter :: maxgrowth = 64.0_qp
           real(qp),parameter :: fudge = 2.0_qp
           integer(ilp),parameter :: maxtry = 6
           integer(ilp),parameter :: allrng = 1
           integer(ilp),parameter :: indrng = 2
           integer(ilp),parameter :: valrng = 3

           ! Local Scalars
           logical(lk) :: forceb,norep,usedqd
           integer(ilp) :: cnt,cnt1,cnt2,i,ibegin,idum,iend,iinfo,in,indl,indu,irange, &
                     j,jblk,mb,mm,wbegin,wend
           real(qp) :: avgap,bsrtol,clwdth,dmax,dpivot,eabs,emax,eold,eps,gl,gu,isleft, &
                      isrght,rtl,rtol,s1,s2,safmin,sgndef,sigma,spdiam,tau,tmp,tmp1
           ! Local Arrays
           integer(ilp) :: iseed(4)
           ! Intrinsic Functions
           intrinsic :: abs,max,min
           ! Executable Statements
           info = 0
           ! quick return if possible
           if (n <= 0) then
              return
           end if
           ! decode range
           if (la_lsame(range,'A')) then
              irange = allrng
           else if (la_lsame(range,'V')) then
              irange = valrng
           else if (la_lsame(range,'I')) then
              irange = indrng
           end if
           m = 0
           ! get machine constants
           safmin = la_qlamch('S')
           eps = la_qlamch('P')
           ! set parameters
           rtl = sqrt(eps)
           bsrtol = sqrt(eps)
           ! treat case of 1x1 matrix for quick return
           if (n == 1) then
              if ((irange == allrng) .or. ((irange == valrng) .and. (d(1) > vl) .and. (d(1) <= vu)) .or. (( &
                        irange == indrng) .and. (il == 1) .and. (iu == 1))) then
                 m = 1
                 w(1) = d(1)
                 ! the computation error of the eigenvalue is zero
                 werr(1) = zero
                 wgap(1) = zero
                 iblock(1) = 1
                 indexw(1) = 1
                 gers(1) = d(1)
                 gers(2) = d(1)
              end if
              ! store the shift for the initial rrr, which is zero in this case
              e(1) = zero
              return
           end if
           ! general case: tridiagonal matrix of order > 1
           ! init werr, wgap. compute gerschgorin intervals and spectral diameter.
           ! compute maximum off-diagonal entry and pivmin.
           gl = d(1)
           gu = d(1)
           eold = zero
           emax = zero
           e(n) = zero
           do i = 1,n
              werr(i) = zero
              wgap(i) = zero
              eabs = abs(e(i))
              if (eabs >= emax) then
                 emax = eabs
              end if
              tmp1 = eabs + eold
              gers(2*i - 1) = d(i) - tmp1
              gl = min(gl,gers(2*i - 1))
              gers(2*i) = d(i) + tmp1
              gu = max(gu,gers(2*i))
              eold = eabs
           end do
           ! the minimum pivot allowed in the sturm sequence for t
           pivmin = safmin*max(one,emax**2)
           ! compute spectral diameter. the gerschgorin bounds give an
           ! estimate that is wrong by at most a factor of sqrt(2)
           spdiam = gu - gl
           ! compute splitting points
           call la_qlarra(n,d,e,e2,spltol,spdiam,nsplit,isplit,iinfo)
           ! can force use of bisection instead of faster dqds.
           ! option left in the code for future multisection work.
           forceb = .false.
           ! initialize usedqd, dqds should be used for allrng unless someone
           ! explicitly wants bisection.
           usedqd = ((irange == allrng) .and. (.not. forceb))
           if ((irange == allrng) .and. (.not. forceb)) then
              ! set interval [vl,vu] that contains all eigenvalues
              vl = gl
              vu = gu
           else
              ! we call la_qlarrd to find crude approximations to the eigenvalues
              ! in the desired range. in case irange = indrng, we also obtain the
              ! interval (vl,vu] that contains all the wanted eigenvalues.
              ! an interval [left,right] has converged if
              ! right-left<rtol*max(abs(left),abs(right))
              ! la_qlarrd needs a work of size 4*n, iwork of size 3*n
              call la_qlarrd(range,'B',n,vl,vu,il,iu,gers,bsrtol,d,e,e2,pivmin, &
                        nsplit,isplit,mm,w,werr,vl,vu,iblock,indexw,work,iwork,iinfo)
              if (iinfo /= 0) then
                 info = -1
                 return
              end if
              ! make sure that the entries m+1 to n in w, werr, iblock, indexw are 0
              do i = mm + 1,n
                 w(i) = zero
                 werr(i) = zero
                 iblock(i) = 0
                 indexw(i) = 0
              end do
           end if
      ! **
           ! loop over unreduced blocks
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,nsplit
              iend = isplit(jblk)
              in = iend - ibegin + 1
              ! 1 x 1 block
              if (in == 1) then
                 if ((irange == allrng) .or. ((irange == valrng) .and. (d(ibegin) > vl) .and. (d( &
                           ibegin) <= vu)) .or. ((irange == indrng) .and. (iblock(wbegin) == jblk))) then
                    m = m + 1
                    w(m) = d(ibegin)
                    werr(m) = zero
                    ! the gap for a single block doesn't matter for the later
                    ! algorithm and is assigned an arbitrary large value
                    wgap(m) = zero
                    iblock(m) = jblk
                    indexw(m) = 1
                    wbegin = wbegin + 1
                 end if
                 ! e( iend ) holds the shift for the initial rrr
                 e(iend) = zero
                 ibegin = iend + 1
                 cycle loop_170
              end if
              ! blocks of size larger than 1x1
              ! e( iend ) will hold the shift for the initial rrr, for now set it =0
              e(iend) = zero
              ! find local outer bounds gl,gu for the block
              gl = d(ibegin)
              gu = d(ibegin)
              do i = ibegin,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              if (.not. ((irange == allrng) .and. (.not. forceb))) then
                 ! count the number of eigenvalues in the current block.
                 mb = 0
                 do i = wbegin,mm
                    if (iblock(i) == jblk) then
                       mb = mb + 1
                    else
                       goto 21
                    end if
                 end do
                 21 continue
                 if (mb == 0) then
                    ! no eigenvalue in the current block lies in the desired range
                    ! e( iend ) holds the shift for the initial rrr
                    e(iend) = zero
                    ibegin = iend + 1
                    cycle loop_170
                 else
                    ! decide whether dqds or bisection is more efficient
                    usedqd = ((mb > fac*in) .and. (.not. forceb))
                    wend = wbegin + mb - 1
                    ! calculate gaps for the current block
                    ! in later stages, when representations for individual
                    ! eigenvalues are different, we use sigma = e( iend ).
                    sigma = zero
                    do i = wbegin,wend - 1
                       wgap(i) = max(zero,w(i + 1) - werr(i + 1) - (w(i) + werr(i)))
                    end do
                    wgap(wend) = max(zero,vu - sigma - (w(wend) + werr(wend)))
                    ! find local index of the first and last desired evalue.
                    indl = indexw(wbegin)
                    indu = indexw(wend)
                 end if
              end if
              if (((irange == allrng) .and. (.not. forceb)) .or. usedqd) then
                 ! case of dqds
                 ! find approximations to the extremal eigenvalues of the block
                 call la_qlarrk(in,1,gl,gu,d(ibegin),e2(ibegin),pivmin,rtl,tmp,tmp1, &
                           iinfo)
                 if (iinfo /= 0) then
                    info = -1
                    return
                 end if
                 isleft = max(gl,tmp - tmp1 - hndrd*eps*abs(tmp - tmp1))
                 call la_qlarrk(in,in,gl,gu,d(ibegin),e2(ibegin),pivmin,rtl,tmp,tmp1, &
                            iinfo)
                 if (iinfo /= 0) then
                    info = -1
                    return
                 end if
                 isrght = min(gu,tmp + tmp1 + hndrd*eps*abs(tmp + tmp1))
                 ! improve the estimate of the spectral diameter
                 spdiam = isrght - isleft
              else
                 ! case of bisection
                 ! find approximations to the wanted extremal eigenvalues
                 isleft = max(gl,w(wbegin) - werr(wbegin) - hndrd*eps*abs(w(wbegin) - werr( &
                           wbegin)))
                 isrght = min(gu,w(wend) + werr(wend) + hndrd*eps*abs(w(wend) + werr(wend)))

              end if
              ! decide whether the base representation for the current block
              ! l_jblk d_jblk l_jblk^t = t_jblk - sigma_jblk i
              ! should be on the left or the right end of the current block.
              ! the strategy is to shift to the end which is "more populated"
              ! furthermore, decide whether to use dqds for the computation of
              ! dqds is chosen if all eigenvalues are desired or the number of
              ! eigenvalues to be computed is large compared to the blocksize.
              if ((irange == allrng) .and. (.not. forceb)) then
                 ! if all the eigenvalues have to be computed, we use dqd
                 usedqd = .true.
                 ! indl is the local index of the first eigenvalue to compute
                 indl = 1
                 indu = in
                 ! mb =  number of eigenvalues to compute
                 mb = in
                 wend = wbegin + mb - 1
                 ! define 1/4 and 3/4 points of the spectrum
                 s1 = isleft + fourth*spdiam
                 s2 = isrght - fourth*spdiam
              else
                 ! la_qlarrd has computed iblock and indexw for each eigenvalue
                 ! approximation.
                 ! choose sigma
                 if (usedqd) then
                    s1 = isleft + fourth*spdiam
                    s2 = isrght - fourth*spdiam
                 else
                    tmp = min(isrght,vu) - max(isleft,vl)
                    s1 = max(isleft,vl) + fourth*tmp
                    s2 = min(isrght,vu) - fourth*tmp
                 end if
              end if
              ! compute the negcount at the 1/4 and 3/4 points
              if (mb > 1) then
                 call la_qlarrc('T',in,s1,s2,d(ibegin),e(ibegin),pivmin,cnt,cnt1, &
                           cnt2,iinfo)
              end if
              if (mb == 1) then
                 sigma = gl
                 sgndef = one
              elseif (cnt1 - indl >= indu - cnt2) then
                 if ((irange == allrng) .and. (.not. forceb)) then
                    sigma = max(isleft,gl)
                 elseif (usedqd) then
                    ! use gerschgorin bound as shift to get pos def matrix
                    ! for dqds
                    sigma = isleft
                 else
                    ! use approximation of the first desired eigenvalue of the
                    ! block as shift
                    sigma = max(isleft,vl)
                 end if
                 sgndef = one
              else
                 if ((irange == allrng) .and. (.not. forceb)) then
                    sigma = min(isrght,gu)
                 elseif (usedqd) then
                    ! use gerschgorin bound as shift to get neg def matrix
                    ! for dqds
                    sigma = isrght
                 else
                    ! use approximation of the first desired eigenvalue of the
                    ! block as shift
                    sigma = min(isrght,vu)
                 end if
                 sgndef = -one
              end if
              ! an initial sigma has been chosen that will be used for computing
              ! t - sigma i = l d l^t
              ! define the increment tau of the shift in case the initial shift
              ! needs to be refined to obtain a factorization with not too much
              ! element growth.
              if (usedqd) then
                 ! the initial sigma was to the outer end of the spectrum
                 ! the matrix is definite and we need not retreat.
                 tau = spdiam*eps*n + two*pivmin
                 tau = max(tau,two*eps*abs(sigma))
              else
                 if (mb > 1) then
                    clwdth = w(wend) + werr(wend) - w(wbegin) - werr(wbegin)
                    avgap = abs(clwdth/real(wend - wbegin,KIND=qp))
                    if (sgndef == one) then
                       tau = half*max(wgap(wbegin),avgap)
                       tau = max(tau,werr(wbegin))
                    else
                       tau = half*max(wgap(wend - 1),avgap)
                       tau = max(tau,werr(wend))
                    end if
                 else
                    tau = werr(wbegin)
                 end if
              end if
              loop_80: do idum = 1,maxtry
                 ! compute l d l^t factorization of tridiagonal matrix t - sigma i.
                 ! store d in work(1:in), l in work(in+1:2*in), and reciprocals of
                 ! pivots in work(2*in+1:3*in)
                 dpivot = d(ibegin) - sigma
                 work(1) = dpivot
                 dmax = abs(work(1))
                 j = ibegin
                 do i = 1,in - 1
                    work(2*in + i) = one/work(i)
                    tmp = e(j)*work(2*in + i)
                    work(in + i) = tmp
                    dpivot = (d(j + 1) - sigma) - tmp*e(j)
                    work(i + 1) = dpivot
                    dmax = max(dmax,abs(dpivot))
                    j = j + 1
                 end do
                 ! check for element growth
                 if (dmax > maxgrowth*spdiam) then
                    norep = .true.
                 else
                    norep = .false.
                 end if
                 if (usedqd .and. .not. norep) then
                    ! ensure the definiteness of the representation
                    ! all entries of d (of l d l^t) must have the same sign
                    do i = 1,in
                       tmp = sgndef*work(i)
                       if (tmp < zero) norep = .true.
                    end do
                 end if
                 if (norep) then
                    ! note that in the case of irange=allrng, we use the gerschgorin
                    ! shift which makes the matrix definite. so we should end up
                    ! here really only in the case of irange = valrng or indrng.
                    if (idum == maxtry - 1) then
                       if (sgndef == one) then
                          ! the fudged gerschgorin shift should succeed
                          sigma = gl - fudge*spdiam*eps*n - fudge*two*pivmin
                       else
                          sigma = gu + fudge*spdiam*eps*n + fudge*two*pivmin
                       end if
                    else
                       sigma = sigma - sgndef*tau
                       tau = two*tau
                    end if
                 else
                    ! an initial rrr is found
                    go to 83
                 end if
              end do loop_80
              ! if the program reaches this point, no base representation could be
              ! found in maxtry iterations.
              info = 2
              return
              83 continue
              ! at this point, we have found an initial base representation
              ! t - sigma i = l d l^t with not too much element growth.
              ! store the shift.
              e(iend) = sigma
              ! store d and l.
              call la_qcopy(in,work,1,d(ibegin),1)
              call la_qcopy(in - 1,work(in + 1),1,e(ibegin),1)
              if (mb > 1) then
                 ! perturb each entry of the base representation by a small
                 ! (but random) relative amount to overcome difficulties with
                 ! glued matrices.
                 do i = 1,4
                    iseed(i) = 1
                 end do
                 call la_qlarnv(2,iseed,2*in - 1,work(1))
                 do i = 1,in - 1
                    d(ibegin + i - 1) = d(ibegin + i - 1)*(one + eps*pert*work(i))
                    e(ibegin + i - 1) = e(ibegin + i - 1)*(one + eps*pert*work(in + i))
                 end do
                 d(iend) = d(iend)*(one + eps*four*work(in))
              end if
              ! don't update the gerschgorin intervals because keeping track
              ! of the updates would be too much work in la_qlarrv.
              ! we update w instead and use it to locate the proper gerschgorin
              ! intervals.
              ! compute the required eigenvalues of l d l' by bisection or dqds
              if (.not. usedqd) then
                 ! if la_qlarrd has been used, shift the eigenvalue approximations
                 ! according to their representation. this is necessary for
                 ! a uniform la_qlarrv since dqds computes eigenvalues of the
                 ! shifted representation. in la_qlarrv, w will always hold the
                 ! unshifted eigenvalue approximation.
                 do j = wbegin,wend
                    w(j) = w(j) - sigma
                    werr(j) = werr(j) + abs(w(j))*eps
                 end do
                 ! call la_qlarrb to reduce eigenvalue error of the approximations
                 ! from la_qlarrd
                 do i = ibegin,iend - 1
                    work(i) = d(i)*e(i)**2
                 end do
                 ! use bisection to find ev from indl to indu
                 call la_qlarrb(in,d(ibegin),work(ibegin),indl,indu,rtol1,rtol2,indl - 1, &
                 w(wbegin),wgap(wbegin),werr(wbegin),work(2*n + 1),iwork,pivmin,spdiam,in, &
                           iinfo)
                 if (iinfo /= 0) then
                    info = -4
                    return
                 end if
                 ! la_qlarrb computes all gaps correctly except for the last one
                 ! record distance to vu/gu
                 wgap(wend) = max(zero, (vu - sigma) - (w(wend) + werr(wend)))
                 do i = indl,indu
                    m = m + 1
                    iblock(m) = jblk
                    indexw(m) = i
                 end do
              else
                 ! call dqds to get all eigs (and then possibly delete unwanted
                 ! eigenvalues).
                 ! note that dqds finds the eigenvalues of the l d l^t representation
                 ! of t to high relative accuracy. high relative accuracy
                 ! might be lost when the shift of the rrr is subtracted to obtain
                 ! the eigenvalues of t. however, t is not guaranteed to define its
                 ! eigenvalues to high relative accuracy anyway.
                 ! set rtol to the order of the tolerance used in la_qlasq2
                 ! this is an estimated error, the worst case bound is 4*n*eps
                 ! which is usually too large and requires unnecessary work to be
                 ! done by bisection when computing the eigenvectors
                 rtol = log(real(in,KIND=qp))*four*eps
                 j = ibegin
                 do i = 1,in - 1
                    work(2*i - 1) = abs(d(j))
                    work(2*i) = e(j)*e(j)*work(2*i - 1)
                    j = j + 1
                 end do
                 work(2*in - 1) = abs(d(iend))
                 work(2*in) = zero
                 call la_qlasq2(in,work,iinfo)
                 if (iinfo /= 0) then
                    ! if iinfo = -5 then an index is part of a tight cluster
                    ! and should be changed. the index is in iwork(1) and the
                    ! gap is in work(n+1)
                    info = -5
                    return
                 else
                    ! test that all eigenvalues are positive as expected
                    do i = 1,in
                       if (work(i) < zero) then
                          info = -6
                          return
                       end if
                    end do
                 end if
                 if (sgndef > zero) then
                    do i = indl,indu
                       m = m + 1
                       w(m) = work(in - i + 1)
                       iblock(m) = jblk
                       indexw(m) = i
                    end do
                 else
                    do i = indl,indu
                       m = m + 1
                       w(m) = -work(i)
                       iblock(m) = jblk
                       indexw(m) = i
                    end do
                 end if
                 do i = m - mb + 1,m
                    ! the value of rtol below should be the tolerance in la_qlasq2
                    werr(i) = rtol*abs(w(i))
                 end do
                 do i = m - mb + 1,m - 1
                    ! compute the right gap between the intervals
                    wgap(i) = max(zero,w(i + 1) - werr(i + 1) - (w(i) + werr(i)))
                 end do
                 wgap(m) = max(zero, (vu - sigma) - (w(m) + werr(m)))
              end if
              ! proceed with next block
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_qlarre

     !> CLAR1V: computes the (scaled) r-th column of the inverse of
     !> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
     !> L D L**T - sigma I. When sigma is close to an eigenvalue, the
     !> computed vector is an accurate eigenvector. Usually, r corresponds
     !> to the index where the eigenvector is largest in magnitude.
     !> The following steps accomplish this computation :
     !> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
     !> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
     !> (c) Computation of the diagonal elements of the inverse of
     !> L D L**T - sigma I by combining the above transforms, and choosing
     !> r as the index where the diagonal of the inverse is (one of the)
     !> largest in magnitude.
     !> (d) Computation of the (scaled) r-th column of the inverse using the
     !> twisted factorization obtained by combining the top part of the
     !> the stationary and the bottom part of the progressive transform.

     pure subroutine la_clar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc, &
               negcnt,ztz,mingma,r,isuppz,nrminv,resid,rqcorr,work)
        use la_constants_sp,only:zero,one,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantnc
           integer(ilp),intent(in) :: b1,bn,n
           integer(ilp),intent(out) :: negcnt
           integer(ilp),intent(inout) :: r
           real(sp),intent(in) :: gaptol,lambda,pivmin
           real(sp),intent(out) :: mingma,nrminv,resid,rqcorr,ztz
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*)
           real(sp),intent(in) :: d(*),l(*),ld(*),lld(*)
           real(sp),intent(out) :: work(*)
           complex(sp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: sawnan1,sawnan2
           integer(ilp) :: i,indlpl,indp,inds,indumn,neg1,neg2,r1,r2
           real(sp) :: dminus,dplus,eps,s,tmp
           ! Intrinsic Functions
           intrinsic :: abs,real
           ! Executable Statements
           eps = la_slamch('PRECISION')
           if (r == 0) then
              r1 = b1
              r2 = bn
           else
              r1 = r
              r2 = r
           end if
           ! storage for lplus
           indlpl = 0
           ! storage for uminus
           indumn = n
           inds = 2*n + 1
           indp = 3*n + 1
           if (b1 == 1) then
              work(inds) = zero
           else
              work(inds + b1 - 1) = lld(b1 - 1)
           end if
           ! compute the stationary transform (using the differential form)
           ! until the index r2.
           sawnan1 = .false.
           neg1 = 0
           s = work(inds + b1 - 1) - lambda
           do i = b1,r1 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              if (dplus < zero) neg1 = neg1 + 1
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_sisnan(s)
           if (sawnan1) goto 60
           do i = r1,r2 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_sisnan(s)
           60 continue
           if (sawnan1) then
              ! runs a slower version of the above loop if a nan is detected
              neg1 = 0
              s = work(inds + b1 - 1) - lambda
              do i = b1,r1 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 if (dplus < zero) neg1 = neg1 + 1
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
              do i = r1,r2 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
           end if
           ! compute the progressive transform (using the differential form)
           ! until the index r1
           sawnan2 = .false.
           neg2 = 0
           work(indp + bn - 1) = d(bn) - lambda
           do i = bn - 1,r1,-1
              dminus = lld(i) + work(indp + i)
              tmp = d(i)/dminus
              if (dminus < zero) neg2 = neg2 + 1
              work(indumn + i) = l(i)*tmp
              work(indp + i - 1) = work(indp + i)*tmp - lambda
           end do
           tmp = work(indp + r1 - 1)
           sawnan2 = la_sisnan(tmp)
           if (sawnan2) then
              ! runs a slower version of the above loop if a nan is detected
              neg2 = 0
              do i = bn - 1,r1,-1
                 dminus = lld(i) + work(indp + i)
                 if (abs(dminus) < pivmin) dminus = -pivmin
                 tmp = d(i)/dminus
                 if (dminus < zero) neg2 = neg2 + 1
                 work(indumn + i) = l(i)*tmp
                 work(indp + i - 1) = work(indp + i)*tmp - lambda
                 if (tmp == zero) work(indp + i - 1) = d(i) - lambda
              end do
           end if
           ! find the index (from r1 to r2) of the largest (in magnitude)
           ! diagonal element of the inverse
           mingma = work(inds + r1 - 1) + work(indp + r1 - 1)
           if (mingma < zero) neg1 = neg1 + 1
           if (wantnc) then
              negcnt = neg1 + neg2
           else
              negcnt = -1
           end if
           if (abs(mingma) == zero) mingma = eps*work(inds + r1 - 1)
           r = r1
           do i = r1,r2 - 1
              tmp = work(inds + i) + work(indp + i)
              if (tmp == zero) tmp = eps*work(inds + i)
              if (abs(tmp) <= abs(mingma)) then
                 mingma = tmp
                 r = i + 1
              end if
           end do
           ! compute the fp vector: solve n^t v = e_r
           isuppz(1) = b1
           isuppz(2) = bn
           z(r) = cone
           ztz = one
           ! compute the fp vector upwards from r
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r - 1,b1,-1
                 z(i) = -(work(indlpl + i)*z(i + 1))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    goto 220
                 end if
                 ztz = ztz + real(z(i)*z(i),KIND=sp)
              end do
              220 continue
           else
              ! run slower loop if nan occurred.
              do i = r - 1,b1,-1
                 if (z(i + 1) == zero) then
                    z(i) = -(ld(i + 1)/ld(i))*z(i + 2)
                 else
                    z(i) = -(work(indlpl + i)*z(i + 1))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    go to 240
                 end if
                 ztz = ztz + real(z(i)*z(i),KIND=sp)
              end do
              240 continue
           end if
           ! compute the fp vector downwards from r in blocks of size blksiz
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r,bn - 1
                 z(i + 1) = -(work(indumn + i)*z(i))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 260
                 end if
                 ztz = ztz + real(z(i + 1)*z(i + 1),KIND=sp)
              end do
              260 continue
           else
              ! run slower loop if nan occurred.
              do i = r,bn - 1
                 if (z(i) == zero) then
                    z(i + 1) = -(ld(i - 1)/ld(i))*z(i - 1)
                 else
                    z(i + 1) = -(work(indumn + i)*z(i))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 280
                 end if
                 ztz = ztz + real(z(i + 1)*z(i + 1),KIND=sp)
              end do
              280 continue
           end if
           ! compute quantities for convergence test
           tmp = one/ztz
           nrminv = sqrt(tmp)
           resid = abs(mingma)*nrminv
           rqcorr = mingma*tmp
           return
     end subroutine la_clar1v
     !> ZLAR1V: computes the (scaled) r-th column of the inverse of
     !> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
     !> L D L**T - sigma I. When sigma is close to an eigenvalue, the
     !> computed vector is an accurate eigenvector. Usually, r corresponds
     !> to the index where the eigenvector is largest in magnitude.
     !> The following steps accomplish this computation :
     !> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
     !> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
     !> (c) Computation of the diagonal elements of the inverse of
     !> L D L**T - sigma I by combining the above transforms, and choosing
     !> r as the index where the diagonal of the inverse is (one of the)
     !> largest in magnitude.
     !> (d) Computation of the (scaled) r-th column of the inverse using the
     !> twisted factorization obtained by combining the top part of the
     !> the stationary and the bottom part of the progressive transform.

     pure subroutine la_zlar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc, &
               negcnt,ztz,mingma,r,isuppz,nrminv,resid,rqcorr,work)
        use la_constants_dp,only:zero,one,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantnc
           integer(ilp),intent(in) :: b1,bn,n
           integer(ilp),intent(out) :: negcnt
           integer(ilp),intent(inout) :: r
           real(dp),intent(in) :: gaptol,lambda,pivmin
           real(dp),intent(out) :: mingma,nrminv,resid,rqcorr,ztz
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*)
           real(dp),intent(in) :: d(*),l(*),ld(*),lld(*)
           real(dp),intent(out) :: work(*)
           complex(dp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: sawnan1,sawnan2
           integer(ilp) :: i,indlpl,indp,inds,indumn,neg1,neg2,r1,r2
           real(dp) :: dminus,dplus,eps,s,tmp
           ! Intrinsic Functions
           intrinsic :: abs,real
           ! Executable Statements
           eps = la_dlamch('PRECISION')
           if (r == 0) then
              r1 = b1
              r2 = bn
           else
              r1 = r
              r2 = r
           end if
           ! storage for lplus
           indlpl = 0
           ! storage for uminus
           indumn = n
           inds = 2*n + 1
           indp = 3*n + 1
           if (b1 == 1) then
              work(inds) = zero
           else
              work(inds + b1 - 1) = lld(b1 - 1)
           end if
           ! compute the stationary transform (using the differential form)
           ! until the index r2.
           sawnan1 = .false.
           neg1 = 0
           s = work(inds + b1 - 1) - lambda
           do i = b1,r1 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              if (dplus < zero) neg1 = neg1 + 1
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_disnan(s)
           if (sawnan1) goto 60
           do i = r1,r2 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_disnan(s)
           60 continue
           if (sawnan1) then
              ! runs a slower version of the above loop if a nan is detected
              neg1 = 0
              s = work(inds + b1 - 1) - lambda
              do i = b1,r1 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 if (dplus < zero) neg1 = neg1 + 1
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
              do i = r1,r2 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
           end if
           ! compute the progressive transform (using the differential form)
           ! until the index r1
           sawnan2 = .false.
           neg2 = 0
           work(indp + bn - 1) = d(bn) - lambda
           do i = bn - 1,r1,-1
              dminus = lld(i) + work(indp + i)
              tmp = d(i)/dminus
              if (dminus < zero) neg2 = neg2 + 1
              work(indumn + i) = l(i)*tmp
              work(indp + i - 1) = work(indp + i)*tmp - lambda
           end do
           tmp = work(indp + r1 - 1)
           sawnan2 = la_disnan(tmp)
           if (sawnan2) then
              ! runs a slower version of the above loop if a nan is detected
              neg2 = 0
              do i = bn - 1,r1,-1
                 dminus = lld(i) + work(indp + i)
                 if (abs(dminus) < pivmin) dminus = -pivmin
                 tmp = d(i)/dminus
                 if (dminus < zero) neg2 = neg2 + 1
                 work(indumn + i) = l(i)*tmp
                 work(indp + i - 1) = work(indp + i)*tmp - lambda
                 if (tmp == zero) work(indp + i - 1) = d(i) - lambda
              end do
           end if
           ! find the index (from r1 to r2) of the largest (in magnitude)
           ! diagonal element of the inverse
           mingma = work(inds + r1 - 1) + work(indp + r1 - 1)
           if (mingma < zero) neg1 = neg1 + 1
           if (wantnc) then
              negcnt = neg1 + neg2
           else
              negcnt = -1
           end if
           if (abs(mingma) == zero) mingma = eps*work(inds + r1 - 1)
           r = r1
           do i = r1,r2 - 1
              tmp = work(inds + i) + work(indp + i)
              if (tmp == zero) tmp = eps*work(inds + i)
              if (abs(tmp) <= abs(mingma)) then
                 mingma = tmp
                 r = i + 1
              end if
           end do
           ! compute the fp vector: solve n^t v = e_r
           isuppz(1) = b1
           isuppz(2) = bn
           z(r) = cone
           ztz = one
           ! compute the fp vector upwards from r
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r - 1,b1,-1
                 z(i) = -(work(indlpl + i)*z(i + 1))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    goto 220
                 end if
                 ztz = ztz + real(z(i)*z(i),KIND=dp)
              end do
              220 continue
           else
              ! run slower loop if nan occurred.
              do i = r - 1,b1,-1
                 if (z(i + 1) == zero) then
                    z(i) = -(ld(i + 1)/ld(i))*z(i + 2)
                 else
                    z(i) = -(work(indlpl + i)*z(i + 1))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    go to 240
                 end if
                 ztz = ztz + real(z(i)*z(i),KIND=dp)
              end do
              240 continue
           end if
           ! compute the fp vector downwards from r in blocks of size blksiz
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r,bn - 1
                 z(i + 1) = -(work(indumn + i)*z(i))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 260
                 end if
                 ztz = ztz + real(z(i + 1)*z(i + 1),KIND=dp)
              end do
              260 continue
           else
              ! run slower loop if nan occurred.
              do i = r,bn - 1
                 if (z(i) == zero) then
                    z(i + 1) = -(ld(i - 1)/ld(i))*z(i - 1)
                 else
                    z(i + 1) = -(work(indumn + i)*z(i))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 280
                 end if
                 ztz = ztz + real(z(i + 1)*z(i + 1),KIND=dp)
              end do
              280 continue
           end if
           ! compute quantities for convergence test
           tmp = one/ztz
           nrminv = sqrt(tmp)
           resid = abs(mingma)*nrminv
           rqcorr = mingma*tmp
           return
     end subroutine la_zlar1v
     !> WLAR1V: computes the (scaled) r-th column of the inverse of
     !> the sumbmatrix in rows B1 through BN of the tridiagonal matrix
     !> L D L**T - sigma I. When sigma is close to an eigenvalue, the
     !> computed vector is an accurate eigenvector. Usually, r corresponds
     !> to the index where the eigenvector is largest in magnitude.
     !> The following steps accomplish this computation :
     !> (a) Stationary qd transform,  L D L**T - sigma I = L(+) D(+) L(+)**T,
     !> (b) Progressive qd transform, L D L**T - sigma I = U(-) D(-) U(-)**T,
     !> (c) Computation of the diagonal elements of the inverse of
     !> L D L**T - sigma I by combining the above transforms, and choosing
     !> r as the index where the diagonal of the inverse is (one of the)
     !> largest in magnitude.
     !> (d) Computation of the (scaled) r-th column of the inverse using the
     !> twisted factorization obtained by combining the top part of the
     !> the stationary and the bottom part of the progressive transform.

     pure subroutine la_wlar1v(n,b1,bn,lambda,d,l,ld,lld,pivmin,gaptol,z,wantnc, &
               negcnt,ztz,mingma,r,isuppz,nrminv,resid,rqcorr,work)
        use la_constants_qp,only:zero,one,cone
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           logical(lk),intent(in) :: wantnc
           integer(ilp),intent(in) :: b1,bn,n
           integer(ilp),intent(out) :: negcnt
           integer(ilp),intent(inout) :: r
           real(qp),intent(in) :: gaptol,lambda,pivmin
           real(qp),intent(out) :: mingma,nrminv,resid,rqcorr,ztz
           ! Array Arguments
           integer(ilp),intent(out) :: isuppz(*)
           real(qp),intent(in) :: d(*),l(*),ld(*),lld(*)
           real(qp),intent(out) :: work(*)
           complex(qp),intent(inout) :: z(*)
        ! =====================================================================

           ! Local Scalars
           logical(lk) :: sawnan1,sawnan2
           integer(ilp) :: i,indlpl,indp,inds,indumn,neg1,neg2,r1,r2
           real(qp) :: dminus,dplus,eps,s,tmp
           ! Intrinsic Functions
           intrinsic :: abs,real
           ! Executable Statements
           eps = la_qlamch('PRECISION')
           if (r == 0) then
              r1 = b1
              r2 = bn
           else
              r1 = r
              r2 = r
           end if
           ! storage for lplus
           indlpl = 0
           ! storage for uminus
           indumn = n
           inds = 2*n + 1
           indp = 3*n + 1
           if (b1 == 1) then
              work(inds) = zero
           else
              work(inds + b1 - 1) = lld(b1 - 1)
           end if
           ! compute the stationary transform (using the differential form)
           ! until the index r2.
           sawnan1 = .false.
           neg1 = 0
           s = work(inds + b1 - 1) - lambda
           do i = b1,r1 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              if (dplus < zero) neg1 = neg1 + 1
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_qisnan(s)
           if (sawnan1) goto 60
           do i = r1,r2 - 1
              dplus = d(i) + s
              work(indlpl + i) = ld(i)/dplus
              work(inds + i) = s*work(indlpl + i)*l(i)
              s = work(inds + i) - lambda
           end do
           sawnan1 = la_qisnan(s)
           60 continue
           if (sawnan1) then
              ! runs a slower version of the above loop if a nan is detected
              neg1 = 0
              s = work(inds + b1 - 1) - lambda
              do i = b1,r1 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 if (dplus < zero) neg1 = neg1 + 1
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
              do i = r1,r2 - 1
                 dplus = d(i) + s
                 if (abs(dplus) < pivmin) dplus = -pivmin
                 work(indlpl + i) = ld(i)/dplus
                 work(inds + i) = s*work(indlpl + i)*l(i)
                 if (work(indlpl + i) == zero) work(inds + i) = lld(i)
                 s = work(inds + i) - lambda
              end do
           end if
           ! compute the progressive transform (using the differential form)
           ! until the index r1
           sawnan2 = .false.
           neg2 = 0
           work(indp + bn - 1) = d(bn) - lambda
           do i = bn - 1,r1,-1
              dminus = lld(i) + work(indp + i)
              tmp = d(i)/dminus
              if (dminus < zero) neg2 = neg2 + 1
              work(indumn + i) = l(i)*tmp
              work(indp + i - 1) = work(indp + i)*tmp - lambda
           end do
           tmp = work(indp + r1 - 1)
           sawnan2 = la_qisnan(tmp)
           if (sawnan2) then
              ! runs a slower version of the above loop if a nan is detected
              neg2 = 0
              do i = bn - 1,r1,-1
                 dminus = lld(i) + work(indp + i)
                 if (abs(dminus) < pivmin) dminus = -pivmin
                 tmp = d(i)/dminus
                 if (dminus < zero) neg2 = neg2 + 1
                 work(indumn + i) = l(i)*tmp
                 work(indp + i - 1) = work(indp + i)*tmp - lambda
                 if (tmp == zero) work(indp + i - 1) = d(i) - lambda
              end do
           end if
           ! find the index (from r1 to r2) of the largest (in magnitude)
           ! diagonal element of the inverse
           mingma = work(inds + r1 - 1) + work(indp + r1 - 1)
           if (mingma < zero) neg1 = neg1 + 1
           if (wantnc) then
              negcnt = neg1 + neg2
           else
              negcnt = -1
           end if
           if (abs(mingma) == zero) mingma = eps*work(inds + r1 - 1)
           r = r1
           do i = r1,r2 - 1
              tmp = work(inds + i) + work(indp + i)
              if (tmp == zero) tmp = eps*work(inds + i)
              if (abs(tmp) <= abs(mingma)) then
                 mingma = tmp
                 r = i + 1
              end if
           end do
           ! compute the fp vector: solve n^t v = e_r
           isuppz(1) = b1
           isuppz(2) = bn
           z(r) = cone
           ztz = one
           ! compute the fp vector upwards from r
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r - 1,b1,-1
                 z(i) = -(work(indlpl + i)*z(i + 1))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    goto 220
                 end if
                 ztz = ztz + real(z(i)*z(i),KIND=qp)
              end do
              220 continue
           else
              ! run slower loop if nan occurred.
              do i = r - 1,b1,-1
                 if (z(i + 1) == zero) then
                    z(i) = -(ld(i + 1)/ld(i))*z(i + 2)
                 else
                    z(i) = -(work(indlpl + i)*z(i + 1))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i) = zero
                    isuppz(1) = i + 1
                    go to 240
                 end if
                 ztz = ztz + real(z(i)*z(i),KIND=qp)
              end do
              240 continue
           end if
           ! compute the fp vector downwards from r in blocks of size blksiz
           if (.not. sawnan1 .and. .not. sawnan2) then
              do i = r,bn - 1
                 z(i + 1) = -(work(indumn + i)*z(i))
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 260
                 end if
                 ztz = ztz + real(z(i + 1)*z(i + 1),KIND=qp)
              end do
              260 continue
           else
              ! run slower loop if nan occurred.
              do i = r,bn - 1
                 if (z(i) == zero) then
                    z(i + 1) = -(ld(i - 1)/ld(i))*z(i - 1)
                 else
                    z(i + 1) = -(work(indumn + i)*z(i))
                 end if
                 if ((abs(z(i)) + abs(z(i + 1)))*abs(ld(i)) < gaptol) then
                    z(i + 1) = zero
                    isuppz(2) = i
                    go to 280
                 end if
                 ztz = ztz + real(z(i + 1)*z(i + 1),KIND=qp)
              end do
              280 continue
           end if
           ! compute quantities for convergence test
           tmp = one/ztz
           nrminv = sqrt(tmp)
           resid = abs(mingma)*nrminv
           rqcorr = mingma*tmp
           return
     end subroutine la_wlar1v

     !> CLARRV: computes the eigenvectors of the tridiagonal matrix
     !> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
     !> The input eigenvalues should have been computed by SLARRE.

     pure subroutine la_clarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1, &
               rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz,isuppz,work,iwork,info)
        use la_constants_sp,only:zero,half,one,two,three,four,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dol,dou,ldz,m,n
           integer(ilp),intent(out) :: info
           real(sp),intent(in) :: minrgp,pivmin,vl,vu
           real(sp),intent(inout) :: rtol1,rtol2
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),indexw(*),isplit(*)
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(sp),intent(inout) :: d(*),l(*),w(*),werr(*),wgap(*)
           real(sp),intent(in) :: gers(*)
           real(sp),intent(out) :: work(*)
           complex(sp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 10

           ! Local Scalars
           logical(lk) :: eskip,needbs,stp2ii,tryrqc,usedbs,usedrq
           integer(ilp) :: done,i,ibegin,idone,iend,ii,iindc1,iindc2,iindr,iindwk,iinfo, &
            im,in,indeig,indld,indlld,indwrk,isupmn,isupmx,iter,itmp1,j,jblk,k, &
            miniwsize,minwsize,nclus,ndepth,negcnt,newcls,newfst,newftt,newlst,newsiz, &
            offset,oldcls,oldfst,oldien,oldlst,oldncl,p,parity,q,wbegin,wend,windex, &
                      windmn,windpl,zfrom,zto,zusedl,zusedu,zusedw
           integer(ilp) :: indin1,indin2
           real(sp) :: bstres,bstw,eps,fudge,gap,gaptol,gl,gu,lambda,left,lgap,mingma, &
           nrminv,resid,rgap,right,rqcorr,rqtol,savgap,sgndef,sigma,spdiam,ssigma,tau, &
                     tmp,tol,ztz
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           intrinsic :: cmplx
           ! Executable Statements
           info = 0
           ! quick return if possible
           if ((n <= 0) .or. (m <= 0)) then
              return
           end if
           ! the first n entries of work are reserved for the eigenvalues
           indld = n + 1
           indlld = 2*n + 1
           indin1 = 3*n + 1
           indin2 = 4*n + 1
           indwrk = 5*n + 1
           minwsize = 12*n
           do i = 1,minwsize
              work(i) = zero
           end do
           ! iwork(iindr+1:iindr+n) hold the twist indices r for the
           ! factorization used to compute the fp vector
           iindr = 0
           ! iwork(iindc1+1:iinc2+n) are used to store the clusters of the current
           ! layer and the one above.
           iindc1 = n
           iindc2 = 2*n
           iindwk = 3*n + 1
           miniwsize = 7*n
           do i = 1,miniwsize
              iwork(i) = 0
           end do
           zusedl = 1
           if (dol > 1) then
              ! set lower bound for use of z
              zusedl = dol - 1
           end if
           zusedu = m
           if (dou < m) then
              ! set lower bound for use of z
              zusedu = dou + 1
           end if
           ! the width of the part of z that is used
           zusedw = zusedu - zusedl + 1
           call la_claset('FULL',n,zusedw,czero,czero,z(1,zusedl),ldz)
           eps = la_slamch('PRECISION')
           rqtol = two*eps
           ! set expert flags for standard code.
           tryrqc = .true.
           if ((dol == 1) .and. (dou == m)) then
           else
              ! only selected eigenpairs are computed. since the other evalues
              ! are not refined by rq iteration, bisection has to compute to full
              ! accuracy.
              rtol1 = four*eps
              rtol2 = four*eps
           end if
           ! the entries wbegin:wend in w, werr, wgap correspond to the
           ! desired eigenvalues. the support of the nonzero eigenvector
           ! entries is contained in the interval ibegin:iend.
           ! remark that if k eigenpairs are desired, then the eigenvectors
           ! are stored in k contiguous columns of z.
           ! done is the number of eigenvectors already computed
           done = 0
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,iblock(m)
              iend = isplit(jblk)
              sigma = l(iend)
              ! find the eigenvectors of the submatrix indexed ibegin
              ! through iend.
              wend = wbegin - 1
              15 continue
              if (wend < m) then
                 if (iblock(wend + 1) == jblk) then
                    wend = wend + 1
                    go to 15
                 end if
              end if
              if (wend < wbegin) then
                 ibegin = iend + 1
                 cycle loop_170
              elseif ((wend < dol) .or. (wbegin > dou)) then
                 ibegin = iend + 1
                 wbegin = wend + 1
                 cycle loop_170
              end if
              ! find local spectral diameter of the block
              gl = gers(2*ibegin - 1)
              gu = gers(2*ibegin)
              do i = ibegin + 1,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              ! oldien is the last index of the previous block
              oldien = ibegin - 1
              ! calculate the size of the current block
              in = iend - ibegin + 1
              ! the number of eigenvalues in the current block
              im = wend - wbegin + 1
              ! this is for a 1x1 block
              if (ibegin == iend) then
                 done = done + 1
                 z(ibegin,wbegin) = cmplx(one,zero,KIND=sp)
                 isuppz(2*wbegin - 1) = ibegin
                 isuppz(2*wbegin) = ibegin
                 w(wbegin) = w(wbegin) + sigma
                 work(wbegin) = w(wbegin)
                 ibegin = iend + 1
                 wbegin = wbegin + 1
                 cycle loop_170
              end if
              ! the desired (shifted) eigenvalues are stored in w(wbegin:wend)
              ! note that these can be approximations, in this case, the corresp.
              ! entries of werr give the size of the uncertainty interval.
              ! the eigenvalue approximations will be refined when necessary as
              ! high relative accuracy is required for the computation of the
              ! corresponding eigenvectors.
              call la_scopy(im,w(wbegin),1,work(wbegin),1)
              ! we store in w the eigenvalue approximations w.r.t. the original
              ! matrix t.
              do i = 1,im
                 w(wbegin + i - 1) = w(wbegin + i - 1) + sigma
              end do
              ! ndepth is the current depth of the representation tree
              ndepth = 0
              ! parity is either 1 or 0
              parity = 1
              ! nclus is the number of clusters for the next level of the
              ! representation tree, we start with nclus = 1 for the root
              nclus = 1
              iwork(iindc1 + 1) = 1
              iwork(iindc1 + 2) = im
              ! idone is the number of eigenvectors already computed in the current
              ! block
              idone = 0
              ! loop while( idone<im )
              ! generate the representation tree for the current block and
              ! compute the eigenvectors
              40 continue
              if (idone < im) then
                 ! this is a crude protection against infinitely deep trees
                 if (ndepth > m) then
                    info = -2
                    return
                 end if
                 ! breadth first processing of the current level of the representation
                 ! tree: oldncl = number of clusters on current level
                 oldncl = nclus
                 ! reset nclus to count the number of child clusters
                 nclus = 0
                 parity = 1 - parity
                 if (parity == 0) then
                    oldcls = iindc1
                    newcls = iindc2
                 else
                    oldcls = iindc2
                    newcls = iindc1
                 end if
                 ! process the clusters on the current level
                 loop_150: do i = 1,oldncl
                    j = oldcls + 2*i
                    ! oldfst, oldlst = first, last index of current cluster.
                                     ! cluster indices start with 1 and are relative
                                     ! to wbegin when accessing w, wgap, werr, z
                    oldfst = iwork(j - 1)
                    oldlst = iwork(j)
                    if (ndepth > 0) then
                       ! retrieve relatively robust representation (rrr) of cluster
                       ! that has been computed at the previous level
                       ! the rrr is stored in z and overwritten once the eigenvectors
                       ! have been computed or when the cluster is refined
                       if ((dol == 1) .and. (dou == m)) then
                          ! get representation from location of the leftmost evalue
                          ! of the cluster
                          j = wbegin + oldfst - 1
                       else
                          if (wbegin + oldfst - 1 < dol) then
                             ! get representation from the left end of z array
                             j = dol - 1
                          elseif (wbegin + oldfst - 1 > dou) then
                             ! get representation from the right end of z array
                             j = dou
                          else
                             j = wbegin + oldfst - 1
                          end if
                       end if
                       do k = 1,in - 1
                          d(ibegin + k - 1) = real(z(ibegin + k - 1,j),KIND=sp)
                          l(ibegin + k - 1) = real(z(ibegin + k - 1,j + 1),KIND=sp)
                       end do
                       d(iend) = real(z(iend,j),KIND=sp)
                       sigma = real(z(iend,j + 1),KIND=sp)
                       ! set the corresponding entries in z to zero
                       call la_claset('FULL',in,2,czero,czero,z(ibegin,j),ldz)

                    end if
                    ! compute dl and dll of current rrr
                    do j = ibegin,iend - 1
                       tmp = d(j)*l(j)
                       work(indld - 1 + j) = tmp
                       work(indlld - 1 + j) = tmp*l(j)
                    end do
                    if (ndepth > 0) then
                       ! p and q are index of the first and last eigenvalue to compute
                       ! within the current block
                       p = indexw(wbegin - 1 + oldfst)
                       q = indexw(wbegin - 1 + oldlst)
                       ! offset for the arrays work, wgap and werr, i.e., the p-offset
                       ! through the q-offset elements of these arrays are to be used.
                        ! offset = p-oldfst
                       offset = indexw(wbegin) - 1
                       ! perform limited bisection (if necessary) to get approximate
                       ! eigenvalues to the precision needed.
                       call la_slarrb(in,d(ibegin),work(indlld + ibegin - 1),p,q,rtol1, &
                       rtol2,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk),iwork( &
                                  iindwk),pivmin,spdiam,in,iinfo)
                       if (iinfo /= 0) then
                          info = -1
                          return
                       end if
                       ! we also recompute the extremal gaps. w holds all eigenvalues
                       ! of the unshifted matrix and must be used for computation
                       ! of wgap, the entries of work might stem from rrrs with
                       ! different shifts. the gaps from wbegin-1+oldfst to
                       ! wbegin-1+oldlst are correctly computed in la_slarrb.
                       ! however, we only allow the gaps to become greater since
                       ! this is what should happen when we decrease werr
                       if (oldfst > 1) then
                          wgap(wbegin + oldfst - 2) = max(wgap(wbegin + oldfst - 2),w(wbegin + oldfst - 1) - &
                          werr(wbegin + oldfst - 1) - w(wbegin + oldfst - 2) - werr(wbegin + oldfst - 2))

                       end if
                       if (wbegin + oldlst - 1 < wend) then
                          wgap(wbegin + oldlst - 1) = max(wgap(wbegin + oldlst - 1),w(wbegin + oldlst) - &
                                    werr(wbegin + oldlst) - w(wbegin + oldlst - 1) - werr(wbegin + oldlst - 1))
                       end if
                       ! each time the eigenvalues in work get refined, we store
                       ! the newly found approximation with all shifts applied in w
                       do j = oldfst,oldlst
                          w(wbegin + j - 1) = work(wbegin + j - 1) + sigma
                       end do
                    end if
                    ! process the current node.
                    newfst = oldfst
                    loop_140: do j = oldfst,oldlst
                       if (j == oldlst) then
                          ! we are at the right end of the cluster, this is also the
                          ! boundary of the child cluster
                          newlst = j
                       else if (wgap(wbegin + j - 1) >= minrgp*abs(work(wbegin + j - 1))) &
                                 then
                          ! the right relative gap is big enough, the child cluster
                          ! (newfst,..,newlst) is well separated from the following
                          newlst = j
                        else
                          ! inside a child cluster, the relative gap is not
                          ! big enough.
                          cycle loop_140
                       end if
                       ! compute size of child cluster found
                       newsiz = newlst - newfst + 1
                       ! newftt is the place in z where the new rrr or the computed
                       ! eigenvector is to be stored
                       if ((dol == 1) .and. (dou == m)) then
                          ! store representation at location of the leftmost evalue
                          ! of the cluster
                          newftt = wbegin + newfst - 1
                       else
                          if (wbegin + newfst - 1 < dol) then
                             ! store representation at the left end of z array
                             newftt = dol - 1
                          elseif (wbegin + newfst - 1 > dou) then
                             ! store representation at the right end of z array
                             newftt = dou
                          else
                             newftt = wbegin + newfst - 1
                          end if
                       end if
                       if (newsiz > 1) then
                          ! current child is not a singleton but a cluster.
                          ! compute and store new representation of child.
                          ! compute left and right cluster gap.
                          ! lgap and rgap are not computed from work because
                          ! the eigenvalue approximations may stem from rrrs
                          ! different shifts. however, w hold all eigenvalues
                          ! of the unshifted matrix. still, the entries in wgap
                          ! have to be computed from work since the entries
                          ! in w might be of the same order so that gaps are not
                          ! exhibited correctly for very close eigenvalues.
                          if (newfst == 1) then
                             lgap = max(zero,w(wbegin) - werr(wbegin) - vl)
                         else
                             lgap = wgap(wbegin + newfst - 2)
                          end if
                          rgap = wgap(wbegin + newlst - 1)
                          ! compute left- and rightmost eigenvalue of child
                          ! to high precision in order to shift as close
                          ! as possible and obtain as large relative gaps
                          ! as possible
                          do k = 1,2
                             if (k == 1) then
                                p = indexw(wbegin - 1 + newfst)
                             else
                                p = indexw(wbegin - 1 + newlst)
                             end if
                             offset = indexw(wbegin) - 1
                             call la_slarrb(in,d(ibegin),work(indlld + ibegin - 1),p,p,rqtol, &
                             rqtol,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk), &
                                       iwork(iindwk),pivmin,spdiam,in,iinfo)
                          end do
                          if ((wbegin + newlst - 1 < dol) .or. (wbegin + newfst - 1 > dou)) then
                             ! if the cluster contains no desired eigenvalues
                             ! skip the computation of that branch of the rep. tree
                             ! we could skip before the refinement of the extremal
                             ! eigenvalues of the child, but then the representation
                             ! tree could be different from the one when nothing is
                             ! skipped. for this reason we skip at this place.
                             idone = idone + newlst - newfst + 1
                             goto 139
                          end if
                          ! compute rrr of child cluster.
                          ! note that the new rrr is stored in z
                          ! la_slarrf needs lwork = 2*n
                          call la_slarrf(in,d(ibegin),l(ibegin),work(indld + ibegin - 1), &
                          newfst,newlst,work(wbegin),wgap(wbegin),werr(wbegin),spdiam,lgap, &
                          rgap,pivmin,tau,work(indin1),work(indin2),work(indwrk),iinfo)

                          ! in the complex case, la_slarrf cannot write
                          ! the new rrr directly into z and needs an intermediate
                          ! workspace
                          do k = 1,in - 1
                             z(ibegin + k - 1,newftt) = cmplx(work(indin1 + k - 1),zero,KIND=sp)

                             z(ibegin + k - 1,newftt + 1) = cmplx(work(indin2 + k - 1),zero,KIND=sp)

                          end do
                          z(iend,newftt) = cmplx(work(indin1 + in - 1),zero,KIND=sp)
                          if (iinfo == 0) then
                             ! a new rrr for the cluster was found by la_slarrf
                             ! update shift and store it
                             ssigma = sigma + tau
                             z(iend,newftt + 1) = cmplx(ssigma,zero,KIND=sp)
                             ! work() are the midpoints and werr() the semi-width
                             ! note that the entries in w are unchanged.
                             do k = newfst,newlst
                                fudge = three*eps*abs(work(wbegin + k - 1))
                                work(wbegin + k - 1) = work(wbegin + k - 1) - tau
                                fudge = fudge + four*eps*abs(work(wbegin + k - 1))
                                ! fudge errors
                                werr(wbegin + k - 1) = werr(wbegin + k - 1) + fudge
                                ! gaps are not fudged. provided that werr is small
                                ! when eigenvalues are close, a zero gap indicates
                                ! that a new representation is needed for resolving
                                ! the cluster. a fudge could lead to a wrong decision
                                ! of judging eigenvalues 'separated' which in
                                ! reality are not. this could have a negative impact
                                ! on the orthogonality of the computed eigenvectors.
                             end do
                             nclus = nclus + 1
                             k = newcls + 2*nclus
                             iwork(k - 1) = newfst
                             iwork(k) = newlst
                          else
                             info = -2
                             return
                          end if
                       else
                          ! compute eigenvector of singleton
                          iter = 0
                          tol = four*log(real(in,KIND=sp))*eps
                          k = newfst
                          windex = wbegin + k - 1
                          windmn = max(windex - 1,1)
                          windpl = min(windex + 1,m)
                          lambda = work(windex)
                          done = done + 1
                          ! check if eigenvector computation is to be skipped
                          if ((windex < dol) .or. (windex > dou)) then
                             eskip = .true.
                             goto 125
                          else
                             eskip = .false.
                          end if
                          left = work(windex) - werr(windex)
                          right = work(windex) + werr(windex)
                          indeig = indexw(windex)
                          ! note that since we compute the eigenpairs for a child,
                          ! all eigenvalue approximations are w.r.t the same shift.
                          ! in this case, the entries in work should be used for
                          ! computing the gaps since they exhibit even very small
                          ! differences in the eigenvalues, as opposed to the
                          ! entries in w which might "look" the same.
                          if (k == 1) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vl, the formula
                             ! lgap = max( zero, (sigma - vl) + lambda )
                             ! can lead to an overestimation of the left gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small left gap.
                             lgap = eps*max(abs(left),abs(right))
                          else
                             lgap = wgap(windmn)
                          end if
                          if (k == im) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vu, the formula
                             ! can lead to an overestimation of the right gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small right gap.
                             rgap = eps*max(abs(left),abs(right))
                          else
                             rgap = wgap(windex)
                          end if
                          gap = min(lgap,rgap)
                          if ((k == 1) .or. (k == im)) then
                             ! the eigenvector support can become wrong
                             ! because significant entries could be cut off due to a
                             ! large gaptol parameter in lar1v. prevent this.
                             gaptol = zero
                          else
                             gaptol = gap*eps
                          end if
                          isupmn = in
                          isupmx = 1
                          ! update wgap so that it holds the minimum gap
                          ! to the left or the right. this is crucial in the
                          ! case where bisection is used to ensure that the
                          ! eigenvalue is refined up to the required precision.
                          ! the correct value is restored afterwards.
                          savgap = wgap(windex)
                          wgap(windex) = gap
                          ! we want to use the rayleigh quotient correction
                          ! as often as possible since it converges quadratically
                          ! when we are close enough to the desired eigenvalue.
                          ! however, the rayleigh quotient can have the wrong sign
                          ! and lead us away from the desired eigenvalue. in this
                          ! case, the best we can do is to use bisection.
                          usedbs = .false.
                          usedrq = .false.
                          ! bisection is initially turned off unless it is forced
                          needbs = .not. tryrqc
                          120 continue
                          ! check if bisection should be used to refine eigenvalue
                          if (needbs) then
                             ! take the bisection as new iterate
                             usedbs = .true.
                             itmp1 = iwork(iindr + windex)
                             offset = indexw(wbegin) - 1
                             call la_slarrb(in,d(ibegin),work(indlld + ibegin - 1),indeig, &
                             indeig,zero,two*eps,offset,work(wbegin),wgap(wbegin),werr(wbegin), &
                                       work(indwrk),iwork(iindwk),pivmin,spdiam,itmp1,iinfo)
                             if (iinfo /= 0) then
                                info = -3
                                return
                             end if
                             lambda = work(windex)
                             ! reset twist index from inaccurate lambda to
                             ! force computation of true mingma
                             iwork(iindr + windex) = 0
                          end if
                          ! given lambda, compute the eigenvector.
                          call la_clar1v(in,1,in,lambda,d(ibegin),l(ibegin),work( &
                          indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z(ibegin,windex &
                          ),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + windex),isuppz( &
                                    2*windex - 1),nrminv,resid,rqcorr,work(indwrk))
                          if (iter == 0) then
                             bstres = resid
                             bstw = lambda
                          elseif (resid < bstres) then
                             bstres = resid
                             bstw = lambda
                          end if
                          isupmn = min(isupmn,isuppz(2*windex - 1))
                          isupmx = max(isupmx,isuppz(2*windex))
                          iter = iter + 1
                          ! sin alpha <= |resid|/gap
                          ! note that both the residual and the gap are
                          ! proportional to the matrix, so ||t|| doesn't play
                          ! a role in the quotient
                          ! convergence test for rayleigh-quotient iteration
                          ! (omitted when bisection has been used)
                          if (resid > tol*gap .and. abs(rqcorr) > rqtol*abs(lambda) .and. .not. &
                                    usedbs) then
                             ! we need to check that the rqcorr update doesn't
                             ! move the eigenvalue away from the desired one and
                             ! towards a neighbor. -> protection with bisection
                             if (indeig <= negcnt) then
                                ! the wanted eigenvalue lies to the left
                                sgndef = -one
                             else
                                ! the wanted eigenvalue lies to the right
                                sgndef = one
                             end if
                             ! we only use the rqcorr if it improves the
                             ! the iterate reasonably.
                             if ((rqcorr*sgndef >= zero) .and. (lambda + rqcorr <= right) .and. ( &
                                       lambda + rqcorr >= left)) then
                                usedrq = .true.
                                ! store new midpoint of bisection interval in work
                                if (sgndef == one) then
                                   ! the current lambda is on the left of the true
                                   ! eigenvalue
                                   left = lambda
                                   ! we prefer to assume that the error estimate
                                   ! is correct. we could make the interval not
                                   ! as a bracket but to be modified if the rqcorr
                                   ! chooses to. in this case, the right side should
                                   ! be modified as follows:
                                    ! right = max(right, lambda + rqcorr)
                                else
                                   ! the current lambda is on the right of the true
                                   ! eigenvalue
                                   right = lambda
                                   ! see comment about assuming the error estimate is
                                   ! correct above.
                                    ! left = min(left, lambda + rqcorr)
                                end if
                                work(windex) = half*(right + left)
                                ! take rqcorr since it has the correct sign and
                                ! improves the iterate reasonably
                                lambda = lambda + rqcorr
                                ! update width of error interval
                                werr(windex) = half*(right - left)
                             else
                                needbs = .true.
                             end if
                             if (right - left < rqtol*abs(lambda)) then
                                   ! the eigenvalue is computed to bisection accuracy
                                   ! compute eigenvector and stop
                                usedbs = .true.
                                goto 120
                             elseif (iter < maxitr) then
                                goto 120
                             elseif (iter == maxitr) then
                                needbs = .true.
                                goto 120
                             else
                                info = 5
                                return
                             end if
                          else
                             stp2ii = .false.
             if (usedrq .and. usedbs .and. bstres <= resid) then
                                lambda = bstw
                                stp2ii = .true.
                             end if
                             if (stp2ii) then
                                ! improve error angle by second step
                                call la_clar1v(in,1,in,lambda,d(ibegin),l(ibegin), &
                                work(indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z( &
                                ibegin,windex),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + &
                                windex),isuppz(2*windex - 1),nrminv,resid,rqcorr,work(indwrk &
                                          ))
                             end if
                             work(windex) = lambda
                          end if
                          ! compute fp-vector support w.r.t. whole matrix
                          isuppz(2*windex - 1) = isuppz(2*windex - 1) + oldien
                          isuppz(2*windex) = isuppz(2*windex) + oldien
                          zfrom = isuppz(2*windex - 1)
                          zto = isuppz(2*windex)
                          isupmn = isupmn + oldien
                          isupmx = isupmx + oldien
                          ! ensure vector is ok if support in the rqi has changed
                          if (isupmn < zfrom) then
                             do ii = isupmn,zfrom - 1
                                z(ii,windex) = zero
                             end do
                          end if
                          if (isupmx > zto) then
                             do ii = zto + 1,isupmx
                                z(ii,windex) = zero
                             end do
                          end if
                          call la_csscal(zto - zfrom + 1,nrminv,z(zfrom,windex),1)
                          125 continue
                          ! update w
                          w(windex) = lambda + sigma
                          ! recompute the gaps on the left and right
                          ! but only allow them to become larger and not
                          ! smaller (which can only happen through "bad"
                          ! cancellation and doesn't reflect the theory
                          ! where the initial gaps are underestimated due
                          ! to werr being too crude.)
                          if (.not. eskip) then
                             if (k > 1) then
                                wgap(windmn) = max(wgap(windmn),w(windex) - werr(windex) - w( &
                                          windmn) - werr(windmn))
                             end if
                             if (windex < wend) then
                                wgap(windex) = max(savgap,w(windpl) - werr(windpl) - w( &
                                          windex) - werr(windex))
                             end if
                          end if
                          idone = idone + 1
                       end if
                       ! here ends the code for the current child
                       139 continue
                       ! proceed to any remaining child nodes
                       newfst = j + 1
                    end do loop_140
                 end do loop_150
                 ndepth = ndepth + 1
                 go to 40
              end if
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_clarrv
     !> ZLARRV: computes the eigenvectors of the tridiagonal matrix
     !> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
     !> The input eigenvalues should have been computed by DLARRE.

     pure subroutine la_zlarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1, &
               rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz,isuppz,work,iwork,info)
        use la_constants_dp,only:zero,half,one,two,three,four,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dol,dou,ldz,m,n
           integer(ilp),intent(out) :: info
           real(dp),intent(in) :: minrgp,pivmin,vl,vu
           real(dp),intent(inout) :: rtol1,rtol2
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),indexw(*),isplit(*)
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(dp),intent(inout) :: d(*),l(*),w(*),werr(*),wgap(*)
           real(dp),intent(in) :: gers(*)
           real(dp),intent(out) :: work(*)
           complex(dp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 10

           ! Local Scalars
           logical(lk) :: eskip,needbs,stp2ii,tryrqc,usedbs,usedrq
           integer(ilp) :: done,i,ibegin,idone,iend,ii,iindc1,iindc2,iindr,iindwk,iinfo, &
            im,in,indeig,indld,indlld,indwrk,isupmn,isupmx,iter,itmp1,j,jblk,k, &
            miniwsize,minwsize,nclus,ndepth,negcnt,newcls,newfst,newftt,newlst,newsiz, &
            offset,oldcls,oldfst,oldien,oldlst,oldncl,p,parity,q,wbegin,wend,windex, &
                      windmn,windpl,zfrom,zto,zusedl,zusedu,zusedw
           integer(ilp) :: indin1,indin2
           real(dp) :: bstres,bstw,eps,fudge,gap,gaptol,gl,gu,lambda,left,lgap,mingma, &
           nrminv,resid,rgap,right,rqcorr,rqtol,savgap,sgndef,sigma,spdiam,ssigma,tau, &
                     tmp,tol,ztz
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           intrinsic :: cmplx
           ! Executable Statements
           info = 0
           ! quick return if possible
           if ((n <= 0) .or. (m <= 0)) then
              return
           end if
           ! the first n entries of work are reserved for the eigenvalues
           indld = n + 1
           indlld = 2*n + 1
           indin1 = 3*n + 1
           indin2 = 4*n + 1
           indwrk = 5*n + 1
           minwsize = 12*n
           do i = 1,minwsize
              work(i) = zero
           end do
           ! iwork(iindr+1:iindr+n) hold the twist indices r for the
           ! factorization used to compute the fp vector
           iindr = 0
           ! iwork(iindc1+1:iinc2+n) are used to store the clusters of the current
           ! layer and the one above.
           iindc1 = n
           iindc2 = 2*n
           iindwk = 3*n + 1
           miniwsize = 7*n
           do i = 1,miniwsize
              iwork(i) = 0
           end do
           zusedl = 1
           if (dol > 1) then
              ! set lower bound for use of z
              zusedl = dol - 1
           end if
           zusedu = m
           if (dou < m) then
              ! set lower bound for use of z
              zusedu = dou + 1
           end if
           ! the width of the part of z that is used
           zusedw = zusedu - zusedl + 1
           call la_zlaset('FULL',n,zusedw,czero,czero,z(1,zusedl),ldz)
           eps = la_dlamch('PRECISION')
           rqtol = two*eps
           ! set expert flags for standard code.
           tryrqc = .true.
           if ((dol == 1) .and. (dou == m)) then
           else
              ! only selected eigenpairs are computed. since the other evalues
              ! are not refined by rq iteration, bisection has to compute to full
              ! accuracy.
              rtol1 = four*eps
              rtol2 = four*eps
           end if
           ! the entries wbegin:wend in w, werr, wgap correspond to the
           ! desired eigenvalues. the support of the nonzero eigenvector
           ! entries is contained in the interval ibegin:iend.
           ! remark that if k eigenpairs are desired, then the eigenvectors
           ! are stored in k contiguous columns of z.
           ! done is the number of eigenvectors already computed
           done = 0
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,iblock(m)
              iend = isplit(jblk)
              sigma = l(iend)
              ! find the eigenvectors of the submatrix indexed ibegin
              ! through iend.
              wend = wbegin - 1
              15 continue
              if (wend < m) then
                 if (iblock(wend + 1) == jblk) then
                    wend = wend + 1
                    go to 15
                 end if
              end if
              if (wend < wbegin) then
                 ibegin = iend + 1
                 cycle loop_170
              elseif ((wend < dol) .or. (wbegin > dou)) then
                 ibegin = iend + 1
                 wbegin = wend + 1
                 cycle loop_170
              end if
              ! find local spectral diameter of the block
              gl = gers(2*ibegin - 1)
              gu = gers(2*ibegin)
              do i = ibegin + 1,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              ! oldien is the last index of the previous block
              oldien = ibegin - 1
              ! calculate the size of the current block
              in = iend - ibegin + 1
              ! the number of eigenvalues in the current block
              im = wend - wbegin + 1
              ! this is for a 1x1 block
              if (ibegin == iend) then
                 done = done + 1
                 z(ibegin,wbegin) = cmplx(one,zero,KIND=dp)
                 isuppz(2*wbegin - 1) = ibegin
                 isuppz(2*wbegin) = ibegin
                 w(wbegin) = w(wbegin) + sigma
                 work(wbegin) = w(wbegin)
                 ibegin = iend + 1
                 wbegin = wbegin + 1
                 cycle loop_170
              end if
              ! the desired (shifted) eigenvalues are stored in w(wbegin:wend)
              ! note that these can be approximations, in this case, the corresp.
              ! entries of werr give the size of the uncertainty interval.
              ! the eigenvalue approximations will be refined when necessary as
              ! high relative accuracy is required for the computation of the
              ! corresponding eigenvectors.
              call la_dcopy(im,w(wbegin),1,work(wbegin),1)
              ! we store in w the eigenvalue approximations w.r.t. the original
              ! matrix t.
              do i = 1,im
                 w(wbegin + i - 1) = w(wbegin + i - 1) + sigma
              end do
              ! ndepth is the current depth of the representation tree
              ndepth = 0
              ! parity is either 1 or 0
              parity = 1
              ! nclus is the number of clusters for the next level of the
              ! representation tree, we start with nclus = 1 for the root
              nclus = 1
              iwork(iindc1 + 1) = 1
              iwork(iindc1 + 2) = im
              ! idone is the number of eigenvectors already computed in the current
              ! block
              idone = 0
              ! loop while( idone<im )
              ! generate the representation tree for the current block and
              ! compute the eigenvectors
              40 continue
              if (idone < im) then
                 ! this is a crude protection against infinitely deep trees
                 if (ndepth > m) then
                    info = -2
                    return
                 end if
                 ! breadth first processing of the current level of the representation
                 ! tree: oldncl = number of clusters on current level
                 oldncl = nclus
                 ! reset nclus to count the number of child clusters
                 nclus = 0
                 parity = 1 - parity
                 if (parity == 0) then
                    oldcls = iindc1
                    newcls = iindc2
                 else
                    oldcls = iindc2
                    newcls = iindc1
                 end if
                 ! process the clusters on the current level
                 loop_150: do i = 1,oldncl
                    j = oldcls + 2*i
                    ! oldfst, oldlst = first, last index of current cluster.
                                     ! cluster indices start with 1 and are relative
                                     ! to wbegin when accessing w, wgap, werr, z
                    oldfst = iwork(j - 1)
                    oldlst = iwork(j)
                    if (ndepth > 0) then
                       ! retrieve relatively robust representation (rrr) of cluster
                       ! that has been computed at the previous level
                       ! the rrr is stored in z and overwritten once the eigenvectors
                       ! have been computed or when the cluster is refined
                       if ((dol == 1) .and. (dou == m)) then
                          ! get representation from location of the leftmost evalue
                          ! of the cluster
                          j = wbegin + oldfst - 1
                       else
                          if (wbegin + oldfst - 1 < dol) then
                             ! get representation from the left end of z array
                             j = dol - 1
                          elseif (wbegin + oldfst - 1 > dou) then
                             ! get representation from the right end of z array
                             j = dou
                          else
                             j = wbegin + oldfst - 1
                          end if
                       end if
                       do k = 1,in - 1
                          d(ibegin + k - 1) = real(z(ibegin + k - 1,j),KIND=dp)
                          l(ibegin + k - 1) = real(z(ibegin + k - 1,j + 1),KIND=dp)
                       end do
                       d(iend) = real(z(iend,j),KIND=dp)
                       sigma = real(z(iend,j + 1),KIND=dp)
                       ! set the corresponding entries in z to zero
                       call la_zlaset('FULL',in,2,czero,czero,z(ibegin,j),ldz)

                    end if
                    ! compute dl and dll of current rrr
                    do j = ibegin,iend - 1
                       tmp = d(j)*l(j)
                       work(indld - 1 + j) = tmp
                       work(indlld - 1 + j) = tmp*l(j)
                    end do
                    if (ndepth > 0) then
                       ! p and q are index of the first and last eigenvalue to compute
                       ! within the current block
                       p = indexw(wbegin - 1 + oldfst)
                       q = indexw(wbegin - 1 + oldlst)
                       ! offset for the arrays work, wgap and werr, i.e., the p-offset
                       ! through the q-offset elements of these arrays are to be used.
                        ! offset = p-oldfst
                       offset = indexw(wbegin) - 1
                       ! perform limited bisection (if necessary) to get approximate
                       ! eigenvalues to the precision needed.
                       call la_dlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,q,rtol1, &
                       rtol2,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk),iwork( &
                                  iindwk),pivmin,spdiam,in,iinfo)
                       if (iinfo /= 0) then
                          info = -1
                          return
                       end if
                       ! we also recompute the extremal gaps. w holds all eigenvalues
                       ! of the unshifted matrix and must be used for computation
                       ! of wgap, the entries of work might stem from rrrs with
                       ! different shifts. the gaps from wbegin-1+oldfst to
                       ! wbegin-1+oldlst are correctly computed in la_dlarrb.
                       ! however, we only allow the gaps to become greater since
                       ! this is what should happen when we decrease werr
                       if (oldfst > 1) then
                          wgap(wbegin + oldfst - 2) = max(wgap(wbegin + oldfst - 2),w(wbegin + oldfst - 1) - &
                          werr(wbegin + oldfst - 1) - w(wbegin + oldfst - 2) - werr(wbegin + oldfst - 2))

                       end if
                       if (wbegin + oldlst - 1 < wend) then
                          wgap(wbegin + oldlst - 1) = max(wgap(wbegin + oldlst - 1),w(wbegin + oldlst) - &
                                    werr(wbegin + oldlst) - w(wbegin + oldlst - 1) - werr(wbegin + oldlst - 1))
                       end if
                       ! each time the eigenvalues in work get refined, we store
                       ! the newly found approximation with all shifts applied in w
                       do j = oldfst,oldlst
                          w(wbegin + j - 1) = work(wbegin + j - 1) + sigma
                       end do
                    end if
                    ! process the current node.
                    newfst = oldfst
                    loop_140: do j = oldfst,oldlst
                       if (j == oldlst) then
                          ! we are at the right end of the cluster, this is also the
                          ! boundary of the child cluster
                          newlst = j
                       else if (wgap(wbegin + j - 1) >= minrgp*abs(work(wbegin + j - 1))) &
                                 then
                          ! the right relative gap is big enough, the child cluster
                          ! (newfst,..,newlst) is well separated from the following
                          newlst = j
                        else
                          ! inside a child cluster, the relative gap is not
                          ! big enough.
                          cycle loop_140
                       end if
                       ! compute size of child cluster found
                       newsiz = newlst - newfst + 1
                       ! newftt is the place in z where the new rrr or the computed
                       ! eigenvector is to be stored
                       if ((dol == 1) .and. (dou == m)) then
                          ! store representation at location of the leftmost evalue
                          ! of the cluster
                          newftt = wbegin + newfst - 1
                       else
                          if (wbegin + newfst - 1 < dol) then
                             ! store representation at the left end of z array
                             newftt = dol - 1
                          elseif (wbegin + newfst - 1 > dou) then
                             ! store representation at the right end of z array
                             newftt = dou
                          else
                             newftt = wbegin + newfst - 1
                          end if
                       end if
                       if (newsiz > 1) then
                          ! current child is not a singleton but a cluster.
                          ! compute and store new representation of child.
                          ! compute left and right cluster gap.
                          ! lgap and rgap are not computed from work because
                          ! the eigenvalue approximations may stem from rrrs
                          ! different shifts. however, w hold all eigenvalues
                          ! of the unshifted matrix. still, the entries in wgap
                          ! have to be computed from work since the entries
                          ! in w might be of the same order so that gaps are not
                          ! exhibited correctly for very close eigenvalues.
                          if (newfst == 1) then
                             lgap = max(zero,w(wbegin) - werr(wbegin) - vl)
                         else
                             lgap = wgap(wbegin + newfst - 2)
                          end if
                          rgap = wgap(wbegin + newlst - 1)
                          ! compute left- and rightmost eigenvalue of child
                          ! to high precision in order to shift as close
                          ! as possible and obtain as large relative gaps
                          ! as possible
                          do k = 1,2
                             if (k == 1) then
                                p = indexw(wbegin - 1 + newfst)
                             else
                                p = indexw(wbegin - 1 + newlst)
                             end if
                             offset = indexw(wbegin) - 1
                             call la_dlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,p,rqtol, &
                             rqtol,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk), &
                                       iwork(iindwk),pivmin,spdiam,in,iinfo)
                          end do
                          if ((wbegin + newlst - 1 < dol) .or. (wbegin + newfst - 1 > dou)) then
                             ! if the cluster contains no desired eigenvalues
                             ! skip the computation of that branch of the rep. tree
                             ! we could skip before the refinement of the extremal
                             ! eigenvalues of the child, but then the representation
                             ! tree could be different from the one when nothing is
                             ! skipped. for this reason we skip at this place.
                             idone = idone + newlst - newfst + 1
                             goto 139
                          end if
                          ! compute rrr of child cluster.
                          ! note that the new rrr is stored in z
                          ! la_dlarrf needs lwork = 2*n
                          call la_dlarrf(in,d(ibegin),l(ibegin),work(indld + ibegin - 1), &
                          newfst,newlst,work(wbegin),wgap(wbegin),werr(wbegin),spdiam,lgap, &
                          rgap,pivmin,tau,work(indin1),work(indin2),work(indwrk),iinfo)

                          ! in the complex case, la_dlarrf cannot write
                          ! the new rrr directly into z and needs an intermediate
                          ! workspace
                          do k = 1,in - 1
                             z(ibegin + k - 1,newftt) = cmplx(work(indin1 + k - 1),zero,KIND=dp)

                             z(ibegin + k - 1,newftt + 1) = cmplx(work(indin2 + k - 1),zero,KIND=dp)

                          end do
                          z(iend,newftt) = cmplx(work(indin1 + in - 1),zero,KIND=dp)
                          if (iinfo == 0) then
                             ! a new rrr for the cluster was found by la_dlarrf
                             ! update shift and store it
                             ssigma = sigma + tau
                             z(iend,newftt + 1) = cmplx(ssigma,zero,KIND=dp)
                             ! work() are the midpoints and werr() the semi-width
                             ! note that the entries in w are unchanged.
                             do k = newfst,newlst
                                fudge = three*eps*abs(work(wbegin + k - 1))
                                work(wbegin + k - 1) = work(wbegin + k - 1) - tau
                                fudge = fudge + four*eps*abs(work(wbegin + k - 1))
                                ! fudge errors
                                werr(wbegin + k - 1) = werr(wbegin + k - 1) + fudge
                                ! gaps are not fudged. provided that werr is small
                                ! when eigenvalues are close, a zero gap indicates
                                ! that a new representation is needed for resolving
                                ! the cluster. a fudge could lead to a wrong decision
                                ! of judging eigenvalues 'separated' which in
                                ! reality are not. this could have a negative impact
                                ! on the orthogonality of the computed eigenvectors.
                             end do
                             nclus = nclus + 1
                             k = newcls + 2*nclus
                             iwork(k - 1) = newfst
                             iwork(k) = newlst
                          else
                             info = -2
                             return
                          end if
                       else
                          ! compute eigenvector of singleton
                          iter = 0
                          tol = four*log(real(in,KIND=dp))*eps
                          k = newfst
                          windex = wbegin + k - 1
                          windmn = max(windex - 1,1)
                          windpl = min(windex + 1,m)
                          lambda = work(windex)
                          done = done + 1
                          ! check if eigenvector computation is to be skipped
                          if ((windex < dol) .or. (windex > dou)) then
                             eskip = .true.
                             goto 125
                          else
                             eskip = .false.
                          end if
                          left = work(windex) - werr(windex)
                          right = work(windex) + werr(windex)
                          indeig = indexw(windex)
                          ! note that since we compute the eigenpairs for a child,
                          ! all eigenvalue approximations are w.r.t the same shift.
                          ! in this case, the entries in work should be used for
                          ! computing the gaps since they exhibit even very small
                          ! differences in the eigenvalues, as opposed to the
                          ! entries in w which might "look" the same.
                          if (k == 1) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vl, the formula
                             ! lgap = max( zero, (sigma - vl) + lambda )
                             ! can lead to an overestimation of the left gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small left gap.
                             lgap = eps*max(abs(left),abs(right))
                          else
                             lgap = wgap(windmn)
                          end if
                          if (k == im) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vu, the formula
                             ! can lead to an overestimation of the right gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small right gap.
                             rgap = eps*max(abs(left),abs(right))
                          else
                             rgap = wgap(windex)
                          end if
                          gap = min(lgap,rgap)
                          if ((k == 1) .or. (k == im)) then
                             ! the eigenvector support can become wrong
                             ! because significant entries could be cut off due to a
                             ! large gaptol parameter in lar1v. prevent this.
                             gaptol = zero
                          else
                             gaptol = gap*eps
                          end if
                          isupmn = in
                          isupmx = 1
                          ! update wgap so that it holds the minimum gap
                          ! to the left or the right. this is crucial in the
                          ! case where bisection is used to ensure that the
                          ! eigenvalue is refined up to the required precision.
                          ! the correct value is restored afterwards.
                          savgap = wgap(windex)
                          wgap(windex) = gap
                          ! we want to use the rayleigh quotient correction
                          ! as often as possible since it converges quadratically
                          ! when we are close enough to the desired eigenvalue.
                          ! however, the rayleigh quotient can have the wrong sign
                          ! and lead us away from the desired eigenvalue. in this
                          ! case, the best we can do is to use bisection.
                          usedbs = .false.
                          usedrq = .false.
                          ! bisection is initially turned off unless it is forced
                          needbs = .not. tryrqc
                          120 continue
                          ! check if bisection should be used to refine eigenvalue
                          if (needbs) then
                             ! take the bisection as new iterate
                             usedbs = .true.
                             itmp1 = iwork(iindr + windex)
                             offset = indexw(wbegin) - 1
                             call la_dlarrb(in,d(ibegin),work(indlld + ibegin - 1),indeig, &
                             indeig,zero,two*eps,offset,work(wbegin),wgap(wbegin),werr(wbegin), &
                                       work(indwrk),iwork(iindwk),pivmin,spdiam,itmp1,iinfo)
                             if (iinfo /= 0) then
                                info = -3
                                return
                             end if
                             lambda = work(windex)
                             ! reset twist index from inaccurate lambda to
                             ! force computation of true mingma
                             iwork(iindr + windex) = 0
                          end if
                          ! given lambda, compute the eigenvector.
                          call la_zlar1v(in,1,in,lambda,d(ibegin),l(ibegin),work( &
                          indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z(ibegin,windex &
                          ),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + windex),isuppz( &
                                    2*windex - 1),nrminv,resid,rqcorr,work(indwrk))
                          if (iter == 0) then
                             bstres = resid
                             bstw = lambda
                          elseif (resid < bstres) then
                             bstres = resid
                             bstw = lambda
                          end if
                          isupmn = min(isupmn,isuppz(2*windex - 1))
                          isupmx = max(isupmx,isuppz(2*windex))
                          iter = iter + 1
                          ! sin alpha <= |resid|/gap
                          ! note that both the residual and the gap are
                          ! proportional to the matrix, so ||t|| doesn't play
                          ! a role in the quotient
                          ! convergence test for rayleigh-quotient iteration
                          ! (omitted when bisection has been used)
                          if (resid > tol*gap .and. abs(rqcorr) > rqtol*abs(lambda) .and. .not. &
                                    usedbs) then
                             ! we need to check that the rqcorr update doesn't
                             ! move the eigenvalue away from the desired one and
                             ! towards a neighbor. -> protection with bisection
                             if (indeig <= negcnt) then
                                ! the wanted eigenvalue lies to the left
                                sgndef = -one
                             else
                                ! the wanted eigenvalue lies to the right
                                sgndef = one
                             end if
                             ! we only use the rqcorr if it improves the
                             ! the iterate reasonably.
                             if ((rqcorr*sgndef >= zero) .and. (lambda + rqcorr <= right) .and. ( &
                                       lambda + rqcorr >= left)) then
                                usedrq = .true.
                                ! store new midpoint of bisection interval in work
                                if (sgndef == one) then
                                   ! the current lambda is on the left of the true
                                   ! eigenvalue
                                   left = lambda
                                   ! we prefer to assume that the error estimate
                                   ! is correct. we could make the interval not
                                   ! as a bracket but to be modified if the rqcorr
                                   ! chooses to. in this case, the right side should
                                   ! be modified as follows:
                                    ! right = max(right, lambda + rqcorr)
                                else
                                   ! the current lambda is on the right of the true
                                   ! eigenvalue
                                   right = lambda
                                   ! see comment about assuming the error estimate is
                                   ! correct above.
                                    ! left = min(left, lambda + rqcorr)
                                end if
                                work(windex) = half*(right + left)
                                ! take rqcorr since it has the correct sign and
                                ! improves the iterate reasonably
                                lambda = lambda + rqcorr
                                ! update width of error interval
                                werr(windex) = half*(right - left)
                             else
                                needbs = .true.
                             end if
                             if (right - left < rqtol*abs(lambda)) then
                                   ! the eigenvalue is computed to bisection accuracy
                                   ! compute eigenvector and stop
                                usedbs = .true.
                                goto 120
                             elseif (iter < maxitr) then
                                goto 120
                             elseif (iter == maxitr) then
                                needbs = .true.
                                goto 120
                             else
                                info = 5
                                return
                             end if
                          else
                             stp2ii = .false.
             if (usedrq .and. usedbs .and. bstres <= resid) then
                                lambda = bstw
                                stp2ii = .true.
                             end if
                             if (stp2ii) then
                                ! improve error angle by second step
                                call la_zlar1v(in,1,in,lambda,d(ibegin),l(ibegin), &
                                work(indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z( &
                                ibegin,windex),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + &
                                windex),isuppz(2*windex - 1),nrminv,resid,rqcorr,work(indwrk &
                                          ))
                             end if
                             work(windex) = lambda
                          end if
                          ! compute fp-vector support w.r.t. whole matrix
                          isuppz(2*windex - 1) = isuppz(2*windex - 1) + oldien
                          isuppz(2*windex) = isuppz(2*windex) + oldien
                          zfrom = isuppz(2*windex - 1)
                          zto = isuppz(2*windex)
                          isupmn = isupmn + oldien
                          isupmx = isupmx + oldien
                          ! ensure vector is ok if support in the rqi has changed
                          if (isupmn < zfrom) then
                             do ii = isupmn,zfrom - 1
                                z(ii,windex) = zero
                             end do
                          end if
                          if (isupmx > zto) then
                             do ii = zto + 1,isupmx
                                z(ii,windex) = zero
                             end do
                          end if
                          call la_zdscal(zto - zfrom + 1,nrminv,z(zfrom,windex),1)
                          125 continue
                          ! update w
                          w(windex) = lambda + sigma
                          ! recompute the gaps on the left and right
                          ! but only allow them to become larger and not
                          ! smaller (which can only happen through "bad"
                          ! cancellation and doesn't reflect the theory
                          ! where the initial gaps are underestimated due
                          ! to werr being too crude.)
                          if (.not. eskip) then
                             if (k > 1) then
                                wgap(windmn) = max(wgap(windmn),w(windex) - werr(windex) - w( &
                                          windmn) - werr(windmn))
                             end if
                             if (windex < wend) then
                                wgap(windex) = max(savgap,w(windpl) - werr(windpl) - w( &
                                          windex) - werr(windex))
                             end if
                          end if
                          idone = idone + 1
                       end if
                       ! here ends the code for the current child
                       139 continue
                       ! proceed to any remaining child nodes
                       newfst = j + 1
                    end do loop_140
                 end do loop_150
                 ndepth = ndepth + 1
                 go to 40
              end if
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_zlarrv
     !> WLARRV: computes the eigenvectors of the tridiagonal matrix
     !> T = L D L**T given L, D and APPROXIMATIONS to the eigenvalues of L D L**T.
     !> The input eigenvalues should have been computed by QLARRE.

     pure subroutine la_wlarrv(n,vl,vu,d,l,pivmin,isplit,m,dol,dou,minrgp,rtol1, &
               rtol2,w,werr,wgap,iblock,indexw,gers,z,ldz,isuppz,work,iwork,info)
        use la_constants_qp,only:zero,half,one,two,three,four,czero
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: dol,dou,ldz,m,n
           integer(ilp),intent(out) :: info
           real(qp),intent(in) :: minrgp,pivmin,vl,vu
           real(qp),intent(inout) :: rtol1,rtol2
           ! Array Arguments
           integer(ilp),intent(in) :: iblock(*),indexw(*),isplit(*)
           integer(ilp),intent(out) :: isuppz(*),iwork(*)
           real(qp),intent(inout) :: d(*),l(*),w(*),werr(*),wgap(*)
           real(qp),intent(in) :: gers(*)
           real(qp),intent(out) :: work(*)
           complex(qp),intent(out) :: z(ldz,*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: maxitr = 10

           ! Local Scalars
           logical(lk) :: eskip,needbs,stp2ii,tryrqc,usedbs,usedrq
           integer(ilp) :: done,i,ibegin,idone,iend,ii,iindc1,iindc2,iindr,iindwk,iinfo, &
            im,in,indeig,indld,indlld,indwrk,isupmn,isupmx,iter,itmp1,j,jblk,k, &
            miniwsize,minwsize,nclus,ndepth,negcnt,newcls,newfst,newftt,newlst,newsiz, &
            offset,oldcls,oldfst,oldien,oldlst,oldncl,p,parity,q,wbegin,wend,windex, &
                      windmn,windpl,zfrom,zto,zusedl,zusedu,zusedw
           integer(ilp) :: indin1,indin2
           real(qp) :: bstres,bstw,eps,fudge,gap,gaptol,gl,gu,lambda,left,lgap,mingma, &
           nrminv,resid,rgap,right,rqcorr,rqtol,savgap,sgndef,sigma,spdiam,ssigma,tau, &
                     tmp,tol,ztz
           ! Intrinsic Functions
           intrinsic :: abs,real,max,min
           intrinsic :: cmplx
           ! Executable Statements
           info = 0
           ! quick return if possible
           if ((n <= 0) .or. (m <= 0)) then
              return
           end if
           ! the first n entries of work are reserved for the eigenvalues
           indld = n + 1
           indlld = 2*n + 1
           indin1 = 3*n + 1
           indin2 = 4*n + 1
           indwrk = 5*n + 1
           minwsize = 12*n
           do i = 1,minwsize
              work(i) = zero
           end do
           ! iwork(iindr+1:iindr+n) hold the twist indices r for the
           ! factorization used to compute the fp vector
           iindr = 0
           ! iwork(iindc1+1:iinc2+n) are used to store the clusters of the current
           ! layer and the one above.
           iindc1 = n
           iindc2 = 2*n
           iindwk = 3*n + 1
           miniwsize = 7*n
           do i = 1,miniwsize
              iwork(i) = 0
           end do
           zusedl = 1
           if (dol > 1) then
              ! set lower bound for use of z
              zusedl = dol - 1
           end if
           zusedu = m
           if (dou < m) then
              ! set lower bound for use of z
              zusedu = dou + 1
           end if
           ! the width of the part of z that is used
           zusedw = zusedu - zusedl + 1
           call la_wlaset('FULL',n,zusedw,czero,czero,z(1,zusedl),ldz)
           eps = la_qlamch('PRECISION')
           rqtol = two*eps
           ! set expert flags for standard code.
           tryrqc = .true.
           if ((dol == 1) .and. (dou == m)) then
           else
              ! only selected eigenpairs are computed. since the other evalues
              ! are not refined by rq iteration, bisection has to compute to full
              ! accuracy.
              rtol1 = four*eps
              rtol2 = four*eps
           end if
           ! the entries wbegin:wend in w, werr, wgap correspond to the
           ! desired eigenvalues. the support of the nonzero eigenvector
           ! entries is contained in the interval ibegin:iend.
           ! remark that if k eigenpairs are desired, then the eigenvectors
           ! are stored in k contiguous columns of z.
           ! done is the number of eigenvectors already computed
           done = 0
           ibegin = 1
           wbegin = 1
           loop_170: do jblk = 1,iblock(m)
              iend = isplit(jblk)
              sigma = l(iend)
              ! find the eigenvectors of the submatrix indexed ibegin
              ! through iend.
              wend = wbegin - 1
              15 continue
              if (wend < m) then
                 if (iblock(wend + 1) == jblk) then
                    wend = wend + 1
                    go to 15
                 end if
              end if
              if (wend < wbegin) then
                 ibegin = iend + 1
                 cycle loop_170
              elseif ((wend < dol) .or. (wbegin > dou)) then
                 ibegin = iend + 1
                 wbegin = wend + 1
                 cycle loop_170
              end if
              ! find local spectral diameter of the block
              gl = gers(2*ibegin - 1)
              gu = gers(2*ibegin)
              do i = ibegin + 1,iend
                 gl = min(gers(2*i - 1),gl)
                 gu = max(gers(2*i),gu)
              end do
              spdiam = gu - gl
              ! oldien is the last index of the previous block
              oldien = ibegin - 1
              ! calculate the size of the current block
              in = iend - ibegin + 1
              ! the number of eigenvalues in the current block
              im = wend - wbegin + 1
              ! this is for a 1x1 block
              if (ibegin == iend) then
                 done = done + 1
                 z(ibegin,wbegin) = cmplx(one,zero,KIND=qp)
                 isuppz(2*wbegin - 1) = ibegin
                 isuppz(2*wbegin) = ibegin
                 w(wbegin) = w(wbegin) + sigma
                 work(wbegin) = w(wbegin)
                 ibegin = iend + 1
                 wbegin = wbegin + 1
                 cycle loop_170
              end if
              ! the desired (shifted) eigenvalues are stored in w(wbegin:wend)
              ! note that these can be approximations, in this case, the corresp.
              ! entries of werr give the size of the uncertainty interval.
              ! the eigenvalue approximations will be refined when necessary as
              ! high relative accuracy is required for the computation of the
              ! corresponding eigenvectors.
              call la_qcopy(im,w(wbegin),1,work(wbegin),1)
              ! we store in w the eigenvalue approximations w.r.t. the original
              ! matrix t.
              do i = 1,im
                 w(wbegin + i - 1) = w(wbegin + i - 1) + sigma
              end do
              ! ndepth is the current depth of the representation tree
              ndepth = 0
              ! parity is either 1 or 0
              parity = 1
              ! nclus is the number of clusters for the next level of the
              ! representation tree, we start with nclus = 1 for the root
              nclus = 1
              iwork(iindc1 + 1) = 1
              iwork(iindc1 + 2) = im
              ! idone is the number of eigenvectors already computed in the current
              ! block
              idone = 0
              ! loop while( idone<im )
              ! generate the representation tree for the current block and
              ! compute the eigenvectors
              40 continue
              if (idone < im) then
                 ! this is a crude protection against infinitely deep trees
                 if (ndepth > m) then
                    info = -2
                    return
                 end if
                 ! breadth first processing of the current level of the representation
                 ! tree: oldncl = number of clusters on current level
                 oldncl = nclus
                 ! reset nclus to count the number of child clusters
                 nclus = 0
                 parity = 1 - parity
                 if (parity == 0) then
                    oldcls = iindc1
                    newcls = iindc2
                 else
                    oldcls = iindc2
                    newcls = iindc1
                 end if
                 ! process the clusters on the current level
                 loop_150: do i = 1,oldncl
                    j = oldcls + 2*i
                    ! oldfst, oldlst = first, last index of current cluster.
                                     ! cluster indices start with 1 and are relative
                                     ! to wbegin when accessing w, wgap, werr, z
                    oldfst = iwork(j - 1)
                    oldlst = iwork(j)
                    if (ndepth > 0) then
                       ! retrieve relatively robust representation (rrr) of cluster
                       ! that has been computed at the previous level
                       ! the rrr is stored in z and overwritten once the eigenvectors
                       ! have been computed or when the cluster is refined
                       if ((dol == 1) .and. (dou == m)) then
                          ! get representation from location of the leftmost evalue
                          ! of the cluster
                          j = wbegin + oldfst - 1
                       else
                          if (wbegin + oldfst - 1 < dol) then
                             ! get representation from the left end of z array
                             j = dol - 1
                          elseif (wbegin + oldfst - 1 > dou) then
                             ! get representation from the right end of z array
                             j = dou
                          else
                             j = wbegin + oldfst - 1
                          end if
                       end if
                       do k = 1,in - 1
                          d(ibegin + k - 1) = real(z(ibegin + k - 1,j),KIND=qp)
                          l(ibegin + k - 1) = real(z(ibegin + k - 1,j + 1),KIND=qp)
                       end do
                       d(iend) = real(z(iend,j),KIND=qp)
                       sigma = real(z(iend,j + 1),KIND=qp)
                       ! set the corresponding entries in z to zero
                       call la_wlaset('FULL',in,2,czero,czero,z(ibegin,j),ldz)

                    end if
                    ! compute dl and dll of current rrr
                    do j = ibegin,iend - 1
                       tmp = d(j)*l(j)
                       work(indld - 1 + j) = tmp
                       work(indlld - 1 + j) = tmp*l(j)
                    end do
                    if (ndepth > 0) then
                       ! p and q are index of the first and last eigenvalue to compute
                       ! within the current block
                       p = indexw(wbegin - 1 + oldfst)
                       q = indexw(wbegin - 1 + oldlst)
                       ! offset for the arrays work, wgap and werr, i.e., the p-offset
                       ! through the q-offset elements of these arrays are to be used.
                        ! offset = p-oldfst
                       offset = indexw(wbegin) - 1
                       ! perform limited bisection (if necessary) to get approximate
                       ! eigenvalues to the precision needed.
                       call la_qlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,q,rtol1, &
                       rtol2,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk),iwork( &
                                  iindwk),pivmin,spdiam,in,iinfo)
                       if (iinfo /= 0) then
                          info = -1
                          return
                       end if
                       ! we also recompute the extremal gaps. w holds all eigenvalues
                       ! of the unshifted matrix and must be used for computation
                       ! of wgap, the entries of work might stem from rrrs with
                       ! different shifts. the gaps from wbegin-1+oldfst to
                       ! wbegin-1+oldlst are correctly computed in la_qlarrb.
                       ! however, we only allow the gaps to become greater since
                       ! this is what should happen when we decrease werr
                       if (oldfst > 1) then
                          wgap(wbegin + oldfst - 2) = max(wgap(wbegin + oldfst - 2),w(wbegin + oldfst - 1) - &
                          werr(wbegin + oldfst - 1) - w(wbegin + oldfst - 2) - werr(wbegin + oldfst - 2))

                       end if
                       if (wbegin + oldlst - 1 < wend) then
                          wgap(wbegin + oldlst - 1) = max(wgap(wbegin + oldlst - 1),w(wbegin + oldlst) - &
                                    werr(wbegin + oldlst) - w(wbegin + oldlst - 1) - werr(wbegin + oldlst - 1))
                       end if
                       ! each time the eigenvalues in work get refined, we store
                       ! the newly found approximation with all shifts applied in w
                       do j = oldfst,oldlst
                          w(wbegin + j - 1) = work(wbegin + j - 1) + sigma
                       end do
                    end if
                    ! process the current node.
                    newfst = oldfst
                    loop_140: do j = oldfst,oldlst
                       if (j == oldlst) then
                          ! we are at the right end of the cluster, this is also the
                          ! boundary of the child cluster
                          newlst = j
                       else if (wgap(wbegin + j - 1) >= minrgp*abs(work(wbegin + j - 1))) &
                                 then
                          ! the right relative gap is big enough, the child cluster
                          ! (newfst,..,newlst) is well separated from the following
                          newlst = j
                        else
                          ! inside a child cluster, the relative gap is not
                          ! big enough.
                          cycle loop_140
                       end if
                       ! compute size of child cluster found
                       newsiz = newlst - newfst + 1
                       ! newftt is the place in z where the new rrr or the computed
                       ! eigenvector is to be stored
                       if ((dol == 1) .and. (dou == m)) then
                          ! store representation at location of the leftmost evalue
                          ! of the cluster
                          newftt = wbegin + newfst - 1
                       else
                          if (wbegin + newfst - 1 < dol) then
                             ! store representation at the left end of z array
                             newftt = dol - 1
                          elseif (wbegin + newfst - 1 > dou) then
                             ! store representation at the right end of z array
                             newftt = dou
                          else
                             newftt = wbegin + newfst - 1
                          end if
                       end if
                       if (newsiz > 1) then
                          ! current child is not a singleton but a cluster.
                          ! compute and store new representation of child.
                          ! compute left and right cluster gap.
                          ! lgap and rgap are not computed from work because
                          ! the eigenvalue approximations may stem from rrrs
                          ! different shifts. however, w hold all eigenvalues
                          ! of the unshifted matrix. still, the entries in wgap
                          ! have to be computed from work since the entries
                          ! in w might be of the same order so that gaps are not
                          ! exhibited correctly for very close eigenvalues.
                          if (newfst == 1) then
                             lgap = max(zero,w(wbegin) - werr(wbegin) - vl)
                         else
                             lgap = wgap(wbegin + newfst - 2)
                          end if
                          rgap = wgap(wbegin + newlst - 1)
                          ! compute left- and rightmost eigenvalue of child
                          ! to high precision in order to shift as close
                          ! as possible and obtain as large relative gaps
                          ! as possible
                          do k = 1,2
                             if (k == 1) then
                                p = indexw(wbegin - 1 + newfst)
                             else
                                p = indexw(wbegin - 1 + newlst)
                             end if
                             offset = indexw(wbegin) - 1
                             call la_qlarrb(in,d(ibegin),work(indlld + ibegin - 1),p,p,rqtol, &
                             rqtol,offset,work(wbegin),wgap(wbegin),werr(wbegin),work(indwrk), &
                                       iwork(iindwk),pivmin,spdiam,in,iinfo)
                          end do
                          if ((wbegin + newlst - 1 < dol) .or. (wbegin + newfst - 1 > dou)) then
                             ! if the cluster contains no desired eigenvalues
                             ! skip the computation of that branch of the rep. tree
                             ! we could skip before the refinement of the extremal
                             ! eigenvalues of the child, but then the representation
                             ! tree could be different from the one when nothing is
                             ! skipped. for this reason we skip at this place.
                             idone = idone + newlst - newfst + 1
                             goto 139
                          end if
                          ! compute rrr of child cluster.
                          ! note that the new rrr is stored in z
                          ! la_qlarrf needs lwork = 2*n
                          call la_qlarrf(in,d(ibegin),l(ibegin),work(indld + ibegin - 1), &
                          newfst,newlst,work(wbegin),wgap(wbegin),werr(wbegin),spdiam,lgap, &
                          rgap,pivmin,tau,work(indin1),work(indin2),work(indwrk),iinfo)

                          ! in the complex case, la_qlarrf cannot write
                          ! the new rrr directly into z and needs an intermediate
                          ! workspace
                          do k = 1,in - 1
                             z(ibegin + k - 1,newftt) = cmplx(work(indin1 + k - 1),zero,KIND=qp)

                             z(ibegin + k - 1,newftt + 1) = cmplx(work(indin2 + k - 1),zero,KIND=qp)

                          end do
                          z(iend,newftt) = cmplx(work(indin1 + in - 1),zero,KIND=qp)
                          if (iinfo == 0) then
                             ! a new rrr for the cluster was found by la_qlarrf
                             ! update shift and store it
                             ssigma = sigma + tau
                             z(iend,newftt + 1) = cmplx(ssigma,zero,KIND=qp)
                             ! work() are the midpoints and werr() the semi-width
                             ! note that the entries in w are unchanged.
                             do k = newfst,newlst
                                fudge = three*eps*abs(work(wbegin + k - 1))
                                work(wbegin + k - 1) = work(wbegin + k - 1) - tau
                                fudge = fudge + four*eps*abs(work(wbegin + k - 1))
                                ! fudge errors
                                werr(wbegin + k - 1) = werr(wbegin + k - 1) + fudge
                                ! gaps are not fudged. provided that werr is small
                                ! when eigenvalues are close, a zero gap indicates
                                ! that a new representation is needed for resolving
                                ! the cluster. a fudge could lead to a wrong decision
                                ! of judging eigenvalues 'separated' which in
                                ! reality are not. this could have a negative impact
                                ! on the orthogonality of the computed eigenvectors.
                             end do
                             nclus = nclus + 1
                             k = newcls + 2*nclus
                             iwork(k - 1) = newfst
                             iwork(k) = newlst
                          else
                             info = -2
                             return
                          end if
                       else
                          ! compute eigenvector of singleton
                          iter = 0
                          tol = four*log(real(in,KIND=qp))*eps
                          k = newfst
                          windex = wbegin + k - 1
                          windmn = max(windex - 1,1)
                          windpl = min(windex + 1,m)
                          lambda = work(windex)
                          done = done + 1
                          ! check if eigenvector computation is to be skipped
                          if ((windex < dol) .or. (windex > dou)) then
                             eskip = .true.
                             goto 125
                          else
                             eskip = .false.
                          end if
                          left = work(windex) - werr(windex)
                          right = work(windex) + werr(windex)
                          indeig = indexw(windex)
                          ! note that since we compute the eigenpairs for a child,
                          ! all eigenvalue approximations are w.r.t the same shift.
                          ! in this case, the entries in work should be used for
                          ! computing the gaps since they exhibit even very small
                          ! differences in the eigenvalues, as opposed to the
                          ! entries in w which might "look" the same.
                          if (k == 1) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vl, the formula
                             ! lgap = max( zero, (sigma - vl) + lambda )
                             ! can lead to an overestimation of the left gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small left gap.
                             lgap = eps*max(abs(left),abs(right))
                          else
                             lgap = wgap(windmn)
                          end if
                          if (k == im) then
                             ! in the case range='i' and with not much initial
                             ! accuracy in lambda and vu, the formula
                             ! can lead to an overestimation of the right gap and
                             ! thus to inadequately early rqi 'convergence'.
                             ! prevent this by forcing a small right gap.
                             rgap = eps*max(abs(left),abs(right))
                          else
                             rgap = wgap(windex)
                          end if
                          gap = min(lgap,rgap)
                          if ((k == 1) .or. (k == im)) then
                             ! the eigenvector support can become wrong
                             ! because significant entries could be cut off due to a
                             ! large gaptol parameter in lar1v. prevent this.
                             gaptol = zero
                          else
                             gaptol = gap*eps
                          end if
                          isupmn = in
                          isupmx = 1
                          ! update wgap so that it holds the minimum gap
                          ! to the left or the right. this is crucial in the
                          ! case where bisection is used to ensure that the
                          ! eigenvalue is refined up to the required precision.
                          ! the correct value is restored afterwards.
                          savgap = wgap(windex)
                          wgap(windex) = gap
                          ! we want to use the rayleigh quotient correction
                          ! as often as possible since it converges quadratically
                          ! when we are close enough to the desired eigenvalue.
                          ! however, the rayleigh quotient can have the wrong sign
                          ! and lead us away from the desired eigenvalue. in this
                          ! case, the best we can do is to use bisection.
                          usedbs = .false.
                          usedrq = .false.
                          ! bisection is initially turned off unless it is forced
                          needbs = .not. tryrqc
                          120 continue
                          ! check if bisection should be used to refine eigenvalue
                          if (needbs) then
                             ! take the bisection as new iterate
                             usedbs = .true.
                             itmp1 = iwork(iindr + windex)
                             offset = indexw(wbegin) - 1
                             call la_qlarrb(in,d(ibegin),work(indlld + ibegin - 1),indeig, &
                             indeig,zero,two*eps,offset,work(wbegin),wgap(wbegin),werr(wbegin), &
                                       work(indwrk),iwork(iindwk),pivmin,spdiam,itmp1,iinfo)
                             if (iinfo /= 0) then
                                info = -3
                                return
                             end if
                             lambda = work(windex)
                             ! reset twist index from inaccurate lambda to
                             ! force computation of true mingma
                             iwork(iindr + windex) = 0
                          end if
                          ! given lambda, compute the eigenvector.
                          call la_wlar1v(in,1,in,lambda,d(ibegin),l(ibegin),work( &
                          indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z(ibegin,windex &
                          ),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + windex),isuppz( &
                                    2*windex - 1),nrminv,resid,rqcorr,work(indwrk))
                          if (iter == 0) then
                             bstres = resid
                             bstw = lambda
                          elseif (resid < bstres) then
                             bstres = resid
                             bstw = lambda
                          end if
                          isupmn = min(isupmn,isuppz(2*windex - 1))
                          isupmx = max(isupmx,isuppz(2*windex))
                          iter = iter + 1
                          ! sin alpha <= |resid|/gap
                          ! note that both the residual and the gap are
                          ! proportional to the matrix, so ||t|| doesn't play
                          ! a role in the quotient
                          ! convergence test for rayleigh-quotient iteration
                          ! (omitted when bisection has been used)
                          if (resid > tol*gap .and. abs(rqcorr) > rqtol*abs(lambda) .and. .not. &
                                    usedbs) then
                             ! we need to check that the rqcorr update doesn't
                             ! move the eigenvalue away from the desired one and
                             ! towards a neighbor. -> protection with bisection
                             if (indeig <= negcnt) then
                                ! the wanted eigenvalue lies to the left
                                sgndef = -one
                             else
                                ! the wanted eigenvalue lies to the right
                                sgndef = one
                             end if
                             ! we only use the rqcorr if it improves the
                             ! the iterate reasonably.
                             if ((rqcorr*sgndef >= zero) .and. (lambda + rqcorr <= right) .and. ( &
                                       lambda + rqcorr >= left)) then
                                usedrq = .true.
                                ! store new midpoint of bisection interval in work
                                if (sgndef == one) then
                                   ! the current lambda is on the left of the true
                                   ! eigenvalue
                                   left = lambda
                                   ! we prefer to assume that the error estimate
                                   ! is correct. we could make the interval not
                                   ! as a bracket but to be modified if the rqcorr
                                   ! chooses to. in this case, the right side should
                                   ! be modified as follows:
                                    ! right = max(right, lambda + rqcorr)
                                else
                                   ! the current lambda is on the right of the true
                                   ! eigenvalue
                                   right = lambda
                                   ! see comment about assuming the error estimate is
                                   ! correct above.
                                    ! left = min(left, lambda + rqcorr)
                                end if
                                work(windex) = half*(right + left)
                                ! take rqcorr since it has the correct sign and
                                ! improves the iterate reasonably
                                lambda = lambda + rqcorr
                                ! update width of error interval
                                werr(windex) = half*(right - left)
                             else
                                needbs = .true.
                             end if
                             if (right - left < rqtol*abs(lambda)) then
                                   ! the eigenvalue is computed to bisection accuracy
                                   ! compute eigenvector and stop
                                usedbs = .true.
                                goto 120
                             elseif (iter < maxitr) then
                                goto 120
                             elseif (iter == maxitr) then
                                needbs = .true.
                                goto 120
                             else
                                info = 5
                                return
                             end if
                          else
                             stp2ii = .false.
             if (usedrq .and. usedbs .and. bstres <= resid) then
                                lambda = bstw
                                stp2ii = .true.
                             end if
                             if (stp2ii) then
                                ! improve error angle by second step
                                call la_wlar1v(in,1,in,lambda,d(ibegin),l(ibegin), &
                                work(indld + ibegin - 1),work(indlld + ibegin - 1),pivmin,gaptol,z( &
                                ibegin,windex),.not. usedbs,negcnt,ztz,mingma,iwork(iindr + &
                                windex),isuppz(2*windex - 1),nrminv,resid,rqcorr,work(indwrk &
                                          ))
                             end if
                             work(windex) = lambda
                          end if
                          ! compute fp-vector support w.r.t. whole matrix
                          isuppz(2*windex - 1) = isuppz(2*windex - 1) + oldien
                          isuppz(2*windex) = isuppz(2*windex) + oldien
                          zfrom = isuppz(2*windex - 1)
                          zto = isuppz(2*windex)
                          isupmn = isupmn + oldien
                          isupmx = isupmx + oldien
                          ! ensure vector is ok if support in the rqi has changed
                          if (isupmn < zfrom) then
                             do ii = isupmn,zfrom - 1
                                z(ii,windex) = zero
                             end do
                          end if
                          if (isupmx > zto) then
                             do ii = zto + 1,isupmx
                                z(ii,windex) = zero
                             end do
                          end if
                          call la_wqscal(zto - zfrom + 1,nrminv,z(zfrom,windex),1)
                          125 continue
                          ! update w
                          w(windex) = lambda + sigma
                          ! recompute the gaps on the left and right
                          ! but only allow them to become larger and not
                          ! smaller (which can only happen through "bad"
                          ! cancellation and doesn't reflect the theory
                          ! where the initial gaps are underestimated due
                          ! to werr being too crude.)
                          if (.not. eskip) then
                             if (k > 1) then
                                wgap(windmn) = max(wgap(windmn),w(windex) - werr(windex) - w( &
                                          windmn) - werr(windmn))
                             end if
                             if (windex < wend) then
                                wgap(windex) = max(savgap,w(windpl) - werr(windpl) - w( &
                                          windex) - werr(windex))
                             end if
                          end if
                          idone = idone + 1
                       end if
                       ! here ends the code for the current child
                       139 continue
                       ! proceed to any remaining child nodes
                       newfst = j + 1
                    end do loop_140
                 end do loop_150
                 ndepth = ndepth + 1
                 go to 40
              end if
              ibegin = iend + 1
              wbegin = wend + 1
           end do loop_170
           return
     end subroutine la_wlarrv

end module la_lapack_eigv_tridiag2
