!> BLAS-like level 1: scaling, conjugation, sums of squares, sorting
module la_lapack_blas_like_l1
     use la_constants
     use la_blas_aux
     use la_blas_level1
     use la_lapack_auxiliary
     implicit none(type,external)
     private

     public :: sp,dp,qp,lk,ilp
     public :: la_slasrt
     public :: la_slassq
     public :: la_srscl
     public :: la_dlasrt
     public :: la_dlassq
     public :: la_drscl
#ifdef LA_WITH_XDP
     public :: la_xlasrt
     public :: la_xlassq
     public :: la_xrscl
#endif
#ifdef LA_WITH_QP
     public :: la_qlasrt
     public :: la_qlassq
     public :: la_qrscl
#endif
     public :: la_csrscl
     public :: la_clacgv
     public :: la_classq
     public :: la_zdrscl
     public :: la_zlacgv
     public :: la_zlassq
#ifdef LA_WITH_XDP
     public :: la_yxrscl
     public :: la_ylacgv
     public :: la_ylassq
#endif
#ifdef LA_WITH_QP
     public :: la_wqrscl
     public :: la_wlacgv
     public :: la_wlassq
#endif

     contains

     !> Sort the numbers in D in increasing order (if ID = 'I') or
     !> in decreasing order (if ID = 'D' ).
     !> Use Quick Sort, reverting to Insertion sort on arrays of
     !> size <= 20. Dimension of STACK limits N to about 2**32.

     pure subroutine la_slasrt(id,n,d,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: id
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(sp),intent(inout) :: d(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: select = 20

           ! Local Scalars
           integer(ilp) :: dir,endd,i,j,start,stkpnt
           real(sp) :: d1,d2,d3,dmnmx,tmp
           ! Local Arrays
           integer(ilp) :: stack(2,32)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           dir = -1
           if (la_lsame(id,'D')) then
              dir = 0
           else if (la_lsame(id,'I')) then
              dir = 1
           end if
           if (dir == -1) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('SLASRT',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           stkpnt = 1
           stack(1,1) = 1
           stack(2,1) = n
           10 continue
           start = stack(1,stkpnt)
           endd = stack(2,stkpnt)
           stkpnt = stkpnt - 1
           if (endd - start <= select .and. endd - start > 0) then
              ! do insertion sort on d( start:endd )
              if (dir == 0) then
                 ! sort into decreasing order
                 loop_30: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) > d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_30
                       end if
                    end do
                 end do loop_30
              else
                 ! sort into increasing order
                 loop_50: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) < d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_50
                       end if
                    end do
                 end do loop_50
              end if
           else if (endd - start > select) then
              ! partition d( start:endd ) and stack parts, largest one first
              ! choose partition entry as median of 3
              d1 = d(start)
              d2 = d(endd)
              i = (start + endd)/2
              d3 = d(i)
              if (d1 < d2) then
                 if (d3 < d1) then
                    dmnmx = d1
                 else if (d3 < d2) then
                    dmnmx = d3
                 else
                    dmnmx = d2
                 end if
              else
                 if (d3 < d2) then
                    dmnmx = d2
                 else if (d3 < d1) then
                    dmnmx = d3
                 else
                    dmnmx = d1
                 end if
              end if
              if (dir == 0) then
                 ! sort into decreasing order
                 i = start - 1
                 j = endd + 1
                 60 continue
                 70 continue
                 j = j - 1
                 if (d(j) < dmnmx) go to 70
                 80 continue
                 i = i + 1
                 if (d(i) > dmnmx) go to 80
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 60
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              else
                 ! sort into increasing order
                 i = start - 1
                 j = endd + 1
                 90 continue
                 100 continue
                 j = j - 1
                 if (d(j) > dmnmx) go to 100
                 110 continue
                 i = i + 1
                 if (d(i) < dmnmx) go to 110
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 90
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              end if
           end if
           if (stkpnt > 0) go to 10
           return
     end subroutine la_slasrt
     !> Sort the numbers in D in increasing order (if ID = 'I') or
     !> in decreasing order (if ID = 'D' ).
     !> Use Quick Sort, reverting to Insertion sort on arrays of
     !> size <= 20. Dimension of STACK limits N to about 2**32.

     pure subroutine la_dlasrt(id,n,d,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: id
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(dp),intent(inout) :: d(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: select = 20

           ! Local Scalars
           integer(ilp) :: dir,endd,i,j,start,stkpnt
           real(dp) :: d1,d2,d3,dmnmx,tmp
           ! Local Arrays
           integer(ilp) :: stack(2,32)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           dir = -1
           if (la_lsame(id,'D')) then
              dir = 0
           else if (la_lsame(id,'I')) then
              dir = 1
           end if
           if (dir == -1) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('DLASRT',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           stkpnt = 1
           stack(1,1) = 1
           stack(2,1) = n
           10 continue
           start = stack(1,stkpnt)
           endd = stack(2,stkpnt)
           stkpnt = stkpnt - 1
           if (endd - start <= select .and. endd - start > 0) then
              ! do insertion sort on d( start:endd )
              if (dir == 0) then
                 ! sort into decreasing order
                 loop_30: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) > d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_30
                       end if
                    end do
                 end do loop_30
              else
                 ! sort into increasing order
                 loop_50: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) < d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_50
                       end if
                    end do
                 end do loop_50
              end if
           else if (endd - start > select) then
              ! partition d( start:endd ) and stack parts, largest one first
              ! choose partition entry as median of 3
              d1 = d(start)
              d2 = d(endd)
              i = (start + endd)/2
              d3 = d(i)
              if (d1 < d2) then
                 if (d3 < d1) then
                    dmnmx = d1
                 else if (d3 < d2) then
                    dmnmx = d3
                 else
                    dmnmx = d2
                 end if
              else
                 if (d3 < d2) then
                    dmnmx = d2
                 else if (d3 < d1) then
                    dmnmx = d3
                 else
                    dmnmx = d1
                 end if
              end if
              if (dir == 0) then
                 ! sort into decreasing order
                 i = start - 1
                 j = endd + 1
                 60 continue
                 70 continue
                 j = j - 1
                 if (d(j) < dmnmx) go to 70
                 80 continue
                 i = i + 1
                 if (d(i) > dmnmx) go to 80
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 60
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              else
                 ! sort into increasing order
                 i = start - 1
                 j = endd + 1
                 90 continue
                 100 continue
                 j = j - 1
                 if (d(j) > dmnmx) go to 100
                 110 continue
                 i = i + 1
                 if (d(i) < dmnmx) go to 110
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 90
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              end if
           end if
           if (stkpnt > 0) go to 10
           return
     end subroutine la_dlasrt
#ifdef LA_WITH_XDP
     !> Sort the numbers in D in increasing order (if ID = 'I') or
     !> in decreasing order (if ID = 'D' ).
     !> Use Quick Sort, reverting to Insertion sort on arrays of
     !> size <= 20. Dimension of STACK limits N to about 2**32.

     pure subroutine la_xlasrt(id,n,d,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: id
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(xdp),intent(inout) :: d(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: select = 20

           ! Local Scalars
           integer(ilp) :: dir,endd,i,j,start,stkpnt
           real(xdp) :: d1,d2,d3,dmnmx,tmp
           ! Local Arrays
           integer(ilp) :: stack(2,32)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           dir = -1
           if (la_lsame(id,'D')) then
              dir = 0
           else if (la_lsame(id,'I')) then
              dir = 1
           end if
           if (dir == -1) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('XLASRT',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           stkpnt = 1
           stack(1,1) = 1
           stack(2,1) = n
           10 continue
           start = stack(1,stkpnt)
           endd = stack(2,stkpnt)
           stkpnt = stkpnt - 1
           if (endd - start <= select .and. endd - start > 0) then
              ! do insertion sort on d( start:endd )
              if (dir == 0) then
                 ! sort into decreasing order
                 loop_30: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) > d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_30
                       end if
                    end do
                 end do loop_30
              else
                 ! sort into increasing order
                 loop_50: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) < d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_50
                       end if
                    end do
                 end do loop_50
              end if
           else if (endd - start > select) then
              ! partition d( start:endd ) and stack parts, largest one first
              ! choose partition entry as median of 3
              d1 = d(start)
              d2 = d(endd)
              i = (start + endd)/2
              d3 = d(i)
              if (d1 < d2) then
                 if (d3 < d1) then
                    dmnmx = d1
                 else if (d3 < d2) then
                    dmnmx = d3
                 else
                    dmnmx = d2
                 end if
              else
                 if (d3 < d2) then
                    dmnmx = d2
                 else if (d3 < d1) then
                    dmnmx = d3
                 else
                    dmnmx = d1
                 end if
              end if
              if (dir == 0) then
                 ! sort into decreasing order
                 i = start - 1
                 j = endd + 1
                 60 continue
                 70 continue
                 j = j - 1
                 if (d(j) < dmnmx) go to 70
                 80 continue
                 i = i + 1
                 if (d(i) > dmnmx) go to 80
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 60
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              else
                 ! sort into increasing order
                 i = start - 1
                 j = endd + 1
                 90 continue
                 100 continue
                 j = j - 1
                 if (d(j) > dmnmx) go to 100
                 110 continue
                 i = i + 1
                 if (d(i) < dmnmx) go to 110
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 90
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              end if
           end if
           if (stkpnt > 0) go to 10
           return
     end subroutine la_xlasrt
#endif
#ifdef LA_WITH_QP
     !> Sort the numbers in D in increasing order (if ID = 'I') or
     !> in decreasing order (if ID = 'D' ).
     !> Use Quick Sort, reverting to Insertion sort on arrays of
     !> size <= 20. Dimension of STACK limits N to about 2**32.

     pure subroutine la_qlasrt(id,n,d,info)
        ! -- lapack computational routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           character,intent(in) :: id
           integer(ilp),intent(out) :: info
           integer(ilp),intent(in) :: n
           ! Array Arguments
           real(qp),intent(inout) :: d(*)
        ! =====================================================================
           ! Parameters
           integer(ilp),parameter :: select = 20

           ! Local Scalars
           integer(ilp) :: dir,endd,i,j,start,stkpnt
           real(qp) :: d1,d2,d3,dmnmx,tmp
           ! Local Arrays
           integer(ilp) :: stack(2,32)
           ! Executable Statements
           ! test the input parameters.
           info = 0
           dir = -1
           if (la_lsame(id,'D')) then
              dir = 0
           else if (la_lsame(id,'I')) then
              dir = 1
           end if
           if (dir == -1) then
              info = -1
           else if (n < 0) then
              info = -2
           end if
           if (info /= 0) then
              call la_xerbla('QLASRT',-info)
              return
           end if
           ! quick return if possible
           if (n <= 1) return
           stkpnt = 1
           stack(1,1) = 1
           stack(2,1) = n
           10 continue
           start = stack(1,stkpnt)
           endd = stack(2,stkpnt)
           stkpnt = stkpnt - 1
           if (endd - start <= select .and. endd - start > 0) then
              ! do insertion sort on d( start:endd )
              if (dir == 0) then
                 ! sort into decreasing order
                 loop_30: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) > d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_30
                       end if
                    end do
                 end do loop_30
              else
                 ! sort into increasing order
                 loop_50: do i = start + 1,endd
                    do j = i,start + 1,-1
                       if (d(j) < d(j - 1)) then
                          dmnmx = d(j)
                          d(j) = d(j - 1)
                          d(j - 1) = dmnmx
                       else
                          cycle loop_50
                       end if
                    end do
                 end do loop_50
              end if
           else if (endd - start > select) then
              ! partition d( start:endd ) and stack parts, largest one first
              ! choose partition entry as median of 3
              d1 = d(start)
              d2 = d(endd)
              i = (start + endd)/2
              d3 = d(i)
              if (d1 < d2) then
                 if (d3 < d1) then
                    dmnmx = d1
                 else if (d3 < d2) then
                    dmnmx = d3
                 else
                    dmnmx = d2
                 end if
              else
                 if (d3 < d2) then
                    dmnmx = d2
                 else if (d3 < d1) then
                    dmnmx = d3
                 else
                    dmnmx = d1
                 end if
              end if
              if (dir == 0) then
                 ! sort into decreasing order
                 i = start - 1
                 j = endd + 1
                 60 continue
                 70 continue
                 j = j - 1
                 if (d(j) < dmnmx) go to 70
                 80 continue
                 i = i + 1
                 if (d(i) > dmnmx) go to 80
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 60
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              else
                 ! sort into increasing order
                 i = start - 1
                 j = endd + 1
                 90 continue
                 100 continue
                 j = j - 1
                 if (d(j) > dmnmx) go to 100
                 110 continue
                 i = i + 1
                 if (d(i) < dmnmx) go to 110
                 if (i < j) then
                    tmp = d(i)
                    d(i) = d(j)
                    d(j) = tmp
                    go to 90
                 end if
                 if (j - start > endd - j - 1) then
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                 else
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                 end if
              end if
           end if
           if (stkpnt > 0) go to 10
           return
     end subroutine la_qlasrt
#endif

     !> !
     !>
     !> SLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_slassq(n,x,incx,scl,sumsq)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(sp),intent(inout) :: scl,sumsq
        ! Array Arguments
        real(sp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(sp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_slassq
     !> !
     !>
     !> DLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_dlassq(n,x,incx,scl,sumsq)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(dp),intent(inout) :: scl,sumsq
        ! Array Arguments
        real(dp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(dp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_dlassq
#ifdef LA_WITH_XDP
     !> !
     !>
     !> XLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_xlassq(n,x,incx,scl,sumsq)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(xdp),intent(inout) :: scl,sumsq
        ! Array Arguments
        real(xdp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(xdp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_xlassq
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> QLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_qlassq(n,x,incx,scl,sumsq)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(qp),intent(inout) :: scl,sumsq
        ! Array Arguments
        real(qp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(qp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(x(ix))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_qlassq
#endif

     !> SRSCL: multiplies an n-element real vector x by the real scalar 1/a.
     !> This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_srscl(n,sa,sx,incx)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(sp),intent(in) :: sa
           ! Array Arguments
           real(sp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(sp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_sscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_srscl
     !> DRSCL: multiplies an n-element real vector x by the real scalar 1/a.
     !> This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_drscl(n,sa,sx,incx)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(dp),intent(in) :: sa
           ! Array Arguments
           real(dp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(dp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_dscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_drscl
#ifdef LA_WITH_XDP
     !> XRSCL: multiplies an n-element real vector x by the real scalar 1/a.
     !> This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_xrscl(n,sa,sx,incx)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(xdp),intent(in) :: sa
           ! Array Arguments
           real(xdp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(xdp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_xlamch('S')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_xscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_xrscl
#endif
#ifdef LA_WITH_QP
     !> QRSCL: multiplies an n-element real vector x by the real scalar 1/a.
     !> This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_qrscl(n,sa,sx,incx)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(qp),intent(in) :: sa
           ! Array Arguments
           real(qp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(qp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_qscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_qrscl
#endif

     !> CSRSCL: multiplies an n-element complex vector x by the real scalar
     !> 1/a.  This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_csrscl(n,sa,sx,incx)
        use la_constants_sp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(sp),intent(in) :: sa
           ! Array Arguments
           complex(sp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(sp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_slamch('S')
           bignum = one/smlnum
           call la_slabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_csscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_csrscl
     !> ZDRSCL: multiplies an n-element complex vector x by the real scalar
     !> 1/a.  This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_zdrscl(n,sa,sx,incx)
        use la_constants_dp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(dp),intent(in) :: sa
           ! Array Arguments
           complex(dp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(dp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_dlamch('S')
           bignum = one/smlnum
           call la_dlabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_zdscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_zdrscl
#ifdef LA_WITH_XDP
     !> YXRSCL: multiplies an n-element complex vector x by the real scalar
     !> 1/a.  This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_yxrscl(n,sa,sx,incx)
        use la_constants_xdp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(xdp),intent(in) :: sa
           ! Array Arguments
           complex(xdp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(xdp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_xlamch('S')
           bignum = one/smlnum
           call la_xlabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_yxscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_yxrscl
#endif
#ifdef LA_WITH_QP
     !> WQRSCL: multiplies an n-element complex vector x by the real scalar
     !> 1/a.  This is done without overflow or underflow as long as
     !> the final result x/a does not overflow or underflow.

     pure subroutine la_wqrscl(n,sa,sx,incx)
        use la_constants_qp,only:zero,one
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           real(qp),intent(in) :: sa
           ! Array Arguments
           complex(qp),intent(inout) :: sx(*)
       ! =====================================================================

           ! Local Scalars
           logical(lk) :: done
           real(qp) :: bignum,cden,cden1,cnum,cnum1,mul,smlnum
           ! Intrinsic Functions
           intrinsic :: abs
           ! Executable Statements
           ! quick return if possible
           if (n <= 0) return
           ! get machine parameters
           smlnum = la_qlamch('S')
           bignum = one/smlnum
           call la_qlabad(smlnum,bignum)
           ! initialize the denominator to sa and the numerator to 1.
           cden = sa
           cnum = one
           10 continue
           cden1 = cden*smlnum
           cnum1 = cnum/bignum
           if (abs(cden1) > abs(cnum) .and. cnum /= zero) then
              ! pre-multiply x by smlnum if cden is large compared to cnum.
              mul = smlnum
              done = .false.
              cden = cden1
           else if (abs(cnum1) > abs(cden)) then
              ! pre-multiply x by bignum if cden is small compared to cnum.
              mul = bignum
              done = .false.
              cnum = cnum1
           else
              ! multiply x by cnum / cden and return.
              mul = cnum/cden
              done = .true.
           end if
           ! scale the vector x by mul
           call la_wqscal(n,mul,sx,incx)
           if (.not. done) go to 10
           return
     end subroutine la_wqrscl
#endif

     !> CLACGV: conjugates a complex vector of length N.

     pure subroutine la_clacgv(n,x,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(sp),intent(inout) :: x(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ioff
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (incx == 1) then
              do i = 1,n
                 x(i) = conjg(x(i))
              end do
           else
              ioff = 1
              if (incx < 0) ioff = 1 - (n - 1)*incx
              do i = 1,n
                 x(ioff) = conjg(x(ioff))
                 ioff = ioff + incx
              end do
           end if
           return
     end subroutine la_clacgv
     !> ZLACGV: conjugates a complex vector of length N.

     pure subroutine la_zlacgv(n,x,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(dp),intent(inout) :: x(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ioff
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (incx == 1) then
              do i = 1,n
                 x(i) = conjg(x(i))
              end do
           else
              ioff = 1
              if (incx < 0) ioff = 1 - (n - 1)*incx
              do i = 1,n
                 x(ioff) = conjg(x(ioff))
                 ioff = ioff + incx
              end do
           end if
           return
     end subroutine la_zlacgv
#ifdef LA_WITH_XDP
     !> YLACGV: conjugates a complex vector of length N.

     pure subroutine la_ylacgv(n,x,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(xdp),intent(inout) :: x(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ioff
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (incx == 1) then
              do i = 1,n
                 x(i) = conjg(x(i))
              end do
           else
              ioff = 1
              if (incx < 0) ioff = 1 - (n - 1)*incx
              do i = 1,n
                 x(ioff) = conjg(x(ioff))
                 ioff = ioff + incx
              end do
           end if
           return
     end subroutine la_ylacgv
#endif
#ifdef LA_WITH_QP
     !> WLACGV: conjugates a complex vector of length N.

     pure subroutine la_wlacgv(n,x,incx)
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
           ! Scalar Arguments
           integer(ilp),intent(in) :: incx,n
           ! Array Arguments
           complex(qp),intent(inout) :: x(*)
       ! =====================================================================
           ! Local Scalars
           integer(ilp) :: i,ioff
           ! Intrinsic Functions
           intrinsic :: conjg
           ! Executable Statements
           if (incx == 1) then
              do i = 1,n
                 x(i) = conjg(x(i))
              end do
           else
              ioff = 1
              if (incx < 0) ioff = 1 - (n - 1)*incx
              do i = 1,n
                 x(ioff) = conjg(x(ioff))
                 ioff = ioff + incx
              end do
           end if
           return
     end subroutine la_wlacgv
#endif

     !> !
     !>
     !> CLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_classq(n,x,incx,scl,sumsq)
        use la_constants_sp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(sp),intent(inout) :: scl,sumsq
        ! Array Arguments
        complex(sp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(sp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=sp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_classq
     !> !
     !>
     !> ZLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_zlassq(n,x,incx,scl,sumsq)
        use la_constants_dp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(dp),intent(inout) :: scl,sumsq
        ! Array Arguments
        complex(dp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(dp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=dp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_zlassq
#ifdef LA_WITH_XDP
     !> !
     !>
     !> YLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_ylassq(n,x,incx,scl,sumsq)
        use la_constants_xdp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(xdp),intent(inout) :: scl,sumsq
        ! Array Arguments
        complex(xdp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(xdp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=xdp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_ylassq
#endif
#ifdef LA_WITH_QP
     !> !
     !>
     !> WLASSQ:  returns the values  scl  and  smsq  such that
     !> ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
     !> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
     !> assumed to be non-negative.
     !> scale and sumsq must be supplied in SCALE and SUMSQ and
     !> scl and smsq are overwritten on SCALE and SUMSQ respectively.
     !> If scale * sqrt( sumsq ) > tbig then
     !> we require:   scale >= sqrt( TINY*EPS ) / sbig   on entry,
     !> and if 0 < scale * sqrt( sumsq ) < tsml then
     !> we require:   scale <= sqrt( HUGE ) / ssml       on entry,
     !> where
     !> tbig -- upper threshold for values whose square is representable;
     !> sbig -- scaling constant for big numbers; \see la_constants.f90
     !> tsml -- lower threshold for values whose square is representable;
     !> ssml -- scaling constant for small numbers; \see la_constants.f90
     !> and
     !> TINY*EPS -- tiniest representable number;
     !> HUGE     -- biggest representable number.

     pure subroutine la_wlassq(n,x,incx,scl,sumsq)
        use la_constants_qp
        ! -- lapack auxiliary routine --
        ! -- lapack is a software package provided by univ. of tennessee,    --
        ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
        ! Scalar Arguments
     integer(ilp),intent(in) :: incx,n
        real(qp),intent(inout) :: scl,sumsq
        ! Array Arguments
        complex(qp),intent(in) :: x(*)
        ! Local Scalars
     integer(ilp) :: i,ix
     logical(lk) :: notbig
        real(qp) :: abig,amed,asml,ax,ymax,ymin
        ! quick return if possible
        if (ieee_is_nan(scl) .or. ieee_is_nan(sumsq)) return
        if (sumsq == zero) scl = one
        if (scl == zero) then
           scl = one
           sumsq = zero
        end if
        if (n <= 0) then
           return
        end if
        ! compute the sum of squares in 3 accumulators:
           ! abig -- sums of squares scaled down to avoid overflow
           ! asml -- sums of squares scaled up to avoid underflow
           ! amed -- sums of squares that do not require scaling
        ! the thresholds and multipliers are
           ! tbig -- values bigger than this are scaled down by sbig
           ! tsml -- values smaller than this are scaled up by ssml
        notbig = .true.
        asml = zero
        amed = zero
        abig = zero
        ix = 1
        if (incx < 0) ix = 1 - (n - 1)*incx
        do i = 1,n
           ax = abs(real(x(ix),KIND=qp))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ax = abs(aimag(x(ix)))
           if (ax > tbig) then
              abig = abig + (ax*sbig)**2
              notbig = .false.
           else if (ax < tsml) then
              if (notbig) asml = asml + (ax*ssml)**2
           else
              amed = amed + ax**2
           end if
           ix = ix + incx
        end do
        ! put the existing sum of squares into one of the accumulators
        if (sumsq > zero) then
           ax = scl*sqrt(sumsq)
           if (ax > tbig) then
              ! we assume scl >= sqrt( tiny*eps ) / sbig
              abig = abig + (scl*sbig)**2*sumsq
           else if (ax < tsml) then
              ! we assume scl <= sqrt( huge ) / ssml
              if (notbig) asml = asml + (scl*ssml)**2*sumsq
           else
              amed = amed + scl**2*sumsq
           end if
        end if
        ! combine abig and amed or amed and asml if more than one
        ! accumulator was used.
        if (abig > zero) then
           ! combine abig and amed if abig > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              abig = abig + (amed*sbig)*sbig
           end if
           scl = one/sbig
           sumsq = abig
        else if (asml > zero) then
           ! combine amed and asml if asml > 0.
           if (amed > zero .or. ieee_is_nan(amed)) then
              amed = sqrt(amed)
              asml = sqrt(asml)/ssml
              if (asml > amed) then
                 ymin = amed
                 ymax = asml
              else
                 ymin = asml
                 ymax = amed
              end if
              scl = one
              sumsq = ymax**2*(one + (ymin/ymax)**2)
           else
              scl = one/ssml
              sumsq = asml
           end if
        else
           ! otherwise all values are mid-range or zero
           scl = one
           sumsq = amed
        end if
        return
     end subroutine la_wlassq
#endif

end module la_lapack_blas_like_l1
